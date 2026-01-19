/**
 * @file Integration.cpp
 * @brief Implementation of numerical integration routines
 * 
 * This file implements the Gauss quadrature and integration rules
 * for finite element methods, following the original Fortran implementation.
 */

#include "Integration.h"
#include "CommonConstants.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>

namespace elmer {

// Static member initialization
bool Integration::initialized_ = false;
std::array<std::vector<double>, Integration::MAXN + 1> Integration::gaussPoints1D_;
std::array<std::vector<double>, Integration::MAXN + 1> Integration::gaussWeights1D_;

std::vector<double> Integration::getGaussPoints1D(int n) {
    validateOrder(n, MAXN);
    if (!initialized_) {
        initialize();
    }
    return gaussPoints1D_[n];
}

std::vector<double> Integration::getGaussWeights1D(int n) {
    validateOrder(n, MAXN);
    if (!initialized_) {
        initialize();
    }
    return gaussWeights1D_[n];
}

TriangleRule Integration::getTriangleRule(int order) {
    TriangleRule rule;
    
    switch (order) {
        case 1: // 1-point rule; exact for polynomials of degree <= 1
            rule.u = {0.3333333333333333};
            rule.v = {0.3333333333333333};
            rule.weights = {1.0};
            break;
            
        case 3: // 3-point rule; exact for polynomials of degree <= 2
            rule.u = {0.16666666666666667, 0.66666666666666667, 0.16666666666666667};
            rule.v = {0.16666666666666667, 0.16666666666666667, 0.66666666666666667};
            rule.weights = {0.33333333333333333, 0.33333333333333333, 0.33333333333333333};
            break;
            
        case 4: // 4-point rule; exact for polynomials of degree <= 3
            rule.u = {0.33333333333333333, 0.2000000000000000, 0.60000000000000000, 0.20000000000000000};
            rule.v = {0.33333333333333333, 0.2000000000000000, 0.20000000000000000, 0.60000000000000000};
            rule.weights = {-0.56250000000000000, 0.52083333333333333, 0.52083333333333333, 0.52083333333333333};
            break;
            
        case 6: // 6-point rule; exact for polynomials of degree <= 4
            rule.u = {0.091576213509771, 0.816847572980459, 0.091576213509771, 
                     0.445948490915965, 0.108103018168070, 0.445948490915965};
            rule.v = {0.091576213509771, 0.091576213509771, 0.816847572980459, 
                     0.445948490915965, 0.445948490915965, 0.108103018168070};
            rule.weights = {0.109951743655322, 0.109951743655322, 0.109951743655322, 
                          0.223381589678011, 0.223381589678011, 0.223381589678011};
            break;
            
        case 7: // 7-point rule; exact for polynomials of degree <= 5
            rule.u = {0.333333333333333, 0.101286507323456, 0.797426985353087, 
                     0.101286507323456, 0.470142064105115, 0.059715871789770, 0.470142064105115};
            rule.v = {0.333333333333333, 0.101286507323456, 0.101286507323456, 
                     0.797426985353087, 0.470142064105115, 0.470142064105115, 0.059715871789770};
            rule.weights = {0.225000000000000, 0.125939180544827, 0.125939180544827, 
                          0.125939180544827, 0.132394152788506, 0.132394152788506, 0.132394152788506};
            break;
            
        case 11: // 11-point rule; exact for polynomials of degree <= 6
            rule.u = {0.063089014491502, 0.873821971016996, 0.063089014491502, 
                     0.249286745170910, 0.501426509658179, 0.249286745170910, 
                     0.053145049844817, 0.310352451033784, 0.636502499121399, 
                     0.310352451033784, 0.053145049844817};
            rule.v = {0.063089014491502, 0.063089014491502, 0.873821971016996, 
                     0.249286745170910, 0.249286745170910, 0.501426509658179, 
                     0.310352451033784, 0.636502499121399, 0.053145049844817, 
                     0.053145049844817, 0.310352451033784};
            rule.weights = {0.050844906370207, 0.050844906370207, 0.050844906370207, 
                          0.116786275726379, 0.116786275726379, 0.116786275726379, 
                          0.082851075618374, 0.082851075618374, 0.082851075618374, 
                          0.082851075618374, 0.082851075618374};
            break;
            
        default:
            throw std::invalid_argument("Unsupported triangle integration order: " + std::to_string(order));
    }
    
    return rule;
}

std::pair<std::vector<double>, std::vector<double>> Integration::getQuadrilateralRule(int order) {
    validateOrder(order, MAXN);
    if (!initialized_) {
        initialize();
    }
    
    // For quadrilaterals, use tensor product of 1D Gauss rules
    auto points1D = getGaussPoints1D(order);
    auto weights1D = getGaussWeights1D(order);
    
    std::vector<double> pointsU, pointsV, weights;
    
    for (size_t i = 0; i < points1D.size(); ++i) {
        for (size_t j = 0; j < points1D.size(); ++j) {
            pointsU.push_back(points1D[i]);
            pointsV.push_back(points1D[j]);
            weights.push_back(weights1D[i] * weights1D[j]);
        }
    }
    
    return std::make_pair(pointsU, weights);
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, 
           std::vector<double>> Integration::getBrickRule(int order) {
    validateOrder(order, MAXN);
    if (!initialized_) {
        initialize();
    }
    
    // For bricks, use tensor product of 1D Gauss rules
    auto points1D = getGaussPoints1D(order);
    auto weights1D = getGaussWeights1D(order);
    
    std::vector<double> pointsU, pointsV, pointsW, weights;
    
    for (size_t i = 0; i < points1D.size(); ++i) {
        for (size_t j = 0; j < points1D.size(); ++j) {
            for (size_t k = 0; k < points1D.size(); ++k) {
                pointsU.push_back(points1D[i]);
                pointsV.push_back(points1D[j]);
                pointsW.push_back(points1D[k]);
                weights.push_back(weights1D[i] * weights1D[j] * weights1D[k]);
            }
        }
    }
    
    return {pointsU, pointsV, pointsW, weights};
}

void Integration::computeGaussPoints1D(int n, double a, double b, 
                                      std::vector<double>& points, 
                                      std::vector<double>& weights) {
    if (n <= 0) {
        throw std::invalid_argument("Number of Gauss points must be positive");
    }
    
    points.resize(n);
    weights.resize(n);
    
    if (n == 1) {
        // Single point at the midpoint
        points[0] = (a + b) / 2.0;
        weights[0] = b - a;
        return;
    }
    
    // Use Newton's method to find roots of Legendre polynomials
    const double tolerance = 1e-15;
    const int maxIterations = 100;
    
    // Initial guesses: Chebyshev nodes
    for (int i = 0; i < n; ++i) {
        double x0 = std::cos(M_PI * (i + 0.5) / n);
        
        // Newton iteration
        for (int iter = 0; iter < maxIterations; ++iter) {
            double p_n, p_prime;
            legendrePolynomial(n, x0, p_n, p_prime);
            
            double delta = p_n / p_prime;
            x0 -= delta;
            
            if (std::abs(delta) < tolerance) {
                break;
            }
        }
        
        points[i] = x0;
        
        // Compute weight using formula: w_i = 2 / [(1 - x_i^2) * (P_n'(x_i))^2]
        double p_n, p_prime;
        legendrePolynomial(n, x0, p_n, p_prime);
        weights[i] = 2.0 / ((1.0 - x0 * x0) * p_prime * p_prime);
    }
    
    // Sort points in ascending order
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (points[i] > points[j]) {
                std::swap(points[i], points[j]);
                std::swap(weights[i], weights[j]);
            }
        }
    }
    
    // Transform from [-1, 1] to [a, b]
    double scale = (b - a) / 2.0;
    double shift = (a + b) / 2.0;
    
    for (int i = 0; i < n; ++i) {
        points[i] = scale * points[i] + shift;
        weights[i] *= scale;
    }
}

void Integration::initialize() {
    if (initialized_) {
        return;
    }
    
    initializeGauss1D();
    initialized_ = true;
}

void Integration::initializeGauss1D() {
    // Precompute Gauss points and weights for orders 1 to MAXN
    for (int n = 1; n <= MAXN; ++n) {
        std::vector<double> points, weights;
        computeGaussPoints1D(n, -1.0, 1.0, points, weights);
        gaussPoints1D_[n] = points;
        gaussWeights1D_[n] = weights;
    }
}

void Integration::validateOrder(int n, int max) {
    if (n < 1 || n > max) {
        throw std::invalid_argument("Integration order " + std::to_string(n) + 
                                   " is not supported. Must be between 1 and " + 
                                   std::to_string(max));
    }
}

void Integration::legendrePolynomial(int n, double x, double& p_n, double& p_prime) {
    if (n == 0) {
        p_n = 1.0;
        p_prime = 0.0;
        return;
    }
    
    if (n == 1) {
        p_n = x;
        p_prime = 1.0;
        return;
    }
    
    // Use recurrence relation: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
    double p_prev = 1.0;      // P_0(x)
    double p_curr = x;         // P_1(x)
    double p_prev_prime = 0.0; // P_0'(x)
    double p_curr_prime = 1.0; // P_1'(x)
    
    for (int k = 1; k < n; ++k) {
        double p_next = ((2.0 * k + 1.0) * x * p_curr - k * p_prev) / (k + 1.0);
        double p_next_prime = ((2.0 * k + 1.0) * (p_curr + x * p_curr_prime) - k * p_prev_prime) / (k + 1.0);
        
        p_prev = p_curr;
        p_curr = p_next;
        p_prev_prime = p_curr_prime;
        p_curr_prime = p_next_prime;
    }
    
    p_n = p_curr;
    p_prime = p_curr_prime;
}

} // namespace elmer