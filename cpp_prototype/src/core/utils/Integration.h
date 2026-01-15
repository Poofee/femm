/**
 * @file Integration.h
 * @brief Numerical integration routines for finite element methods
 * 
 * This module provides Gauss quadrature points and weights for various
 * element types, including 1D, triangles, quadrilaterals, and bricks.
 * It implements the same functionality as the original Fortran Integration.F90 module.
 */

#pragma once

#include <vector>
#include <array>
#include <stdexcept>

namespace elmer {

/**
 * @brief Structure to hold triangle integration rule data
 */
struct TriangleRule {
    std::vector<double> u;        ///< First barycentric coordinate
    std::vector<double> v;        ///< Second barycentric coordinate
    std::vector<double> weights;  ///< Integration weights
};

/**
 * @brief Numerical integration utilities
 * 
 * This class provides static methods for computing Gauss quadrature
 * points and weights for various element types and integration orders.
 */
class Integration {
public:
    // Maximum number of integration points supported
    static constexpr int MAXN = 13;
    static constexpr int MAXNPAD = 16; // Padded for 64-byte alignment
    static constexpr int MAX_INTEGRATION_POINTS = MAXN * MAXN * MAXN;
    
    /**
     * @brief Get 1D Gauss quadrature points for given order
     * @param n Number of integration points (1-10 supported)
     * @return Vector of Gauss points in [-1, 1]
     */
    static std::vector<double> getGaussPoints1D(int n);
    
    /**
     * @brief Get 1D Gauss quadrature weights for given order
     * @param n Number of integration points (1-10 supported)
     * @return Vector of Gauss weights
     */
    static std::vector<double> getGaussWeights1D(int n);
    
    /**
     * @brief Get triangle integration rule for given order
     * @param order Integration order (1, 3, 4, 6, 7, 11 supported)
     * @return TriangleRule structure with points and weights
     */
    static TriangleRule getTriangleRule(int order);
    
    /**
     * @brief Get quadrilateral integration rule for given order
     * @param order Integration order (1-10 supported)
     * @return Quadrilateral rule with points and weights
     */
    static std::pair<std::vector<double>, std::vector<double>> getQuadrilateralRule(int order);
    
    /**
     * @brief Get brick (hexahedron) integration rule for given order
     * @param order Integration order (1-10 supported)
     * @return Brick rule with points and weights
     */
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, 
                      std::vector<double>> getBrickRule(int order);
    
    /**
     * @brief Compute Gauss quadrature points and weights on the fly
     * 
     * This method computes Gauss quadrature points and weights using
     * the Golub-Welsch algorithm, which finds the roots of Legendre polynomials.
     * 
     * @param n Number of points
     * @param a Lower limit of integration
     * @param b Upper limit of integration
     * @param points Output vector for points
     * @param weights Output vector for weights
     */
    static void computeGaussPoints1D(int n, double a, double b, 
                                    std::vector<double>& points, 
                                    std::vector<double>& weights);
    
    /**
     * @brief Initialize integration module if needed
     * 
     * This method ensures that the integration points and weights
     * are computed and cached for efficient access.
     */
    static void initialize();
    
private:
    // Private constructor to prevent instantiation
    Integration() = delete;
    
    // Static initialization flag
    static bool initialized_;
    
    // Cache for 1D Gauss points and weights
    static std::array<std::vector<double>, MAXN + 1> gaussPoints1D_;
    static std::array<std::vector<double>, MAXN + 1> gaussWeights1D_;
    
    /**
     * @brief Initialize 1D Gauss quadrature data
     */
    static void initializeGauss1D();
    
    /**
     * @brief Validate integration order
     * @param n Order to validate
     * @param max Maximum allowed order
     */
    static void validateOrder(int n, int max);
    
    /**
     * @brief Compute Legendre polynomial and its derivative
     * @param n Polynomial degree
     * @param x Evaluation point
     * @param p_n Polynomial value
     * @param p_prime Derivative value
     */
    static void legendrePolynomial(int n, double x, double& p_n, double& p_prime);
};

} // namespace elmer