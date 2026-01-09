/**
 * @file test_integration.cpp
 * @brief Test suite for the Integration module
 * 
 * This file contains tests for numerical integration routines,
 * including Gauss quadrature rules for various element types.
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../src/Integration.h"

using namespace elmer;

// Utility function to check if two numbers are close within tolerance
bool isClose(double a, double b, double tolerance = 1e-12) {
    return std::abs(a - b) < tolerance;
}

// Test 1D Gauss quadrature points and weights
void testGauss1DPointsWeights() {
    std::cout << "=== Testing 1D Gauss Quadrature Points and Weights ===" << std::endl;
    
    // Test 1-point Gauss quadrature
    auto points1 = Integration::getGaussPoints1D(1);
    auto weights1 = Integration::getGaussWeights1D(1);
    
    if (points1.size() != 1 || weights1.size() != 1) {
        std::cout << "FAILED: 1-point quadrature size mismatch" << std::endl;
        return;
    }
    
    if (!isClose(points1[0], 0.0) || !isClose(weights1[0], 2.0)) {
        std::cout << "FAILED: 1-point quadrature values incorrect" << std::endl;
        return;
    }
    
    std::cout << "PASSED: 1-point Gauss quadrature" << std::endl;
    
    // Test 2-point Gauss quadrature
    auto points2 = Integration::getGaussPoints1D(2);
    auto weights2 = Integration::getGaussWeights1D(2);
    
    if (points2.size() != 2 || weights2.size() != 2) {
        std::cout << "FAILED: 2-point quadrature size mismatch" << std::endl;
        return;
    }
    
    double expected_points2[] = {-0.5773502691896257, 0.5773502691896257};
    double expected_weights2[] = {1.0, 1.0};
    
    bool passed = true;
    for (int i = 0; i < 2; ++i) {
        if (!isClose(points2[i], expected_points2[i]) || !isClose(weights2[i], expected_weights2[i])) {
            passed = false;
            break;
        }
    }
    
    if (!passed) {
        std::cout << "FAILED: 2-point quadrature values incorrect" << std::endl;
        return;
    }
    
    std::cout << "PASSED: 2-point Gauss quadrature" << std::endl;
}

// Test numerical integration accuracy
void testIntegrationAccuracy() {
    std::cout << "=== Testing Numerical Integration Accuracy ===" << std::endl;
    
    // Test integration of polynomial x^2 over [-1, 1]
    // ∫x^2 dx from -1 to 1 = 2/3 ≈ 0.6666666666666666
    
    auto points = Integration::getGaussPoints1D(3);
    auto weights = Integration::getGaussWeights1D(3);
    
    double integral = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        double x = points[i];
        integral += weights[i] * (x * x);
    }
    
    if (!isClose(integral, 2.0/3.0)) {
        std::cout << "FAILED: Integration of x^2 incorrect. Got: " << integral 
                  << ", Expected: " << (2.0/3.0) << std::endl;
        return;
    }
    
    std::cout << "PASSED: Integration of x^2 polynomial" << std::endl;
}

// Test triangle integration rules
void testTriangleIntegration() {
    std::cout << "=== Testing Triangle Integration Rules ===" << std::endl;
    
    // Test order 1 triangle rule
    auto rule1 = Integration::getTriangleRule(1);
    
    if (rule1.u.size() != 1 || rule1.v.size() != 1 || rule1.weights.size() != 1) {
        std::cout << "FAILED: Triangle order 1 rule size mismatch" << std::endl;
        return;
    }
    
    // For order 1, should have single point at barycenter with weight 1.0
    if (!isClose(rule1.u[0], 1.0/3.0) || !isClose(rule1.v[0], 1.0/3.0) || 
        !isClose(rule1.weights[0], 1.0)) {
        std::cout << "FAILED: Triangle order 1 rule values incorrect" << std::endl;
        return;
    }
    
    std::cout << "PASSED: Triangle order 1 integration rule" << std::endl;
}

// Test error handling for invalid orders
void testErrorHandling() {
    std::cout << "=== Testing Error Handling ===" << std::endl;
    
    try {
        auto points = Integration::getGaussPoints1D(0);
        std::cout << "FAILED: Should have thrown exception for order 0" << std::endl;
        return;
    } catch (const std::exception& e) {
        std::cout << "PASSED: Correctly threw exception for invalid order" << std::endl;
    }
    
    try {
        auto points = Integration::getGaussPoints1D(Integration::MAXN + 1);
        std::cout << "FAILED: Should have thrown exception for order > MAXN" << std::endl;
        return;
    } catch (const std::exception& e) {
        std::cout << "PASSED: Correctly threw exception for order > MAXN" << std::endl;
    }
}

int main() {
    std::cout << "Starting Integration Module Tests..." << std::endl;
    std::cout << "====================================" << std::endl;
    
    testGauss1DPointsWeights();
    testIntegrationAccuracy();
    testTriangleIntegration();
    testErrorHandling();
    
    std::cout << "====================================" << std::endl;
    std::cout << "All tests completed!" << std::endl;
    
    return 0;
}