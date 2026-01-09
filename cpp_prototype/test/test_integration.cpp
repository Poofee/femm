#include "../src/Integration.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace elmer;

void testGaussPoints1D() {
    std::cout << "=== Testing 1D Gauss Points ===" << std::endl;
    
    // Test order 1
    auto points1 = Integration::getGaussPoints1D(1);
    std::cout << "Order 1: " << points1.size() << " points" << std::endl;
    for (const auto& p : points1) {
        std::cout << "  Point: (" << p.u << "), Weight: " << p.weight << std::endl;
    }
    
    // Test order 2
    auto points2 = Integration::getGaussPoints1D(2);
    std::cout << "Order 2: " << points2.size() << " points" << std::endl;
    for (const auto& p : points2) {
        std::cout << "  Point: (" << p.u << "), Weight: " << p.weight << std::endl;
    }
    
    // Test order 3
    auto points3 = Integration::getGaussPoints1D(3);
    std::cout << "Order 3: " << points3.size() << " points" << std::endl;
    for (const auto& p : points3) {
        std::cout << "  Point: (" << p.u << "), Weight: " << p.weight << std::endl;
    }
}

void testGaussPoints2D() {
    std::cout << "\n=== Testing 2D Gauss Points ===" << std::endl;
    
    auto points = Integration::getGaussPoints2D(2, 2);
    std::cout << "2x2 Gauss points: " << points.size() << " points" << std::endl;
    
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        std::cout << "  Point " << i << ": (" << p.u << ", " << p.v 
                  << "), Weight: " << p.weight << std::endl;
    }
}

void testGaussPoints3D() {
    std::cout << "\n=== Testing 3D Gauss Points ===" << std::endl;
    
    auto points = Integration::getGaussPoints3D(2, 2, 2);
    std::cout << "2x2x2 Gauss points: " << points.size() << " points" << std::endl;
    
    for (size_t i = 0; i < std::min(points.size(), size_t(4)); ++i) {
        const auto& p = points[i];
        std::cout << "  Point " << i << ": (" << p.u << ", " << p.v 
                  << ", " << p.w << "), Weight: " << p.weight << std::endl;
    }
    std::cout << "  ... and " << (points.size() - 4) << " more points" << std::endl;
}

void testTrianglePoints() {
    std::cout << "\n=== Testing Triangle Integration Points ===" << std::endl;
    
    // Test order 1
    auto points1 = Integration::getTrianglePoints(1);
    std::cout << "Triangle order 1: " << points1.size() << " points" << std::endl;
    for (const auto& p : points1) {
        std::cout << "  Point: (" << p.u << ", " << p.v 
                  << "), Weight: " << p.weight << std::endl;
    }
    
    // Test order 2
    auto points2 = Integration::getTrianglePoints(2);
    std::cout << "Triangle order 2: " << points2.size() << " points" << std::endl;
    for (const auto& p : points2) {
        std::cout << "  Point: (" << p.u << ", " << p.v 
                  << "), Weight: " << p.weight << std::endl;
    }
}

void testTetrahedronPoints() {
    std::cout << "\n=== Testing Tetrahedron Integration Points ===" << std::endl;
    
    // Test order 1
    auto points1 = Integration::getTetrahedronPoints(1);
    std::cout << "Tetrahedron order 1: " << points1.size() << " points" << std::endl;
    for (const auto& p : points1) {
        std::cout << "  Point: (" << p.u << ", " << p.v 
                  << ", " << p.w << "), Weight: " << p.weight << std::endl;
    }
    
    // Test order 2
    auto points2 = Integration::getTetrahedronPoints(2);
    std::cout << "Tetrahedron order 2: " << points2.size() << " points" << std::endl;
    for (const auto& p : points2) {
        std::cout << "  Point: (" << p.u << ", " << p.v 
                  << ", " << p.w << "), Weight: " << p.weight << std::endl;
    }
}

void testShapeFunctions() {
    std::cout << "\n=== Testing Shape Functions ===" << std::endl;
    
    // Test 1D linear shape functions
    std::cout << "1D Linear Shape Functions at xi=0.5:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        double N = Integration::shapeFunction1DLinear(0.5, i);
        std::cout << "  N" << i << " = " << N << std::endl;
    }
    
    // Test 2D linear shape functions
    std::cout << "2D Linear Shape Functions at (0.5, 0.5):" << std::endl;
    for (int i = 0; i < 4; ++i) {
        double N = Integration::shapeFunction2DLinear(0.5, 0.5, i);
        std::cout << "  N" << i << " = " << N << std::endl;
    }
    
    // Test 3D linear shape functions
    std::cout << "3D Linear Shape Functions at (0.5, 0.5, 0.5):" << std::endl;
    for (int i = 0; i < 8; ++i) {
        double N = Integration::shapeFunction3DLinear(0.5, 0.5, 0.5, i);
        std::cout << "  N" << i << " = " << N << std::endl;
    }
}

void testShapeFunctionDerivatives() {
    std::cout << "\n=== Testing Shape Function Derivatives ===" << std::endl;
    
    // Test 1D linear shape function derivatives
    std::cout << "1D Linear Shape Function Derivatives:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        double dN = Integration::shapeFunctionDerivative1DLinear(0.5, i);
        std::cout << "  dN" << i << "/dxi = " << dN << std::endl;
    }
    
    // Test 2D linear shape function derivatives
    std::cout << "2D Linear Shape Function Derivatives at (0.5, 0.5):" << std::endl;
    for (int i = 0; i < 4; ++i) {
        double dNdxi = Integration::shapeFunctionDerivative2DLinear(0.5, 0.5, i, 0);
        double dNdeta = Integration::shapeFunctionDerivative2DLinear(0.5, 0.5, i, 1);
        std::cout << "  Node " << i << ": dN/dxi = " << dNdxi 
                  << ", dN/deta = " << dNdeta << std::endl;
    }
}

void testJacobian() {
    std::cout << "\n=== Testing Jacobian Computation ===" << std::endl;
    
    // Create a simple 1D element with nodes at x=0 and x=2
    std::vector<Node> nodes1D;
    nodes1D.push_back(Node(0.0, 0.0, 0.0));
    nodes1D.push_back(Node(2.0, 0.0, 0.0));
    
    double jac1D = Integration::computeJacobian1D(nodes1D, 0.0);
    std::cout << "1D Jacobian at xi=0.0: " << jac1D << " (expected: 1.0)" << std::endl;
    
    // Create a simple 2D quadrilateral element
    std::vector<Node> nodes2D;
    nodes2D.push_back(Node(0.0, 0.0, 0.0));
    nodes2D.push_back(Node(2.0, 0.0, 0.0));
    nodes2D.push_back(Node(2.0, 1.0, 0.0));
    nodes2D.push_back(Node(0.0, 1.0, 0.0));
    
    double jac2D = Integration::computeJacobian2D(nodes2D, 0.0, 0.0);
    std::cout << "2D Jacobian at (0,0): " << jac2D << " (expected: 0.5)" << std::endl;
    
    auto jacMatrix2D = Integration::computeJacobianMatrix2D(nodes2D, 0.0, 0.0);
    std::cout << "2D Jacobian Matrix:" << std::endl;
    std::cout << "  [" << jacMatrix2D[0][0] << ", " << jacMatrix2D[0][1] << "]" << std::endl;
    std::cout << "  [" << jacMatrix2D[1][0] << ", " << jacMatrix2D[1][1] << "]" << std::endl;
}

void testIntegrationAccuracy() {
    std::cout << "\n=== Testing Integration Accuracy ===" << std::endl;
    
    // Test integrating x^2 from -1 to 1 using Gauss quadrature
    // Exact result: ∫x^2 dx from -1 to 1 = 2/3 ≈ 0.6666667
    
    auto points = Integration::getGaussPoints1D(3);
    double integral = 0.0;
    
    for (const auto& p : points) {
        double x = p.u;
        integral += p.weight * (x * x);
    }
    
    double exact = 2.0 / 3.0;
    double error = std::abs(integral - exact);
    
    std::cout << "Integral of x^2 from -1 to 1:" << std::endl;
    std::cout << "  Numerical: " << integral << std::endl;
    std::cout << "  Exact:     " << exact << std::endl;
    std::cout << "  Error:     " << error << std::endl;
    std::cout << "  Order 3 Gauss quadrature should be exact for polynomials up to degree 5" << std::endl;
}

int main() {
    std::cout << "Integration Module Test Suite" << std::endl;
    std::cout << "==============================" << std::endl;
    
    try {
        testGaussPoints1D();
        testGaussPoints2D();
        testGaussPoints3D();
        testTrianglePoints();
        testTetrahedronPoints();
        testShapeFunctions();
        testShapeFunctionDerivatives();
        testJacobian();
        testIntegrationAccuracy();
        
        std::cout << "\n=== All tests completed successfully! ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}