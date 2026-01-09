#include "../src/ShapeFunctions.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace elmer;

void test1DLinearElement() {
    std::cout << "=== Testing 1D Linear Element ===" << std::endl;
    
    // Create a simple 1D element with nodes at x=0 and x=2
    std::vector<Node> nodes = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
    
    // Test at center (xi = 0)
    auto result = ShapeFunctions::computeShapeFunctions(ElementType::LINEAR, nodes, 0.0, 0.0, 0.0);
    
    std::cout << "Nodes: (" << nodes[0].x << ", " << nodes[0].y << ", " << nodes[0].z << "), "
              << "(" << nodes[1].x << ", " << nodes[1].y << ", " << nodes[1].z << ")" << std::endl;
    std::cout << "Natural coordinate: xi = 0.0" << std::endl;
    std::cout << "Shape function values: N0 = " << result.values[0] << ", N1 = " << result.values[1] << std::endl;
    std::cout << "Sum of shape functions: " << result.values[0] + result.values[1] << " (should be 1.0)" << std::endl;
    std::cout << "Jacobian determinant: " << result.detJ << " (should be 1.0)" << std::endl;
    std::cout << "dN/dx: dN0/dx = " << result.dNdx[0] << ", dN1/dx = " << result.dNdx[1] << std::endl;
    
    // Verify partition of unity
    double sum = result.values[0] + result.values[1];
    if (std::abs(sum - 1.0) < 1e-10) {
        std::cout << "✓ Partition of unity test PASSED" << std::endl;
    } else {
        std::cout << "✗ Partition of unity test FAILED" << std::endl;
    }
    
    std::cout << std::endl;
}

void test2DQuadrilateralElement() {
    std::cout << "=== Testing 2D Quadrilateral Element ===" << std::endl;
    
    // Create a simple quadrilateral element
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0},  // Node 0
        {2.0, 0.0, 0.0},  // Node 1
        {2.0, 1.0, 0.0},  // Node 2
        {0.0, 1.0, 0.0}   // Node 3
    };
    
    // Test at center (xi = 0, eta = 0)
    auto result = ShapeFunctions::computeShapeFunctions(ElementType::LINEAR, nodes, 0.0, 0.0, 0.0);
    
    std::cout << "Natural coordinates: xi = 0.0, eta = 0.0" << std::endl;
    std::cout << "Shape function values: ";
    for (size_t i = 0; i < result.values.size(); ++i) {
        std::cout << "N" << i << " = " << result.values[i] << " ";
    }
    std::cout << std::endl;
    
    double sum = 0.0;
    for (double val : result.values) {
        sum += val;
    }
    std::cout << "Sum of shape functions: " << sum << " (should be 1.0)" << std::endl;
    std::cout << "Jacobian determinant: " << result.detJ << std::endl;
    
    // Verify partition of unity
    if (std::abs(sum - 1.0) < 1e-10) {
        std::cout << "✓ Partition of unity test PASSED" << std::endl;
    } else {
        std::cout << "✗ Partition of unity test FAILED" << std::endl;
    }
    
    // Verify that Jacobian is positive (valid element)
    if (result.detJ > 0.0) {
        std::cout << "✓ Jacobian determinant test PASSED" << std::endl;
    } else {
        std::cout << "✗ Jacobian determinant test FAILED" << std::endl;
    }
    
    std::cout << std::endl;
}

void test3DHexahedronElement() {
    std::cout << "=== Testing 3D Hexahedron Element ===" << std::endl;
    
    // Create a simple hexahedral element
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0},  // Node 0
        {1.0, 0.0, 0.0},  // Node 1
        {1.0, 1.0, 0.0},  // Node 2
        {0.0, 1.0, 0.0},  // Node 3
        {0.0, 0.0, 1.0},  // Node 4
        {1.0, 0.0, 1.0},  // Node 5
        {1.0, 1.0, 1.0},  // Node 6
        {0.0, 1.0, 1.0}   // Node 7
    };
    
    // Test at center (xi = 0, eta = 0, zeta = 0)
    auto result = ShapeFunctions::computeShapeFunctions(ElementType::LINEAR, nodes, 0.0, 0.0, 0.0);
    
    std::cout << "Natural coordinates: xi = 0.0, eta = 0.0, zeta = 0.0" << std::endl;
    std::cout << "Shape function values: ";
    for (size_t i = 0; i < result.values.size(); ++i) {
        std::cout << "N" << i << " = " << result.values[i] << " ";
    }
    std::cout << std::endl;
    
    double sum = 0.0;
    for (double val : result.values) {
        sum += val;
    }
    std::cout << "Sum of shape functions: " << sum << " (should be 1.0)" << std::endl;
    std::cout << "Jacobian determinant: " << result.detJ << std::endl;
    
    // Verify partition of unity
    if (std::abs(sum - 1.0) < 1e-10) {
        std::cout << "✓ Partition of unity test PASSED" << std::endl;
    } else {
        std::cout << "✗ Partition of unity test FAILED" << std::endl;
    }
    
    // Verify that Jacobian is positive (valid element)
    if (result.detJ > 0.0) {
        std::cout << "✓ Jacobian determinant test PASSED" << std::endl;
    } else {
        std::cout << "✗ Jacobian determinant test FAILED" << std::endl;
    }
    
    std::cout << std::endl;
}

void testTriangleElement() {
    std::cout << "=== Testing Triangle Element ===" << std::endl;
    
    // Create a simple triangle element
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0},  // Node 0
        {2.0, 0.0, 0.0},  // Node 1
        {0.0, 1.0, 0.0}   // Node 2
    };
    
    // Test at center (xi = 1/3, eta = 1/3)
    auto result = ShapeFunctions::computeShapeFunctions(ElementType::LINEAR, nodes, 1.0/3.0, 1.0/3.0, 0.0);
    
    std::cout << "Natural coordinates: xi = 1/3, eta = 1/3" << std::endl;
    std::cout << "Shape function values: ";
    for (size_t i = 0; i < result.values.size(); ++i) {
        std::cout << "N" << i << " = " << result.values[i] << " ";
    }
    std::cout << std::endl;
    
    double sum = 0.0;
    for (double val : result.values) {
        sum += val;
    }
    std::cout << "Sum of shape functions: " << sum << " (should be 1.0)" << std::endl;
    std::cout << "Jacobian determinant (area element): " << result.detJ << std::endl;
    
    // Verify partition of unity
    if (std::abs(sum - 1.0) < 1e-10) {
        std::cout << "✓ Partition of unity test PASSED" << std::endl;
    } else {
        std::cout << "✗ Partition of unity test FAILED" << std::endl;
    }
    
    std::cout << std::endl;
}

void testTetrahedronElement() {
    std::cout << "=== Testing Tetrahedron Element ===" << std::endl;
    
    // Create a simple tetrahedral element
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0},  // Node 0
        {1.0, 0.0, 0.0},  // Node 1
        {0.0, 1.0, 0.0},  // Node 2
        {0.0, 0.0, 1.0}   // Node 3
    };
    
    // Test at center (xi = 1/4, eta = 1/4, zeta = 1/4)
    auto result = ShapeFunctions::computeShapeFunctions(ElementType::TETRAHEDRON, nodes, 0.25, 0.25, 0.25);
    
    std::cout << "Natural coordinates: xi = 0.25, eta = 0.25, zeta = 0.25" << std::endl;
    std::cout << "Shape function values: ";
    for (size_t i = 0; i < result.values.size(); ++i) {
        std::cout << "N" << i << " = " << result.values[i] << " ";
    }
    std::cout << std::endl;
    
    double sum = 0.0;
    for (double val : result.values) {
        sum += val;
    }
    std::cout << "Sum of shape functions: " << sum << " (should be 1.0)" << std::endl;
    std::cout << "Jacobian determinant (volume element): " << result.detJ << std::endl;
    
    // Verify partition of unity
    if (std::abs(sum - 1.0) < 1e-10) {
        std::cout << "✓ Partition of unity test PASSED" << std::endl;
    } else {
        std::cout << "✗ Partition of unity test FAILED" << std::endl;
    }
    
    std::cout << std::endl;
}

void testShapeFunctionProperties() {
    std::cout << "=== Testing Shape Function Properties ===" << std::endl;
    
    // Test 1D quadratic shape functions
    auto quad1D = ShapeFunctions::linear1D(0.5, 3);
    std::cout << "1D Quadratic shape functions at xi=0.5: ";
    for (double val : quad1D) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // Test 2D quadrilateral quadratic shape functions
    auto quad2D = ShapeFunctions::linearQuadrilateral(0.0, 0.0, 8);
    std::cout << "2D Quadratic shape functions at center: ";
    for (double val : quad2D) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // Test triangle quadratic shape functions
    auto quadTri = ShapeFunctions::linearTriangle(1.0/3.0, 1.0/3.0, 6);
    std::cout << "Triangle Quadratic shape functions at center: ";
    for (double val : quadTri) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << std::endl;
}

int main() {
    std::cout << "Elmer FEM C++ Shape Functions Test Program" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
    
    try {
        test1DLinearElement();
        test2DQuadrilateralElement();
        test3DHexahedronElement();
        testTriangleElement();
        testTetrahedronElement();
        testShapeFunctionProperties();
        
        std::cout << "All tests completed successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}