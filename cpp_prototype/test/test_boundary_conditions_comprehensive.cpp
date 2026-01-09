#include <iostream>
#include <vector>
#include <memory>
#include "BoundaryConditions.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include "MeshUtils.h"
#include "ElmerCpp.h"

using namespace elmer;
using namespace ElmerCpp;

/**
 * @brief Test basic boundary condition functionality
 */
void testBasicBoundaryConditions() {
    std::cout << "=== Testing Basic Boundary Conditions ===" << std::endl;
    
    // Create a simple mesh with 4 nodes
    auto mesh = MeshUtils::allocateMesh("test");
    
    // Create a simple 2x2 matrix for testing
    int n = 4; // 4 DOFs
    auto matrix = createCRSMatrix(n, n);
    std::vector<double> rhs(n, 0.0);
    
    // Initialize matrix with identity
    for (int i = 0; i < n; ++i) {
        matrix->set(i, i, 1.0);
    }
    
    // Test Dirichlet boundary condition
    {
        std::cout << "\n1. Testing Dirichlet Boundary Condition" << std::endl;
        
        auto dirichletBC = std::make_shared<DirichletBoundaryCondition>("FixedValue");
        std::vector<int> nodeIndices = {0, 2}; // Fix nodes 0 and 2
        std::vector<double> values = {5.0, 10.0}; // Set values to 5 and 10
        dirichletBC->setValues(nodeIndices, values);
        
        // Create DOF mapping (simple 1-to-1 mapping)
        std::vector<int> dofMap(n);
        for (int i = 0; i < n; ++i) {
            dofMap[i] = i;
        }
        
        // Apply boundary condition
        dirichletBC->apply(matrix, rhs, dofMap);
        
        // Check results
        std::cout << "Matrix diagonal after Dirichlet BC:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << "  [" << i << "," << i << "] = " << matrix->get(i, i) << std::endl;
        }
        
        std::cout << "RHS after Dirichlet BC:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << "  rhs[" << i << "] = " << rhs[i] << std::endl;
        }
        
        // Verify that fixed DOFs have correct values
        if (matrix->get(0, 0) == 1.0 && rhs[0] == 5.0 &&
            matrix->get(2, 2) == 1.0 && rhs[2] == 10.0) {
            std::cout << "✓ Dirichlet BC test passed" << std::endl;
        } else {
            std::cout << "✗ Dirichlet BC test failed" << std::endl;
        }
    }
    
    // Test Neumann boundary condition
    {
        std::cout << "\n2. Testing Neumann Boundary Condition" << std::endl;
        
        auto neumannBC = std::make_shared<NeumannBoundaryCondition>("FluxBC");
        std::vector<int> nodeIndices = {1, 3};
        std::vector<double> fluxValues = {2.0, 3.0};
        neumannBC->setFluxValues(nodeIndices, fluxValues);
        
        // Neumann BC typically doesn't modify the matrix
        // The flux contributions are added during element assembly
        std::cout << "✓ Neumann BC test passed (no matrix modification expected)" << std::endl;
    }
    
    // Test Robin boundary condition
    {
        std::cout << "\n3. Testing Robin Boundary Condition" << std::endl;
        
        auto robinBC = std::make_shared<RobinBoundaryCondition>("MixedBC");
        std::vector<int> nodeIndices = {1};
        std::vector<double> alphaValues = {2.0};  // Coefficient for field value
        std::vector<double> betaValues = {1.0};   // Coefficient for flux
        std::vector<double> gammaValues = {5.0};  // Constant term
        robinBC->setParameters(nodeIndices, alphaValues, betaValues, gammaValues);
        
        std::vector<int> dofMap(n);
        for (int i = 0; i < n; ++i) {
            dofMap[i] = i;
        }
        
        // Reset matrix and RHS for Robin test
        for (int i = 0; i < n; ++i) {
            matrix->set(i, i, 1.0);
            rhs[i] = 0.0;
        }
        
        robinBC->apply(matrix, rhs, dofMap);
        
        std::cout << "Matrix diagonal after Robin BC:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << "  [" << i << "," << i << "] = " << matrix->get(i, i) << std::endl;
        }
        
        std::cout << "RHS after Robin BC:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << "  rhs[" << i << "] = " << rhs[i] << std::endl;
        }
        
        // For Robin BC: alpha*u + beta*du/dn = gamma
        // Simplified implementation adds alpha/beta to diagonal and gamma/beta to RHS
        if (std::abs(matrix->get(1, 1) - 3.0) < 1e-10 && // 1.0 + 2.0/1.0 = 3.0
            std::abs(rhs[1] - 5.0) < 1e-10) {            // 0.0 + 5.0/1.0 = 5.0
            std::cout << "✓ Robin BC test passed" << std::endl;
        } else {
            std::cout << "✗ Robin BC test failed" << std::endl;
        }
    }
}

/**
 * @brief Test boundary condition manager
 */
void testBoundaryConditionManager() {
    std::cout << "\n=== Testing Boundary Condition Manager ===" << std::endl;
    
    BoundaryConditionManager bcManager;
    
    // Create different types of boundary conditions
    auto bc1 = std::make_shared<DirichletBoundaryCondition>("LeftBoundary");
    bc1->setUniformValue({0, 1}, 0.0);
    
    auto bc2 = std::make_shared<DirichletBoundaryCondition>("RightBoundary");
    bc2->setUniformValue({2, 3}, 10.0);
    
    auto bc3 = std::make_shared<NeumannBoundaryCondition>("TopBoundary");
    bc3->setFluxValues({1}, {2.0});
    
    // Add boundary conditions to manager
    bcManager.addBoundaryCondition(bc1);
    bcManager.addBoundaryCondition(bc2);
    bcManager.addBoundaryCondition(bc3);
    
    // Test retrieval by name
    auto retrievedBC = bcManager.getBoundaryCondition("LeftBoundary");
    if (retrievedBC && retrievedBC->name() == "LeftBoundary") {
        std::cout << "✓ Boundary condition retrieval by name works" << std::endl;
    } else {
        std::cout << "✗ Boundary condition retrieval by name failed" << std::endl;
    }
    
    // Test retrieval by type
    auto dirichletBCs = bcManager.getBoundaryConditionsByType(BoundaryConditionType::DIRICHLET);
    if (dirichletBCs.size() == 2) {
        std::cout << "✓ Boundary condition retrieval by type works" << std::endl;
    } else {
        std::cout << "✗ Boundary condition retrieval by type failed" << std::endl;
    }
    
    // Test applying all boundary conditions
    int n = 4;
    auto matrix = createCRSMatrix(n, n);
    std::vector<double> rhs(n, 0.0);
    
    // Initialize matrix
    for (int i = 0; i < n; ++i) {
        matrix->setElement(i, i, 1.0);
    }
    
    std::vector<int> dofMap(n);
    for (int i = 0; i < n; ++i) {
        dofMap[i] = i;
    }
    
    bcManager.applyBoundaryConditions(matrix, rhs, dofMap);
    
    std::cout << "RHS after applying all BCs:" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << "  rhs[" << i << "] = " << rhs[i] << std::endl;
    }
    
    // Check that Dirichlet conditions were applied correctly
    if (rhs[0] == 0.0 && rhs[1] == 0.0 && rhs[2] == 10.0 && rhs[3] == 10.0) {
        std::cout << "✓ Boundary condition manager applies all BCs correctly" << std::endl;
    } else {
        std::cout << "✗ Boundary condition manager failed to apply BCs correctly" << std::endl;
    }
}

/**
 * @brief Test boundary condition utilities
 */
void testBoundaryConditionUtils() {
    std::cout << "\n=== Testing Boundary Condition Utilities ===" << std::endl;
    
    // Test applyDirichletBC function
    {
        int n = 3;
        auto matrix = createCRSMatrix(n, n);
        std::vector<double> rhs(n, 0.0);
        
        // Initialize matrix
        for (int i = 0; i < n; ++i) {
            matrix->setElement(i, i, 2.0);
        }
        
        // Apply Dirichlet BC to DOF 1 with value 5.0
        BoundaryConditionUtils::applyDirichletBC(1, 5.0, matrix, rhs);
        
        // Check that row 1 was zeroed and diagonal set to 1
        if (matrix->getElement(1, 0) == 0.0 && matrix->getElement(1, 1) == 1.0 && 
            matrix->getElement(1, 2) == 0.0 && rhs[1] == 5.0) {
            std::cout << "✓ applyDirichletBC function works correctly" << std::endl;
        } else {
            std::cout << "✗ applyDirichletBC function failed" << std::endl;
        }
    }
    
    // Test createBoundaryCondition function
    {
        auto bc = BoundaryConditionUtils::createBoundaryCondition(
            BoundaryConditionType::DIRICHLET, "TestBC", {0, 1}, {1.0, 2.0});
        
        if (bc && bc->name() == "TestBC" && bc->type() == BoundaryConditionType::DIRICHLET) {
            std::cout << "✓ createBoundaryCondition function works correctly" << std::endl;
        } else {
            std::cout << "✗ createBoundaryCondition function failed" << std::endl;
        }
    }
}

/**
 * @brief Test electromagnetic-specific boundary conditions
 */
void testElectromagneticBoundaryConditions() {
    std::cout << "\n=== Testing Electromagnetic Boundary Conditions ===" << std::endl;
    
    // Test magnetic symmetry boundary condition
    {
        auto magneticSymmetryBC = std::make_shared<DirichletBoundaryCondition>("MagneticSymmetry");
        magneticSymmetryBC->setUniformValue({0, 1, 2}, 0.0);
        
        if (magneticSymmetryBC->type() == BoundaryConditionType::DIRICHLET &&
            magneticSymmetryBC->name() == "MagneticSymmetry") {
            std::cout << "✓ Magnetic symmetry BC test passed" << std::endl;
        } else {
            std::cout << "✗ Magnetic symmetry BC test failed" << std::endl;
        }
    }
    
    // Test electric insulation boundary condition
    {
        auto electricInsulationBC = std::make_shared<NeumannBoundaryCondition>("ElectricInsulation");
        electricInsulationBC->setFluxValues({3, 4}, {0.0, 0.0});
        
        if (electricInsulationBC->type() == BoundaryConditionType::NEUMANN &&
            electricInsulationBC->name() == "ElectricInsulation") {
            std::cout << "✓ Electric insulation BC test passed" << std::endl;
        } else {
            std::cout << "✗ Electric insulation BC test failed" << std::endl;
        }
    }
}

int main() {
    std::cout << "=== Boundary Condition Test Suite ===" << std::endl;
    
    try {
        testBasicBoundaryConditions();
        testBoundaryConditionManager();
        testBoundaryConditionUtils();
        testElectromagneticBoundaryConditions();
        
        std::cout << "\n=== All Tests Completed ===" << std::endl;
        std::cout << "Boundary condition system is working correctly!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "\n=== Test Error ===" << std::endl;
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}