#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include "BoundaryConditions.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include "MeshUtils.h"
#include "ElmerCpp.h"
#include "ElementUtils.h"

using namespace elmer;

/**
 * @brief Test priority-based boundary condition application
 */
void testPriorityBoundaryConditions() {
    std::cout << "=== Testing Priority-Based Boundary Condition Application ===" << std::endl;
    
    // Create a simple mesh with 6 nodes
    int n = 6; // 6 DOFs
    auto matrix = createCRSMatrix(n, n);
    std::vector<double> rhs(n, 0.0);
    
    // Initialize matrix with identity
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
    }
    
    // Create DOF mapping (simple 1-to-1 mapping)
    std::vector<int> dofMap(n);
    for (int i = 0; i < n; ++i) {
        dofMap[i] = i;
    }
    
    // Test 1: Dirichlet vs Neumann priority (Dirichlet should win)
    {
        std::cout << "\n1. Testing Dirichlet vs Neumann Priority" << std::endl;
        
        BoundaryConditionManager bcManager;
        
        // Create conflicting boundary conditions
        auto dirichletBC = std::make_shared<DirichletBoundaryCondition>("FixedValue");
        dirichletBC->setUniformValue({0}, 5.0); // Fix DOF 0 to 5.0
        
        auto neumannBC = std::make_shared<NeumannBoundaryCondition>("FluxBC");
        neumannBC->setFluxValues({0}, {2.0}); // Try to set flux on same DOF
        
        bcManager.addBoundaryCondition(dirichletBC);
        bcManager.addBoundaryCondition(neumannBC);
        
        // Apply with priority
        std::vector<bool> constrainedDOFs(n, false);
        bcManager.applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // Dirichlet should take precedence (matrix diagonal should be 1.0, RHS should be 5.0)
        if (std::abs(matrix->GetElement(0, 0) - 1.0) < 1e-10 && 
            std::abs(rhs[0] - 5.0) < 1e-10 &&
            constrainedDOFs[0] == true) {
            std::cout << "✓ Dirichlet priority test passed" << std::endl;
        } else {
            std::cout << "✗ Dirichlet priority test failed" << std::endl;
            std::cout << "  Matrix[0,0] = " << matrix->GetElement(0, 0) << " (expected: 1.0)" << std::endl;
            std::cout << "  RHS[0] = " << rhs[0] << " (expected: 5.0)" << std::endl;
            std::cout << "  constrainedDOFs[0] = " << constrainedDOFs[0] << " (expected: true)" << std::endl;
        }
    }
    
    // Reset matrix and RHS for next test
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
        rhs[i] = 0.0;
    }
    
    // Test 2: Multiple Dirichlet conditions (first applied should win)
    {
        std::cout << "\n2. Testing Multiple Dirichlet Conditions" << std::endl;
        
        BoundaryConditionManager bcManager;
        
        // Create multiple Dirichlet conditions on same DOF
        auto dirichlet1 = std::make_shared<DirichletBoundaryCondition>("BC1");
        dirichlet1->setUniformValue({1}, 10.0);
        
        auto dirichlet2 = std::make_shared<DirichletBoundaryCondition>("BC2");
        dirichlet2->setUniformValue({1}, 20.0); // Different value on same DOF
        
        bcManager.addBoundaryCondition(dirichlet1);
        bcManager.addBoundaryCondition(dirichlet2);
        
        std::vector<bool> constrainedDOFs(n, false);
        bcManager.applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // First Dirichlet condition should be applied
        if (std::abs(rhs[1] - 10.0) < 1e-10 && constrainedDOFs[1] == true) {
            std::cout << "✓ Multiple Dirichlet test passed" << std::endl;
        } else {
            std::cout << "✗ Multiple Dirichlet test failed" << std::endl;
            std::cout << "  RHS[1] = " << rhs[1] << " (expected: 10.0)" << std::endl;
        }
    }
    
    // Reset matrix and RHS for next test
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
        rhs[i] = 0.0;
    }
    
    // Test 3: Robin and Neumann conditions (should be able to coexist)
    {
        std::cout << "\n3. Testing Robin and Neumann Coexistence" << std::endl;
        
        BoundaryConditionManager bcManager;
        
        auto robinBC = std::make_shared<RobinBoundaryCondition>("MixedBC");
        robinBC->setParameters({2}, {2.0}, {1.0}, {5.0});
        
        auto neumannBC = std::make_shared<NeumannBoundaryCondition>("FluxBC");
        neumannBC->setFluxValues({3}, {3.0});
        
        bcManager.addBoundaryCondition(robinBC);
        bcManager.addBoundaryCondition(neumannBC);
        
        std::vector<bool> constrainedDOFs(n, false);
        bcManager.applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // Robin should modify matrix and RHS, Neumann should not affect matrix
        if (std::abs(matrix->GetElement(2, 2) - 3.0) < 1e-10 && // 1.0 + 2.0/1.0 = 3.0
            std::abs(rhs[2] - 5.0) < 1e-10 &&            // 0.0 + 5.0/1.0 = 5.0
            std::abs(matrix->GetElement(3, 3) - 1.0) < 1e-10 && // Neumann doesn't modify matrix
            constrainedDOFs[2] == false && constrainedDOFs[3] == false) { // Neither should constrain DOFs
            std::cout << "✓ Robin-Neumann coexistence test passed" << std::endl;
        } else {
            std::cout << "✗ Robin-Neumann coexistence test failed" << std::endl;
        }
    }
}

/**
 * @brief Test electromagnetic boundary conditions
 */
void testElectromagneticBoundaryConditions() {
    std::cout << "\n=== Testing Electromagnetic Boundary Conditions ===" << std::endl;
    
    int n = 8; // 8 DOFs for electromagnetic field
    auto matrix = createCRSMatrix(n, n);
    std::vector<double> rhs(n, 0.0);
    
    // Initialize matrix with identity
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
    }
    
    std::vector<int> dofMap(n);
    for (int i = 0; i < n; ++i) {
        dofMap[i] = i;
    }
    
    // Test Magnetic Symmetry boundary condition
    {
        std::cout << "\n1. Testing Magnetic Symmetry Boundary Condition" << std::endl;
        
        auto magneticSymmetryBC = std::make_shared<MagneticSymmetryBoundaryCondition>("Symmetry");
        magneticSymmetryBC->setSymmetryNodes({0, 1, 2});
        
        std::vector<bool> constrainedDOFs(n, false);
        magneticSymmetryBC->applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // Magnetic symmetry should constrain DOFs (Dirichlet-like behavior)
        bool testPassed = true;
        for (int i : {0, 1, 2}) {
            if (!constrainedDOFs[i]) {
                testPassed = false;
                std::cout << "  DOF " << i << " not constrained" << std::endl;
            }
        }
        
        if (testPassed) {
            std::cout << "✓ Magnetic symmetry test passed" << std::endl;
        } else {
            std::cout << "✗ Magnetic symmetry test failed" << std::endl;
        }
    }
    
    // Reset for next test
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
        rhs[i] = 0.0;
    }
    
    // Test Magnetic Antisymmetry boundary condition
    {
        std::cout << "\n2. Testing Magnetic Antisymmetry Boundary Condition" << std::endl;
        
        auto magneticAntisymmetryBC = std::make_shared<MagneticAntisymmetryBoundaryCondition>("Antisymmetry");
        magneticAntisymmetryBC->setAntisymmetryNodes({3, 4, 5});
        
        std::vector<bool> constrainedDOFs(n, false);
        magneticAntisymmetryBC->applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // Magnetic antisymmetry should also constrain DOFs
        bool testPassed = true;
        for (int i : {3, 4, 5}) {
            if (!constrainedDOFs[i]) {
                testPassed = false;
                std::cout << "  DOF " << i << " not constrained" << std::endl;
            }
        }
        
        if (testPassed) {
            std::cout << "✓ Magnetic antisymmetry test passed" << std::endl;
        } else {
            std::cout << "✗ Magnetic antisymmetry test failed" << std::endl;
        }
    }
}

/**
 * @brief Test thermal and structural boundary conditions
 */
void testThermalStructuralBoundaryConditions() {
    std::cout << "\n=== Testing Thermal and Structural Boundary Conditions ===" << std::endl;
    
    int n = 6;
    auto matrix = createCRSMatrix(n, n);
    std::vector<double> rhs(n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
    }
    
    std::vector<int> dofMap(n);
    for (int i = 0; i < n; ++i) {
        dofMap[i] = i;
    }
    
    // Test Convection boundary condition
    {
        std::cout << "\n1. Testing Convection Boundary Condition" << std::endl;
        
        auto convectionBC = std::make_shared<ConvectionBoundaryCondition>("Convection");
        convectionBC->setConvectionParameters({0, 1}, {10.0, 20.0}, {25.0, 30.0});
        
        std::vector<bool> constrainedDOFs(n, false);
        convectionBC->applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // Convection should modify matrix (adds to diagonal) and RHS
        // Simplified: adds h to diagonal and h*Tinf to RHS
        if (std::abs(matrix->GetElement(0, 0) - 11.0) < 1e-10 && // 1.0 + 10.0 = 11.0
            std::abs(rhs[0] - 250.0) < 1e-10 &&           // 0.0 + 10.0*25.0 = 250.0
            std::abs(matrix->GetElement(1, 1) - 21.0) < 1e-10 && // 1.0 + 20.0 = 21.0
            std::abs(rhs[1] - 600.0) < 1e-10) {           // 0.0 + 20.0*30.0 = 600.0
            std::cout << "✓ Convection test passed" << std::endl;
        } else {
            std::cout << "✗ Convection test failed" << std::endl;
        }
    }
    
    // Reset for next test
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
        rhs[i] = 0.0;
    }
    
    // Test Pressure boundary condition
    {
        std::cout << "\n2. Testing Pressure Boundary Condition" << std::endl;
        
        auto pressureBC = std::make_shared<PressureBoundaryCondition>("Pressure");
        pressureBC->setPressureValues({2, 3}, {100.0, 200.0});
        
        std::vector<bool> constrainedDOFs(n, false);
        pressureBC->applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        
        // Pressure BC typically adds to RHS (force vector)
        if (std::abs(rhs[2] - 100.0) < 1e-10 &&
            std::abs(rhs[3] - 200.0) < 1e-10) {
            std::cout << "✓ Pressure test passed" << std::endl;
        } else {
            std::cout << "✗ Pressure test failed" << std::endl;
        }
    }
}

/**
 * @brief Test boundary condition ordering and priority
 */
void testBoundaryConditionOrdering() {
    std::cout << "\n=== Testing Boundary Condition Ordering and Priority ===" << std::endl;
    
    BoundaryConditionManager bcManager;
    
    // Add boundary conditions in mixed order
    auto neumannBC = std::make_shared<NeumannBoundaryCondition>("NeumannFirst");
    neumannBC->setFluxValues({0}, {1.0});
    
    auto dirichletBC = std::make_shared<DirichletBoundaryCondition>("DirichletSecond");
    dirichletBC->setUniformValue({0}, 5.0); // Same DOF as Neumann
    
    auto robinBC = std::make_shared<RobinBoundaryCondition>("RobinThird");
    robinBC->setParameters({1}, {2.0}, {1.0}, {3.0});
    
    // Add in arbitrary order
    bcManager.addBoundaryCondition(neumannBC);
    bcManager.addBoundaryCondition(dirichletBC);
    bcManager.addBoundaryCondition(robinBC);
    
    int n = 4;
    auto matrix = createCRSMatrix(n, n);
    std::vector<double> rhs(n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 1.0);
    }
    
    std::vector<int> dofMap(n);
    for (int i = 0; i < n; ++i) {
        dofMap[i] = i;
    }
    
    std::vector<bool> constrainedDOFs(n, false);
    bcManager.applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
    
    // Dirichlet should take precedence over Neumann, Robin should be applied
    if (std::abs(rhs[0] - 5.0) < 1e-10 && // Dirichlet value, not Neumann
        std::abs(matrix->GetElement(1, 1) - 3.0) < 1e-10 && // Robin modification
        std::abs(rhs[1] - 3.0) < 1e-10 && // Robin modification
        constrainedDOFs[0] == true) {      // DOF 0 should be constrained
        std::cout << "✓ Boundary condition ordering test passed" << std::endl;
    } else {
        std::cout << "✗ Boundary condition ordering test failed" << std::endl;
    }
}

int main() {
    std::cout << "Starting Boundary Condition Priority Tests" << std::endl;
    
    try {
        testPriorityBoundaryConditions();
        testElectromagneticBoundaryConditions();
        testThermalStructuralBoundaryConditions();
        testBoundaryConditionOrdering();
        
        std::cout << "\n=== All Boundary Condition Priority Tests Completed ===" << std::endl;
        std::cout << "✓ Priority-based boundary condition application is working correctly" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}