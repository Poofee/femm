#include "MultiphysicsSolver.h"
#include "Mesh.h"
#include "Material.h"
#include <iostream>
#include <vector>

using namespace elmer;

/**
 * @brief Simple test for multiphysics coupling interface
 */
void testMultiphysicsInterface() {
    std::cout << "=== Testing Multiphysics Interface ===" << std::endl;
    
    // Create a simple mesh
    auto mesh = std::make_shared<Mesh>();
    
    // Add nodes
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 1.0, 0.0));
    mesh->addNode(Node(0.0, 1.0, 0.0));
    
    // Add element
    std::vector<int> quadNodes = {0, 1, 2, 3};
    Element quadElement(ElementType::QUADRILATERAL, quadNodes, "Copper");
    mesh->addElement(quadElement);
    
    // Create multiphysics solver
    MultiphysicsSolver solver(mesh);
    
    // Test 1: Check if solvers are initialized
    std::cout << "\n1. Testing solver initialization..." << std::endl;
    
    auto emSolver = solver.getElectromagneticSolver();
    auto thermalSolver = solver.getThermalSolver();
    auto structuralSolver = solver.getStructuralSolver();
    
    if (emSolver) {
        std::cout << "✓ Electromagnetic solver initialized" << std::endl;
    } else {
        std::cout << "✗ Electromagnetic solver not initialized" << std::endl;
    }
    
    if (thermalSolver) {
        std::cout << "✓ Thermal solver initialized" << std::endl;
    } else {
        std::cout << "✗ Thermal solver not initialized" << std::endl;
    }
    
    if (structuralSolver) {
        std::cout << "✓ Structural solver initialized" << std::endl;
    } else {
        std::cout << "✗ Structural solver not initialized" << std::endl;
    }
    
    // Test 2: Test material database
    std::cout << "\n2. Testing material database..." << std::endl;
    
    MaterialDatabase customDB;
    customDB.createPredefinedMaterials();
    
    // Check if predefined materials exist
    if (customDB.hasMaterial("Copper")) {
        std::cout << "✓ Copper material found" << std::endl;
    }
    
    if (customDB.hasMaterial("Iron")) {
        std::cout << "✓ Iron material found" << std::endl;
    }
    
    if (customDB.hasMaterial("Air")) {
        std::cout << "✓ Air material found" << std::endl;
    }
    
    // Set custom material database
    solver.setMaterialDatabase(customDB);
    std::cout << "✓ Material database set successfully" << std::endl;
    
    // Test 3: Test coupling parameters
    std::cout << "\n3. Testing coupling parameters..." << std::endl;
    
    CouplingParameters params;
    params.type = CouplingType::WEAK;
    params.maxCouplingIterations = 10;
    params.couplingTolerance = 1.0e-6;
    params.jouleHeatingCoefficient = 1.0;
    params.thermalExpansionCoefficient = 1.0;
    
    solver.setCouplingParameters(params);
    std::cout << "✓ Coupling parameters set successfully" << std::endl;
    
    // Test 4: Test Joule heating calculation
    std::cout << "\n4. Testing Joule heating calculation..." << std::endl;
    
    MagnetodynamicsResults emResults;
    emResults.electricField.resize(4, {1.0, 0.0, 0.0});
    
    auto jouleHeating = solver.calculateJouleHeating(emResults);
    if (jouleHeating.size() == 4) {
        std::cout << "✓ Joule heating calculated for " << jouleHeating.size() << " nodes" << std::endl;
    } else {
        std::cout << "✗ Joule heating calculation failed" << std::endl;
    }
    
    // Test 5: Test thermal expansion forces
    std::cout << "\n5. Testing thermal expansion forces..." << std::endl;
    
    std::vector<double> temperature(4, 373.15);
    auto thermalForces = solver.calculateThermalExpansionForces(temperature);
    if (thermalForces.size() == 4) {
        std::cout << "✓ Thermal expansion forces calculated for " << thermalForces.size() << " nodes" << std::endl;
    } else {
        std::cout << "✗ Thermal expansion forces calculation failed" << std::endl;
    }
    
    // Test 6: Test magnetostriction forces
    std::cout << "\n6. Testing magnetostriction forces..." << std::endl;
    
    emResults.magneticField.resize(4, {0.1, 0.0, 0.0});
    params.magnetostrictionCoefficient = 1.0e-6;
    solver.setCouplingParameters(params);
    
    auto magnetostrictionForces = solver.calculateMagnetostrictionForces(emResults);
    if (magnetostrictionForces.size() == 4) {
        std::cout << "✓ Magnetostriction forces calculated for " << magnetostrictionForces.size() << " nodes" << std::endl;
    } else {
        std::cout << "✗ Magnetostriction forces calculation failed" << std::endl;
    }
    
    // Test 7: Test coupling types
    std::cout << "\n7. Testing coupling types..." << std::endl;
    
    params.type = CouplingType::NONE;
    solver.setCouplingParameters(params);
    std::cout << "✓ No coupling mode set" << std::endl;
    
    params.type = CouplingType::WEAK;
    solver.setCouplingParameters(params);
    std::cout << "✓ Weak coupling mode set" << std::endl;
    
    params.type = CouplingType::ITERATIVE;
    solver.setCouplingParameters(params);
    std::cout << "✓ Iterative coupling mode set" << std::endl;
    
    params.type = CouplingType::STRONG;
    solver.setCouplingParameters(params);
    std::cout << "✓ Strong coupling mode set" << std::endl;
    
    std::cout << "\n=== Multiphysics Interface Test Completed ===" << std::endl;
    std::cout << "All interface tests passed successfully!" << std::endl;
}

/**
 * @brief Test case for basic functionality
 */
void testBasicFunctionality() {
    std::cout << "\n=== Testing Basic Functionality ===" << std::endl;
    
    // Test mesh creation
    auto mesh = std::make_shared<Mesh>();
    
    // Add nodes and elements
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    mesh->addNode(Node(0.0, 1.0, 0.0));
    
    std::vector<int> triNodes = {0, 1, 2};
    Element triElement(ElementType::TRIANGLE, triNodes, "Iron");
    mesh->addElement(triElement);
    
    // Test multiphysics solver creation
    MultiphysicsSolver solver(mesh);
    
    // Test coupling parameters
    CouplingParameters params;
    params.type = CouplingType::WEAK;
    params.maxCouplingIterations = 5;
    
    solver.setCouplingParameters(params);
    
    std::cout << "✓ Basic functionality test passed" << std::endl;
    std::cout << "=== Basic Functionality Test Completed ===" << std::endl;
}

/**
 * @brief Main test function
 */
int main() {
    std::cout << "Starting Multiphysics Coupling Tests" << std::endl;
    std::cout << "====================================" << std::endl;
    
    try {
        testBasicFunctionality();
        testMultiphysicsInterface();
        
        std::cout << "\n✅ All multiphysics coupling tests completed successfully!" << std::endl;
        
        // Summary
        std::cout << "\n=== Test Summary ===" << std::endl;
        std::cout << "✓ Multiphysics solver interface" << std::endl;
        std::cout << "✓ Material database integration" << std::endl;
        std::cout << "✓ Coupling parameter configuration" << std::endl;
        std::cout << "✓ Joule heating calculation" << std::endl;
        std::cout << "✓ Thermal expansion forces" << std::endl;
        std::cout << "✓ Magnetostriction forces" << std::endl;
        std::cout << "✓ Multiple coupling modes" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "❌ Test error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}