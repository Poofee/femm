#include "MultiphysicsSolver.h"
#include "Mesh.h"
#include "Material.h"
#include <iostream>
#include <vector>

using namespace elmer;

/**
 * @brief Test case for multiphysics coupling functionality
 */
void testMultiphysicsCoupling() {
    std::cout << "=== Testing Multiphysics Coupling ===" << std::endl;
    
    // Create a simple mesh for testing
    auto mesh = std::make_shared<Mesh>();
    
    // Add some nodes (simplified 2D mesh)
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 1.0, 0.0));
    mesh->addNode(Node(0.0, 1.0, 0.0));
    
    // Add a quadrilateral element
    std::vector<int> quadNodes = {0, 1, 2, 3};
    Element quadElement(ElementType::QUADRILATERAL, quadNodes, "Copper");
    mesh->addElement(quadElement);
    
    // Create multiphysics solver
    MultiphysicsSolver multiphysicsSolver(mesh);
    
    // Test different coupling types
    CouplingParameters params;
    
    std::cout << "\n1. Testing uncoupled simulation..." << std::endl;
    params.type = CouplingType::NONE;
    multiphysicsSolver.setCouplingParameters(params);
    
    try {
        multiphysicsSolver.solveCoupledProblem();
        std::cout << "Uncoupled simulation completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Uncoupled simulation error: " << e.what() << std::endl;
    }
    
    std::cout << "\n2. Testing weak coupling..." << std::endl;
    params.type = CouplingType::WEAK;
    params.maxCouplingIterations = 3;
    params.couplingTolerance = 1.0e-4;
    params.jouleHeatingCoefficient = 0.5;
    multiphysicsSolver.setCouplingParameters(params);
    
    try {
        multiphysicsSolver.solveCoupledProblem();
        std::cout << "Weak coupling simulation completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Weak coupling simulation error: " << e.what() << std::endl;
    }
    
    std::cout << "\n3. Testing iterative coupling..." << std::endl;
    params.type = CouplingType::ITERATIVE;
    params.maxCouplingIterations = 5;
    params.couplingTolerance = 1.0e-6;
    multiphysicsSolver.setCouplingParameters(params);
    
    try {
        multiphysicsSolver.solveCoupledProblem();
        std::cout << "Iterative coupling simulation completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Iterative coupling simulation error: " << e.what() << std::endl;
    }
    
    std::cout << "\n4. Testing strong coupling..." << std::endl;
    params.type = CouplingType::STRONG;
    multiphysicsSolver.setCouplingParameters(params);
    
    try {
        multiphysicsSolver.solveCoupledProblem();
        std::cout << "Strong coupling simulation completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Strong coupling simulation error: " << e.what() << std::endl;
    }
    
    // Test individual solver access
    std::cout << "\n5. Testing individual solver access..." << std::endl;
    auto emSolver = multiphysicsSolver.getElectromagneticSolver();
    auto thermalSolver = multiphysicsSolver.getThermalSolver();
    auto structuralSolver = multiphysicsSolver.getStructuralSolver();
    
    if (emSolver) {
        std::cout << "Electromagnetic solver available" << std::endl;
    }
    
    if (thermalSolver) {
        std::cout << "Thermal solver available" << std::endl;
    }
    
    if (structuralSolver) {
        std::cout << "Structural solver available" << std::endl;
    }
    
    // Test material database
    std::cout << "\n6. Testing material database..." << std::endl;
    MaterialDatabase customDB;
    customDB.createPredefinedMaterials();
    
    // Add a custom material
    Material customMaterial;
    customMaterial.name = "CustomMaterial";
    customMaterial.relativePermittivity = 2.0;
    customMaterial.relativePermeability = 1000.0;
    customMaterial.conductivity = 1.0e6;
    customMaterial.thermalConductivity = 50.0;
    customMaterial.density = 5000.0;
    customMaterial.specificHeat = 500.0;
    customMaterial.youngsModulus = 100e9;
    customMaterial.poissonsRatio = 0.3;
    customMaterial.thermalExpansion = 15e-6;
    
    customDB.addMaterial("CustomMaterial", customMaterial);
    multiphysicsSolver.setMaterialDatabase(customDB);
    
    std::cout << "Material database updated successfully" << std::endl;
    
    // Test Joule heating calculation
    std::cout << "\n7. Testing Joule heating calculation..." << std::endl;
    MagnetodynamicsResults emResults;
    emResults.electricField.resize(4, {1.0, 0.0, 0.0}); // Simple electric field
    
    auto jouleHeating = multiphysicsSolver.calculateJouleHeating(emResults);
    std::cout << "Joule heating calculated for " << jouleHeating.size() << " nodes" << std::endl;
    
    // Test thermal expansion forces calculation
    std::cout << "\n8. Testing thermal expansion forces calculation..." << std::endl;
    std::vector<double> temperature(4, 373.15); // 100°C temperature
    
    auto thermalForces = multiphysicsSolver.calculateThermalExpansionForces(temperature);
    std::cout << "Thermal expansion forces calculated for " << thermalForces.size() << " nodes" << std::endl;
    
    // Test magnetostriction forces calculation
    std::cout << "\n9. Testing magnetostriction forces calculation..." << std::endl;
    emResults.magneticField.resize(4, {0.1, 0.0, 0.0}); // Simple magnetic field
    params.magnetostrictionCoefficient = 1.0e-6;
    multiphysicsSolver.setCouplingParameters(params);
    
    auto magnetostrictionForces = multiphysicsSolver.calculateMagnetostrictionForces(emResults);
    std::cout << "Magnetostriction forces calculated for " << magnetostrictionForces.size() << " nodes" << std::endl;
    
    std::cout << "\n=== Multiphysics Coupling Test Completed ===" << std::endl;
}

/**
 * @brief Test case for electromagnetic-thermal coupling
 */
void testElectromagneticThermalCoupling() {
    std::cout << "\n=== Testing Electromagnetic-Thermal Coupling ===" << std::endl;
    
    // Create a simple mesh for testing
    auto mesh = std::make_shared<Mesh>();
    
    // Add nodes for a simple conductor
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 1.0, 0.0));
    mesh->addNode(Node(0.0, 1.0, 0.0));
    
    // Add a quadrilateral element representing a conductor
    std::vector<int> quadNodes = {0, 1, 2, 3};
    Element quadElement(ElementType::QUADRILATERAL, quadNodes, "Copper");
    mesh->addElement(quadElement);
    
    // Create multiphysics solver with weak coupling
    MultiphysicsSolver multiphysicsSolver(mesh);
    
    CouplingParameters params;
    params.type = CouplingType::WEAK;
    params.jouleHeatingCoefficient = 1.0; // Full Joule heating effect
    multiphysicsSolver.setCouplingParameters(params);
    
    std::cout << "Simulating electromagnetic-thermal coupling..." << std::endl;
    
    try {
        multiphysicsSolver.solveCoupledProblem();
        std::cout << "Electromagnetic-thermal coupling simulation completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Electromagnetic-thermal coupling error: " << e.what() << std::endl;
    }
    
    std::cout << "=== Electromagnetic-Thermal Coupling Test Completed ===" << std::endl;
}

/**
 * @brief Test case for thermal-structural coupling
 */
void testThermalStructuralCoupling() {
    std::cout << "\n=== Testing Thermal-Structural Coupling ===" << std::endl;
    
    // Create a simple mesh for testing
    auto mesh = std::make_shared<Mesh>();
    
    // Add nodes for a simple structural element
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 1.0, 0.0));
    mesh->addNode(Node(0.0, 1.0, 0.0));
    
    // Add a quadrilateral element
    std::vector<int> quadNodes = {0, 1, 2, 3};
    Element quadElement(ElementType::QUADRILATERAL, quadNodes, "Iron");
    mesh->addElement(quadElement);
    
    // Create multiphysics solver
    MultiphysicsSolver multiphysicsSolver(mesh);
    
    CouplingParameters params;
    params.type = CouplingType::WEAK;
    params.thermalExpansionCoefficient = 1.0; // Full thermal expansion effect
    multiphysicsSolver.setCouplingParameters(params);
    
    std::cout << "Simulating thermal-structural coupling..." << std::endl;
    
    try {
        multiphysicsSolver.solveCoupledProblem();
        std::cout << "Thermal-structural coupling simulation completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Thermal-structural coupling error: " << e.what() << std::endl;
    }
    
    std::cout << "=== Thermal-Structural Coupling Test Completed ===" << std::endl;
}

/**
 * @brief Main test function
 */
int main() {
    std::cout << "Starting Multiphysics Coupling Tests" << std::endl;
    std::cout << "====================================" << std::endl;
    
    try {
        testMultiphysicsCoupling();
        testElectromagneticThermalCoupling();
        testThermalStructuralCoupling();
        
        std::cout << "\nAll multiphysics coupling tests completed successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Test error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}