#include "MultiphysicsSolver.h"
#include "HeatTransferSolver.h"
#include "StructuralSolver.h"
#include <iostream>

namespace elmer {

/**
 * @brief Initialize individual solvers
 */
void MultiphysicsSolver::initializeSolvers() {
    emSolver = std::make_shared<MagnetodynamicsSolver>(mesh);
    thermalSolver = std::make_shared<HeatTransferSolver>(mesh);
    structuralSolver = std::make_shared<StructuralSolver>(mesh);
    
    // Set material database for all solvers
    if (emSolver) emSolver->setMaterialDatabase(materialDB);
    if (thermalSolver) thermalSolver->setMaterialDatabase(materialDB);
    if (structuralSolver) structuralSolver->setMaterialDatabase(materialDB);
}

/**
 * @brief Solve uncoupled problems
 */
void MultiphysicsSolver::solveUncoupled() {
    std::cout << "Solving uncoupled multiphysics problems..." << std::endl;
    
    // Solve electromagnetic problem
    if (emSolver) {
        auto emResults = emSolver->solve();
        potentialField = emResults.potentialReal;
        std::cout << "Electromagnetic solver completed with " 
                  << emResults.iterations << " iterations" << std::endl;
    }
    
    // Solve thermal problem
    if (thermalSolver) {
        auto thermalResults = thermalSolver->solve();
        temperatureField = thermalResults.temperature;
        std::cout << "Thermal solver completed with " 
                  << thermalResults.iterations << " iterations" << std::endl;
    }
    
    // Solve structural problem
    if (structuralSolver) {
        auto structuralResults = structuralSolver->solve();
        displacementField = structuralResults.displacement;
        std::cout << "Structural solver completed with " 
                  << structuralResults.iterations << " iterations" << std::endl;
    }
}

/**
 * @brief Solve with weak coupling (sequential)
 */
void MultiphysicsSolver::solveWeakCoupling() {
    std::cout << "Solving with weak coupling..." << std::endl;
    
    // Step 1: Solve electromagnetic problem
    auto emResults = emSolver->solve();
    potentialField = emResults.potentialReal;
    
    // Step 2: Calculate coupling effects
    auto jouleHeating = calculateJouleHeating(emResults);
    auto lorentzForces = emResults.lorentzForce;
    
    // Step 3: Solve thermal problem with Joule heating
    thermalSolver->setHeatSource(jouleHeating);
    auto thermalResults = thermalSolver->solve();
    temperatureField = thermalResults.temperature;
    
    // Step 4: Calculate thermal expansion forces
    auto thermalExpansionForces = calculateThermalExpansionForces(temperatureField);
    
    // Step 5: Solve structural problem with thermal and electromagnetic forces
    structuralSolver->setBodyForces(lorentzForces);
    structuralSolver->setThermalForces(thermalExpansionForces);
    structuralSolver->setTemperatureField(temperatureField);
    auto structuralResults = structuralSolver->solve();
    displacementField = structuralResults.displacement;
    
    std::cout << "Weak coupling completed:" << std::endl;
    std::cout << "- EM iterations: " << emResults.iterations << std::endl;
    std::cout << "- Thermal iterations: " << thermalResults.iterations << std::endl;
    std::cout << "- Structural iterations: " << structuralResults.iterations << std::endl;
}

/**
 * @brief Solve with strong coupling (simultaneous)
 */
void MultiphysicsSolver::solveStrongCoupling() {
    std::cout << "Solving with strong coupling (monolithic system)..." << std::endl;
    
    // Assemble monolithic system matrix
    assembleMonolithicSystem();
    
    // Solve the coupled system
    // This would involve a larger linear system with all fields
    
    std::cout << "Strong coupling completed (simplified implementation)" << std::endl;
}

/**
 * @brief Solve with iterative coupling
 */
void MultiphysicsSolver::solveIterativeCoupling() {
    std::cout << "Solving with iterative coupling..." << std::endl;
    
    double residual = 1.0;
    int iteration = 0;
    
    // Initialize solution fields
    std::vector<double> prevTemperature = temperatureField;
    std::vector<double> prevDisplacement = displacementField;
    
    while (residual > couplingParams.couplingTolerance && 
           iteration < couplingParams.maxCouplingIterations) {
        
        std::cout << "Iteration " << iteration + 1 << ":" << std::endl;
        
        // Solve electromagnetic problem
        auto emResults = emSolver->solve();
        potentialField = emResults.potentialReal;
        
        // Calculate coupling terms
        auto jouleHeating = calculateJouleHeating(emResults);
        auto lorentzForces = emResults.lorentzForce;
        
        // Solve thermal problem with Joule heating
        thermalSolver->setHeatSource(jouleHeating);
        auto thermalResults = thermalSolver->solve();
        temperatureField = thermalResults.temperature;
        
        // Calculate thermal expansion forces
        auto thermalExpansionForces = calculateThermalExpansionForces(temperatureField);
        
        // Solve structural problem with thermal and electromagnetic forces
        structuralSolver->setBodyForces(lorentzForces);
        structuralSolver->setThermalForces(thermalExpansionForces);
        structuralSolver->setTemperatureField(temperatureField);
        auto structuralResults = structuralSolver->solve();
        displacementField = structuralResults.displacement;
        
        // Calculate coupling residual
        residual = calculateCouplingResidual(prevTemperature, prevDisplacement,
                                           temperatureField, displacementField);
        
        std::cout << "  Residual: " << residual << std::endl;
        
        // Update previous solutions for next iteration
        prevTemperature = temperatureField;
        prevDisplacement = displacementField;
        
        iteration++;
    }
    
    std::cout << "Iterative coupling completed after " << iteration << " iterations" << std::endl;
}

/**
 * @brief Assemble monolithic system matrix for strong coupling
 */
void MultiphysicsSolver::assembleMonolithicSystem() {
    // This would assemble a large system matrix that includes
    // all coupled physics in a single matrix
    
    // Implementation depends on the specific coupling
    // For now, this is a placeholder implementation
    
    std::cout << "Assembling monolithic system matrix (simplified)..." << std::endl;
}

/**
 * @brief Calculate coupling residual for iterative method
 */
double MultiphysicsSolver::calculateCouplingResidual(
    const std::vector<double>& prevTemp,
    const std::vector<double>& prevDisp,
    const std::vector<double>& currentTemp,
    const std::vector<double>& currentDisp) {
    
    if (prevTemp.empty() || prevDisp.empty()) {
        return 1.0; // First iteration
    }
    
    double tempResidual = 0.0;
    double dispResidual = 0.0;
    
    // Calculate temperature residual
    for (size_t i = 0; i < currentTemp.size(); ++i) {
        double diff = currentTemp[i] - prevTemp[i];
        tempResidual += diff * diff;
    }
    tempResidual = std::sqrt(tempResidual) / currentTemp.size();
    
    // Calculate displacement residual
    for (size_t i = 0; i < currentDisp.size(); ++i) {
        double diff = currentDisp[i] - prevDisp[i];
        dispResidual += diff * diff;
    }
    dispResidual = std::sqrt(dispResidual) / currentDisp.size();
    
    // Return maximum residual
    return std::max(tempResidual, dispResidual);
}

/**
 * @brief Calculate thermal expansion forces (enhanced implementation)
 */
std::vector<std::array<double, 3>> MultiphysicsSolver::calculateThermalExpansionForces(
    const std::vector<double>& temperature) {
    
    int nNodes = mesh->getNodes().size();
    std::vector<std::array<double, 3>> thermalForces(nNodes, {0.0, 0.0, 0.0});
    
    ShapeFunctions shapeFunc;
    
    for (const auto& element : mesh->getElements()) {
        auto material = materialDB.getMaterial(element.getMaterialName());
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate at element center
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), {0.0, 0.0, 0.0});
        
        if (!shapeResult.shapeFunctionDerivatives.empty()) {
            // Calculate thermal strain: ε_thermal = α * ΔT
            double avgTemp = 0.0;
            for (int nodeIdx : nodeIndices) {
                avgTemp += temperature[nodeIdx];
            }
            avgTemp /= nodeIndices.size();
            
            double thermalStrain = material.thermalExpansion * (avgTemp - 293.15); // 20°C reference
            
            // Convert to thermal stress: σ_thermal = E * ε_thermal
            double thermalStress = material.youngsModulus * thermalStrain;
            
            // Distribute thermal forces to nodes
            for (size_t i = 0; i < nodeIndices.size(); ++i) {
                for (int dim = 0; dim < 3; ++dim) {
                    thermalForces[nodeIndices[i]][dim] += thermalStress * 
                                                         shapeResult.shapeFunctions[i] * 0.1; // Simplified
                }
            }
        }
    }
    
    return thermalForces;
}

/**
 * @brief Calculate magnetostriction forces (enhanced implementation)
 */
std::vector<std::array<double, 3>> MultiphysicsSolver::calculateMagnetostrictionForces(
    const MagnetodynamicsResults& emResults) {
    
    int nNodes = mesh->getNodes().size();
    std::vector<std::array<double, 3>> magnetostrictionForces(nNodes, {0.0, 0.0, 0.0});
    
    // Simplified implementation based on magnetic field strength
    for (size_t i = 0; i < nNodes; ++i) {
        double B_squared = 0.0;
        for (int dim = 0; dim < 3; ++dim) {
            B_squared += emResults.magneticField[i][dim] * emResults.magneticField[i][dim];
        }
        
        // Magnetostriction force proportional to B²
        double magnetostrictionForce = couplingParams.magnetostrictionCoefficient * B_squared;
        
        // Distribute force equally in all directions (simplified)
        for (int dim = 0; dim < 3; ++dim) {
            magnetostrictionForces[i][dim] = magnetostrictionForce * 0.01; // Scaling factor
        }
    }
    
    return magnetostrictionForces;
}

} // namespace elmer