#pragma once

#include "MagnetodynamicsSolver.h"
#include "HeatTransferSolver.h"
#include "StructuralSolver.h"
#include "Mesh.h"
#include "Material.h"
#include <memory>
#include <vector>
#include <array>
#include <iostream>

// Use explicit namespace to avoid conflicts
#include "Material.h"

namespace elmer {

/**
 * @brief Types of multiphysics coupling
 */
enum class CouplingType {
    NONE,                   ///< No coupling
    WEAK,                   ///< Weak coupling (sequential)
    STRONG,                 ///< Strong coupling (simultaneous)
    ITERATIVE              ///< Iterative coupling
};

/**
 * @brief Coupling parameters for multiphysics simulations
 */
struct CouplingParameters {
    CouplingType type = CouplingType::WEAK;
    int maxCouplingIterations = 10;          ///< Maximum iterations for iterative coupling
    double couplingTolerance = 1.0e-6;       ///< Tolerance for coupling convergence
    bool updateGeometry = false;             ///< Update mesh geometry during coupling
    
    // Physical coupling parameters
    double thermalExpansionCoefficient = 0.0; ///< For thermal-structural coupling
    double jouleHeatingCoefficient = 1.0;     ///< For electromagnetic-thermal coupling
    double magnetostrictionCoefficient = 0.0; ///< For electromagnetic-structural coupling
    
    // Coupling flags
    bool enableElectroThermalCoupling = true;    ///< Enable electromagnetic-thermal coupling
    bool enableElectroStructuralCoupling = true; ///< Enable electromagnetic-structural coupling
    bool enableThermoStructuralCoupling = true;  ///< Enable thermal-structural coupling
    
    CouplingParameters() = default;
};

/**
 * @brief Multiphysics solver for coupled problems
 * 
 * This solver handles coupling between different physical fields:
 * - Electromagnetic-thermal coupling (Joule heating)
 * - Electromagnetic-structural coupling (Lorentz forces, magnetostriction)
 * - Thermal-structural coupling (thermal expansion)
 */
class MultiphysicsSolver {
private:
    std::shared_ptr<Mesh> mesh;
    MaterialDatabase materialDB;
    CouplingParameters couplingParams;
    
    // Individual solvers
    std::shared_ptr<MagnetodynamicsSolver> emSolver;
    std::shared_ptr<HeatTransferSolver> thermalSolver;
    std::shared_ptr<StructuralSolver> structuralSolver;
    
    // Coupling matrices and vectors
    std::shared_ptr<CRSMatrix> couplingMatrix;
    std::shared_ptr<Vector> couplingRHS;
    
    // Solution fields
    std::vector<double> temperatureField;     ///< Temperature field
    std::vector<double> displacementField;    ///< Structural displacement field
    std::vector<double> potentialField;       ///< Electromagnetic potential field
    
public:
    /**
     * @brief Constructor
     */
    MultiphysicsSolver(std::shared_ptr<Mesh> m = nullptr)
        : mesh(m) {
        materialDB.createPredefinedMaterials();
        if (mesh) {
            initializeSolvers();
        }
    }
    
    /**
     * @brief Set the mesh for the simulation
     */
    void setMesh(std::shared_ptr<Mesh> m) {
        mesh = m;
        if (mesh) {
            initializeSolvers();
        }
    }
    
    /**
     * @brief Set coupling parameters
     */
    void setCouplingParameters(const CouplingParameters& params) {
        couplingParams = params;
    }
    
    /**
     * @brief Set coupling type
     */
    void setCouplingType(CouplingType type) {
        couplingParams.type = type;
    }
    
    /**
     * @brief Solve coupled multiphysics problem
     */
    void solveCoupledProblem() {
        if (!mesh) {
            throw std::runtime_error("Mesh not set for multiphysics simulation");
        }
        
        switch (couplingParams.type) {
            case CouplingType::WEAK:
                solveWeakCoupling();
                break;
            case CouplingType::STRONG:
                solveStrongCoupling();
                break;
            case CouplingType::ITERATIVE:
                solveIterativeCoupling();
                break;
            default:
                solveUncoupled();
                break;
        }
    }
    
    /**
     * @brief Get electromagnetic solver
     */
    std::shared_ptr<MagnetodynamicsSolver> getElectromagneticSolver() {
        return emSolver;
    }
    
    /**
     * @brief Get thermal solver
     */
    std::shared_ptr<HeatTransferSolver> getThermalSolver() {
        return thermalSolver;
    }
    
    /**
     * @brief Get structural solver
     */
    std::shared_ptr<StructuralSolver> getStructuralSolver() {
        return structuralSolver;
    }
    
    /**
     * @brief Get potential field (electromagnetic)
     */
    const std::vector<double>& getPotentialField() const {
        return potentialField;
    }
    
    /**
     * @brief Get temperature field
     */
    const std::vector<double>& getTemperatureField() const {
        return temperatureField;
    }
    
    /**
     * @brief Get displacement field
     */
    const std::vector<double>& getDisplacementField() const {
        return displacementField;
    }
    
    /**
     * @brief Set material database for all solvers
     */
    void setMaterialDatabase(const MaterialDatabase& db) {
        materialDB = db;
        if (emSolver) emSolver->setMaterialDatabase(materialDB);
        if (thermalSolver) thermalSolver->setMaterialDatabase(materialDB);
        if (structuralSolver) structuralSolver->setMaterialDatabase(materialDB);
    }
    
    /**
     * @brief Calculate Joule heating from electromagnetic solution
     */
    std::vector<double> calculateJouleHeating(const MagnetodynamicsResults& emResults) {
        size_t nNodes = mesh->getNodes().numberOfNodes();
        std::vector<double> jouleHeating(nNodes, 0.0);
        
        // Q = J·E = σ|E|² for conductive heating
        // Simplified implementation using a default material
        auto material = materialDB.getMaterial("Copper");
        
        for (const auto& element : mesh->getBulkElements()) {
            auto nodeIndices = element.getNodeIndices();
            
            for (size_t nodeIdx : nodeIndices) {
                double E_squared = 0.0;
                for (int i = 0; i < 3; ++i) {
                    E_squared += emResults.electricField[nodeIdx][i] * 
                                emResults.electricField[nodeIdx][i];
                }
                jouleHeating[nodeIdx] += material.conductivity * E_squared * 
                                       couplingParams.jouleHeatingCoefficient;
            }
        }
        
        return jouleHeating;
    }
    
    /**
     * @brief Calculate thermal expansion forces
     */
    std::vector<std::array<double, 3>> calculateThermalExpansionForces(
        const std::vector<double>& temperature) {
        
        size_t nNodes = mesh->getNodes().numberOfNodes();
        std::vector<std::array<double, 3>> thermalForces(nNodes, {0.0, 0.0, 0.0});
        
        // Simplified implementation
        // In practice, this would involve stress-strain relationships
        
        return thermalForces;
    }
    
    /**
     * @brief Calculate magnetostriction forces
     */
    std::vector<std::array<double, 3>> calculateMagnetostrictionForces(
        const MagnetodynamicsResults& emResults) {
        
        size_t nNodes = mesh->getNodes().numberOfNodes();
        std::vector<std::array<double, 3>> magnetostrictionForces(nNodes, {0.0, 0.0, 0.0});
        
        // Simplified implementation
        // In practice, this would involve magnetostrictive coefficients
        
        return magnetostrictionForces;
    }
    
public:
    /**
     * @brief Initialize individual solvers
     */
    void initializeSolvers() {
        emSolver = std::make_shared<MagnetodynamicsSolver>(mesh);
        // Initialize other solvers here
    }
    
    /**
     * @brief Solve uncoupled problems
     */
    void solveUncoupled() {
        // Solve each physics independently
        if (emSolver) {
            auto emResults = emSolver->solve();
            // Store results
        }
        // Solve other physics here
    }
    
    /**
     * @brief Solve with weak coupling (sequential)
     */
    void solveWeakCoupling() {
        // Sequential solution: solve one physics, then the other
        
        // Step 1: Solve electromagnetic problem
        auto emResults = emSolver->solve();
        
        // Step 2: Calculate coupling effects
        auto jouleHeating = calculateJouleHeating(emResults);
        auto lorentzForces = emResults.lorentzForce;
        
        // Step 3: Solve thermal problem with Joule heating
        // thermalSolver->setHeatSource(jouleHeating);
        // auto thermalResults = thermalSolver->solve();
        
        // Step 4: Solve structural problem with thermal and electromagnetic forces
        // structuralSolver->setBodyForces(lorentzForces);
        // structuralSolver->setThermalForces(thermalExpansionForces);
        // auto structuralResults = structuralSolver->solve();
        
        // Store results
        potentialField = emResults.potentialReal;
    }
    
    /**
     * @brief Solve with strong coupling (simultaneous)
     */
    void solveStrongCoupling() {
        // Assemble monolithic system matrix
        assembleMonolithicSystem();
        
        // Solve the coupled system
        // This would involve a larger linear system with all fields
    }
    
    /**
     * @brief Solve with iterative coupling
     */
    void solveIterativeCoupling() {
        double residual = 1.0;
        int iteration = 0;
        
        while (residual > couplingParams.couplingTolerance && 
               iteration < couplingParams.maxCouplingIterations) {
            
            // Solve electromagnetic problem
            auto emResults = emSolver->solve();
            
            // Calculate coupling terms
            auto jouleHeating = calculateJouleHeating(emResults);
            auto lorentzForces = emResults.lorentzForce;
            
            // Solve other physics with coupling terms
            // thermalSolver->setHeatSource(jouleHeating);
            // auto thermalResults = thermalSolver->solve();
            
            // structuralSolver->setBodyForces(lorentzForces);
            // auto structuralResults = structuralSolver->solve();
            
            // Calculate coupling residual
            residual = calculateCouplingResidual();
            
            iteration++;
        }
    }

private:
    
    /**
     * @brief Assemble monolithic system matrix for strong coupling
     */
    void assembleMonolithicSystem() {
        // This would assemble a large system matrix that includes
        // all coupled physics in a single matrix
        
        // Implementation depends on the specific coupling
    }
    
    /**
     * @brief Calculate coupling residual for iterative method
     */
    double calculateCouplingResidual() {
        // Calculate the residual between coupled solutions
        // This measures how much the solutions have changed between iterations
        
        return 0.0; // Placeholder
    }
};

} // namespace elmer