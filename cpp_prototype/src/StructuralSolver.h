#pragma once

#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include "ElementMatrix.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include "Material.h"
#include <memory>
#include <vector>
#include <array>

namespace elmer {

/**
 * @brief Boundary condition types for structural simulations
 */
enum class Structural_BoundaryType {
    DISPLACEMENT,       ///< Fixed displacement
    FORCE,              ///< Applied force
    PRESSURE,           ///< Applied pressure
    SPRING,             ///< Spring support
    PERIODIC,           ///< Periodic boundary condition
    SYMMETRY            ///< Symmetry boundary condition
};

/**
 * @brief Boundary condition for structural simulations
 */
struct StructuralBoundaryCondition {
    Structural_BoundaryType type;
    std::vector<int> nodeIndices;     ///< Nodes affected by this BC
    std::vector<std::array<double, 3>> values; ///< Boundary values (displacement or force)
    std::string name;                 ///< Boundary condition name
    
    // Additional parameters for specific boundary types
    double springConstant = 0.0;      ///< For spring BC
    double pressureValue = 0.0;       ///< For pressure BC
    
    StructuralBoundaryCondition(Structural_BoundaryType t = Structural_BoundaryType::DISPLACEMENT, 
                               const std::string& n = "")
        : type(t), name(n) {}
};

/**
 * @brief Solver parameters for structural simulations
 */
struct StructuralParameters {
    // Solver control
    double tolerance = 1.0e-8;        ///< Convergence tolerance
    int maxIterations = 1000;         ///< Maximum iterations
    
    // Analysis type
    bool isTransient = false;         ///< Transient analysis flag
    double timeStep = 1.0;            ///< Time step for transient analysis
    int timeSteps = 100;              ///< Number of time steps
    
    // Physical parameters
    bool includeThermalExpansion = false; ///< Include thermal expansion effects
    bool includeBodyForces = false;       ///< Include gravity and other body forces
    bool includeInertia = false;          ///< Include inertia effects (transient)
    
    // Body forces
    std::array<double, 3> gravity = {0.0, 0.0, -9.81}; ///< Gravity vector [m/s²]
    
    // Output control
    bool calculateStresses = true;    ///< Calculate stresses
    bool calculateStrains = true;     ///< Calculate strains
    bool calculateReactionForces = true; ///< Calculate reaction forces
    
    // Preconditioning
    std::string preconditioner = "ILU0"; ///< Preconditioner type
    
    StructuralParameters() = default;
};

/**
 * @brief Results from structural simulation
 */
struct StructuralResults {
    // Primary solution
    std::vector<std::array<double, 3>> displacement;   ///< Displacement field [m]
    
    // Derived fields
    std::vector<std::array<double, 6>> stress;         ///< Stress tensor [Pa] (σ_xx, σ_yy, σ_zz, τ_xy, τ_xz, τ_yz)
    std::vector<std::array<double, 6>> strain;         ///< Strain tensor (ε_xx, ε_yy, ε_zz, γ_xy, γ_xz, γ_yz)
    std::vector<std::array<double, 3>> reactionForce;  ///< Reaction forces [N]
    
    // Time-dependent results (for transient analysis)
    std::vector<std::vector<std::array<double, 3>>> displacementHistory; ///< Displacement history
    std::vector<double> timePoints;                    ///< Time points
    
    // Convergence information
    int iterations = 0;                      ///< Number of iterations
    double residual = 0.0;                   ///< Final residual
    bool converged = false;                  ///< Convergence status
    
    // Global quantities
    double totalStrainEnergy = 0.0;          ///< Total strain energy [J]
    double maxVonMisesStress = 0.0;          ///< Maximum von Mises stress [Pa]
    double maxDisplacement = 0.0;            ///< Maximum displacement [m]
    
    StructuralResults() = default;
};

/**
 * @brief Structural solver for mechanical analysis
 * 
 * This solver handles linear elastic structural mechanics problems:
 * - Static and dynamic analysis
 * - Thermal expansion effects
 * - Various boundary conditions (displacement, force, pressure)
 * - Stress and strain calculation
 */
class StructuralSolver {
private:
    std::shared_ptr<Mesh> mesh;
    MaterialDatabase materialDB;
    StructuralParameters parameters;
    
    // Boundary conditions
    std::vector<StructuralBoundaryCondition> boundaryConditions;
    
    // System matrices and vectors
    std::shared_ptr<CRSMatrix> stiffnessMatrix;       ///< Stiffness matrix
    std::shared_ptr<CRSMatrix> massMatrix;            ///< Mass matrix (transient)
    std::shared_ptr<Vector> forceVector;              ///< Force vector
    std::shared_ptr<Vector> boundaryVector;           ///< Boundary condition vector
    
    // Solution vectors
    std::shared_ptr<Vector> displacementVector;       ///< Displacement solution
    std::shared_ptr<Vector> previousDisplacement;     ///< Previous time step displacement
    std::shared_ptr<Vector> velocityVector;           ///< Velocity (for transient)
    std::shared_ptr<Vector> accelerationVector;       ///< Acceleration (for transient)
    
    // Temperature field (for thermal expansion)
    std::vector<double> temperatureField;
    
    // Iterative solver
    std::unique_ptr<IterativeSolver> solver;
    
public:
    /**
     * @brief Constructor
     */
    StructuralSolver(std::shared_ptr<Mesh> meshPtr = nullptr);
    
    /**
     * @brief Set the mesh for the simulation
     */
    void setMesh(std::shared_ptr<Mesh> meshPtr);
    
    /**
     * @brief Set solver parameters
     */
    void setParameters(const StructuralParameters& params);
    
    /**
     * @brief Add a boundary condition
     */
    void addBoundaryCondition(const StructuralBoundaryCondition& bc);
    
    /**
     * @brief Clear all boundary conditions
     */
    void clearBoundaryConditions();
    
    /**
     * @brief Set temperature field for thermal expansion
     */
    void setTemperatureField(const std::vector<double>& temperature);
    
    /**
     * @brief Assemble the system matrices
     */
    void assembleSystem();
    
    /**
     * @brief Solve the structural problem
     */
    StructuralResults solve();
    
    /**
     * @brief Solve transient structural problem
     */
    StructuralResults solveTransient();
    
    /**
     * @brief Calculate stresses from displacement field
     */
    std::vector<std::array<double, 6>> calculateStresses(const std::vector<std::array<double, 3>>& displacement);
    
    /**
     * @brief Set material database
     */
    void setMaterialDatabase(const MaterialDatabase& db) {
        materialDB = db;
    }
    
    /**
     * @brief Set body forces for the simulation
     */
    void setBodyForces(const std::vector<std::array<double, 3>>& bodyForces) {
        // This would set the body force vector
        // For now, this is a placeholder implementation
        std::cout << "Setting body forces for structural simulation" << std::endl;
    }
    
    /**
     * @brief Set thermal forces for thermal expansion
     */
    void setThermalForces(const std::vector<std::array<double, 3>>& thermalForces) {
        // This would set the thermal expansion forces
        // For now, this is a placeholder implementation
        std::cout << "Setting thermal forces for structural simulation" << std::endl;
    }
    
    /**
     * @brief Calculate strains from displacement field
     */
    std::vector<std::array<double, 6>> calculateStrains(const std::vector<std::array<double, 3>>& displacement);
    
    /**
     * @brief Calculate reaction forces
     */
    std::vector<std::array<double, 3>> calculateReactionForces(const std::vector<std::array<double, 3>>& displacement);
    
    /**
     * @brief Apply boundary conditions to the system
     */
    void applyBoundaryConditions();
    
    /**
     * @brief Calculate element stiffness matrix
     */
    std::vector<std::vector<double>> calculateElementStiffnessMatrix(const Element& element);
    
    /**
     * @brief Calculate element mass matrix
     */
    std::vector<std::vector<double>> calculateElementMassMatrix(const Element& element);
    
    /**
     * @brief Calculate element force vector
     */
    std::vector<double> calculateElementForceVector(const Element& element);
    
    /**
     * @brief Calculate element thermal expansion vector
     */
    std::vector<double> calculateElementThermalExpansion(const Element& element);
    
private:
    /**
     * @brief Initialize the solver
     */
    void initialize();
    
    /**
     * @brief Check if the system is properly set up
     */
    bool isSystemReady() const;
    
    /**
     * @brief Update time step for transient analysis
     */
    void updateTimeStep(int step, double time);
    
    /**
     * @brief Calculate constitutive matrix (D matrix) for linear elasticity
     */
    std::vector<std::vector<double>> calculateConstitutiveMatrix(const Material& material);
    
    /**
     * @brief Calculate strain-displacement matrix (B matrix)
     */
    std::vector<std::vector<double>> calculateStrainDisplacementMatrix(
        const std::vector<std::vector<double>>& shapeDerivatives);
    
    /**
     * @brief Calculate von Mises stress
     */
    double calculateVonMisesStress(const std::array<double, 6>& stress) const;
};

} // namespace elmer