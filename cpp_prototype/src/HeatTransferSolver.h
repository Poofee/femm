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
 * @brief Boundary condition types for heat transfer simulations
 */
enum class Thermal_BoundaryType {
    DIRICHLET,          ///< Fixed temperature
    NEUMANN,            ///< Fixed heat flux
    ROBIN,              ///< Convective boundary condition
    RADIATION,          ///< Radiation boundary condition
    PERIODIC,           ///< Periodic boundary condition
    SYMMETRY            ///< Symmetry boundary condition
};

/**
 * @brief Boundary condition for heat transfer simulations
 */
struct ThermalBoundaryCondition {
    Thermal_BoundaryType type;
    std::vector<int> nodeIndices;     ///< Nodes affected by this BC
    std::vector<double> values;       ///< Boundary values
    std::string name;                 ///< Boundary condition name
    
    // Additional parameters for specific boundary types
    double heatTransferCoefficient = 0.0;  ///< For Robin BC
    double ambientTemperature = 0.0;      ///< For Robin/Radiation BC
    double emissivity = 0.0;              ///< For Radiation BC
    
    ThermalBoundaryCondition(Thermal_BoundaryType t = Thermal_BoundaryType::DIRICHLET, 
                           const std::string& n = "")
        : type(t), name(n) {}
};

/**
 * @brief Solver parameters for heat transfer simulations
 */
struct HeatTransferParameters {
    // Solver control
    double tolerance = 1.0e-8;        ///< Convergence tolerance
    int maxIterations = 1000;         ///< Maximum iterations
    
    // Analysis type
    bool isTransient = false;         ///< Transient analysis flag
    double timeStep = 1.0;            ///< Time step for transient analysis
    int timeSteps = 100;              ///< Number of time steps
    
    // Physical parameters
    bool includeConvection = false;   ///< Include natural convection
    bool includeRadiation = false;    ///< Include radiation effects
    bool includePhaseChange = false;  ///< Include phase change effects
    
    // Initial conditions
    double initialTemperature = 293.15; ///< Initial temperature [K]
    
    // Output control
    bool calculateHeatFlux = true;   ///< Calculate heat flux
    bool calculateTemperatureGradient = true; ///< Calculate temperature gradient
    
    // Preconditioning
    std::string preconditioner = "ILU0"; ///< Preconditioner type
    
    HeatTransferParameters() = default;
};

/**
 * @brief Results from heat transfer simulation
 */
struct HeatTransferResults {
    // Primary solution
    std::vector<double> temperature;         ///< Temperature field [K]
    
    // Derived fields
    std::vector<std::array<double, 3>> heatFlux;          ///< Heat flux [W/m²]
    std::vector<std::array<double, 3>> temperatureGradient; ///< Temperature gradient [K/m]
    
    // Time-dependent results (for transient analysis)
    std::vector<std::vector<double>> temperatureHistory;  ///< Temperature history
    std::vector<double> timePoints;                       ///< Time points
    
    // Convergence information
    int iterations = 0;                      ///< Number of iterations
    double residual = 0.0;                   ///< Final residual
    bool converged = false;                  ///< Convergence status
    
    // Energy quantities
    double totalHeatFlux = 0.0;              ///< Total heat flux [W]
    double averageTemperature = 0.0;          ///< Average temperature [K]
    
    HeatTransferResults() = default;
};

/**
 * @brief Heat transfer solver for thermal analysis
 * 
 * This solver handles both steady-state and transient heat transfer problems:
 * - Conduction (Fourier's law)
 * - Convection (natural and forced)
 * - Radiation (Stefan-Boltzmann law)
 * - Phase change (latent heat)
 */
class HeatTransferSolver {
private:
    std::shared_ptr<Mesh> mesh;
    MaterialDatabase materialDB;
    HeatTransferParameters parameters;
    
    // Boundary conditions
    std::vector<ThermalBoundaryCondition> boundaryConditions;
    
    // System matrices and vectors
    std::shared_ptr<CRSMatrix> conductivityMatrix;    ///< Thermal conductivity matrix
    std::shared_ptr<CRSMatrix> capacityMatrix;       ///< Heat capacity matrix (transient)
    std::shared_ptr<Vector> heatSourceVector;       ///< Heat source vector
    std::shared_ptr<Vector> boundaryVector;          ///< Boundary condition vector
    
    // Solution vectors
    std::shared_ptr<Vector> temperatureVector;       ///< Temperature solution
    std::shared_ptr<Vector> previousTemperature;      ///< Previous time step temperature
    
    // Iterative solver
    std::unique_ptr<IterativeSolver> solver;
    
public:
    /**
     * @brief Constructor
     */
    HeatTransferSolver(std::shared_ptr<Mesh> meshPtr = nullptr);
    
    /**
     * @brief Set the mesh for the simulation
     */
    void setMesh(std::shared_ptr<Mesh> meshPtr);
    
    /**
     * @brief Set solver parameters
     */
    void setParameters(const HeatTransferParameters& params);
    
    /**
     * @brief Add a boundary condition
     */
    void addBoundaryCondition(const ThermalBoundaryCondition& bc);
    
    /**
     * @brief Clear all boundary conditions
     */
    void clearBoundaryConditions();
    
    /**
     * @brief Assemble the system matrices
     */
    void assembleSystem();
    
    /**
     * @brief Solve the heat transfer problem
     */
    HeatTransferResults solve();
    
    /**
     * @brief Solve transient heat transfer problem
     */
    HeatTransferResults solveTransient();
    
    /**
     * @brief Calculate heat flux from temperature field
     */
    std::vector<std::array<double, 3>> calculateHeatFlux(const std::vector<double>& temperature);
    
    /**
     * @brief Calculate temperature gradient
     */
    std::vector<std::array<double, 3>> calculateTemperatureGradient(const std::vector<double>& temperature);
    
    /**
     * @brief Set material database
     */
    void setMaterialDatabase(const MaterialDatabase& db) {
        materialDB = db;
    }
    
    /**
     * @brief Set heat source for the simulation
     */
    void setHeatSource(const std::vector<double>& heatSource) {
        // This would set the internal heat source vector
        // For now, this is a placeholder implementation
        std::cout << "Setting heat source for thermal simulation" << std::endl;
    }
    
    /**
     * @brief Apply boundary conditions to the system
     */
    void applyBoundaryConditions();
    
    /**
     * @brief Calculate element conductivity matrix
     */
    std::vector<std::vector<double>> calculateElementConductivityMatrix(const Element& element);
    
    /**
     * @brief Calculate element capacity matrix
     */
    std::vector<std::vector<double>> calculateElementCapacityMatrix(const Element& element);
    
    /**
     * @brief Calculate element heat source vector
     */
    std::vector<double> calculateElementHeatSource(const Element& element);
    
    /**
     * @brief Calculate element boundary flux vector
     */
    std::vector<double> calculateElementBoundaryFlux(const Element& element, 
                                                    const ThermalBoundaryCondition& bc);
    
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
};

} // namespace elmer