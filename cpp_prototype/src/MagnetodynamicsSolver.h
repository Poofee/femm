#pragma once

#include "ElectromagneticMaterial.h"
#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include "ElementMatrix.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include <memory>
#include <vector>
#include <array>
#include <complex>

namespace elmer {

/**
 * @brief Boundary condition types for electromagnetic simulations
 */
enum class EM_BoundaryType {
    DIRICHLET,          ///< Fixed potential/field value
    NEUMANN,            ///< Fixed flux/current density
    ROBIN,              ///< Mixed boundary condition
    PERIODIC,           ///< Periodic boundary condition
    SYMMETRY,           ///< Symmetry boundary condition
    ANTISYMMETRY        ///< Antisymmetry boundary condition
};

/**
 * @brief Boundary condition for electromagnetic simulations
 */
struct EMBoundaryCondition {
    EM_BoundaryType type;
    std::vector<int> nodeIndices;     ///< Nodes affected by this BC
    std::vector<double> values;       ///< Boundary values
    std::vector<double> valuesImag;   ///< Imaginary values (for harmonic analysis)
    std::string name;                 ///< Boundary condition name
    
    EMBoundaryCondition(EM_BoundaryType t = EM_BoundaryType::DIRICHLET, 
                       const std::string& n = "")
        : type(t), name(n) {}
};

/**
 * @brief Solver parameters for magnetodynamics simulations
 */
struct MagnetodynamicsParameters {
    // Solver control
    double tolerance = 1.0e-8;        ///< Convergence tolerance
    int maxIterations = 1000;         ///< Maximum iterations
    int restart = 100;                ///< Restart parameter for GCR
    
    // Analysis type
    bool isHarmonic = false;          ///< Harmonic analysis flag
    double frequency = 0.0;           ///< Frequency for harmonic analysis
    
    // Physical parameters
    bool includeDisplacementCurrent = false;  ///< Include displacement current
    bool includeConvection = false;           ///< Include convection terms
    
    // Output control
    bool calculateMagneticField = true;       ///< Calculate magnetic field
    bool calculateElectricField = false;      ///< Calculate electric field
    bool calculateCurrentDensity = false;     ///< Calculate current density
    bool calculateLorentzForce = false;       ///< Calculate Lorentz force
    
    // Preconditioning
    std::string preconditioner = "ILU0";      ///< Preconditioner type
    
    MagnetodynamicsParameters() = default;
};

/**
 * @brief Results from magnetodynamics simulation
 */
struct MagnetodynamicsResults {
    // Primary solution
    std::vector<double> potentialReal;        ///< Real part of potential
    std::vector<double> potentialImag;        ///< Imaginary part of potential
    
    // Derived fields
    std::vector<std::array<double, 3>> magneticField;      ///< Magnetic field [T]
    std::vector<std::array<double, 3>> electricField;      ///< Electric field [V/m]
    std::vector<std::array<double, 3>> currentDensity;     ///< Current density [A/m²]
    std::vector<std::array<double, 3>> lorentzForce;       ///< Lorentz force [N/m³]
    
    // Convergence information
    int iterations = 0;                      ///< Number of iterations
    double residual = 0.0;                   ///< Final residual
    bool converged = false;                  ///< Convergence status
    
    // Energy quantities
    double magneticEnergy = 0.0;             ///< Magnetic energy [J]
    double electricEnergy = 0.0;             ///< Electric energy [J]
    
    MagnetodynamicsResults() = default;
};

/**
 * @brief Magnetodynamics solver for low-frequency electromagnetic problems
 * 
 * This solver implements the A-φ formulation for magnetodynamics problems,
 * suitable for low-frequency applications like transformers, motors, and
 * electromagnetic actuators.
 */
class MagnetodynamicsSolver {
private:
    std::shared_ptr<Mesh> mesh;
    MaterialDatabase materialDB;
    MagnetodynamicsParameters parameters;
    
    // System matrices
    std::shared_ptr<CRSMatrix> stiffnessMatrix;
    std::shared_ptr<CRSMatrix> massMatrix;
    std::shared_ptr<Vector> rhsVector;
    
    // Boundary conditions
    std::vector<EMBoundaryCondition> boundaryConditions;
    
public:
    /**
     * @brief Constructor
     */
    MagnetodynamicsSolver(std::shared_ptr<Mesh> m = nullptr)
        : mesh(m) {
        materialDB.createPredefinedMaterials();
    }
    
    /**
     * @brief Set the mesh for the simulation
     */
    void setMesh(std::shared_ptr<Mesh> m) {
        mesh = m;
    }
    
    /**
     * @brief Set solver parameters
     */
    void setParameters(const MagnetodynamicsParameters& params) {
        parameters = params;
    }
    
    /**
     * @brief Add a boundary condition
     */
    void addBoundaryCondition(const EMBoundaryCondition& bc) {
        boundaryConditions.push_back(bc);
    }
    
    /**
     * @brief Assemble the system matrices for magnetodynamics
     */
    void assembleSystem() {
        if (!mesh) {
            throw std::runtime_error("Mesh not set for assembly");
        }
        
        int nNodes = mesh->getNodes().size();
        int dofPerNode = 4; // A_x, A_y, A_z, φ for 3D problems
        
        // Initialize system matrices
        stiffnessMatrix = std::make_shared<CRSMatrix>(nNodes * dofPerNode, nNodes * dofPerNode);
        massMatrix = std::make_shared<CRSMatrix>(nNodes * dofPerNode, nNodes * dofPerNode);
        rhsVector = std::make_shared<Vector>(nNodes * dofPerNode);
        
        // Assemble element contributions
        for (const auto& element : mesh->getElements()) {
            assembleElementContribution(element);
        }
        
        // Apply boundary conditions
        applyBoundaryConditions();
    }
    
    /**
     * @brief Solve the magnetodynamics problem
     */
    MagnetodynamicsResults solve() {
        MagnetodynamicsResults results;
        
        if (!stiffnessMatrix || !rhsVector) {
            assembleSystem();
        }
        
        // Create solution vector
        Vector solution(rhsVector->size());
        
        // Solve the linear system
        ConjugateGradientSolver cgSolver;
        cgSolver.setMaxIterations(parameters.maxIterations);
        cgSolver.setTolerance(parameters.tolerance);
        
        results.converged = cgSolver.solve(*stiffnessMatrix, solution, *rhsVector);
        results.iterations = cgSolver.getIterationCount();
        results.residual = cgSolver.getResidual();
        
        // Extract solution
        extractSolution(results, solution);
        
        // Calculate derived fields if requested
        if (parameters.calculateMagneticField) {
            calculateMagneticField(results);
        }
        if (parameters.calculateElectricField) {
            calculateElectricField(results);
        }
        if (parameters.calculateCurrentDensity) {
            calculateCurrentDensity(results);
        }
        if (parameters.calculateLorentzForce) {
            calculateLorentzForce(results);
        }
        
        return results;
    }
    
    /**
     * @brief Get the material database
     */
    MaterialDatabase& getMaterialDatabase() {
        return materialDB;
    }
    
private:
    /**
     * @brief Assemble element contribution to system matrices
     */
    void assembleElementContribution(const Element& element) {
        // Get element nodes
        auto nodes = element.getNodes();
        
        // Get material properties
        auto material = materialDB.getMaterial(element.getMaterialName());
        
        // Convert material properties for element matrix
        ElementMatrix::MaterialProperties matProps;
        matProps.conductivity = material.conductivity;
        // Additional properties would be set here
        
        // Assemble element matrices based on problem type
        if (parameters.isHarmonic) {
            assembleHarmonicElementMatrices(element, nodes, material, matProps);
        } else {
            assembleTransientElementMatrices(element, nodes, material, matProps);
        }
    }
    
    /**
     * @brief Assemble matrices for harmonic analysis
     */
    void assembleHarmonicElementMatrices(const Element& element, 
                                        const std::vector<Node>& nodes,
                                        const ElectromagneticMaterial& material,
                                        const ElementMatrix::MaterialProperties& matProps) {
        // For harmonic analysis, we need complex matrices
        // This is a simplified implementation
        
        // Assemble conductivity matrix (real part)
        auto conductivityMatrix = ElementMatrix::computeConductivityMatrix(
            element.getType(), nodes, matProps, 2);
        
        // Assemble other matrices as needed
        // In practice, this would involve complex arithmetic
    }
    
    /**
     * @brief Assemble matrices for transient analysis
     */
    void assembleTransientElementMatrices(const Element& element,
                                         const std::vector<Node>& nodes,
                                         const ElectromagneticMaterial& material,
                                         const ElementMatrix::MaterialProperties& matProps) {
        // For transient analysis, assemble standard matrices
        
        // Mass matrix for time derivative terms
        auto massMat = ElementMatrix::computeMassMatrix(
            element.getType(), nodes, matProps, 2);
        
        // Stiffness matrix for diffusion terms
        auto stiffnessMat = ElementMatrix::computeStiffnessMatrix(
            element.getType(), nodes, matProps, 2);
        
        // Add contributions to global matrices
        addElementMatrixToGlobal(massMat, element, *massMatrix);
        addElementMatrixToGlobal(stiffnessMat, element, *stiffnessMatrix);
    }
    
    /**
     * @brief Add element matrix to global matrix
     */
    void addElementMatrixToGlobal(const std::vector<std::vector<double>>& elementMatrix,
                                 const Element& element,
                                 CRSMatrix& globalMatrix) {
        auto nodeIndices = element.getNodeIndices();
        int nNodes = nodeIndices.size();
        
        for (int i = 0; i < nNodes; ++i) {
            for (int j = 0; j < nNodes; ++j) {
                // For simplicity, assuming scalar potential
                // In practice, this would handle vector potentials
                globalMatrix.add(nodeIndices[i], nodeIndices[j], elementMatrix[i][j]);
            }
        }
    }
    
    /**
     * @brief Apply boundary conditions to the system
     */
    void applyBoundaryConditions() {
        for (const auto& bc : boundaryConditions) {
            applyBoundaryCondition(bc);
        }
    }
    
    /**
     * @brief Apply a single boundary condition
     */
    void applyBoundaryCondition(const EMBoundaryCondition& bc) {
        // Implementation depends on boundary condition type
        switch (bc.type) {
            case EM_BoundaryType::DIRICHLET:
                applyDirichletBC(bc);
                break;
            case EM_BoundaryType::NEUMANN:
                applyNeumannBC(bc);
                break;
            case EM_BoundaryType::ROBIN:
                applyRobinBC(bc);
                break;
            default:
                // Handle other boundary types
                break;
        }
    }
    
    /**
     * @brief Apply Dirichlet boundary condition
     */
    void applyDirichletBC(const EMBoundaryCondition& bc) {
        for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
            int nodeIdx = bc.nodeIndices[i];
            double value = (i < bc.values.size()) ? bc.values[i] : 0.0;
            
            // Set diagonal to 1 and RHS to value
            stiffnessMatrix->set(nodeIdx, nodeIdx, 1.0);
            rhsVector->set(nodeIdx, value);
            
            // Zero out row and column
            for (int j = 0; j < stiffnessMatrix->getCols(); ++j) {
                if (j != nodeIdx) {
                    stiffnessMatrix->set(nodeIdx, j, 0.0);
                    stiffnessMatrix->set(j, nodeIdx, 0.0);
                }
            }
        }
    }
    
    /**
     * @brief Apply Neumann boundary condition
     */
    void applyNeumannBC(const EMBoundaryCondition& bc) {
        // Neumann conditions are natural boundary conditions in FEM
        // They are automatically satisfied and don't require modification
        // of the system matrices for standard formulations
    }
    
    /**
     * @brief Apply Robin boundary condition
     */
    void applyRobinBC(const EMBoundaryCondition& bc) {
        // Robin conditions require modification of both stiffness matrix and RHS
        // This is a simplified implementation
        for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
            int nodeIdx = bc.nodeIndices[i];
            double alpha = (i < bc.values.size()) ? bc.values[i] : 0.0;
            double beta = (i < bc.values.size()) ? bc.values[i] : 0.0;
            
            // Modify stiffness matrix and RHS
            double currentDiag = stiffnessMatrix->get(nodeIdx, nodeIdx);
            stiffnessMatrix->set(nodeIdx, nodeIdx, currentDiag + alpha);
            
            double currentRhs = rhsVector->get(nodeIdx);
            rhsVector->set(nodeIdx, currentRhs + beta);
        }
    }
    
    /**
     * @brief Extract solution from the solved system
     */
    void extractSolution(MagnetodynamicsResults& results, const Vector& solution) {
        int nNodes = mesh->getNodes().size();
        results.potentialReal.resize(nNodes);
        
        for (int i = 0; i < nNodes; ++i) {
            results.potentialReal[i] = solution.get(i);
        }
        
        // For harmonic analysis, extract imaginary part as well
        if (parameters.isHarmonic) {
            results.potentialImag.resize(nNodes);
            for (int i = 0; i < nNodes; ++i) {
                results.potentialImag[i] = solution.get(i + nNodes);
            }
        }
    }
    
    /**
     * @brief Calculate magnetic field from potential
     */
    void calculateMagneticField(MagnetodynamicsResults& results) {
        int nNodes = mesh->getNodes().size();
        results.magneticField.resize(nNodes, {0.0, 0.0, 0.0});
        
        // Simplified calculation - in practice, this would use curl of A
        for (const auto& element : mesh->getElements()) {
            calculateElementMagneticField(element, results);
        }
    }
    
    /**
     * @brief Calculate magnetic field for a single element
     */
    void calculateElementMagneticField(const Element& element, MagnetodynamicsResults& results) {
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Get shape functions and derivatives at integration points
        ShapeFunctions shapeFunc;
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), 
                                                          {0.0, 0.0, 0.0}); // At element center
        
        if (!shapeResult.jacobianMatrix.empty()) {
            // Calculate magnetic field B = ∇ × A
            // For 2D problems: B_z = ∂A_y/∂x - ∂A_x/∂y
            // For 3D problems: B = ∇ × A
            
            // This is a simplified implementation
            // In practice, we would compute the curl of the vector potential
            
            // For now, use a simple gradient-based approximation
            if (nodes.size() >= 2) {
                double dx = nodes[1].x - nodes[0].x;
                double dy = nodes[1].y - nodes[0].y;
                
                if (std::abs(dx) > 1e-10 || std::abs(dy) > 1e-10) {
                    // Calculate approximate magnetic field
                    for (int nodeIdx : nodeIndices) {
                        results.magneticField[nodeIdx][2] = (results.potentialReal[nodeIdx] / 
                                                            std::sqrt(dx*dx + dy*dy));
                    }
                }
            }
        }
    }
    
    /**
     * @brief Calculate electric field
     */
    void calculateElectricField(MagnetodynamicsResults& results) {
        int nNodes = mesh->getNodes().size();
        results.electricField.resize(nNodes, {0.0, 0.0, 0.0});
        
        // E = -∇φ - ∂A/∂t for transient problems
        // For harmonic: E = -∇φ - jωA
        
        // Simplified implementation using potential gradients
        for (const auto& element : mesh->getElements()) {
            calculateElementElectricField(element, results);
        }
    }
    
    /**
     * @brief Calculate electric field for a single element
     */
    void calculateElementElectricField(const Element& element, MagnetodynamicsResults& results) {
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate electric field using potential gradient
        if (nodes.size() >= 2) {
            double dx = nodes[1].x - nodes[0].x;
            double dy = nodes[1].y - nodes[0].y;
            
            if (std::abs(dx) > 1e-10) {
                for (int nodeIdx : nodeIndices) {
                    results.electricField[nodeIdx][0] = -results.potentialReal[nodeIdx] / dx;
                }
            }
            if (std::abs(dy) > 1e-10) {
                for (int nodeIdx : nodeIndices) {
                    results.electricField[nodeIdx][1] = -results.potentialReal[nodeIdx] / dy;
                }
            }
        }
    }
    
    /**
     * @brief Calculate current density
     */
    void calculateCurrentDensity(MagnetodynamicsResults& results) {
        int nNodes = mesh->getNodes().size();
        results.currentDensity.resize(nNodes, {0.0, 0.0, 0.0});
        
        // J = σE for conductive materials
        // For harmonic: J = (σ + jωε)E
        
        for (const auto& element : mesh->getElements()) {
            calculateElementCurrentDensity(element, results);
        }
    }
    
    /**
     * @brief Calculate current density for a single element
     */
    void calculateElementCurrentDensity(const Element& element, MagnetodynamicsResults& results) {
        auto material = materialDB.getMaterial(element.getMaterialName());
        auto nodeIndices = element.getNodeIndices();
        
        // J = σE
        for (int nodeIdx : nodeIndices) {
            for (int i = 0; i < 3; ++i) {
                results.currentDensity[nodeIdx][i] = material.conductivity * 
                                                    results.electricField[nodeIdx][i];
            }
        }
    }
    
    /**
     * @brief Calculate Lorentz force
     */
    void calculateLorentzForce(MagnetodynamicsResults& results) {
        int nNodes = mesh->getNodes().size();
        results.lorentzForce.resize(nNodes, {0.0, 0.0, 0.0});
        
        // F = J × B (Lorentz force density)
        
        for (const auto& element : mesh->getElements()) {
            calculateElementLorentzForce(element, results);
        }
    }
    
    /**
     * @brief Calculate Lorentz force for a single element
     */
    void calculateElementLorentzForce(const Element& element, MagnetodynamicsResults& results) {
        auto nodeIndices = element.getNodeIndices();
        
        // F = J × B
        for (int nodeIdx : nodeIndices) {
            // Cross product: F_x = J_y * B_z - J_z * B_y
            results.lorentzForce[nodeIdx][0] = 
                results.currentDensity[nodeIdx][1] * results.magneticField[nodeIdx][2] -
                results.currentDensity[nodeIdx][2] * results.magneticField[nodeIdx][1];
            
            // F_y = J_z * B_x - J_x * B_z
            results.lorentzForce[nodeIdx][1] = 
                results.currentDensity[nodeIdx][2] * results.magneticField[nodeIdx][0] -
                results.currentDensity[nodeIdx][0] * results.magneticField[nodeIdx][2];
            
            // F_z = J_x * B_y - J_y * B_x
            results.lorentzForce[nodeIdx][2] = 
                results.currentDensity[nodeIdx][0] * results.magneticField[nodeIdx][1] -
                results.currentDensity[nodeIdx][1] * results.magneticField[nodeIdx][0];
        }
    }
};

} // namespace elmer