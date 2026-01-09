#include "StructuralSolver.h"
#include "IterativeSolver.h"
#include <iostream>
#include <cmath>

namespace elmer {

StructuralSolver::StructuralSolver(std::shared_ptr<Mesh> meshPtr) 
    : mesh(meshPtr) {
    initialize();
}

void StructuralSolver::setMesh(std::shared_ptr<Mesh> meshPtr) {
    mesh = meshPtr;
    initialize();
}

void StructuralSolver::setParameters(const StructuralParameters& params) {
    parameters = params;
}

void StructuralSolver::addBoundaryCondition(const StructuralBoundaryCondition& bc) {
    boundaryConditions.push_back(bc);
}

void StructuralSolver::clearBoundaryConditions() {
    boundaryConditions.clear();
}

void StructuralSolver::setTemperatureField(const std::vector<double>& temperature) {
    temperatureField = temperature;
}

void StructuralSolver::initialize() {
    if (!mesh) {
        return;
    }
    
    int nNodes = mesh->getNodes().size();
    int nDOF = nNodes * 3; // 3 degrees of freedom per node (x, y, z)
    
    // Initialize system matrices and vectors
    stiffnessMatrix = std::make_shared<CRSMatrix>(nDOF, nDOF);
    massMatrix = std::make_shared<CRSMatrix>(nDOF, nDOF);
    forceVector = std::make_shared<Vector>(nDOF);
    boundaryVector = std::make_shared<Vector>(nDOF);
    displacementVector = std::make_shared<Vector>(nDOF);
    previousDisplacement = std::make_shared<Vector>(nDOF);
    velocityVector = std::make_shared<Vector>(nDOF);
    accelerationVector = std::make_shared<Vector>(nDOF);
    
    // Initialize with zero values
    displacementVector->clear();
    previousDisplacement->clear();
    velocityVector->clear();
    accelerationVector->clear();
    
    // Create iterative solver
    solver = std::make_unique<ConjugateGradientSolver>();
    solver->setTolerance(parameters.tolerance);
    solver->setMaxIterations(parameters.maxIterations);
}

void StructuralSolver::assembleSystem() {
    if (!isSystemReady()) {
        throw std::runtime_error("System not ready for assembly");
    }
    
    // Clear existing matrices and vectors
    stiffnessMatrix->clear();
    massMatrix->clear();
    forceVector->clear();
    boundaryVector->clear();
    
    // Assemble element contributions
    for (const auto& element : mesh->getElements()) {
        auto elementStiffness = calculateElementStiffnessMatrix(element);
        auto elementMass = calculateElementMassMatrix(element);
        auto elementForce = calculateElementForceVector(element);
        
        auto nodeIndices = element.getNodeIndices();
        int nNodes = nodeIndices.size();
        
        // Add to global matrices (3 DOF per node)
        for (int i = 0; i < nNodes; ++i) {
            int nodeI = nodeIndices[i];
            for (int j = 0; j < nNodes; ++j) {
                int nodeJ = nodeIndices[j];
                
                // Add 3x3 submatrices for each node pair
                for (int dofI = 0; dofI < 3; ++dofI) {
                    for (int dofJ = 0; dofJ < 3; ++dofJ) {
                        int globalRow = nodeI * 3 + dofI;
                        int globalCol = nodeJ * 3 + dofJ;
                        int localRow = i * 3 + dofI;
                        int localCol = j * 3 + dofJ;
                        
                        stiffnessMatrix->add(globalRow, globalCol, elementStiffness[localRow][localCol]);
                        massMatrix->add(globalRow, globalCol, elementMass[localRow][localCol]);
                    }
                }
            }
            
            // Add force vector contributions
            for (int dof = 0; dof < 3; ++dof) {
                int globalDOF = nodeI * 3 + dof;
                int localDOF = i * 3 + dof;
                (*forceVector)[globalDOF] += elementForce[localDOF];
            }
        }
    }
    
    // Apply boundary conditions
    applyBoundaryConditions();
}

StructuralResults StructuralSolver::solve() {
    if (!isSystemReady()) {
        throw std::runtime_error("System not ready for solving");
    }
    
    assembleSystem();
    
    StructuralResults results;
    
    if (parameters.isTransient) {
        results = solveTransient();
    } else {
        // Static solution
        auto rhs = *forceVector + *boundaryVector;
        
        // Solve K * u = F
        auto solution = solver->solve(*stiffnessMatrix, rhs);
        
        // Convert solution to displacement field
        int nNodes = mesh->getNodes().size();
        results.displacement.resize(nNodes);
        
        for (int i = 0; i < nNodes; ++i) {
            for (int dof = 0; dof < 3; ++dof) {
                results.displacement[i][dof] = solution[i * 3 + dof];
            }
        }
        
        results.iterations = solver->getIterations();
        results.residual = solver->getResidual();
        results.converged = solver->hasConverged();
    }
    
    // Calculate derived fields
    results.stress = calculateStresses(results.displacement);
    results.strain = calculateStrains(results.displacement);
    results.reactionForce = calculateReactionForces(results.displacement);
    
    // Calculate global quantities
    results.totalStrainEnergy = 0.0;
    results.maxVonMisesStress = 0.0;
    results.maxDisplacement = 0.0;
    
    for (const auto& disp : results.displacement) {
        double mag = std::sqrt(disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2]);
        if (mag > results.maxDisplacement) {
            results.maxDisplacement = mag;
        }
    }
    
    for (const auto& stress : results.stress) {
        double vonMises = calculateVonMisesStress(stress);
        if (vonMises > results.maxVonMisesStress) {
            results.maxVonMisesStress = vonMises;
        }
    }
    
    return results;
}

StructuralResults StructuralSolver::solveTransient() {
    StructuralResults results;
    results.timePoints.resize(parameters.timeSteps + 1);
    results.displacementHistory.resize(parameters.timeSteps + 1);
    
    double time = 0.0;
    results.timePoints[0] = time;
    results.displacementHistory[0] = displacementVector->getData();
    
    // For simplicity, use Newmark-beta method for time integration
    double beta = 0.25;  // Newmark parameter
    double gamma = 0.5;  // Newmark parameter
    
    for (int step = 1; step <= parameters.timeSteps; ++step) {
        time += parameters.timeStep;
        updateTimeStep(step, time);
        
        // Newmark time integration
        // This is a simplified implementation
        // In practice, this would be more complex
        
        auto systemMatrix = *stiffnessMatrix + 
                           (1.0/(beta*parameters.timeStep*parameters.timeStep)) * (*massMatrix);
        
        auto rhs = *forceVector + *boundaryVector + 
                  (*massMatrix) * ((1.0/(beta*parameters.timeStep*parameters.timeStep)) * (*displacementVector) +
                                 (1.0/(beta*parameters.timeStep)) * (*velocityVector) +
                                 (1.0/(2.0*beta) - 1.0) * (*accelerationVector));
        
        auto newDisplacement = solver->solve(systemMatrix, rhs);
        
        // Update velocity and acceleration
        auto newAcceleration = (1.0/(beta*parameters.timeStep*parameters.timeStep)) * 
                              (newDisplacement - *displacementVector) -
                              (1.0/(beta*parameters.timeStep)) * (*velocityVector) -
                              (1.0/(2.0*beta) - 1.0) * (*accelerationVector);
        
        auto newVelocity = *velocityVector + 
                          parameters.timeStep * ((1.0 - gamma) * (*accelerationVector) + 
                                               gamma * newAcceleration);
        
        *displacementVector = newDisplacement;
        *velocityVector = newVelocity;
        *accelerationVector = newAcceleration;
        
        results.timePoints[step] = time;
        results.displacementHistory[step] = displacementVector->getData();
    }
    
    // Convert final displacement to results
    int nNodes = mesh->getNodes().size();
    results.displacement.resize(nNodes);
    
    for (int i = 0; i < nNodes; ++i) {
        for (int dof = 0; dof < 3; ++dof) {
            results.displacement[i][dof] = (*displacementVector)[i * 3 + dof];
        }
    }
    
    results.iterations = solver->getIterations();
    results.residual = solver->getResidual();
    results.converged = solver->hasConverged();
    
    return results;
}

std::vector<std::array<double, 6>> StructuralSolver::calculateStresses(
    const std::vector<std::array<double, 3>>& displacement) {
    
    std::vector<std::array<double, 6>> stresses(displacement.size(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    
    ShapeFunctions shapeFunc;
    
    for (const auto& element : mesh->getElements()) {
        auto material = materialDB.getMaterial(element.getMaterialName());
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate at element center
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), {0.0, 0.0, 0.0});
        
        if (!shapeResult.shapeFunctionDerivatives.empty()) {
            // Calculate strain-displacement matrix (B matrix)
            auto B = calculateStrainDisplacementMatrix(shapeResult.shapeFunctionDerivatives);
            
            // Calculate constitutive matrix (D matrix)
            auto D = calculateConstitutiveMatrix(material);
            
            // Calculate element displacement vector
            std::vector<double> elementDisplacement;
            for (int nodeIdx : nodeIndices) {
                for (int dof = 0; dof < 3; ++dof) {
                    elementDisplacement.push_back(displacement[nodeIdx][dof]);
                }
            }
            
            // Calculate strain: ε = B * u
            std::vector<double> strain(6, 0.0);
            for (int i = 0; i < 6; ++i) {
                for (size_t j = 0; j < elementDisplacement.size(); ++j) {
                    strain[i] += B[i][j] * elementDisplacement[j];
                }
            }
            
            // Calculate stress: σ = D * ε
            std::vector<double> stress(6, 0.0);
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    stress[i] += D[i][j] * strain[j];
                }
            }
            
            // Distribute to nodes (simple averaging)
            for (int nodeIdx : nodeIndices) {
                for (int i = 0; i < 6; ++i) {
                    stresses[nodeIdx][i] += stress[i];
                }
            }
        }
    }
    
    // Average stresses at nodes
    for (auto& stress : stresses) {
        for (int i = 0; i < 6; ++i) {
            stress[i] /= mesh->getNodes().size();
        }
    }
    
    return stresses;
}

std::vector<std::array<double, 6>> StructuralSolver::calculateStrains(
    const std::vector<std::array<double, 3>>& displacement) {
    
    std::vector<std::array<double, 6>> strains(displacement.size(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    
    ShapeFunctions shapeFunc;
    
    for (const auto& element : mesh->getElements()) {
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate at element center
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), {0.0, 0.0, 0.0});
        
        if (!shapeResult.shapeFunctionDerivatives.empty()) {
            // Calculate strain-displacement matrix (B matrix)
            auto B = calculateStrainDisplacementMatrix(shapeResult.shapeFunctionDerivatives);
            
            // Calculate element displacement vector
            std::vector<double> elementDisplacement;
            for (int nodeIdx : nodeIndices) {
                for (int dof = 0; dof < 3; ++dof) {
                    elementDisplacement.push_back(displacement[nodeIdx][dof]);
                }
            }
            
            // Calculate strain: ε = B * u
            std::vector<double> strain(6, 0.0);
            for (int i = 0; i < 6; ++i) {
                for (size_t j = 0; j < elementDisplacement.size(); ++j) {
                    strain[i] += B[i][j] * elementDisplacement[j];
                }
            }
            
            // Distribute to nodes (simple averaging)
            for (int nodeIdx : nodeIndices) {
                for (int i = 0; i < 6; ++i) {
                    strains[nodeIdx][i] += strain[i];
                }
            }
        }
    }
    
    // Average strains at nodes
    for (auto& strain : strains) {
        for (int i = 0; i < 6; ++i) {
            strain[i] /= mesh->getNodes().size();
        }
    }
    
    return strains;
}

std::vector<std::array<double, 3>> StructuralSolver::calculateReactionForces(
    const std::vector<std::array<double, 3>>& displacement) {
    
    std::vector<std::array<double, 3>> reactionForces(displacement.size(), {0.0, 0.0, 0.0});
    
    // Calculate reaction forces as R = K * u - F
    // This is a simplified implementation
    
    for (const auto& element : mesh->getElements()) {
        auto elementStiffness = calculateElementStiffnessMatrix(element);
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate element displacement vector
        std::vector<double> elementDisplacement;
        for (int nodeIdx : nodeIndices) {
            for (int dof = 0; dof < 3; ++dof) {
                elementDisplacement.push_back(displacement[nodeIdx][dof]);
            }
        }
        
        // Calculate element force: F_e = K_e * u_e
        std::vector<double> elementForce(elementDisplacement.size(), 0.0);
        for (size_t i = 0; i < elementDisplacement.size(); ++i) {
            for (size_t j = 0; j < elementDisplacement.size(); ++j) {
                elementForce[i] += elementStiffness[i][j] * elementDisplacement[j];
            }
        }
        
        // Distribute to nodes
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int nodeIdx = nodeIndices[i];
            for (int dof = 0; dof < 3; ++dof) {
                reactionForces[nodeIdx][dof] += elementForce[i * 3 + dof];
            }
        }
    }
    
    return reactionForces;
}

void StructuralSolver::applyBoundaryConditions() {
    for (const auto& bc : boundaryConditions) {
        switch (bc.type) {
            case Structural_BoundaryType::DISPLACEMENT: {
                // Fixed displacement - modify matrix and RHS
                for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
                    int nodeIdx = bc.nodeIndices[i];
                    
                    for (int dof = 0; dof < 3; ++dof) {
                        int globalDOF = nodeIdx * 3 + dof;
                        double dispValue = bc.values[i][dof];
                        
                        // Set diagonal to 1 and RHS to displacement value
                        stiffnessMatrix->set(globalDOF, globalDOF, 1.0);
                        (*boundaryVector)[globalDOF] = dispValue;
                        
                        // Zero out off-diagonal entries in this row
                        for (int j = 0; j < stiffnessMatrix->getCols(); ++j) {
                            if (j != globalDOF) {
                                stiffnessMatrix->set(globalDOF, j, 0.0);
                            }
                        }
                    }
                }
                break;
            }
            case Structural_BoundaryType::FORCE: {
                // Applied force - add to RHS
                for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
                    int nodeIdx = bc.nodeIndices[i];
                    
                    for (int dof = 0; dof < 3; ++dof) {
                        int globalDOF = nodeIdx * 3 + dof;
                        double forceValue = bc.values[i][dof];
                        (*boundaryVector)[globalDOF] += forceValue;
                    }
                }
                break;
            }
            default:
                // Other boundary types not implemented yet
                break;
        }
    }
}

std::vector<std::vector<double>> StructuralSolver::calculateElementStiffnessMatrix(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    int nDOF = nNodes * 3;
    
    std::vector<std::vector<double>> kMatrix(nDOF, std::vector<double>(nDOF, 0.0));
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // Calculate strain-displacement matrix (B matrix)
            auto B = calculateStrainDisplacementMatrix(shapeResult.shapeFunctionDerivatives);
            
            // Calculate constitutive matrix (D matrix)
            auto D = calculateConstitutiveMatrix(material);
            
            // Calculate element stiffness matrix: K_e = ∫ B^T * D * B dΩ
            // K_ij = B_i^T * D * B_j * weight
            
            for (int i = 0; i < nDOF; ++i) {
                for (int j = 0; j < nDOF; ++j) {
                    double sum = 0.0;
                    for (int k = 0; k < 6; ++k) {
                        for (int l = 0; l < 6; ++l) {
                            sum += B[k][i] * D[k][l] * B[l][j];
                        }
                    }
                    kMatrix[i][j] += sum * weight;
                }
            }
        }
    }
    
    return kMatrix;
}

std::vector<std::vector<double>> StructuralSolver::calculateElementMassMatrix(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    int nDOF = nNodes * 3;
    
    std::vector<std::vector<double>> mMatrix(nDOF, std::vector<double>(nDOF, 0.0));
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // Consistent mass matrix: M_ij = ∫ ρ * N_i * N_j dΩ
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    double massTerm = material.density * shapeResult.shapeFunctions[i] * 
                                    shapeResult.shapeFunctions[j] * weight;
                    
                    // Add to diagonal 3x3 block
                    for (int dof = 0; dof < 3; ++dof) {
                        int row = i * 3 + dof;
                        int col = j * 3 + dof;
                        mMatrix[row][col] += massTerm;
                    }
                }
            }
        }
    }
    
    return mMatrix;
}

std::vector<double> StructuralSolver::calculateElementForceVector(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    int nDOF = nNodes * 3;
    
    std::vector<double> forceVector(nDOF, 0.0);
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // Body forces (gravity)
            if (parameters.includeBodyForces) {
                for (int i = 0; i < nNodes; ++i) {
                    double bodyForceTerm = material.density * shapeResult.shapeFunctions[i] * weight;
                    
                    for (int dof = 0; dof < 3; ++dof) {
                        forceVector[i * 3 + dof] += bodyForceTerm * parameters.gravity[dof];
                    }
                }
            }
        }
    }
    
    return forceVector;
}

std::vector<double> StructuralSolver::calculateElementThermalExpansion(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    int nDOF = nNodes * 3;
    
    std::vector<double> thermalVector(nDOF, 0.0);
    
    if (!parameters.includeThermalExpansion || temperatureField.empty()) {
        return thermalVector;
    }
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // Calculate average temperature in element
            double avgTemp = 0.0;
            for (int nodeIdx : nodeIndices) {
                avgTemp += temperatureField[nodeIdx];
            }
            avgTemp /= nNodes;
            
            // Reference temperature (could be parameterized)
            double refTemp = 293.15; // 20°C
            double deltaT = avgTemp - refTemp;
            
            // Thermal strain: ε_thermal = α * ΔT
            // Thermal force: F_thermal = ∫ B^T * D * ε_thermal dΩ
            
            // This is a simplified implementation
            // In practice, this would be more complex
            
            auto D = calculateConstitutiveMatrix(material);
            std::vector<double> thermalStrain = {material.thermalExpansion * deltaT, 
                                                material.thermalExpansion * deltaT,
                                                material.thermalExpansion * deltaT,
                                                0.0, 0.0, 0.0};
            
            auto B = calculateStrainDisplacementMatrix(shapeResult.shapeFunctionDerivatives);
            
            // Calculate B^T * D * ε_thermal
            for (int i = 0; i < nDOF; ++i) {
                for (int j = 0; j < 6; ++j) {
                    for (int k = 0; k < 6; ++k) {
                        thermalVector[i] += B[j][i] * D[j][k] * thermalStrain[k] * weight;
                    }
                }
            }
        }
    }
    
    return thermalVector;
}

bool StructuralSolver::isSystemReady() const {
    return mesh != nullptr && stiffnessMatrix != nullptr && 
           forceVector != nullptr && solver != nullptr;
}

void StructuralSolver::updateTimeStep(int step, double time) {
    // Update boundary conditions or material properties if they are time-dependent
    // For now, this is a placeholder for future time-dependent functionality
}

std::vector<std::vector<double>> StructuralSolver::calculateConstitutiveMatrix(const Material& material) {
    // Constitutive matrix for linear isotropic elasticity (plane stress/3D)
    // This is a simplified implementation for isotropic materials
    
    std::vector<std::vector<double>> D(6, std::vector<double>(6, 0.0));
    
    double E = material.youngsModulus;
    double nu = material.poissonsRatio;
    
    // For 3D isotropic elasticity
    double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    D[0][0] = lambda + 2.0 * mu; D[0][1] = lambda; D[0][2] = lambda;
    D[1][0] = lambda; D[1][1] = lambda + 2.0 * mu; D[1][2] = lambda;
    D[2][0] = lambda; D[2][1] = lambda; D[2][2] = lambda + 2.0 * mu;
    D[3][3] = mu;
    D[4][4] = mu;
    D[5][5] = mu;
    
    return D;
}

std::vector<std::vector<double>> StructuralSolver::calculateStrainDisplacementMatrix(
    const std::vector<std::vector<double>>& shapeDerivatives) {
    
    int nNodes = shapeDerivatives.size();
    int nDOF = nNodes * 3;
    
    std::vector<std::vector<double>> B(6, std::vector<double>(nDOF, 0.0));
    
    // Strain-displacement matrix for 3D elasticity
    // ε_xx = ∂u/∂x, ε_yy = ∂v/∂y, ε_zz = ∂w/∂z
    // γ_xy = ∂u/∂y + ∂v/∂x, γ_xz = ∂u/∂z + ∂w/∂x, γ_yz = ∂v/∂z + ∂w/∂y
    
    for (int i = 0; i < nNodes; ++i) {
        double dN_dx = shapeDerivatives[i][0];
        double dN_dy = shapeDerivatives[i][1];
        double dN_dz = shapeDerivatives[i][2];
        
        // ε_xx
        B[0][i * 3 + 0] = dN_dx;
        
        // ε_yy
        B[1][i * 3 + 1] = dN_dy;
        
        // ε_zz
        B[2][i * 3 + 2] = dN_dz;
        
        // γ_xy
        B[3][i * 3 + 0] = dN_dy;
        B[3][i * 3 + 1] = dN_dx;
        
        // γ_xz
        B[4][i * 3 + 0] = dN_dz;
        B[4][i * 3 + 2] = dN_dx;
        
        // γ_yz
        B[5][i * 3 + 1] = dN_dz;
        B[5][i * 3 + 2] = dN_dy;
    }
    
    return B;
}

double StructuralSolver::calculateVonMisesStress(const std::array<double, 6>& stress) const {
    // von Mises stress for 3D stress state
    double sxx = stress[0], syy = stress[1], szz = stress[2];
    double txy = stress[3], txz = stress[4], tyz = stress[5];
    
    double vonMises = std::sqrt(0.5 * ((sxx - syy) * (sxx - syy) +
                                     (syy - szz) * (syy - szz) +
                                     (szz - sxx) * (szz - sxx) +
                                     6.0 * (txy * txy + txz * txz + tyz * tyz)));
    
    return vonMises;
}

} // namespace elmer