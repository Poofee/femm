#include "HeatTransferSolver.h"
#include "IterativeSolver.h"
#include "Material.h"
#include <iostream>
#include <cmath>

namespace elmer {

HeatTransferSolver::HeatTransferSolver(std::shared_ptr<Mesh> meshPtr) 
    : mesh(meshPtr) {
    initialize();
}

void HeatTransferSolver::setMesh(std::shared_ptr<Mesh> meshPtr) {
    mesh = meshPtr;
    initialize();
}

void HeatTransferSolver::setParameters(const HeatTransferParameters& params) {
    parameters = params;
}

void HeatTransferSolver::addBoundaryCondition(const ThermalBoundaryCondition& bc) {
    boundaryConditions.push_back(bc);
}

void HeatTransferSolver::clearBoundaryConditions() {
    boundaryConditions.clear();
}

void HeatTransferSolver::initialize() {
    if (!mesh) {
        return;
    }
    
    int nNodes = mesh->getNodes().size();
    
    // Initialize system matrices and vectors
    conductivityMatrix = std::make_shared<CRSMatrix>(nNodes, nNodes);
    capacityMatrix = std::make_shared<CRSMatrix>(nNodes, nNodes);
    heatSourceVector = std::make_shared<Vector>(nNodes);
    boundaryVector = std::make_shared<Vector>(nNodes);
    temperatureVector = std::make_shared<Vector>(nNodes);
    previousTemperature = std::make_shared<Vector>(nNodes);
    
    // Initialize with initial temperature
    for (int i = 0; i < nNodes; ++i) {
        (*temperatureVector)[i] = parameters.initialTemperature;
        (*previousTemperature)[i] = parameters.initialTemperature;
    }
    
    // Create iterative solver
    solver = std::make_unique<ConjugateGradientSolver>();
    solver->setTolerance(parameters.tolerance);
    solver->setMaxIterations(parameters.maxIterations);
}

void HeatTransferSolver::assembleSystem() {
    if (!isSystemReady()) {
        throw std::runtime_error("System not ready for assembly");
    }
    
    // Clear existing matrices and vectors
    conductivityMatrix->clear();
    capacityMatrix->clear();
    heatSourceVector->clear();
    boundaryVector->clear();
    
    // Assemble element contributions
    for (const auto& element : mesh->getElements()) {
        auto elementConductivity = calculateElementConductivityMatrix(element);
        auto elementCapacity = calculateElementCapacityMatrix(element);
        auto elementHeatSource = calculateElementHeatSource(element);
        
        auto nodeIndices = element.getNodeIndices();
        
        // Add to global matrices
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            for (size_t j = 0; j < nodeIndices.size(); ++j) {
                conductivityMatrix->add(nodeIndices[i], nodeIndices[j], elementConductivity[i][j]);
                capacityMatrix->add(nodeIndices[i], nodeIndices[j], elementCapacity[i][j]);
            }
            (*heatSourceVector)[nodeIndices[i]] += elementHeatSource[i];
        }
    }
    
    // Apply boundary conditions
    applyBoundaryConditions();
}

HeatTransferResults HeatTransferSolver::solve() {
    if (!isSystemReady()) {
        throw std::runtime_error("System not ready for solving");
    }
    
    assembleSystem();
    
    HeatTransferResults results;
    
    if (parameters.isTransient) {
        results = solveTransient();
    } else {
        // Steady-state solution
        auto rhs = *heatSourceVector + *boundaryVector;
        
        // Solve K * T = Q
        auto solution = solver->solve(*conductivityMatrix, rhs);
        
        // Store results
        results.temperature = solution.getData();
        results.iterations = solver->getIterations();
        results.residual = solver->getResidual();
        results.converged = solver->hasConverged();
    }
    
    // Calculate derived fields
    results.heatFlux = calculateHeatFlux(results.temperature);
    results.temperatureGradient = calculateTemperatureGradient(results.temperature);
    
    // Calculate global quantities
    results.averageTemperature = 0.0;
    for (double temp : results.temperature) {
        results.averageTemperature += temp;
    }
    results.averageTemperature /= results.temperature.size();
    
    return results;
}

HeatTransferResults HeatTransferSolver::solveTransient() {
    HeatTransferResults results;
    results.timePoints.resize(parameters.timeSteps + 1);
    results.temperatureHistory.resize(parameters.timeSteps + 1);
    
    double time = 0.0;
    results.timePoints[0] = time;
    results.temperatureHistory[0] = temperatureVector->getData();
    
    for (int step = 1; step <= parameters.timeSteps; ++step) {
        time += parameters.timeStep;
        updateTimeStep(step, time);
        
        // For simplicity, use backward Euler method
        // (C + Δt * K) * T^{n+1} = C * T^n + Δt * Q
        
        auto systemMatrix = *capacityMatrix + parameters.timeStep * (*conductivityMatrix);
        auto rhs = (*capacityMatrix) * (*temperatureVector) + 
                   parameters.timeStep * (*heatSourceVector + *boundaryVector);
        
        auto newTemperature = solver->solve(systemMatrix, rhs);
        *temperatureVector = newTemperature;
        
        results.timePoints[step] = time;
        results.temperatureHistory[step] = temperatureVector->getData();
    }
    
    results.temperature = temperatureVector->getData();
    results.iterations = solver->getIterations();
    results.residual = solver->getResidual();
    results.converged = solver->hasConverged();
    
    return results;
}

std::vector<std::array<double, 3>> HeatTransferSolver::calculateHeatFlux(const std::vector<double>& temperature) {
    std::vector<std::array<double, 3>> heatFlux(temperature.size(), {0.0, 0.0, 0.0});
    
    ShapeFunctions shapeFunc;
    
    for (const auto& element : mesh->getElements()) {
        auto material = materialDB.getMaterial(element.getMaterialName());
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate at element center
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), {0.0, 0.0, 0.0});
        
        if (!shapeResult.shapeFunctionDerivatives.empty()) {
            // Calculate temperature gradient
            std::array<double, 3> gradT = {0.0, 0.0, 0.0};
            for (size_t i = 0; i < nodeIndices.size(); ++i) {
                for (int dim = 0; dim < 3; ++dim) {
                    gradT[dim] += temperature[nodeIndices[i]] * 
                                 shapeResult.shapeFunctionDerivatives[i][dim];
                }
            }
            
            // Calculate heat flux: q = -k * ∇T
            for (int nodeIdx : nodeIndices) {
                for (int dim = 0; dim < 3; ++dim) {
                    heatFlux[nodeIdx][dim] = -material.thermalConductivity * gradT[dim];
                }
            }
        }
    }
    
    return heatFlux;
}

std::vector<std::array<double, 3>> HeatTransferSolver::calculateTemperatureGradient(const std::vector<double>& temperature) {
    std::vector<std::array<double, 3>> gradT(temperature.size(), {0.0, 0.0, 0.0});
    
    ShapeFunctions shapeFunc;
    
    for (const auto& element : mesh->getElements()) {
        auto nodes = element.getNodes();
        auto nodeIndices = element.getNodeIndices();
        
        // Calculate at element center
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), {0.0, 0.0, 0.0});
        
        if (!shapeResult.shapeFunctionDerivatives.empty()) {
            // Calculate temperature gradient
            std::array<double, 3> elementGradT = {0.0, 0.0, 0.0};
            for (size_t i = 0; i < nodeIndices.size(); ++i) {
                for (int dim = 0; dim < 3; ++dim) {
                    elementGradT[dim] += temperature[nodeIndices[i]] * 
                                        shapeResult.shapeFunctionDerivatives[i][dim];
                }
            }
            
            // Distribute to nodes
            for (int nodeIdx : nodeIndices) {
                for (int dim = 0; dim < 3; ++dim) {
                    gradT[nodeIdx][dim] += elementGradT[dim];
                }
            }
        }
    }
    
    // Average gradients at nodes
    for (auto& grad : gradT) {
        // Simple averaging - could be improved with proper weighting
        for (int dim = 0; dim < 3; ++dim) {
            grad[dim] /= mesh->getNodes().size();
        }
    }
    
    return gradT;
}

void HeatTransferSolver::applyBoundaryConditions() {
    for (const auto& bc : boundaryConditions) {
        switch (bc.type) {
            case Thermal_BoundaryType::DIRICHLET: {
                // Fixed temperature - modify matrix and RHS
                for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
                    int nodeIdx = bc.nodeIndices[i];
                    double tempValue = bc.values[i];
                    
                    // Set diagonal to 1 and RHS to temperature value
                    conductivityMatrix->set(nodeIdx, nodeIdx, 1.0);
                    (*boundaryVector)[nodeIdx] = tempValue;
                    
                    // Zero out off-diagonal entries in this row
                    for (int j = 0; j < conductivityMatrix->getCols(); ++j) {
                        if (j != nodeIdx) {
                            conductivityMatrix->set(nodeIdx, j, 0.0);
                        }
                    }
                }
                break;
            }
            case Thermal_BoundaryType::NEUMANN: {
                // Fixed heat flux - add to RHS
                for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
                    int nodeIdx = bc.nodeIndices[i];
                    double fluxValue = bc.values[i];
                    (*boundaryVector)[nodeIdx] += fluxValue;
                }
                break;
            }
            case Thermal_BoundaryType::ROBIN: {
                // Convective boundary: q = h * (T - T_inf)
                // This adds h * T_inf to RHS and h to matrix diagonal
                for (size_t i = 0; i < bc.nodeIndices.size(); ++i) {
                    int nodeIdx = bc.nodeIndices[i];
                    double h = bc.heatTransferCoefficient;
                    double T_inf = bc.ambientTemperature;
                    
                    conductivityMatrix->add(nodeIdx, nodeIdx, h);
                    (*boundaryVector)[nodeIdx] += h * T_inf;
                }
                break;
            }
            default:
                // Other boundary types not implemented yet
                break;
        }
    }
}

std::vector<std::vector<double>> HeatTransferSolver::calculateElementConductivityMatrix(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    
    std::vector<std::vector<double>> kMatrix(nNodes, std::vector<double>(nNodes, 0.0));
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // K_ij = ∫ k * ∇N_i · ∇N_j dΩ
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    double dotProduct = 0.0;
                    for (int dim = 0; dim < 3; ++dim) {
                        dotProduct += shapeResult.shapeFunctionDerivatives[i][dim] * 
                                     shapeResult.shapeFunctionDerivatives[j][dim];
                    }
                    kMatrix[i][j] += material.thermalConductivity * dotProduct * weight;
                }
            }
        }
    }
    
    return kMatrix;
}

std::vector<std::vector<double>> HeatTransferSolver::calculateElementCapacityMatrix(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    
    std::vector<std::vector<double>> cMatrix(nNodes, std::vector<double>(nNodes, 0.0));
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // C_ij = ∫ ρ * c_p * N_i * N_j dΩ
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    cMatrix[i][j] += material.density * material.specificHeat * 
                                   shapeResult.shapeFunctions[i] * 
                                   shapeResult.shapeFunctions[j] * weight;
                }
            }
        }
    }
    
    return cMatrix;
}

std::vector<double> HeatTransferSolver::calculateElementHeatSource(const Element& element) {
    auto material = materialDB.getMaterial(element.getMaterialName());
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    
    std::vector<double> heatSource(nNodes, 0.0);
    
    // Use Gaussian integration
    GaussIntegration gauss;
    auto integrationPoints = gauss.getIntegrationPoints(element.getType());
    
    ShapeFunctions shapeFunc;
    
    for (const auto& point : integrationPoints) {
        auto shapeResult = shapeFunc.computeShapeFunctions(element.getType(), point.coords);
        
        if (!shapeResult.jacobianDeterminant.empty()) {
            double detJ = shapeResult.jacobianDeterminant[0];
            double weight = point.weight * detJ;
            
            // Q_i = ∫ q * N_i dΩ
            for (int i = 0; i < nNodes; ++i) {
                heatSource[i] += material.heatSource * shapeResult.shapeFunctions[i] * weight;
            }
        }
    }
    
    return heatSource;
}

std::vector<double> HeatTransferSolver::calculateElementBoundaryFlux(const Element& element, 
                                                                    const ThermalBoundaryCondition& bc) {
    auto nodes = element.getNodes();
    auto nodeIndices = element.getNodeIndices();
    int nNodes = nodeIndices.size();
    
    std::vector<double> flux(nNodes, 0.0);
    
    // This would require boundary element integration
    // For now, return zero vector
    return flux;
}

bool HeatTransferSolver::isSystemReady() const {
    return mesh != nullptr && conductivityMatrix != nullptr && 
           heatSourceVector != nullptr && solver != nullptr;
}

void HeatTransferSolver::updateTimeStep(int step, double time) {
    // Update boundary conditions or material properties if they are time-dependent
    // For now, this is a placeholder for future time-dependent functionality
}

} // namespace elmer