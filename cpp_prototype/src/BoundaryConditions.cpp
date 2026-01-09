#include "BoundaryConditions.h"
#include "LinearAlgebra.h"
#include <algorithm>
#include <stdexcept>

namespace elmer {

// DirichletBoundaryCondition implementation
void DirichletBoundaryCondition::apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                                      std::vector<double>& rhs,
                                      const std::vector<int>& dofMap) const {
    if (nodeIndices_.size() != values_.size()) {
        throw std::invalid_argument("Number of node indices and values must match for Dirichlet BC");
    }
    
    for (size_t i = 0; i < nodeIndices_.size(); ++i) {
        int nodeIndex = nodeIndices_[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue; // Skip invalid node indices
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0) { // Only apply to active DOFs
            BoundaryConditionUtils::applyDirichletBC(dofIndex, values_[i], matrix, rhs);
        }
    }
}

bool DirichletBoundaryCondition::isActiveForElement(const Element& element) const {
    // Check if any of the element's nodes are in the boundary condition
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double DirichletBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    // For Dirichlet BC, we can return the value if the parameter name matches
    if (name == "Value" && !values_.empty()) {
        return values_[0]; // Return first value as default
    }
    return defaultValue;
}

bool DirichletBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    // Dirichlet BC doesn't typically have boolean parameters
    return defaultValue;
}

std::string DirichletBoundaryCondition::getStringParameter(const std::string& name, 
                                                          const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

// NeumannBoundaryCondition implementation
void NeumannBoundaryCondition::apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                                     std::vector<double>& rhs,
                                     const std::vector<int>& dofMap) const {
    // Neumann boundary conditions are natural boundary conditions in FEM
    // They are typically handled during element assembly by adding flux contributions
    // to the right-hand side vector
    
    if (nodeIndices_.size() != fluxValues_.size()) {
        throw std::invalid_argument("Number of node indices and flux values must match for Neumann BC");
    }
    
    // For standard formulations, Neumann BCs don't require modification
    // of the system matrices. The flux contributions are added during
    // element assembly in the solver-specific code.
}

bool NeumannBoundaryCondition::isActiveForElement(const Element& element) const {
    // Check if any of the element's nodes are in the boundary condition
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double NeumannBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    if (name == "Flux" && !fluxValues_.empty()) {
        return fluxValues_[0]; // Return first flux value as default
    }
    return defaultValue;
}

bool NeumannBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    return defaultValue;
}

std::string NeumannBoundaryCondition::getStringParameter(const std::string& name, 
                                                        const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

// RobinBoundaryCondition implementation
void RobinBoundaryCondition::apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                                    std::vector<double>& rhs,
                                    const std::vector<int>& dofMap) const {
    if (nodeIndices_.size() != alphaValues_.size() || 
        nodeIndices_.size() != betaValues_.size() ||
        nodeIndices_.size() != gammaValues_.size()) {
        throw std::invalid_argument("Number of node indices and parameter values must match for Robin BC");
    }
    
    // Robin boundary condition: alpha*u + beta*du/dn = gamma
    // This requires modification of both the stiffness matrix and right-hand side
    
    for (size_t i = 0; i < nodeIndices_.size(); ++i) {
        int nodeIndex = nodeIndices_[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue; // Skip invalid node indices
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0) { // Only apply to active DOFs
            double alpha = alphaValues_[i];
            double beta = betaValues_[i];
            double gamma = gammaValues_[i];
            
            // For Robin BC, we need to modify both the matrix and RHS
            // This is typically done during element assembly
            // For now, we'll implement a simplified version
            
            if (beta != 0.0) {
                // Add contribution to matrix diagonal
                double currentValue = matrix->get(dofIndex, dofIndex);
                matrix->set(dofIndex, dofIndex, currentValue + alpha / beta);
                
                // Add contribution to RHS
                rhs[dofIndex] += gamma / beta;
            }
        }
    }
}

bool RobinBoundaryCondition::isActiveForElement(const Element& element) const {
    // Check if any of the element's nodes are in the boundary condition
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double RobinBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    if (name == "Alpha" && !alphaValues_.empty()) {
        return alphaValues_[0];
    }
    if (name == "Beta" && !betaValues_.empty()) {
        return betaValues_[0];
    }
    if (name == "Gamma" && !gammaValues_.empty()) {
        return gammaValues_[0];
    }
    return defaultValue;
}

bool RobinBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    return defaultValue;
}

std::string RobinBoundaryCondition::getStringParameter(const std::string& name, 
                                                      const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

// BoundaryConditionUtils implementation
void BoundaryConditionUtils::applyDirichletBC(int dofIndex, double value, 
                                             std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                                             std::vector<double>& rhs) {
    if (dofIndex < 0 || dofIndex >= static_cast<int>(rhs.size())) {
        return; // Invalid DOF index
    }
    
    // Zero out the row in the matrix
    matrix->zeroRow(dofIndex);
    
    // Set the diagonal element to 1
    matrix->set(dofIndex, dofIndex, 1.0);
    
    // Set the corresponding RHS value
    rhs[dofIndex] = value;
}

void BoundaryConditionUtils::applyDirichletBCs(const std::vector<int>& dofIndices,
                                              const std::vector<double>& values,
                                              std::shared_ptr<ElmerCpp::IMatrix> matrix,
                                              std::vector<double>& rhs) {
    if (dofIndices.size() != values.size()) {
        throw std::invalid_argument("Number of DOF indices and values must match");
    }
    
    for (size_t i = 0; i < dofIndices.size(); ++i) {
        applyDirichletBC(dofIndices[i], values[i], matrix, rhs);
    }
}

bool BoundaryConditionUtils::hasBoundaryCondition(const Element& element, 
                                                 const std::vector<std::shared_ptr<BoundaryCondition>>& bcs) {
    for (const auto& bc : bcs) {
        if (bc->isActiveForElement(element)) {
            return true;
        }
    }
    return false;
}

std::vector<std::shared_ptr<BoundaryCondition>> BoundaryConditionUtils::getBoundaryConditionsForElement(
    const Element& element, 
    const std::vector<std::shared_ptr<BoundaryCondition>>& bcs) {
    
    std::vector<std::shared_ptr<BoundaryCondition>> result;
    
    for (const auto& bc : bcs) {
        if (bc->isActiveForElement(element)) {
            result.push_back(bc);
        }
    }
    
    return result;
}

std::shared_ptr<BoundaryCondition> BoundaryConditionUtils::createBoundaryCondition(
    BoundaryConditionType type, 
    const std::string& name,
    const std::vector<int>& nodeIndices,
    const std::vector<double>& values) {
    
    switch (type) {
        case BoundaryConditionType::DIRICHLET: {
            auto bc = std::make_shared<DirichletBoundaryCondition>(name);
            bc->setValues(nodeIndices, values);
            return bc;
        }
        
        case BoundaryConditionType::NEUMANN: {
            auto bc = std::make_shared<NeumannBoundaryCondition>(name);
            bc->setFluxValues(nodeIndices, values);
            return bc;
        }
        
        case BoundaryConditionType::ROBIN: {
            // For Robin BC, we need three sets of parameters
            // For simplicity, we'll use the values for alpha and set beta=1, gamma=0
            auto bc = std::make_shared<RobinBoundaryCondition>(name);
            std::vector<double> betaValues(nodeIndices.size(), 1.0);
            std::vector<double> gammaValues(nodeIndices.size(), 0.0);
            bc->setParameters(nodeIndices, values, betaValues, gammaValues);
            return bc;
        }
        
        default:
            throw std::invalid_argument("Unsupported boundary condition type");
    }
}

} // namespace elmer