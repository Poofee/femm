#include "BoundaryConditions.h"
#include "LinearAlgebra.h"
#include <algorithm>
#include <stdexcept>

using namespace elmer;

namespace elmer {

// DirichletBoundaryCondition implementation
void DirichletBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
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
            elmerBoundaryConditionUtils::applyDirichletBC(dofIndex, values_[i], matrix, rhs);
        }
    }
}

void PressureBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                                   std::vector<double>& rhs,
                                                   const std::vector<int>& dofMap,
                                                   std::vector<bool>& constrainedDOFs) const {
    // Pressure boundary condition is a Neumann-like condition
    // It can be applied even if DOFs are already constrained by Dirichlet BCs
    apply(matrix, rhs, dofMap);
}
void DirichletBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                                  std::vector<double>& rhs,
                                                  const std::vector<int>& dofMap,
                                                  std::vector<bool>& constrainedDOFs) const {
    if (nodeIndices_.size() != values_.size()) {
        throw std::invalid_argument("Number of node indices and values must match for Dirichlet BC");
    }
    
    for (size_t i = 0; i < nodeIndices_.size(); ++i) {
        int nodeIndex = nodeIndices_[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue; // Skip invalid node indices
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0 && !constrainedDOFs[dofIndex]) { // Only apply to unconstrained DOFs
            elmerBoundaryConditionUtils::applyDirichletBC(dofIndex, values_[i], matrix, rhs);
            constrainedDOFs[dofIndex] = true; // Mark this DOF as constrained
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
void NeumannBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
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

void NeumannBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                                 std::vector<double>& rhs,
                                                 const std::vector<int>& dofMap,
                                                 std::vector<bool>& constrainedDOFs) const {
    // Neumann boundary conditions are natural conditions and don't constrain DOFs
    // They can be applied even if DOFs are already constrained by Dirichlet BCs
    apply(matrix, rhs, dofMap);
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
void RobinBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
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
                double currentValue = matrix->GetElement(dofIndex, dofIndex);
                matrix->SetElement(dofIndex, dofIndex, currentValue + alpha / beta);
                
                // Add contribution to RHS
                rhs[dofIndex] += gamma / beta;
            }
        }
    }
}

void RobinBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                               std::vector<double>& rhs,
                                               const std::vector<int>& dofMap,
                                               std::vector<bool>& constrainedDOFs) const {
    // Robin boundary conditions can be applied even if DOFs are already constrained
    // They modify both matrix and RHS, so they should be applied after Dirichlet BCs
    // but before Neumann BCs in the priority order
    apply(matrix, rhs, dofMap);
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



// MagneticSymmetryBoundaryCondition implementation
void MagneticSymmetryBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
                                              std::vector<double>& rhs,
                                              const std::vector<int>& dofMap) const {
    // Magnetic symmetry: B_normal = 0
    // This is typically implemented as a Dirichlet condition for the normal component
    // of the magnetic vector potential or magnetic field
    for (int nodeIndex : nodeIndices_) {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue;
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0) {
            // Set the normal component to zero
            // This is a simplified implementation; actual implementation depends on formulation
            elmerBoundaryConditionUtils::applyDirichletBC(dofIndex, 0.0, matrix, rhs);
        }
    }
}

void MagneticSymmetryBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                                          std::vector<double>& rhs,
                                                          const std::vector<int>& dofMap,
                                                          std::vector<bool>& constrainedDOFs) const {
    // Magnetic symmetry is a Dirichlet-like condition for the normal component
    // It should be applied with priority handling similar to Dirichlet BCs
    for (int nodeIndex : nodeIndices_) {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue;
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0 && !constrainedDOFs[dofIndex]) {
            elmerBoundaryConditionUtils::applyDirichletBC(dofIndex, 0.0, matrix, rhs);
            constrainedDOFs[dofIndex] = true;
        }
    }
}

bool MagneticSymmetryBoundaryCondition::isActiveForElement(const Element& element) const {
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double MagneticSymmetryBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    return defaultValue;
}

bool MagneticSymmetryBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    return defaultValue;
}

std::string MagneticSymmetryBoundaryCondition::getStringParameter(const std::string& name, 
                                                                  const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

// MagneticAntisymmetryBoundaryCondition implementation
void MagneticAntisymmetryBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
                                                  std::vector<double>& rhs,
                                                  const std::vector<int>& dofMap) const {
    // Magnetic antisymmetry: B_tangential = 0
    // This is typically implemented as a Dirichlet condition for the tangential component
    // of the magnetic vector potential or magnetic field
    for (int nodeIndex : nodeIndices_) {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue;
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0) {
            // Set the tangential component to zero
            // This is a simplified implementation; actual implementation depends on formulation
            elmerBoundaryConditionUtils::applyDirichletBC(dofIndex, 0.0, matrix, rhs);
        }
    }
}

void MagneticAntisymmetryBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                                              std::vector<double>& rhs,
                                                              const std::vector<int>& dofMap,
                                                              std::vector<bool>& constrainedDOFs) const {
    // Magnetic antisymmetry is a Dirichlet-like condition for the tangential component
    // It should be applied with priority handling similar to Dirichlet BCs
    for (int nodeIndex : nodeIndices_) {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue;
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0 && !constrainedDOFs[dofIndex]) {
            elmerBoundaryConditionUtils::applyDirichletBC(dofIndex, 0.0, matrix, rhs);
            constrainedDOFs[dofIndex] = true;
        }
    }
}

bool MagneticAntisymmetryBoundaryCondition::isActiveForElement(const Element& element) const {
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double MagneticAntisymmetryBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    return defaultValue;
}

bool MagneticAntisymmetryBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    return defaultValue;
}

std::string MagneticAntisymmetryBoundaryCondition::getStringParameter(const std::string& name, 
                                                                      const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

// ConvectionBoundaryCondition implementation
void ConvectionBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
                                        std::vector<double>& rhs,
                                        const std::vector<int>& dofMap) const {
    if (nodeIndices_.size() != hValues_.size() || 
        nodeIndices_.size() != TinfValues_.size()) {
        throw std::invalid_argument("Parameter arrays must have same size for Convection BC");
    }
    
    // Convection boundary condition: h(T - T_inf) = -k∂T/∂n
    // In FEM, this modifies both the matrix and RHS
    for (size_t i = 0; i < nodeIndices_.size(); ++i) {
        int nodeIndex = nodeIndices_[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue;
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0) {
            // Add convection contributions to matrix and RHS
            double currentValue = matrix->GetElement(dofIndex, dofIndex);
            matrix->SetElement(dofIndex, dofIndex, currentValue + hValues_[i]);
            rhs[dofIndex] += hValues_[i] * TinfValues_[i];
        }
    }
}

void ConvectionBoundaryCondition::applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                                    std::vector<double>& rhs,
                                                    const std::vector<int>& dofMap,
                                                    std::vector<bool>& constrainedDOFs) const {
    // Convection boundary condition is a Robin-like condition
    // It can be applied even if DOFs are already constrained by Dirichlet BCs
    apply(matrix, rhs, dofMap);
}

bool ConvectionBoundaryCondition::isActiveForElement(const Element& element) const {
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double ConvectionBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    if (name == "HeatTransferCoefficient" && !hValues_.empty()) {
        return hValues_[0];
    }
    if (name == "AmbientTemperature" && !TinfValues_.empty()) {
        return TinfValues_[0];
    }
    return defaultValue;
}

bool ConvectionBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    return defaultValue;
}

std::string ConvectionBoundaryCondition::getStringParameter(const std::string& name, 
                                                            const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

// PressureBoundaryCondition implementation
void PressureBoundaryCondition::apply(std::shared_ptr<Matrix> matrix, 
                                      std::vector<double>& rhs,
                                      const std::vector<int>& dofMap) const {
    // Pressure boundary condition for structural problems
    // This adds pressure forces to the right-hand side vector
    for (size_t i = 0; i < nodeIndices_.size(); ++i) {
        int nodeIndex = nodeIndices_[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(dofMap.size())) {
            continue;
        }
        
        int dofIndex = dofMap[nodeIndex];
        if (dofIndex >= 0) {
            // Add pressure force to RHS
            // Note: The actual implementation depends on the element formulation
            // and the direction of the pressure
            rhs[dofIndex] += pressureValues_[i];
        }
    }
}

bool PressureBoundaryCondition::isActiveForElement(const Element& element) const {
    for (int nodeIndex : nodeIndices_) {
        for (size_t elementNode : element.getNodeIndices()) {
            if (nodeIndex == static_cast<int>(elementNode)) {
                return true;
            }
        }
    }
    return false;
}

double PressureBoundaryCondition::getParameter(const std::string& name, double defaultValue) const {
    if (name == "Pressure" && !pressureValues_.empty()) {
        return pressureValues_[0];
    }
    return defaultValue;
}

bool PressureBoundaryCondition::getBooleanParameter(const std::string& name, bool defaultValue) const {
    return defaultValue;
}

std::string PressureBoundaryCondition::getStringParameter(const std::string& name, 
                                                          const std::string& defaultValue) const {
    if (name == "Name") {
        return name_;
    }
    return defaultValue;
}

} // namespace elmer

// elmerBoundaryConditionUtils implementation
void elmerBoundaryConditionUtils::applyDirichletBC(int dofIndex, double value, 
                                                   std::shared_ptr<Matrix> matrix, 
                                                   std::vector<double>& rhs) {
    if (dofIndex < 0 || dofIndex >= static_cast<int>(rhs.size())) {
        return; // Invalid DOF index
    }
    
    // Zero out the row in the matrix by setting all elements to 0
    for (int j = 0; j < matrix->GetNumCols(); ++j) {
        matrix->SetElement(dofIndex, j, 0.0);
    }
    
    // Set the diagonal element to 1
    matrix->SetElement(dofIndex, dofIndex, 1.0);
    
    // Set the corresponding RHS value
    rhs[dofIndex] = value;
}

void elmerBoundaryConditionUtils::applyDirichletBCs(const std::vector<int>& dofIndices,
                                                    const std::vector<double>& values,
                                                    std::shared_ptr<Matrix> matrix,
                                                    std::vector<double>& rhs) {
    if (dofIndices.size() != values.size()) {
        throw std::invalid_argument("Number of DOF indices and values must match");
    }
    
    for (size_t i = 0; i < dofIndices.size(); ++i) {
        applyDirichletBC(dofIndices[i], values[i], matrix, rhs);
    }
}

bool elmerBoundaryConditionUtils::hasBoundaryCondition(const Element& element, 
                                                       const std::vector<std::shared_ptr<BoundaryCondition>>& bcs) {
    for (const auto& bc : bcs) {
        if (bc->isActiveForElement(element)) {
            return true;
        }
    }
    return false;
}

std::vector<std::shared_ptr<BoundaryCondition>> elmerBoundaryConditionUtils::getBoundaryConditionsForElement(
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

std::shared_ptr<BoundaryCondition> elmerBoundaryConditionUtils::createBoundaryCondition(
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
        
        case BoundaryConditionType::MAGNETIC_SYMMETRY: {
            auto bc = std::make_shared<MagneticSymmetryBoundaryCondition>(name);
            bc->setBoundaryNodes(nodeIndices);
            return bc;
        }
        
        case BoundaryConditionType::MAGNETIC_ANTISYMMETRY: {
            auto bc = std::make_shared<MagneticAntisymmetryBoundaryCondition>(name);
            bc->setBoundaryNodes(nodeIndices);
            return bc;
        }
        
        case BoundaryConditionType::CONVECTION: {
            auto bc = std::make_shared<ConvectionBoundaryCondition>(name);
            // For convection BC, we need heat transfer coefficients and ambient temperatures
            // For simplicity, use values for h and set Tinf=0
            std::vector<double> ambientTemps(nodeIndices.size(), 0.0);
            bc->setConvectionParameters(nodeIndices, values, ambientTemps);
            return bc;
        }
        
        case BoundaryConditionType::PRESSURE: {
            auto bc = std::make_shared<PressureBoundaryCondition>(name);
            bc->setPressureValues(nodeIndices, values);
            return bc;
        }
        
        default:
            throw std::invalid_argument("Unsupported boundary condition type");
    }
}

void elmerBoundaryConditionUtils::GetElementMeshEdgeInfo(const Mesh& mesh, const Element& element,
                                                        std::vector<int>& edgeDegree,
                                                        std::vector<std::vector<int>>& edgeDirection,
                                                        int& edgeMaxDegree) {
    // 获取元素的边信息
    auto elementEdges = element.getEdges();
    edgeDegree.clear();
    edgeDirection.clear();
    edgeMaxDegree = 0;
    
    // 简化实现：假设所有边的度数为1
    for (const auto& edgeIndex : elementEdges) {
        // 简化实现：固定度数
        int degree = 1;
        edgeDegree.push_back(degree);
        
        // 更新最大度数
        if (degree > edgeMaxDegree) {
            edgeMaxDegree = degree;
        }
        
        // 设置边方向（简化实现）
        std::vector<int> direction = {1, 1}; // 默认方向
        edgeDirection.push_back(direction);
    }
    
    // 如果没有边，设置默认值
    if (edgeDegree.empty()) {
        edgeMaxDegree = 1;
    }
}

void elmerBoundaryConditionUtils::GetElementMeshFaceInfo(const Mesh& mesh, const Element& element,
                                                        std::vector<int>& faceDegree,
                                                        std::vector<std::vector<int>>& faceDirection,
                                                        int& faceMaxDegree) {
    // 获取元素的面信息
    auto elementFaces = element.getFaces();
    faceDegree.clear();
    faceDirection.clear();
    faceMaxDegree = 0;
    
    // 简化实现：假设所有面的度数为1
    for (const auto& faceIndex : elementFaces) {
        // 简化实现：固定度数
        int degree = 1;
        faceDegree.push_back(degree);
        
        // 更新最大度数
        if (degree > faceMaxDegree) {
            faceMaxDegree = degree;
        }
        
        // 设置面方向（简化实现）
        std::vector<int> direction = {1, 1, 1}; // 默认方向
        faceDirection.push_back(direction);
    }
    
    // 如果没有面，设置默认值
    if (faceDegree.empty()) {
        faceMaxDegree = 1;
    }
}

void elmerBoundaryConditionUtils::FaceElementBasisOrdering(const Element& element,
                                                          std::vector<std::vector<int>>& fDofMap,
                                                          int faceNumber,
                                                          std::vector<bool>* reverseSign) {
    // 获取元素的面信息
    auto elementFaces = element.getFaces();
    
    // 清空输出参数
    fDofMap.clear();
    if (reverseSign) {
        reverseSign->clear();
    }
    
    // 根据面编号处理
    if (faceNumber >= 0 && faceNumber < static_cast<int>(elementFaces.size())) {
        // 处理特定面
        auto faceIndex = elementFaces[faceNumber];
        
        // 简化的自由度映射
        // 简化实现：固定2个自由度
        std::vector<int> faceDofs(2);
        for (int i = 0; i < face.bDofs; ++i) {
            faceDofs[i] = i;
        }
        fDofMap.push_back(faceDofs);
        
        // 设置反向符号标志
        if (reverseSign) {
            reverseSign->push_back(false);
        }
    } else {
        // 处理所有面
        for (const auto& face : elementFaces) {
            std::vector<int> faceDofs(face.bDofs);
            for (int i = 0; i < face.bDofs; ++i) {
                faceDofs[i] = i;
            }
            fDofMap.push_back(faceDofs);
            
            if (reverseSign) {
                reverseSign->push_back(false);
            }
        }
    }
}

void elmerBoundaryConditionUtils::GetEdgeBasis(const Element& element,
                                              std::vector<double>& wBasis,
                                              std::vector<double>& rotWBasis,
                                              std::vector<double>& basis,
                                              std::vector<std::vector<double>>& dBasisdx) {
    // 获取元素的边信息
    auto elementEdges = element.getEdges();
    
    // 清空输出参数
    wBasis.clear();
    rotWBasis.clear();
    basis.clear();
    dBasisdx.clear();
    
    // 简化的边基函数计算
    // 简化实现：假设每个边有2个自由度
    for (const auto& edgeIndex : elementEdges) {
        // 边基函数值（简化实现）
        for (int i = 0; i < 2; ++i) { // 固定2个自由度
            wBasis.push_back(1.0); // 常数基函数
            rotWBasis.push_back(0.0); // 旋转基函数
            basis.push_back(1.0); // 基函数值
            
            // 基函数导数（简化实现）
            std::vector<double> deriv = {0.0, 0.0, 0.0};
            dBasisdx.push_back(deriv);
        }
    }
    
    // 如果没有边，设置默认值
    if (wBasis.empty()) {
        wBasis.push_back(1.0);
        rotWBasis.push_back(0.0);
        basis.push_back(1.0);
        std::vector<double> deriv = {0.0, 0.0, 0.0};
        dBasisdx.push_back(deriv);
    }
}