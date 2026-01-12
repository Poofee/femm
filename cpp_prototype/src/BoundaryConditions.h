#pragma once

#include "LinearAlgebra.h"
#include "Mesh.h"
#include "ElmerCpp.h"
#include <algorithm>
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace elmer;

namespace elmer {

/**
 * @brief Base boundary condition type
 */
enum class BoundaryConditionType {
    DIRICHLET,          ///< Fixed value boundary condition
    NEUMANN,            ///< Fixed flux/gradient boundary condition  
    ROBIN,              ///< Mixed boundary condition
    PERIODIC,           ///< Periodic boundary condition
    SYMMETRY,           ///< Symmetry boundary condition
    ANTISYMMETRY,       ///< Antisymmetry boundary condition
    
    // Electromagnetic specific boundary conditions
    MAGNETIC_SYMMETRY,   ///< Magnetic symmetry (B normal = 0)
    MAGNETIC_ANTISYMMETRY, ///< Magnetic antisymmetry (B tangential = 0)
    ELECTRIC_INSULATION, ///< Electric insulation (J normal = 0)
    INFINITY_BC,        ///< Infinity boundary condition
    
    // Harmonic analysis specific boundary conditions
    HARMONIC_EXCITATION, ///< Harmonic excitation boundary condition
    HARMONIC_CURRENT,    ///< Harmonic current source
    HARMONIC_VOLTAGE,    ///< Harmonic voltage source
    HARMONIC_FLUX,       ///< Harmonic magnetic flux source
    
    // Thermal specific boundary conditions
    CONVECTION,         ///< Convective boundary condition
    RADIATION,          ///< Radiative boundary condition
    
    // Structural specific boundary conditions
    DISPLACEMENT,       ///< Fixed displacement
    TRACTION,           ///< Surface traction
    PRESSURE,           ///< Pressure boundary condition
    SPRING              ///< Spring boundary condition
};

/**
 * @brief Base class for all boundary conditions
 */
class BoundaryCondition {
public:
    BoundaryCondition(BoundaryConditionType type, const std::string& name = "")
        : type_(type), name_(name) {}
    
    virtual ~BoundaryCondition() = default;
    
    /**
     * @brief Get boundary condition type
     */
    BoundaryConditionType type() const { return type_; }
    
    /**
     * @brief Get boundary condition name
     */
    const std::string& name() const { return name_; }
    
    /**
     * @brief Get nodes affected by this boundary condition
     */
    virtual const std::vector<int>& getNodeIndices() const = 0;
    
    /**
     * @brief Get boundary condition values
     */
    virtual const std::vector<double>& getValues() const = 0;
    
    /**
     * @brief Apply boundary condition to system matrix and right-hand side
     */
    virtual void apply(std::shared_ptr<Matrix> matrix, 
                      std::vector<double>& rhs,
                      const std::vector<int>& dofMap) const = 0;
    
    /**
     * @brief Apply boundary condition with priority handling
     */
    virtual void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                                  std::vector<double>& rhs,
                                  const std::vector<int>& dofMap,
                                  std::vector<bool>& constrainedDOFs) const = 0;
    
    /**
     * @brief Check if boundary condition is active for given element
     */
    virtual bool isActiveForElement(const Element& element) const = 0;
    
    /**
     * @brief Check if boundary condition is relevant for specific physics type
     */
    virtual bool isRelevantForPhysics(const std::string& physicsType) const {
        // Default implementation: check if physics type is in the name or parameters
        std::string physicsLower = physicsType;
        std::transform(physicsLower.begin(), physicsLower.end(), physicsLower.begin(), ::tolower);
        
        std::string nameLower = name_;
        std::transform(nameLower.begin(), nameLower.end(), nameLower.begin(), ::tolower);
        
        return nameLower.find(physicsLower) != std::string::npos;
    }
    
    /**
     * @brief Check if boundary condition couples multiple physics types
     */
    virtual bool isCoupledForPhysics(const std::vector<std::string>& physicsTypes) const {
        // Default implementation: check if boundary condition is relevant for all physics types
        for (const auto& physicsType : physicsTypes) {
            if (!isRelevantForPhysics(physicsType)) {
                return false;
            }
        }
        return !physicsTypes.empty();
    }
    
    /**
     * @brief Get parameter value by name
     */
    virtual double getParameter(const std::string& name, double defaultValue = 0.0) const = 0;
    
    /**
     * @brief Get boolean parameter value by name
     */
    virtual bool getBooleanParameter(const std::string& name, bool defaultValue = false) const = 0;
    
    /**
     * @brief Get string parameter value by name
     */
    virtual std::string getStringParameter(const std::string& name, 
                                          const std::string& defaultValue = "") const = 0;

protected:
    BoundaryConditionType type_;
    std::string name_;
};

/**
 * @brief Dirichlet boundary condition implementation
 */
class DirichletBoundaryCondition : public BoundaryCondition {
public:
    DirichletBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::DIRICHLET, name) {}
    
    /**
     * @brief Set boundary condition values for specific nodes
     */
    void setValues(const std::vector<int>& nodeIndices, const std::vector<double>& values) {
        nodeIndices_ = nodeIndices;
        values_ = values;
    }
    
    /**
     * @brief Set uniform value for all nodes
     */
    void setUniformValue(const std::vector<int>& nodeIndices, double value) {
        nodeIndices_ = nodeIndices;
        values_.assign(nodeIndices.size(), value);
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return values_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> values_;
};

/**
 * @brief Neumann boundary condition implementation
 */
class NeumannBoundaryCondition : public BoundaryCondition {
public:
    NeumannBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::NEUMANN, name) {}
    
    /**
     * @brief Set flux values for boundary nodes
     */
    void setFluxValues(const std::vector<int>& nodeIndices, const std::vector<double>& fluxValues) {
        nodeIndices_ = nodeIndices;
        fluxValues_ = fluxValues;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return fluxValues_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> fluxValues_; ///< Flux values for Neumann BC
};

/**
 * @brief Robin boundary condition implementation
 */
class RobinBoundaryCondition : public BoundaryCondition {
public:
    RobinBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::ROBIN, name) {}
    
    /**
     * @brief Set Robin boundary condition parameters
     */
    void setParameters(const std::vector<int>& nodeIndices, 
                      const std::vector<double>& alphaValues,
                      const std::vector<double>& betaValues,
                      const std::vector<double>& gammaValues) {
        nodeIndices_ = nodeIndices;
        alphaValues_ = alphaValues;
        betaValues_ = betaValues;
        gammaValues_ = gammaValues;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return alphaValues_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> alphaValues_; ///< Coefficient for field value
    std::vector<double> betaValues_;  ///< Coefficient for flux
    std::vector<double> gammaValues_; ///< Constant term
};

/**
 * @brief Magnetic symmetry boundary condition for electromagnetic problems
 * B_normal = 0 (magnetic field normal component is zero)
 */
class MagneticSymmetryBoundaryCondition : public BoundaryCondition {
public:
    MagneticSymmetryBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::MAGNETIC_SYMMETRY, name) {}
    
    /**
     * @brief Set boundary nodes for magnetic symmetry
     */
    void setBoundaryNodes(const std::vector<int>& nodeIndices) {
        nodeIndices_ = nodeIndices;
    }
    
    /**
     * @brief Set symmetry nodes (alias for setBoundaryNodes)
     */
    void setSymmetryNodes(const std::vector<int>& nodeIndices) {
        nodeIndices_ = nodeIndices;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return dummyValues_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> dummyValues_;
};

/**
 * @brief Magnetic antisymmetry boundary condition for electromagnetic problems
 * B_tangential = 0 (magnetic field tangential component is zero)
 */
class MagneticAntisymmetryBoundaryCondition : public BoundaryCondition {
public:
    MagneticAntisymmetryBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::MAGNETIC_ANTISYMMETRY, name) {}
    
    /**
     * @brief Set boundary nodes for magnetic antisymmetry
     */
    void setBoundaryNodes(const std::vector<int>& nodeIndices) {
        nodeIndices_ = nodeIndices;
    }
    
    /**
     * @brief Set antisymmetry nodes (alias for setBoundaryNodes)
     */
    void setAntisymmetryNodes(const std::vector<int>& nodeIndices) {
        nodeIndices_ = nodeIndices;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return dummyValues_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> dummyValues_;
};

/**
 * @brief Convection boundary condition for thermal problems
 */
class ConvectionBoundaryCondition : public BoundaryCondition {
public:
    ConvectionBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::CONVECTION, name) {}
    
    /**
     * @brief Set convection parameters
     */
    void setConvectionParameters(const std::vector<int>& nodeIndices,
                                const std::vector<double>& heatTransferCoeffs,
                                const std::vector<double>& ambientTemps) {
        nodeIndices_ = nodeIndices;
        hValues_ = heatTransferCoeffs;
        TinfValues_ = ambientTemps;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return hValues_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> hValues_;      ///< Heat transfer coefficients
    std::vector<double> TinfValues_;   ///< Ambient temperatures
};

/**
 * @brief Pressure boundary condition for structural problems
 */
class PressureBoundaryCondition : public BoundaryCondition {
public:
    PressureBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::PRESSURE, name) {}
    
    /**
     * @brief Set pressure values
     */
    void setPressureValues(const std::vector<int>& nodeIndices, const std::vector<double>& pressures) {
        nodeIndices_ = nodeIndices;
        pressureValues_ = pressures;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { return pressureValues_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override;
    
    bool isActiveForElement(const Element& element) const override;
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override;
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override;
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override;

private:
    std::vector<int> nodeIndices_;
    std::vector<double> pressureValues_;
};

/**
 * @brief Boundary condition manager for handling multiple boundary conditions
 */
class BoundaryConditionManager {
public:
    BoundaryConditionManager() = default;
    
    /**
     * @brief Add a boundary condition
     */
    void addBoundaryCondition(std::shared_ptr<BoundaryCondition> bc) {
        boundaryConditions_.push_back(bc);
        bcMap_[bc->name()] = bc;
    }
    
    /**
     * @brief Get boundary condition by name
     */
    std::shared_ptr<BoundaryCondition> getBoundaryCondition(const std::string& name) const {
        auto it = bcMap_.find(name);
        if (it != bcMap_.end()) {
            return it->second;
        }
        return nullptr;
    }
    
    /**
     * @brief Get all boundary conditions
     */
    const std::vector<std::shared_ptr<BoundaryCondition>>& getAllBoundaryConditions() const {
        return boundaryConditions_;
    }
    
    /**
     * @brief Apply all boundary conditions to system with priority handling
     */
    void applyBoundaryConditions(std::shared_ptr<Matrix> matrix, 
                                std::vector<double>& rhs,
                                const std::vector<int>& dofMap) const {
        // Track which DOFs have been constrained
        std::vector<bool> constrainedDOFs(dofMap.size(), false);
        
        // Apply with priority
        applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
    }
    
    /**
     * @brief Apply boundary conditions with explicit DOF constraint tracking
     */
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const {
        // Apply boundary conditions in priority order
        // Higher priority conditions override lower priority ones
        std::vector<std::shared_ptr<BoundaryCondition>> sortedBCs = boundaryConditions_;
        
        // Sort by priority: Dirichlet > Robin > Neumann > Others
        std::sort(sortedBCs.begin(), sortedBCs.end(), 
                  [](const std::shared_ptr<BoundaryCondition>& a, 
                     const std::shared_ptr<BoundaryCondition>& b) {
                      
            // Priority order: Dirichlet (highest) > Robin > Neumann > Others
            auto getPriority = [](BoundaryConditionType type) {
                switch (type) {
                    case BoundaryConditionType::DIRICHLET: return 0;
                    case BoundaryConditionType::ROBIN: return 1;
                    case BoundaryConditionType::NEUMANN: return 2;
                    default: return 3;
                }
            };
            
            return getPriority(a->type()) < getPriority(b->type());
        });
        
        for (const auto& bc : sortedBCs) {
            bc->applyWithPriority(matrix, rhs, dofMap, constrainedDOFs);
        }
    }
    
    /**
     * @brief Apply boundary conditions for specific physics type
     */
    void applyPhysicsBoundaryConditions(const std::string& physicsType,
                                       std::shared_ptr<Matrix> matrix, 
                                       std::vector<double>& rhs,
                                       const std::vector<int>& dofMap) const {
        // Apply only boundary conditions relevant to specific physics
        for (const auto& bc : boundaryConditions_) {
            if (bc->isRelevantForPhysics(physicsType)) {
                bc->apply(matrix, rhs, dofMap);
            }
        }
    }
    
    /**
     * @brief Apply coupled boundary conditions for multi-physics problems
     */
    void applyCoupledBoundaryConditions(std::shared_ptr<Matrix> matrix, 
                                       std::vector<double>& rhs,
                                       const std::vector<int>& dofMap,
                                       const std::vector<std::string>& physicsTypes) const {
        // Apply boundary conditions that couple multiple physics
        for (const auto& bc : boundaryConditions_) {
            if (bc->isCoupledForPhysics(physicsTypes)) {
                bc->apply(matrix, rhs, dofMap);
            }
        }
    }
    
    /**
     * @brief Clear all boundary conditions
     */
    void clear() {
        boundaryConditions_.clear();
        bcMap_.clear();
    }
    
    /**
     * @brief Get boundary conditions of specific type
     */
    std::vector<std::shared_ptr<BoundaryCondition>> getBoundaryConditionsByType(BoundaryConditionType type) const {
        std::vector<std::shared_ptr<BoundaryCondition>> result;
        for (const auto& bc : boundaryConditions_) {
            if (bc->type() == type) {
                result.push_back(bc);
            }
        }
        return result;
    }

private:
    std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions_;
    std::unordered_map<std::string, std::shared_ptr<BoundaryCondition>> bcMap_;
};

/**
 * @brief Utility functions for boundary condition handling
 */

/**
 * @brief 谐波激励边界条件
 * 
 * 用于谐波分析的激励源边界条件，支持复数激励
 */
class HarmonicExcitationBoundaryCondition : public BoundaryCondition {
public:
    HarmonicExcitationBoundaryCondition(const std::string& name = "")
        : BoundaryCondition(BoundaryConditionType::HARMONIC_EXCITATION, name) {}
    
    /**
     * @brief 设置谐波激励参数
     */
    void setHarmonicParameters(const std::vector<int>& nodeIndices,
                              const std::vector<double>& amplitudes,
                              const std::vector<double>& phases,
                              double frequency) {
        nodeIndices_ = nodeIndices;
        amplitudes_ = amplitudes;
        phases_ = phases;
        frequency_ = frequency;
    }
    
    /**
     * @brief 设置复数激励值
     */
    void setComplexExcitation(const std::vector<int>& nodeIndices,
                             const std::vector<std::complex<double>>& complexValues) {
        nodeIndices_ = nodeIndices;
        complexValues_ = complexValues;
        useComplexValues_ = true;
    }
    
    const std::vector<int>& getNodeIndices() const override { return nodeIndices_; }
    const std::vector<double>& getValues() const override { 
        // 返回幅值作为实数值
        return amplitudes_; 
    }
    
    /**
     * @brief 获取复数激励值
     */
    const std::vector<std::complex<double>>& getComplexValues() const {
        if (useComplexValues_) {
            return complexValues_;
        } else {
            // 从幅值和相位计算复数值
            static std::vector<std::complex<double>> computedValues;
            computedValues.resize(amplitudes_.size());
            for (size_t i = 0; i < amplitudes_.size(); ++i) {
                computedValues[i] = std::polar(amplitudes_[i], phases_[i]);
            }
            return computedValues;
        }
    }
    
    /**
     * @brief 获取频率
     */
    double getFrequency() const { return frequency_; }
    
    void apply(std::shared_ptr<Matrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override {
        // 谐波边界条件通常应用于复数系统
        // 对于实数系统，使用幅值作为激励
        applyRealExcitation(rhs, dofMap);
    }
    
    void applyWithPriority(std::shared_ptr<Matrix> matrix, 
                          std::vector<double>& rhs,
                          const std::vector<int>& dofMap,
                          std::vector<bool>& constrainedDOFs) const override {
        applyRealExcitation(rhs, dofMap);
    }
    
    /**
     * @brief 应用复数激励到复数系统
     */
    void applyComplexExcitation(std::shared_ptr<Matrix> matrix,
                               std::vector<std::complex<double>>& complexRhs,
                               const std::vector<int>& dofMap) const {
        auto complexValues = getComplexValues();
        
        for (size_t i = 0; i < nodeIndices_.size(); ++i) {
            int nodeIdx = nodeIndices_[i];
            if (nodeIdx < dofMap.size()) {
                int dofIdx = dofMap[nodeIdx];
                if (dofIdx >= 0 && dofIdx < complexRhs.size()) {
                    complexRhs[dofIdx] += complexValues[i];
                }
            }
        }
    }
    
    bool isActiveForElement(const Element& element) const override {
        // 检查元素是否包含边界节点
        auto elemNodeIndices = element.getNodeIndices();
        for (int nodeIdx : nodeIndices_) {
            if (std::find(elemNodeIndices.begin(), elemNodeIndices.end(), nodeIdx) != elemNodeIndices.end()) {
                return true;
            }
        }
        return false;
    }
    
    bool isRelevantForPhysics(const std::string& physicsType) const override {
        // 谐波激励适用于电磁问题
        return physicsType == "MagnetoDynamics" || 
               physicsType == "ElectroMagnetics" ||
               physicsType == "HarmonicAnalysis";
    }
    
    double getParameter(const std::string& name, double defaultValue = 0.0) const override {
        if (name == "Frequency") return frequency_;
        if (name == "Amplitude" && !amplitudes_.empty()) return amplitudes_[0];
        if (name == "Phase" && !phases_.empty()) return phases_[0];
        return defaultValue;
    }
    
    bool getBooleanParameter(const std::string& name, bool defaultValue = false) const override {
        if (name == "UseComplex") return useComplexValues_;
        return defaultValue;
    }
    
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override {
        if (name == "ExcitationType") return "Harmonic";
        return defaultValue;
    }

private:
    std::vector<int> nodeIndices_;
    std::vector<double> amplitudes_;  ///< 激励幅值
    std::vector<double> phases_;      ///< 激励相位 [rad]
    double frequency_ = 0.0;          ///< 激励频率 [Hz]
    
    std::vector<std::complex<double>> complexValues_; ///< 复数激励值
    bool useComplexValues_ = false;   ///< 是否使用复数激励值
    
    /**
     * @brief 应用实数激励到右端向量
     */
    void applyRealExcitation(std::vector<double>& rhs, const std::vector<int>& dofMap) const {
        for (size_t i = 0; i < nodeIndices_.size(); ++i) {
            int nodeIdx = nodeIndices_[i];
            if (nodeIdx < dofMap.size()) {
                int dofIdx = dofMap[nodeIdx];
                if (dofIdx >= 0 && dofIdx < rhs.size()) {
                    rhs[dofIdx] += amplitudes_[i];
                }
            }
        }
    }
};

/**
 * @brief 谐波电流源边界条件
 */
class HarmonicCurrentBoundaryCondition : public HarmonicExcitationBoundaryCondition {
public:
    HarmonicCurrentBoundaryCondition(const std::string& name = "")
        : HarmonicExcitationBoundaryCondition(name) {
        // 设置边界条件类型
        type_ = BoundaryConditionType::HARMONIC_CURRENT;
    }
    
    bool isRelevantForPhysics(const std::string& physicsType) const override {
        return physicsType == "MagnetoDynamics" || 
               physicsType == "ElectroMagnetics" ||
               physicsType == "HarmonicAnalysis";
    }
    
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override {
        if (name == "ExcitationType") return "Current";
        return defaultValue;
    }
};

/**
 * @brief 谐波电压源边界条件
 */
class HarmonicVoltageBoundaryCondition : public HarmonicExcitationBoundaryCondition {
public:
    HarmonicVoltageBoundaryCondition(const std::string& name = "")
        : HarmonicExcitationBoundaryCondition(name) {
        // 设置边界条件类型
        type_ = BoundaryConditionType::HARMONIC_VOLTAGE;
    }
    
    bool isRelevantForPhysics(const std::string& physicsType) const override {
        return physicsType == "ElectroMagnetics" || 
               physicsType == "HarmonicAnalysis";
    }
    
    std::string getStringParameter(const std::string& name, 
                                  const std::string& defaultValue = "") const override {
        if (name == "ExcitationType") return "Voltage";
        return defaultValue;
    }
};

} // namespace elmer

namespace elmerBoundaryConditionUtils {
    
    /**
     * @brief Apply Dirichlet boundary condition to system matrix and RHS
     */
    void applyDirichletBC(int dofIndex, double value, 
                         std::shared_ptr<Matrix> matrix, 
                         std::vector<double>& rhs);
    
    /**
     * @brief Apply multiple Dirichlet boundary conditions
     */
    void applyDirichletBCs(const std::vector<int>& dofIndices,
                          const std::vector<double>& values,
                          std::shared_ptr<Matrix> matrix,
                          std::vector<double>& rhs);
    
    /**
     * @brief Check if element has boundary condition
     */
    bool hasBoundaryCondition(const Element& element, 
                             const std::vector<std::shared_ptr<BoundaryCondition>>& bcs);
    
    /**
     * @brief Get boundary conditions for specific element
     */
    std::vector<std::shared_ptr<BoundaryCondition>> getBoundaryConditionsForElement(
        const Element& element, 
        const std::vector<std::shared_ptr<BoundaryCondition>>& bcs);
    
    /**
     * @brief Create boundary condition from parameters
     */
    std::shared_ptr<BoundaryCondition> createBoundaryCondition(
        BoundaryConditionType type, 
        const std::string& name,
        const std::vector<int>& nodeIndices,
        const std::vector<double>& values);
} // namespace elmerBoundaryConditionUtils