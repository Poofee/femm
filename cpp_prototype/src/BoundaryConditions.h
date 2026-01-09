#pragma once

#include "LinearAlgebra.h"
#include "Mesh.h"
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>

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
    virtual void apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                       std::vector<double>& rhs,
                       const std::vector<int>& dofMap) const = 0;
    
    /**
     * @brief Check if boundary condition is active for given element
     */
    virtual bool isActiveForElement(const Element& element) const = 0;
    
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
    
    void apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
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
    
    void apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
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
    
    void apply(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
               std::vector<double>& rhs,
               const std::vector<int>& dofMap) const override;
    
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
     * @brief Apply all boundary conditions to system
     */
    void applyBoundaryConditions(std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                                std::vector<double>& rhs,
                                const std::vector<int>& dofMap) const {
        for (const auto& bc : boundaryConditions_) {
            bc->apply(matrix, rhs, dofMap);
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
namespace BoundaryConditionUtils {
    
    /**
     * @brief Apply Dirichlet boundary condition to system matrix and RHS
     */
    void applyDirichletBC(int dofIndex, double value, 
                         std::shared_ptr<ElmerCpp::IMatrix> matrix, 
                         std::vector<double>& rhs);
    
    /**
     * @brief Apply multiple Dirichlet boundary conditions
     */
    void applyDirichletBCs(const std::vector<int>& dofIndices,
                          const std::vector<double>& values,
                          std::shared_ptr<ElmerCpp::IMatrix> matrix,
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
}

} // namespace elmer