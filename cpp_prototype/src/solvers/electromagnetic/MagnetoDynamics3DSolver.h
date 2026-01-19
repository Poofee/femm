#pragma once

#include "MagnetoDynamicsSolverBase.h"
#include "ElectromagneticMaterial.h"
#include <memory>
#include <vector>
#include <array>
#include <complex>

namespace elmer {

/**
 * @brief 3D磁动力学求解器
 * 
 * 基于Elmer Fortran WhitneyAVSolver求解器的C++实现
 * 支持三维低频电磁场分析，使用Whitney边元方法
 */
class MagnetoDynamics3DSolver : public MagnetoDynamicsSolverBase {
public:
    MagnetoDynamics3DSolver();
    virtual ~MagnetoDynamics3DSolver() = default;
    
    // 重写基类方法
    virtual bool initialize() override;
    virtual bool assembleSystem() override;
    virtual bool solve() override;
    
    // 3D特定方法
    void setUseWhitneyElements(bool useWhitney);
    bool getUseWhitneyElements() const;
    
    void setUsePiolaTransformation(bool usePiola);
    bool getUsePiolaTransformation() const;
    
    void setUseSecondOrderElements(bool useSecondOrder);
    bool getUseSecondOrderElements() const;
    
    // 结果获取
    std::vector<std::array<double, 3>> getVectorPotential() const;
    std::vector<std::array<std::complex<double>, 3>> getComplexVectorPotential() const;
    
    // 3D特定场量计算
    std::vector<std::array<double, 3>> calculateMagneticFluxDensity3D() const;
    std::vector<std::array<double, 3>> calculateMagneticFieldStrength3D() const;
    std::vector<std::array<double, 3>> calculateCurrentDensity3D() const;
    
    // 复数场量计算（谐波分析）
    std::vector<std::array<std::complex<double>, 3>> calculateComplexMagneticFluxDensity3D() const;
    std::vector<std::array<std::complex<double>, 3>> calculateComplexMagneticFieldStrength3D() const;
    std::vector<std::array<std::complex<double>, 3>> calculateComplexCurrentDensity3D() const;
    
    // 集总参数计算
    double calculateTorque3D() const;
    double calculateMagneticEnergy3D() const;
    double calculateInductance3D() const;
    
    // 复数集总参数（谐波分析）
    std::complex<double> calculateComplexTorque3D() const;
    std::complex<double> calculateComplexMagneticEnergy3D() const;
    std::complex<double> calculateComplexInductance3D() const;
    
protected:
    // 重写保护方法
    virtual bool assembleElementMatrix(int elementId, ElementMatrix& elementMatrix) override;
    virtual bool assembleElementRHS(int elementId, std::vector<double>& elementRHS) override;
    virtual bool applyBoundaryConditions() override;
    virtual bool updateMaterialProperties() override;
    virtual bool calculateDerivedFields() override;
    
    // 维度特定方法
    virtual int getDegreesOfFreedom() const override;
    virtual std::string getVariableName() const override;
    
    // 3D特定工具方法
    bool assembleWhitneyElementMatrix(int elementId, ElementMatrix& elementMatrix);
    bool assembleLagrangeElementMatrix(int elementId, ElementMatrix& elementMatrix);
    
    double calculateElementVolume3D(int elementId) const;
    std::array<double, 3> calculateElementCentroid3D(int elementId) const;
    double calculateElementConductivity3D(int elementId) const;
    double calculateElementPermeability3D(int elementId) const;
    
    // 3D特定场量计算
    std::array<double, 3> calculateElementMagneticFluxDensity(int elementId) const;
    std::array<double, 3> calculateElementMagneticFieldStrength(int elementId) const;
    std::array<double, 3> calculateElementCurrentDensity(int elementId) const;
    std::array<double, 3> getElementCurrentDensity3D(int elementId) const;
    
    // 复数场量计算（谐波分析）
    std::array<std::complex<double>, 3> calculateElementComplexMagneticFluxDensity(int elementId) const;
    std::array<std::complex<double>, 3> calculateElementComplexMagneticFieldStrength(int elementId) const;
    std::array<std::complex<double>, 3> calculateElementComplexCurrentDensity(int elementId) const;
    
    // Whitney边元特定方法
    bool calculateWhitneyShapeFunctions(int elementId, std::vector<std::array<double, 3>>& shapeFuncs) const;
    bool calculateWhitneyCurlShapeFunctions(int elementId, std::vector<std::array<double, 3>>& curlShapeFuncs) const;
    std::vector<std::array<double, 3>> calculateWhitneyShapeFunctions(int elementId) const;
    
    // Lagrange单元特定方法
    std::vector<double> calculateLagrangeShapeFunctions(int elementId) const;
    
private:
    // 辅助函数声明
    bool isTransientAnalysis() const;
    bool getElementVectorPotential(int elementId, std::vector<double>& vectorPotential) const;
    bool hasExternalCurrentSource(int elementId) const;
    std::array<double, 3> getExternalCurrentDensity(int elementId) const;
    
    // 几何计算辅助函数
    std::array<double, 3> calculateEdgeVector(const std::array<double, 3>& node1, 
                                             const std::array<double, 3>& node2) const;
    std::vector<std::array<double, 3>> getElementNodeCoordinates(int elementId) const;
    
    // 3D特定参数
    bool useWhitneyElements_ = true;          ///< 使用Whitney边元
    bool usePiolaTransformation_ = false;     ///< 使用Piola变换
    bool useSecondOrderElements_ = false;     ///< 使用二阶单元
    bool useLagrangeGauge_ = false;           ///< 使用拉格朗日规范
    
    // 3D特定场量存储
    std::vector<std::array<double, 3>> magneticFluxDensity3D_;      ///< 磁通密度 (Bx, By, Bz)
    std::vector<std::array<double, 3>> magneticFieldStrength3D_;    ///< 磁场强度 (Hx, Hy, Hz)
    std::vector<std::array<double, 3>> currentDensity3D_;           ///< 电流密度 (Jx, Jy, Jz)
    
    // 复数场量存储（谐波分析）
    std::vector<std::array<std::complex<double>, 3>> complexMagneticFluxDensity3D_;
    std::vector<std::array<std::complex<double>, 3>> complexMagneticFieldStrength3D_;
    std::vector<std::array<std::complex<double>, 3>> complexCurrentDensity3D_;
    
    // 集总参数
    double torque3D_ = 0.0;
    double magneticEnergy3D_ = 0.0;
    double inductance3D_ = 0.0;
    
    // 复数集总参数（谐波分析）
    std::complex<double> complexTorque3D_ = {0.0, 0.0};
    
    // 系统矩阵
    std::shared_ptr<elmer::Matrix> systemMatrix_;          ///< 系统刚度矩阵
    std::shared_ptr<elmer::Matrix> massMatrix_;            ///< 质量矩阵
    std::shared_ptr<elmer::Matrix> dampingMatrix_;         ///< 阻尼矩阵
    std::vector<double> rhsVector_;                        ///< 右端向量
    
    // 复数系统矩阵（谐波分析）
    std::shared_ptr<elmer::Matrix> complexSystemMatrix_;   ///< 复数系统矩阵
    std::shared_ptr<elmer::Matrix> complexMassMatrix_;     ///< 复数质量矩阵
    std::vector<std::complex<double>> complexRhsVector_;   ///< 复数右端向量
    std::complex<double> complexMagneticEnergy3D_ = {0.0, 0.0};
    std::complex<double> complexInductance3D_ = {0.0, 0.0};
    
    // Whitney边元相关数据
    std::vector<int> edgeDegreesOfFreedom_;   ///< 边自由度映射
    std::vector<std::array<int, 2>> edgeConnectivity_; ///< 边连接关系
    
    // 非线性材料模型支持
    bool useNonlinearMaterialModel_ = false;  ///< 使用非线性材料模型
    std::vector<double> nonlinearPermeability_; ///< 非线性磁导率存储
    
    // 谐波分析支持
    bool useHarmonicAnalysis_ = false;        ///< 使用谐波分析
    double frequency_ = 0.0;                  ///< 频率 (Hz)
    
private:
    // 非线性材料模型支持方法
    bool updateNonlinearMaterialProperties();
    double calculateNonlinearPermeability(double H_magnitude, int elementId) const;
    
    // 谐波分析方法
    bool solveHarmonic();
    bool solveSteadyState();
    bool solveTransient();
    bool assembleComplexSystem();
    bool solveComplexLinearSystem();
    bool calculateComplexFields();
    
    // 边界条件处理方法
    bool applyDirichletBoundaryConditions();
    bool applyNeumannBoundaryConditions();
    bool applyPeriodicBoundaryConditions();
    bool applyCircuitCouplingBoundaryConditions();
    
    // 边界条件辅助方法
    std::vector<int> getBoundaryConditions() const;
    bool isDirichletBoundary(const elmer::Element& element) const;
    bool isNeumannBoundary(const elmer::Element& element) const;
    double getDirichletBoundaryValue(const elmer::Element& element) const;
    double getNeumannBoundaryValue(const elmer::Element& element) const;
    bool applyDirichletToSystem(const elmer::Element& element, double value);
    bool applyNeumannToRHS(const elmer::Element& element, double fluxValue);
    std::vector<std::pair<int, int>> getPeriodicBoundaryPairs() const;
    bool applyPeriodicConstraint(int masterDOF, int slaveDOF);
    std::vector<double> getCircuitCouplingData() const;
    bool applyCircuitCouplingToSystem(const std::vector<double>& data);
    bool applyCircuitCouplingToRHS(const std::vector<double>& data);
    
    // 系统组装辅助方法
    bool initializeSystemMatrices();
    bool assembleElementContributions(int elementId);
    bool assembleWhitneyElementContributions(int elementId);
    bool assembleLagrangeElementContributions(int elementId);
    size_t calculateTotalDegreesOfFreedom() const;
    size_t applyDOFConstraints(size_t totalDOFs) const;
    
    // 全局组装辅助方法
    std::vector<int> getElementDOFMapping(int elementId) const;
    bool assembleMatrixToGlobal(const std::vector<std::vector<double>>& elementMatrix,
                               const std::vector<int>& dofMapping,
                               std::shared_ptr<elmer::Matrix>& globalMatrix);
    bool assembleVectorToGlobal(const std::vector<double>& elementVector,
                               const std::vector<int>& dofMapping,
                               std::vector<double>& globalVector);
    
    // 求解器方法
    bool solveLinearSystem();
    bool assembleTransientSystem(double currentTime);
    bool solveTransientSystem();
    double getTimeStep() const;
    double getEndTime() const;
    int getMaxTimeSteps() const;
    bool shouldOutputTimeStep(int step) const;
    void outputTimeStepResults(int step, double currentTime);
    
    // Whitney边元组装方法
    std::vector<std::vector<double>> calculateWhitneyElementStiffnessMatrix(int elementId, double conductivity, double permeability) const;
    std::vector<std::vector<double>> calculateWhitneyElementMassMatrix(int elementId, double conductivity) const;
    std::vector<double> calculateWhitneyElementLoadVector(int elementId) const;
    bool assembleToGlobalSystem(int elementId, 
                               const std::vector<std::vector<double>>& stiffnessMatrix,
                               const std::vector<std::vector<double>>& massMatrix,
                               const std::vector<double>& loadVector);
    int getNumberOfElementEdges(int elementId) const;
    
    // Lagrange单元组装方法
    std::vector<std::vector<double>> calculateLagrangeElementStiffnessMatrix(int elementId, double conductivity, double permeability) const;
    std::vector<std::vector<double>> calculateLagrangeElementMassMatrix(int elementId, double conductivity) const;
    std::vector<double> calculateLagrangeElementLoadVector(int elementId) const;
    int getNumberOfElementNodes(int elementId) const;
};

} // namespace elmer