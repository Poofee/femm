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
    double calculateElementConductivity3D(int elementId) const;
    double calculateElementPermeability3D(int elementId) const;
    
    // 3D特定场量计算
    std::array<double, 3> calculateElementMagneticFluxDensity(int elementId) const;
    std::array<double, 3> calculateElementMagneticFieldStrength(int elementId) const;
    std::array<double, 3> calculateElementCurrentDensity(int elementId) const;
    
    // 复数场量计算（谐波分析）
    std::array<std::complex<double>, 3> calculateElementComplexMagneticFluxDensity(int elementId) const;
    std::array<std::complex<double>, 3> calculateElementComplexMagneticFieldStrength(int elementId) const;
    std::array<std::complex<double>, 3> calculateElementComplexCurrentDensity(int elementId) const;
    
    // Whitney边元特定方法
    bool calculateWhitneyShapeFunctions(int elementId, std::vector<std::array<double, 3>>& shapeFuncs) const;
    bool calculateWhitneyCurlShapeFunctions(int elementId, std::vector<std::array<double, 3>>& curlShapeFuncs) const;
    
private:
    // 辅助函数声明
    bool isTransientAnalysis() const;
    bool getElementVectorPotential(int elementId, std::vector<double>& vectorPotential) const;
    bool hasExternalCurrentSource(int elementId) const;
    std::array<double, 3> getExternalCurrentDensity(int elementId) const;
    
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
    std::complex<double> complexMagneticEnergy3D_ = {0.0, 0.0};
    std::complex<double> complexInductance3D_ = {0.0, 0.0};
    
    // Whitney边元相关数据
    std::vector<int> edgeDegreesOfFreedom_;   ///< 边自由度映射
    std::vector<std::array<int, 2>> edgeConnectivity_; ///< 边连接关系
};

} // namespace elmer