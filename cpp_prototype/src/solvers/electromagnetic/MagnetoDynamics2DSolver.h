#pragma once

#include "MagnetoDynamicsSolverBase.h"
#include "ElectromagneticMaterial.h"
#include <memory>
#include <vector>
#include <array>
#include <complex>

namespace elmer {

/**
 * @brief 2D磁动力学求解器
 * 
 * 基于Elmer Fortran MagnetoDynamics2D求解器的C++实现
 * 支持笛卡尔和轴对称坐标系下的低频电磁场分析
 */
class MagnetoDynamics2DSolver : public MagnetoDynamicsSolverBase {
public:
    MagnetoDynamics2DSolver(CoordinateSystemType coordSystem = CoordinateSystemType::CARTESIAN);
    virtual ~MagnetoDynamics2DSolver() = default;
    
    // 重写基类方法
    virtual bool initialize() override;
    virtual bool assembleSystem() override;
    virtual bool solve() override;
    
    // 2D特定方法
    void setAxisymmetricRadius(double radius);
    double getAxisymmetricRadius() const;
    
    // 结果获取
    std::vector<double> getVectorPotential() const;
    std::vector<std::complex<double>> getComplexVectorPotential() const;
    
    // 2D特定场量计算
    std::vector<std::array<double, 2>> calculateMagneticFluxDensity2D() const;
    std::vector<std::array<double, 2>> calculateMagneticFieldStrength2D() const;
    std::vector<double> calculateCurrentDensity2D() const;
    
    // 复数场量计算（谐波分析）
    std::vector<std::array<std::complex<double>, 2>> calculateComplexMagneticFluxDensity2D() const;
    std::vector<std::array<std::complex<double>, 2>> calculateComplexMagneticFieldStrength2D() const;
    std::vector<std::complex<double>> calculateComplexCurrentDensity2D() const;
    
    // 集总参数计算
    double calculateTorque2D() const;
    double calculateMagneticEnergy2D() const;
    double calculateInductance2D() const;
    
    // 复数集总参数（谐波分析）
    std::complex<double> calculateComplexTorque2D() const;
    std::complex<double> calculateComplexMagneticEnergy2D() const;
    std::complex<double> calculateComplexInductance2D() const;
    
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
    
    // 2D特定工具方法
    bool assembleCartesianElementMatrix(int elementId, ElementMatrix& elementMatrix);
    bool assembleAxisymmetricElementMatrix(int elementId, ElementMatrix& elementMatrix);
    bool assembleCylindricElementMatrix(int elementId, ElementMatrix& elementMatrix);
    
    double calculateElementArea2D(int elementId) const;
    double calculateElementConductivity2D(int elementId) const;
    double calculateElementPermeability2D(int elementId) const;
    
    // 2D特定场量计算
    std::array<double, 2> calculateElementMagneticFluxDensity(int elementId) const;
    std::array<double, 2> calculateElementMagneticFieldStrength(int elementId) const;
    double calculateElementCurrentDensity(int elementId) const;
    
    // 复数场量计算（谐波分析）
    std::array<std::complex<double>, 2> calculateElementComplexMagneticFluxDensity(int elementId) const;
    std::array<std::complex<double>, 2> calculateElementComplexMagneticFieldStrength(int elementId) const;
    std::complex<double> calculateElementComplexCurrentDensity(int elementId) const;
    
private:
    // 2D特定参数
    double axisymmetricRadius_ = 1.0;  ///< 轴对称半径（轴对称坐标系）
    
    // 2D特定场量存储
    std::vector<std::array<double, 2>> magneticFluxDensity2D_;      ///< 磁通密度 (Bx, By)
    std::vector<std::array<double, 2>> magneticFieldStrength2D_;    ///< 磁场强度 (Hx, Hy)
    std::vector<double> currentDensity2D_;                          ///< 电流密度 Jz
    
    // 复数场量存储（谐波分析）
    std::vector<std::array<std::complex<double>, 2>> complexMagneticFluxDensity2D_;
    std::vector<std::array<std::complex<double>, 2>> complexMagneticFieldStrength2D_;
    std::vector<std::complex<double>> complexCurrentDensity2D_;
    
    // 集总参数
    double torque2D_ = 0.0;
    double magneticEnergy2D_ = 0.0;
    double inductance2D_ = 0.0;
    
    // 复数集总参数（谐波分析）
    std::complex<double> complexTorque2D_ = {0.0, 0.0};
    std::complex<double> complexMagneticEnergy2D_ = {0.0, 0.0};
    std::complex<double> complexInductance2D_ = {0.0, 0.0};
};

} // namespace elmer