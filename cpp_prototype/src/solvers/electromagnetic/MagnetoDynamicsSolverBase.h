#pragma once

#include "ElectromagneticMaterial.h"
#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include "ElementMatrix.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "MaterialDatabase.h"
#include "IterativeSolver.h"
#include "NonlinearSolver.h"
#include "LoggerFactory.h"
#include "CommonConstants.h"
#include <memory>
#include <vector>
#include <array>
#include <complex>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cmath>

namespace elmer {

/**
 * @brief 低频电磁求解器维度类型
 */
enum class MagnetoDynamicsDimension {
    DIM_2D,           ///< 二维求解（笛卡尔或轴对称）
    DIM_3D            ///< 三维求解
};

/**
 * @brief 低频电磁分析类型
 */
enum class MagnetoDynamicsAnalysisType {
    STEADY_STATE,     ///< 稳态分析
    TRANSIENT,        ///< 瞬态分析
    HARMONIC          ///< 谐波分析
};

/**
 * @brief 坐标系类型
 */
enum class CoordinateSystemType {
    CARTESIAN,           ///< 笛卡尔坐标系
    AXISYMMETRIC,        ///< 轴对称坐标系
    CYLINDRIC_SYMMETRIC  ///< 柱对称坐标系
};

/**
 * @brief 低频电磁求解器通用参数
 */
struct MagnetoDynamicsParameters {
    // 基本求解器参数
    double tolerance = 1.0e-8;        ///< 收敛容差
    int maxIterations = 1000;         ///< 最大迭代次数
    int maxNonlinearIterations = 10;  ///< 最大非线性迭代次数
    
    // 分析类型
    MagnetoDynamicsAnalysisType analysisType = MagnetoDynamicsAnalysisType::STEADY_STATE;
    bool isTransient = false;         ///< 瞬态分析标志
    bool isHarmonic = false;          ///< 谐波分析标志
    
    // 时间相关参数
    double timeStep = 0.0;            ///< 时间步长（瞬态分析）
    double endTime = 0.0;             ///< 结束时间（瞬态分析）
    double frequency = 0.0;           ///< 谐波分析频率 [Hz]
    
    // 物理模型选项
    bool includeEddyCurrents = true;  ///< 包含涡流效应
    bool includeDisplacementCurrent = false;  ///< 包含位移电流
    bool includeConvection = false;           ///< 包含对流项
    bool useNewtonRaphson = false;            ///< 使用牛顿-拉夫逊迭代
    bool useLagrangeGauge = false;            ///< 使用拉格朗日规范
    
    // 坐标系
    CoordinateSystemType coordinateSystem = CoordinateSystemType::CARTESIAN;
    
    // 输出选项
    bool calculateMagneticField = true;       ///< 计算磁场
    bool calculateElectricField = false;      ///< 计算电场
    bool calculateCurrentDensity = false;     ///< 计算电流密度
    bool calculateTorque = false;             ///< 计算转矩
    bool calculateLumpedParameters = false;   ///< 计算集总参数
    
    // 边界条件选项
    bool useInfinityBC = false;               ///< 使用无限远边界条件
    bool useAirGapBC = false;                 ///< 使用气隙边界条件
    bool useSkinBC = false;                   ///< 使用集肤效应边界条件
    
    MagnetoDynamicsParameters() = default;
};

/**
 * @brief 低频电磁求解器通用结果
 */
struct MagnetoDynamicsResults {
    // 求解状态
    bool success = false;                    ///< 求解是否成功
    double executionTime = 0.0;              ///< 执行时间
    int iterations = 0;                      ///< 迭代次数
    double finalResidual = 0.0;              ///< 最终残差
    std::string message;                     ///< 结果消息
    
    // 集总参数
    double torque = 0.0;                     ///< 转矩 [N·m]
    double magneticEnergy = 0.0;             ///< 磁能 [J]
    double inductance = 0.0;                 ///< 电感 [H]
    double resistance = 0.0;                 ///< 电阻 [Ω]
    
    // 复数集总参数（谐波分析）
    std::complex<double> complexTorque = {0.0, 0.0};
    std::complex<double> complexMagneticEnergy = {0.0, 0.0};
    std::complex<double> complexInductance = {0.0, 0.0};
    std::complex<double> complexResistance = {0.0, 0.0};
    
    MagnetoDynamicsResults() = default;
};

/**
 * @brief 低频电磁求解器基类
 * 
 * 提供统一的接口和基础功能，整合2D和3D求解器
 */
class MagnetoDynamicsSolverBase {
public:
    MagnetoDynamicsSolverBase(MagnetoDynamicsDimension dimension, 
                             CoordinateSystemType coordSystem = CoordinateSystemType::CARTESIAN);
    virtual ~MagnetoDynamicsSolverBase() = default;
    
    // 求解器配置
    virtual void setParameters(const MagnetoDynamicsParameters& params);
    virtual MagnetoDynamicsParameters getParameters() const;
    
    // 网格和材料管理
    virtual bool loadMesh(std::shared_ptr<Mesh> mesh);
    virtual bool setMaterialProperties(const std::string& materialName, 
                                      const ElectromagneticMaterial& material);
    
    // 边界条件管理
    virtual void addBoundaryCondition(std::shared_ptr<BoundaryCondition> bc);
    virtual void setBoundaryCondition(const std::string& bcName, 
                                     const std::vector<double>& values);
    
    // 求解过程
    virtual bool initialize();
    virtual bool assembleSystem();
    virtual bool solve();
    virtual MagnetoDynamicsResults execute();
    
    // 结果获取
    virtual std::vector<double> getSolution() const;
    virtual std::vector<std::complex<double>> getComplexSolution() const;
    virtual MagnetoDynamicsResults getResults() const;
    
    // 导出场量计算
    virtual std::vector<std::array<double, 3>> calculateMagneticFluxDensity() const;
    virtual std::vector<std::array<double, 3>> calculateMagneticFieldStrength() const;
    virtual std::vector<double> calculateCurrentDensity() const;
    
    // 复数场量计算（谐波分析）
    virtual std::vector<std::array<std::complex<double>, 3>> calculateComplexMagneticFluxDensity() const;
    virtual std::vector<std::array<std::complex<double>, 3>> calculateComplexMagneticFieldStrength() const;
    virtual std::vector<std::complex<double>> calculateComplexCurrentDensity() const;
    
    // 状态查询
    virtual bool isInitialized() const;
    virtual bool isMeshLoaded() const;
    virtual std::string getStatus() const;
    
    // 维度信息
    MagnetoDynamicsDimension getDimension() const { return dimension_; }
    CoordinateSystemType getCoordinateSystem() const { return coordinateSystem_; }
    
protected:
    // 维度相关参数
    MagnetoDynamicsDimension dimension_;
    CoordinateSystemType coordinateSystem_;
    
    // 求解器状态
    bool initialized_ = false;
    bool meshLoaded_ = false;
    bool systemAssembled_ = false;
    
    // 求解器参数
    MagnetoDynamicsParameters parameters_;
    MagnetoDynamicsResults results_;
    
    // 网格和材料
    std::shared_ptr<elmer::Mesh> mesh_;
    std::shared_ptr<MaterialDatabase> materialDB_;
    std::shared_ptr<BoundaryConditionManager> bcManager_;
    
    // 线性代数系统
    std::shared_ptr<Matrix> systemMatrix_;
    std::vector<double> rhsVector_;
    std::vector<double> solutionVector_;
    std::vector<std::complex<double>> complexSolutionVector_;
    
    // 求解器
    // TODO: 修复抽象类实例化问题，暂时使用具体实现类
    std::unique_ptr<ConjugateGradientSolver> linearSolver_;
    std::unique_ptr<NewtonRaphsonSolver> nonlinearSolver_;
    
    // 保护方法 - 由具体实现类重写
    virtual bool assembleElementMatrix(int elementId, ElementMatrix& elementMatrix);
    virtual bool assembleElementRHS(int elementId, std::vector<double>& elementRHS);
    virtual bool applyBoundaryConditions();
    virtual bool updateMaterialProperties();
    virtual bool calculateDerivedFields();
    
    // 工具方法
    virtual double calculateElementConductivity(int elementId) const;
    virtual double calculateElementPermeability(int elementId) const;
    virtual std::array<double, 3> calculateElementMagneticField(int elementId) const;
    
    // 维度特定的虚方法
    virtual int getDegreesOfFreedom() const = 0;
    virtual std::string getVariableName() const = 0;
    virtual bool isComplexAnalysisRequired() const;
};

} // namespace elmer