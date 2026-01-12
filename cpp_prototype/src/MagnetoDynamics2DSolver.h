#pragma once

#include "ElectromagneticMaterial.h"
#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include "ElementMatrix.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "IterativeSolver.h"
#include "NonlinearSolver.h"
#include <memory>
#include <vector>
#include <array>
#include <complex>
#include <unordered_map>
#include <omp.h>

namespace elmer {

/**
 * @brief 2D磁动力学求解器参数
 * 
 * 基于Fortran MagnetoDynamics2D求解器的参数配置
 */
struct MagnetoDynamics2DParameters {
    // 求解器控制
    double tolerance = 1.0e-8;        ///< 收敛容差
    int maxIterations = 1000;         ///< 最大迭代次数
    int maxNonlinearIterations = 10;  ///< 最大非线性迭代次数
    
    // 分析类型
    bool isTransient = false;         ///< 瞬态分析标志
    bool isHarmonic = false;          ///< 谐波分析标志
    double frequency = 0.0;           ///< 谐波分析频率 [Hz]
    double timeStep = 0.0;            ///< 时间步长（瞬态分析）
    
    // 谐波分析特定参数
    bool useComplexMatrices = false;  ///< 使用复数矩阵（谐波分析）
    double harmonicPhase = 0.0;       ///< 谐波相位 [rad]
    bool includeEddyCurrents = true;  ///< 包含涡流效应
    bool includeDisplacementCurrent = false;  ///< 包含位移电流（高频）
    
    // 物理参数
    bool includeConvection = false;           ///< 包含对流项
    bool useNewtonRaphson = false;            ///< 使用牛顿-拉夫逊迭代
    
    // 坐标系
    enum CoordinateSystem {
        CARTESIAN,           ///< 笛卡尔坐标系
        AXISYMMETRIC,        ///< 轴对称坐标系
        CYLINDRIC_SYMMETRIC  ///< 柱对称坐标系
    } coordinateSystem = CARTESIAN;
    
    // 输出控制
    bool calculateMagneticField = true;       ///< 计算磁场
    bool calculateElectricField = false;      ///< 计算电场
    bool calculateCurrentDensity = false;     ///< 计算电流密度
    bool calculateTorque = false;             ///< 计算转矩
    bool calculateLumpedParameters = false;   ///< 计算集总参数
    
    // 边界条件类型
    bool useInfinityBC = false;               ///< 使用无限远边界条件
    bool useAirGapBC = false;                 ///< 使用气隙边界条件
    bool useSkinBC = false;                   ///< 使用集肤效应边界条件
    
    MagnetoDynamics2DParameters() = default;
};

/**
 * @brief 2D磁动力学求解结果
 */
struct MagnetoDynamics2DResults {
    // 主解（磁矢势A_z）
    std::vector<double> vectorPotential;      ///< 磁矢势 [Wb/m]（实数解）
    std::vector<std::complex<double>> complexVectorPotential; ///< 复数磁矢势（谐波分析）
    
    // 导出场量
    std::vector<std::array<double, 2>> magneticFluxDensity;  ///< 磁通密度 [T] (B_x, B_y)
    std::vector<std::array<double, 2>> magneticFieldStrength; ///< 磁场强度 [A/m] (H_x, H_y)
    std::vector<double> currentDensity;       ///< 电流密度 [A/m²] (J_z)
    
    // 复数导出场量（谐波分析）
    std::vector<std::array<std::complex<double>, 2>> complexMagneticFluxDensity;  ///< 复数磁通密度
    std::vector<std::array<std::complex<double>, 2>> complexMagneticFieldStrength; ///< 复数磁场强度
    std::vector<std::complex<double>> complexCurrentDensity; ///< 复数电流密度
    
    // 集总参数
    double torque = 0.0;                      ///< 转矩 [N·m]
    double magneticEnergy = 0.0;              ///< 磁能 [J]
    double inductance = 0.0;                  ///< 电感 [H]
    
    // 复数集总参数（谐波分析）
    std::complex<double> complexImpedance;    ///< 复数阻抗 [Ω]
    std::complex<double> complexInductance;   ///< 复数电感 [H]
    double powerLoss = 0.0;                   ///< 功率损耗 [W]
    
    // 收敛信息
    int iterations = 0;                       ///< 迭代次数
    int nonlinearIterations = 0;              ///< 非线性迭代次数
    double residual = 0.0;                    ///< 最终残差
    bool converged = false;                   ///< 收敛状态
    
    MagnetoDynamics2DResults() = default;
};

/**
 * @brief 2D磁动力学求解器
 * 
 * 实现基于A-φ公式的二维磁动力学求解，支持笛卡尔和柱对称坐标系。
 * 对应Fortran的MagnetoDynamics2D求解器。
 */
class MagnetoDynamics2DSolver {
private:
    std::shared_ptr<Mesh> mesh;
    elmer::MaterialDatabase materialDB;
    MagnetoDynamics2DParameters parameters;
    
    // 系统矩阵
    std::shared_ptr<Matrix> stiffnessMatrix;      ///< 刚度矩阵
    std::shared_ptr<Matrix> massMatrix;           ///< 质量矩阵
    std::shared_ptr<Matrix> dampingMatrix;        ///< 阻尼矩阵
    std::shared_ptr<Vector> rhsVector;            ///< 右端向量
    
    // 复数系统矩阵（谐波分析）
    std::shared_ptr<Matrix> complexStiffnessMatrix;      ///< 复数刚度矩阵
    std::shared_ptr<Matrix> complexMassMatrix;           ///< 复数质量矩阵
    std::shared_ptr<Matrix> complexDampingMatrix;        ///< 复数阻尼矩阵
    std::shared_ptr<Vector> complexRhsVector;            ///< 复数右端向量
    
    // 边界条件管理器
    BoundaryConditionManager bcManager;
    
    // 内部状态
    bool systemAssembled = false;
    bool useBasisFunctionsCache = false;
    
    // 非线性求解器
    std::unique_ptr<NonlinearSolver> nonlinearSolver;
    NonlinearMaterialDatabase nonlinearMaterialDB;
    
    // 非线性迭代状态
    std::vector<double> currentPotential;      ///< 当前向量势解
    std::vector<double> magneticFluxDensity;   ///< 当前磁通密度
    std::vector<double> magneticFieldStrength; ///< 当前磁场强度
    
    // 基函数缓存（性能优化）
    struct BasisFunctionCache {
        std::vector<double> basis;         ///< 基函数值
        std::vector<std::array<double, 3>> dBasisdx; ///< 基函数导数
        double weight;                     ///< 积分权重
        double detJ;                       ///< 雅可比行列式
        std::array<double, 3> coords;     ///< 积分点坐标
    };
    
    // 优化的缓存结构
    struct ElementCache {
        std::vector<BasisFunctionCache> integrationPoints; ///< 积分点缓存
        std::vector<int> nodeIndices;                      ///< 节点索引
        int elementType;                                   ///< 单元类型
        double area;                                       ///< 单元面积
    };
    
    std::vector<ElementCache> elementCache;               ///< 单元级缓存
    std::unordered_map<int, std::vector<double>> shapeFunctionCache; ///< 形状函数缓存
    
    // 矩阵组装优化
    struct MatrixAssemblyCache {
        std::vector<double> elementStiffness;     ///< 单元刚度矩阵缓存
        std::vector<double> elementMass;          ///< 单元质量矩阵缓存
        std::vector<double> elementDamping;        ///< 单元阻尼矩阵缓存
        std::vector<double> elementRHS;           ///< 单元右端向量缓存
        std::vector<int> elementDOFs;             ///< 单元自由度索引
    };
    
    MatrixAssemblyCache assemblyCache;            ///< 矩阵组装缓存

public:
    /**
     * @brief 构造函数
     */
    MagnetoDynamics2DSolver(std::shared_ptr<Mesh> m = nullptr)
        : mesh(m) {
        materialDB.createPredefinedMaterials();
        // 创建默认的非线性求解器
        nonlinearSolver = std::make_unique<NewtonRaphsonSolver>();
    }
    
    /**
     * @brief 设置网格
     */
    void setMesh(std::shared_ptr<Mesh> m) {
        mesh = m;
        systemAssembled = false;
    }
    
    /**
     * @brief 设置材料数据库
     */
    void setMaterialDatabase(const elmer::MaterialDatabase& db) {
        materialDB = db;
    }
    
    /**
     * @brief 设置求解器参数
     */
    void setParameters(const MagnetoDynamics2DParameters& params) {
        parameters = params;
        systemAssembled = false;
    }
    
    /**
     * @brief 添加边界条件
     */
    void addBoundaryCondition(std::shared_ptr<BoundaryCondition> bc) {
        bcManager.addBoundaryCondition(bc);
        systemAssembled = false;
    }
    
    /**
     * @brief 获取边界条件管理器
     */
    BoundaryConditionManager& getBoundaryConditionManager() {
        return bcManager;
    }
    
    /**
     * @brief 获取边界条件管理器（常量版本）
     */
    const BoundaryConditionManager& getBoundaryConditionManager() const {
        return bcManager;
    }
    
    // 性能优化方法
    /**
     * @brief 启用基函数缓存
     */
    void enableBasisFunctionCache(bool enable = true) {
        useBasisFunctionsCache = enable;
    }
    
    /**
     * @brief 预计算基函数缓存
     */
    void precomputeBasisFunctionCache();
    
    /**
     * @brief 优化的矩阵组装方法
     */
    void assembleSystemOptimized();
    
    /**
     * @brief 并行矩阵组装
     */
    void assembleSystemParallel();
    
    /**
     * @brief 增量式矩阵更新（非线性迭代优化）
     */
    void updateSystemIncremental();
    
    /**
     * @brief 清除缓存
     */
    void clearCache();
    
    /**
     * @brief 获取缓存统计信息
     */
    void getCacheStatistics() const;
    
    /**
     * @brief 初始化求解器"
     * 
     * 对应Fortran的MagnetoDynamics2D_Init子程序
     */
    void initialize() {
        if (!mesh) {
            throw std::runtime_error("Mesh not set for initialization");
        }
        
        // 检查坐标系
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            // 轴对称情况下的特殊处理
            if (parameters.useInfinityBC) {
                std::cout << "Warning: Infinity BC not yet available in axisymmetric case!" << std::endl;
            }
        }
        
        // 初始化基函数缓存（如果启用）
        if (useBasisFunctionsCache) {
            initializeBasisFunctionCache();
        }
        
        std::cout << "MagnetoDynamics2D solver initialized" << std::endl;
    }
    
    /**
     * @brief 组装系统矩阵
     * 
     * 对应Fortran的LocalMatrix子程序
     */
    void assembleSystem() {
        if (!mesh) {
            throw std::runtime_error("Mesh not set for assembly");
        }
        
        size_t nNodes = mesh->getNodes().numberOfNodes();
        
        // 初始化系统矩阵（每个节点1个自由度：A_z）
        stiffnessMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
        rhsVector = std::shared_ptr<Vector>(Vector::Create(nNodes));
        
        if (parameters.isTransient) {
            massMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
            dampingMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
        }
        
        // 组装单元贡献
        assembleElementContributions();
        
        // 组装边界条件贡献
        assembleBoundaryContributions();
        
        // 应用边界条件
        applyBoundaryConditions();
        
        systemAssembled = true;
        
        std::cout << "System assembly completed" << std::endl;
    }
    
    /**
     * @brief 设置非线性求解器
     */
    void setNonlinearSolver(std::unique_ptr<NonlinearSolver> solver) {
        nonlinearSolver = std::move(solver);
    }
    
    /**
     * @brief 添加非线性材料
     */
    void addNonlinearMaterial(const std::string& name, const ElectromagneticMaterial& material) {
        nonlinearMaterialDB.addMaterial(name, material);
    }
    
    /**
     * @brief 求解磁动力学问题
     * 
     * 对应Fortran的MagnetoDynamics2D主求解器
     */
    MagnetoDynamics2DResults solve() {
        if (!systemAssembled) {
            assembleSystem();
        }
        
        MagnetoDynamics2DResults results;
        
        // 检查是否需要非线性求解
        bool hasNonlinearMaterials = !nonlinearMaterialDB.getMaterialNames().empty();
        
        if (hasNonlinearMaterials && parameters.useNewtonRaphson) {
            // 使用非线性求解器
            results = solveNonlinear();
        } else if (hasNonlinearMaterials) {
            // 使用简单的非线性迭代
            results = solveSimpleNonlinear();
        } else {
            // 线性求解
            results = solveLinear();
        }
        
        return results;
    }
    
    /**
     * @brief 线性求解
     */
    MagnetoDynamics2DResults solveLinear() {
        MagnetoDynamics2DResults results;
        
        auto solution = solveLinearSystem();
        results.vectorPotential = solution;
        results.converged = true;
        results.iterations = 1;
        
        // 计算导出场量
        calculateDerivedFields(results);
        
        // 计算集总参数
        if (parameters.calculateLumpedParameters) {
            calculateLumpedParameters(results);
        }
        
        return results;
    }
    
    /**
     * @brief 简单非线性迭代求解
     */
    MagnetoDynamics2DResults solveSimpleNonlinear() {
        MagnetoDynamics2DResults results;
        
        // 初始解
        auto solution = solveLinearSystem();
        currentPotential = solution;
        
        // 非线性迭代循环
        for (int nonlinearIter = 0; nonlinearIter < parameters.maxNonlinearIterations; ++nonlinearIter) {
            std::cout << "Performing nonlinear iteration: " << nonlinearIter + 1 << std::endl;
            
            // 更新材料参数
            updateMaterialParameters();
            
            // 重新组装系统
            reassembleNonlinearSystem();
            
            // 求解线性系统
            solution = solveLinearSystem();
            
            // 更新结果
            results.vectorPotential = solution;
            results.nonlinearIterations = nonlinearIter + 1;
            currentPotential = solution;
            
            // 检查收敛
            if (checkConvergence(results)) {
                results.converged = true;
                break;
            }
        }
        
        // 计算导出场量
        calculateDerivedFields(results);
        
        // 计算集总参数
        if (parameters.calculateLumpedParameters) {
            calculateLumpedParameters(results);
        }
        
        return results;
    }
    
    /**
     * @brief 非线性求解器求解
     */
    MagnetoDynamics2DResults solveNonlinear() {
        MagnetoDynamics2DResults results;
        
        if (!nonlinearSolver) {
            throw std::runtime_error("Nonlinear solver not set");
        }
        
        // 初始猜测
        auto initialGuess = std::shared_ptr<Vector>(Vector::Create(mesh->getNodes().numberOfNodes()));
        for (size_t i = 0; i < initialGuess->Size(); ++i) {
            (*initialGuess)[i] = 0.0; // 零初始猜测
        }
        
        // 定义残差函数
        auto residualFunction = [this](const std::shared_ptr<Vector>& x) -> std::shared_ptr<Vector> {
            // 更新当前解
            currentPotential.resize(x->Size());
            for (size_t i = 0; i < x->Size(); ++i) {
                currentPotential[i] = (*x)[i];
            }
            
            // 更新材料参数
            updateMaterialParameters();
            
            // 重新组装系统
            reassembleNonlinearSystem();
            
            // 计算残差：r = Kx - f
            auto residual = std::shared_ptr<Vector>(Vector::Create(x->Size()));
            auto Kx = std::shared_ptr<Vector>(Vector::Create(x->Size()));
            stiffnessMatrix->Multiply(*x, *Kx);
            
            for (size_t i = 0; i < x->Size(); ++i) {
                (*residual)[i] = (*Kx)[i] - (*rhsVector)[i];
            }
            
            return residual;
        };
        
        // 定义雅可比函数
        auto jacobianFunction = [this](const std::shared_ptr<Vector>& x) -> std::shared_ptr<Matrix> {
            // 更新当前解
            currentPotential.resize(x->Size());
            for (size_t i = 0; i < x->Size(); ++i) {
                currentPotential[i] = (*x)[i];
            }
            
            // 更新材料参数
            updateMaterialParameters();
            
            // 重新组装系统（包含雅可比项）
            reassembleNonlinearSystemWithJacobian();
            
            return stiffnessMatrix; // 雅可比矩阵就是刚度矩阵
        };
        
        // 求解非线性系统
        auto nonlinearResults = nonlinearSolver->solve(initialGuess, residualFunction, jacobianFunction);
        
        // 更新结果
        results.vectorPotential = currentPotential;
        results.converged = nonlinearResults.converged;
        results.iterations = nonlinearResults.iterations;
        results.residual = nonlinearResults.residualNorm;
        
        // 计算导出场量
        calculateDerivedFields(results);
        
        // 计算集总参数
        if (parameters.calculateLumpedParameters) {
            calculateLumpedParameters(results);
        }
        
        return results;
    }
    
    /**
     * @brief 设置并行线程数
     */
    void setParallelThreads(int numThreads);
    
    /**
     * @brief 获取并行线程数
     */
    int getParallelThreads() const;
    
    // 谐波分析方法
    /**
     * @brief 谐波分析求解
     * 
     * 实现频率域分析，支持复数矩阵和涡流效应
     */
    MagnetoDynamics2DResults solveHarmonic();
    
    /**
     * @brief 组装复数系统矩阵
     * 
     * 用于谐波分析，包含涡流项和位移电流项
     */
    void assembleComplexSystem();
    
    /**
     * @brief 求解复数线性系统
     */
    std::vector<std::complex<double>> solveComplexLinearSystem();
    
    /**
     * @brief 计算复数导出场量
     */
    void calculateComplexDerivedFields(MagnetoDynamics2DResults& results);
    
    /**
     * @brief 计算复数集总参数
     */
    void calculateComplexLumpedParameters(MagnetoDynamics2DResults& results);
    
    /**
     * @brief 设置谐波激励源
     */
    void setHarmonicExcitation(double amplitude, double frequency, double phase = 0.0);
    
    /**
     * @brief 获取谐波分析结果（时域重构）
     */
    std::vector<double> reconstructTimeDomain(double time) const;
    
private:
    /**
     * @brief 组装单元贡献
     */
    void assembleElementContributions();
    
    /**
     * @brief 组装边界贡献
     */
    void assembleBoundaryContributions();
    
    /**
     * @brief 应用边界条件
     */
    void applyBoundaryConditions();
    
    /**
     * @brief 重新组装非线性系统
     */
    void reassembleNonlinearSystem();
    
    /**
     * @brief 重新组装非线性系统（包含雅可比项）
     */
    void reassembleNonlinearSystemWithJacobian();
    
    /**
     * @brief 更新材料参数
     */
    void updateMaterialParameters();
    
    /**
     * @brief 计算磁场强度
     */
    std::vector<double> computeMagneticFieldStrength();
    
    /**
     * @brief 计算磁通密度
     */
    std::vector<double> computeMagneticFluxDensity();
    
    /**
     * @brief 求解线性系统
     */
    std::vector<double> solveLinearSystem();
    
    /**
     * @brief 检查收敛性
     */
    bool checkConvergence(const MagnetoDynamics2DResults& results);
    
    /**
     * @brief 计算导出场量
     */
    void calculateDerivedFields(MagnetoDynamics2DResults& results);
    
    /**
     * @brief 计算集总参数
     */
    void calculateLumpedParameters(MagnetoDynamics2DResults& results);
    
    /**
     * @brief 初始化基函数缓存
     */
    void initializeBasisFunctionCache();
    
    /**
     * @brief 计算单元局部矩阵
     * 
     * 对应Fortran的LocalMatrix子程序的核心实现
     */
    void computeLocalMatrix(const Element& element, 
                           std::vector<std::vector<double>>& stiffness,
                           std::vector<std::vector<double>>& mass,
                           std::vector<std::vector<double>>& damping,
                           std::vector<double>& force);
    
    /**
     * @brief 计算磁阻率
     * 
     * 对应Fortran的GetReluctivity子程序
     */
    std::array<std::array<double, 2>, 2> computeReluctivity(const std::string& materialName, 
                                                           double magneticFluxDensity,
                                                           const Element& element);
    
    /**
     * @brief 处理无限远边界条件
     * 
     * 对应Fortran的LocalMatrixInfinityBC子程序
     */
    void applyInfinityBoundaryCondition(const Element& element, 
                                       std::vector<std::vector<double>>& stiffness,
                                       std::vector<double>& force);
    
    /**
     * @brief 处理气隙边界条件
     * 
     * 对应Fortran的LocalMatrixAirGapBC子程序
     */
    void applyAirGapBoundaryCondition(const Element& element, 
                                     std::vector<std::vector<double>>& stiffness,
                                     std::vector<double>& force);
    
    // 并行化相关方法
    /**
     * @brief 并行组装单元贡献（OpenMP实现）
     */
    void assembleElementContributionsParallel();
    
    /**
     * @brief 线程安全的矩阵元素添加
     */
    void addToMatrixThreadSafe(std::shared_ptr<Matrix>& matrix, int row, int col, double value);
    
    /**
     * @brief 线程安全的向量元素添加
     */
    void addToVectorThreadSafe(std::shared_ptr<Vector>& vector, int index, double value);
    
    // 并行化控制参数
    int numThreads = 4;  ///< 并行线程数（默认4线程）
    bool useParallelAssembly = false;  ///< 是否使用并行组装
};

} // namespace elmer