/**
 * @file MagneticSolve.h
 * @brief 磁动力学求解器模块
 * 
 * 对应Fortran模块：MagneticSolve.F90
 * 实现MHD Maxwell方程（或感应方程）的求解器。
 */

#pragma once

#include "SolverBase.h"
#include "ElectromagneticMaterial.h"
#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include "ElementMatrix.h"
#include "LinearAlgebra.h"
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "IterativeSolver.h"
#include <memory>
#include <vector>
#include <array>

namespace elmer {

// 引入elmer命名空间中的GaussIntegration
using elmer::GaussIntegration;

/**
 * @brief 磁动力学求解器参数
 */
struct MagneticSolveParameters {
    // 求解器控制
    double tolerance = 1.0e-8;        ///< 收敛容差
    int maxIterations = 1000;         ///< 最大迭代次数
    bool verbose = true;              ///< 详细输出
    
    // 物理参数
    bool includeDisplacementCurrent = false;  ///< 包含位移电流
    bool includeConvection = true;            ///< 包含对流项
    bool stabilize = true;                    ///< 稳定化
    
    // 分析类型
    bool transientSimulation = false;         ///< 瞬态模拟
    double timeStep = 1.0;                    ///< 时间步长
    
    // 输出控制
    bool calculateElectricField = false;      ///< 计算电场
    bool calculateLorentzForce = false;       ///< 计算洛伦兹力
    bool calculateCurrentDensity = false;     ///< 计算电流密度
    
    // 边界条件
    std::string velocitySolutionName = "Velocity";  ///< 速度解名称
    std::string forceVariableName = "Lorentz Force"; ///< 力变量名称
    
    MagneticSolveParameters() = default;
};

/**
 * @brief 磁动力学求解结果
 */
struct MagneticSolveResults {
    // 主解
    std::vector<double> magneticField;        ///< 磁场解
    
    // 导出场
    std::vector<std::array<double, 3>> electricField;      ///< 电场 [V/m]
    std::vector<std::array<double, 3>> currentDensity;     ///< 电流密度 [A/m²]
    std::vector<std::array<double, 3>> lorentzForce;       ///< 洛伦兹力 [N/m³]
    
    // 收敛信息
    int iterations = 0;                      ///< 迭代次数
    double residual = 0.0;                   ///< 最终残差
    bool converged = false;                  ///< 收敛状态
    
    // 能量量
    double magneticEnergy = 0.0;             ///< 磁能 [J]
    
    MagneticSolveResults() = default;
};

/**
 * @brief 元素节点坐标结构
 */
struct ElementNodes {
    std::vector<double> x;  ///< x坐标数组
    std::vector<double> y;  ///< y坐标数组
    std::vector<double> z;  ///< z坐标数组
    
    ElementNodes() = default;
    explicit ElementNodes(size_t size) {
        x.resize(size);
        y.resize(size);
        z.resize(size);
    }
};

/**
 * @brief 磁动力学求解器
 * 
 * 实现MHD Maxwell方程（或感应方程）的求解器，
 * 支持笛卡尔坐标系、轴对称坐标系和一般坐标系。
 */
class MagneticSolve : public SolverBase {
private:
    /**
     * @brief Shape function evaluation result
     */
    struct ShapeResult {
        std::vector<double> values;        ///< Shape function values
        std::vector<double> dNdxi;        ///< Derivatives in xi direction
        std::vector<double> dNdeta;       ///< Derivatives in eta direction
        std::vector<double> dNdzeta;      ///< Derivatives in zeta direction
        std::vector<double> dNdx;         ///< Derivatives in x direction
        std::vector<double> dNdy;         ///< Derivatives in y direction
        std::vector<double> dNdz;         ///< Derivatives in z direction
        double detJ;                      ///< Jacobian determinant
        
        ShapeResult(size_t nNodes) 
            : values(nNodes, 0.0), dNdxi(nNodes, 0.0), dNdeta(nNodes, 0.0),
              dNdzeta(nNodes, 0.0), dNdx(nNodes, 0.0), dNdy(nNodes, 0.0),
              dNdz(nNodes, 0.0), detJ(0.0) {}
    };
    
    MagneticSolveParameters parameters;
    
    // 系统矩阵
    std::shared_ptr<Matrix> massMatrix;
    std::shared_ptr<Matrix> stiffnessMatrix;
    std::shared_ptr<Vector> forceVector;
    
    // 临时存储
    std::vector<double> conductivity;        ///< 电导率
    std::vector<double> permeability;        ///< 磁导率
    std::vector<double> appliedFieldX;       ///< 施加的磁场X分量
    std::vector<double> appliedFieldY;       ///< 施加的磁场Y分量
    std::vector<double> appliedFieldZ;       ///< 施加的磁场Z分量
    std::vector<double> velocityX;           ///< 速度X分量
    std::vector<double> velocityY;           ///< 速度Y分量
    std::vector<double> velocityZ;           ///< 速度Z分量
    std::vector<double> meshVelocityX;       ///< 网格速度X分量
    std::vector<double> meshVelocityY;       ///< 网格速度Y分量
    std::vector<double> meshVelocityZ;       ///< 网格速度Z分量
    
    bool allocationsDone = false;            ///< 内存分配完成标志
    
public:
    /**
     * @brief 构造函数
     */
    MagneticSolve();
    
    /**
     * @brief 析构函数
     */
    virtual ~MagneticSolve();
    
    /**
     * @brief 设置求解器参数
     */
    void setParameters(const SolverBaseParameters& params) override {
        // 将通用求解器参数转换为磁动力学专用参数
        parameters.tolerance = params.tolerance;
        parameters.maxIterations = params.maxIterations;
        parameters.verbose = params.verbose;
    }
    
    /**
     * @brief 获取求解器参数
     */
    SolverBaseParameters getParameters() const override {
        SolverBaseParameters params;
        params.solverName = "MagneticSolve";
        params.tolerance = parameters.tolerance;
        params.maxIterations = parameters.maxIterations;
        params.verbose = parameters.verbose;
        return params;
    }
    
    /**
     * @brief 设置网格数据
     */
    void setMesh(std::shared_ptr<Mesh> mesh) {
        mesh_ = mesh;
    }
    
    /**
     * @brief 执行求解
     */
    virtual bool solve() override;
    
    /**
     * @brief 组装系统矩阵（SolverBase纯虚函数实现）
     */
    virtual bool assemble() override;
    
    /**
     * @brief 获取求解结果（SolverBase纯虚函数实现）
     */
    virtual std::vector<double> getSolution() const override;
    
    /**
     * @brief 组装系统矩阵
     */
    virtual void assembleSystem();
    
    /**
     * @brief 应用边界条件
     */
    virtual void applyBoundaryConditions();
    
    /**
     * @brief 求解线性系统
     */
    virtual void solveLinearSystem();
    
    /**
     * @brief 计算导出场
     */
    void computeDerivedFields(MagneticSolveResults& results);
    
    /**
     * @brief 计算洛伦兹力
     */
    void computeLorentzForce(std::vector<std::array<double, 3>>& lorentzForce);
    
    /**
     * @brief 计算电场
     */
    void computeElectricField(std::vector<std::array<double, 3>>& electricField);
    
    /**
     * @brief 计算电流密度
     */
    void computeCurrentDensity(std::vector<std::array<double, 3>>& currentDensity);
    
    /**
     * @brief 计算磁能
     */
    double computeMagneticEnergy();
    
    /**
     * @brief 检查自由表面
     */
    bool checkFreeSurface();
    
    /**
     * @brief 获取材料参数
     */
    void getMaterialParameters(const Element& element);
    
    /**
     * @brief 获取边界条件参数
     */
    void getBoundaryConditionParameters(const Element& element);
    
    /**
     * @brief 组装单元贡献（笛卡尔坐标系）
     */
    void assembleCartesianElement(const Element& element, const ElementNodes& nodes);
    
    /**
     * @brief 组装单元贡献（轴对称坐标系）
     */
    void assembleAxisymmetricElement(const Element& element, const ElementNodes& nodes);
    
    /**
     * @brief 组装单元贡献（一般坐标系）
     */
    void assembleGeneralElement(const Element& element, const ElementNodes& nodes);
    
    /**
     * @brief 组装边界条件贡献
     */
    void assembleBoundaryCondition(const Element& element, const ElementNodes& nodes);
    
    /**
     * @brief 获取元素的高斯积分点
     */
    // std::vector<elmer::GaussIntegration::IntegrationPoint> getGaussPointsForElement(const Element& element);
    
    /**
     * @brief 计算基函数及其导数
     */
    struct ShapeResult evaluateShapeFunctions(const Element& element, 
                                              const ElementNodes& nodes, 
                                              double xi, double eta, double zeta);
    
    /**
     * @brief 插值场量
     */
    double interpolateField(const std::vector<double>& shapeValues, 
                           const std::vector<double>& nodalValues, 
                           int nNodes);
    
    /**
     * @brief 插值场量导数
     */
    std::array<double, 3> interpolateFieldDerivatives(const std::vector<double>& shapeDerivatives, 
                                                     const std::vector<double>& nodalValues, 
                                                     int nNodes);
    
    /**
     * @brief 将局部矩阵组装到全局系统
     */
    void assembleLocalToGlobal(const Element& element, 
                              const std::vector<std::vector<double>>& localMassMatrix,
                              const std::vector<std::vector<double>>& localStiffMatrix,
                              const std::vector<double>& localForceVector);
    
    /**
     * @brief 初始化临时存储
     */
    void initializeTemporaryStorage();
    
    /**
     * @brief 清理临时存储
     */
    void cleanupTemporaryStorage();
};

/**
 * @brief 磁动力学求解器工厂函数
 */
std::shared_ptr<MagneticSolve> CreateMagneticSolve();

} // namespace elmer