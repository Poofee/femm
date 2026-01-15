#pragma once

#include "../core/math/LinearAlgebra.h"
#include "../core/mesh/Mesh.h"
#include "../core/base/SolverBase.h"
#include "../core/base/Types.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace elmer {

/**
 * @brief 磁场求解器类 - 移植自Fortran版本的MagneticSolve.F90
 * 
 * 实现MHD Maxwell方程（或感应方程）的求解器
 * 支持笛卡尔坐标系、柱对称坐标系和一般坐标系
 */
class MagneticSolver : public SolverBase {
public:
    // ===== 构造函数和析构函数 =====
    MagneticSolver();
    ~MagneticSolver();

    // ===== 求解器接口函数 =====
    std::string getName() const override;
    bool initialize() override;
    bool assemble() override;
    bool solve() override;
    std::vector<double> getSolution() const override;

    // ===== 磁场求解器特定功能 =====
    void setTransientSimulation(bool transient);
    void setTimeStep(double dt);
    
    // 物理量计算
    bool computeLorentzForce();
    bool computeElectricField();
    
    // 结果获取
    std::vector<double> getMagneticField() const;
    std::vector<double> getElectricCurrent() const;
    std::vector<double> getLorentzForce() const;
    std::vector<double> getElectricField() const;

private:
    // ===== 内部数据结构 =====
    
    // 场变量（基于Fortran数据结构）
    std::vector<double> magneticField_;      // MagneticField
    std::vector<double> electricCurrent_;    // ElectricCurrent
    std::vector<double> elecField_;          // ElecField
    std::vector<double> lorentzForce_;       // LrF
    
    // 速度场
    std::vector<double> U_, V_, W_;           // 流体速度
    std::vector<double> MU_, MV_, MW_;        // 网格速度
    
    // 材料参数
    std::vector<double> permeability_;        // Permeability
    std::vector<double> conductivity_;      // Conductivity
    
    // 辅助变量
    std::vector<double> divB_;               // divB
    std::vector<double> B1_, B2_, B3_;       // B场分量
    std::vector<double> ExBx_, ExBy_, ExBz_; // E × B项
    
    // 求解器参数
    bool transientSimulation_ = false;       // TransientSimulation
    double timeStep_ = 0.0;                  // dt
    int maxNonlinearIterations_ = 10;        // NonlinearIter
    double nonlinearTolerance_ = 1e-6;       // NonlinearTol
    bool stabilize_ = true;                  // Stabilize
    
    // 状态标志
    bool allocationsDone_ = false;          // AllocationsDone
    bool userDefinedVelo_ = false;           // UserDefinedVelo
    bool calculateMagneticForce_ = false;    // CalculateMagneticForce

    // ===== 私有辅助函数 =====
    
    // 内存管理
    bool allocateMemory();
    void deallocateMemory();
    
    // 参数获取
    bool getMaterialParameters();
    bool getVelocityField();
    
    // 坐标系检测和组装
    std::string detectCoordinateSystem() const;
    bool assembleCartesian();
    bool assembleAxisymmetric();
    bool assembleGeneral();
    bool assembleBoundaryConditions();
    
    // 非线性迭代求解
    bool solveNonlinearIteration(int maxIterations, double tolerance);
    std::shared_ptr<Vector> computeResidualVector();
    double computeVectorNorm(const Vector& vec) const;
    bool checkConvergence(double prevNorm, double currentNorm) const;
    bool updateJacobianMatrix();
    std::shared_ptr<Vector> solveLinearSystem(const Vector& residual);
    void updateSolutionVector(const Vector& deltaX);
    
    // 物理量更新
    void updateMagneticFieldFromSolution();
    
    // 元素级计算函数（基于Fortran代码）
    bool computeElementMatrices(const Element& element,
                               std::vector<std::vector<double>>& mass,
                               std::vector<std::vector<double>>& stiff,
                               std::vector<double>& load);
    
    // 边界条件处理
    bool applyDirichletBoundaryConditions();
    bool applyNeumannBoundaryConditions();
    bool applyMagneticBoundaryConditions();
    
    // 工具函数
    double getElementConductivity(const Element& element) const;
    double getElementPermeability(const Element& element) const;
    void getElementVelocity(const Element& element, 
                           std::vector<double>& u, 
                           std::vector<double>& v, 
                           std::vector<double>& w) const;
    
    // 洛伦兹力计算（基于Fortran子程序）
    bool computeLorentzForceNodal();
    
    // 高斯积分和形状函数计算
    void getGaussPoint3D(int gaussPoint, double& xi, double& eta, double& zeta, double& weight);
    void computeShapeFunctions3D(double xi, double eta, double zeta,
                                std::vector<double>& shapeFunctions,
                                std::vector<double>& dShapeDXi,
                                std::vector<double>& dShapeDEta,
                                std::vector<double>& dShapeDZeta);
    
    // 内存管理
    void deallocateMemory();
    
    // 边界条件应用函数
    bool applyMagneticForceBoundaryCondition(int bcId, const Element& boundaryElement);
    bool applyDirichletBoundaryCondition(int bcId, const Element& boundaryElement);
    bool applyNeumannBoundaryCondition(int bcId, const Element& boundaryElement);
    
    // 坐标系特定组装函数
    bool assembleCartesian();
    bool assembleAxisymmetric();
    bool assembleGeneral();
    
    // 单元矩阵计算函数
    bool computeAxisymmetricElementMatrices(int elementId, 
                                           std::vector<std::vector<double>>& elementStiffness,
                                           std::vector<double>& elementRHS);
    bool computeGeneralElementMatrices(int elementId, 
                                      std::vector<std::vector<double>>& elementStiffness,
                                      std::vector<double>& elementRHS);
    bool computeElementMatrices(int elementId, 
                               std::vector<std::vector<double>>& elementStiffness,
                               std::vector<double>& elementRHS);
    bool assembleToGlobalSystem(int elementId,
                               const std::vector<std::vector<double>>& elementStiffness,
                               const std::vector<double>& elementRHS);
    
    // 非线性迭代辅助函数
    std::unique_ptr<Vector> computeResidualVector();
    double computeVectorNorm(const Vector& vec);
    bool updateJacobianMatrix();
    std::unique_ptr<Vector> solveLinearSystem(const Vector& residual);
    bool updateSolutionVector(const Vector& deltaX);
    bool updateMagneticFieldFromSolution();
    
    // 坐标系检测
    std::string detectCoordinateSystem();
    
    // 旋度计算
    bool computeCurl(const std::vector<double>& fieldX,
                    const std::vector<double>& fieldY, 
                    const std::vector<double>& fieldZ,
                    std::vector<double>& curlX,
                    std::vector<double>& curlY,
                    std::vector<double>& curlZ);
    
    // 节点场计算
    bool computeNodalField();
    
    // 边界条件检测和周期性边界条件
    std::string detectBoundaryConditionType(int bcId, const Element& boundaryElement);
    bool applyPeriodicBoundaryCondition(int bcId, const Element& boundaryElement);
    
    // 非线性问题检测和雅可比矩阵计算
    bool checkNonlinearity();
    bool computeNonlinearElementJacobian(int elementId, 
                                       std::vector<std::vector<double>>& elementJacobian,
                                       std::vector<double>& elementRHS);
    
    // 旋度计算辅助函数
    bool computeElementCurl(int elementId,
                           const std::vector<double>& nodeCoordsX,
                           const std::vector<double>& nodeCoordsY,
                           const std::vector<double>& nodeCoordsZ,
                           const std::vector<double>& fieldX,
                           const std::vector<double>& fieldY,
                           const std::vector<double>& fieldZ,
                           std::vector<double>& curlX,
                           std::vector<double>& curlY,
                           std::vector<double>& curlZ);
};

} // namespace elmer

} // namespace elmer