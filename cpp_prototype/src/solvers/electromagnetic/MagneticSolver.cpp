/**
 * @file MagneticSolver.cpp
 * @brief 磁场求解器实现 - 移植自Fortran版本的MagneticSolve.F90
 * 
 * 实现MHD Maxwell方程（或感应方程）的求解器
 * 
 * TODO: 需要后续进一步开发的功能
 * - [ ] 实现完整的Maxwell方程求解
 * - [ ] 实现坐标系转换（笛卡尔、柱对称、一般坐标系）
 * - [ ] 实现瞬态仿真支持
 * - [ ] 实现非线性迭代求解
 * - [ ] 实现边界条件处理
 * - [ ] 实现洛伦兹力计算
 * - [ ] 实现电场计算
 * - [ ] 实现材料参数处理
 */

#include "MagneticSolver.h"
#include "Types.h"
#include "LoggerFactory.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>

namespace elmer {

// ===== 构造函数和析构函数 =====

MagneticSolver::MagneticSolver() {
    // 基于Fortran版本的MagneticSolve.F90初始化参数
    ELMER_DEBUG("初始化MagneticSolver...");
    
    // 设置默认参数（基于Fortran代码中的默认值）
    maxNonlinearIterations_ = 10;        // NonlinearIter默认值
    nonlinearTolerance_ = 1e-6;           // Nonlinear System Convergence Tolerance
    stabilize_ = true;                    // Stabilize默认值
    transientSimulation_ = false;         // TransientSimulation默认值
    timeStep_ = 0.0;                      // dt默认值
    
    // 初始化状态标志
    allocationsDone_ = false;
    userDefinedVelo_ = false;
    calculateMagneticForce_ = false;
    
    ELMER_DEBUG("MagneticSolver构造函数完成");
}

MagneticSolver::~MagneticSolver() {
    ELMER_DEBUG("开始清理MagneticSolver资源...");
    
    // 释放内存（基于Fortran的DEALLOCATE）
    deallocateMemory();
    
    ELMER_DEBUG("MagneticSolver析构函数完成");
}

// ===== 基本接口函数 =====

std::string MagneticSolver::getName() const {
    return "MagneticSolver";
}

bool MagneticSolver::initialize() {
    // TODO: 实现磁场求解器的完整初始化流程
    ELMER_INFO("开始初始化MagneticSolver...");
    
    // 1. 检查必要的组件是否已设置
    if (!mesh_ || !materialDB_ || !bc_) {
        ELMER_ERROR("求解器未正确初始化，缺少必要的组件");
        return false;
    }
    
    // 2. 分配内存
    if (!allocateMemory()) {
        ELMER_ERROR("内存分配失败");
        return false;
    }
    
    // 3. 获取材料参数
    if (!getMaterialParameters()) {
        ELMER_ERROR("材料参数获取失败");
        return false;
    }
    
    // 4. 获取速度场（如果存在）
    if (!getVelocityField()) {
        ELMER_WARN("速度场获取失败，将使用默认值");
    }
    
    // 5. 设置求解器状态
    status_ = SolverStatus::INITIALIZED;
    
    ELMER_INFO("MagneticSolver初始化完成");
    return true;
}

bool MagneticSolver::assemble() {
    // TODO: 实现磁场求解器的系统矩阵组装
    ELMER_INFO("开始组装磁场求解器系统矩阵...");
    
    if (status_ != SolverStatus::INITIALIZED) {
        ELMER_ERROR("求解器未初始化，无法组装系统矩阵");
        return false;
    }
    
    if (!globalStiffnessMatrix_ || !solutionVector_) {
        ELMER_ERROR("刚度矩阵或右端向量未初始化");
        return false;
    }
    
    // 清零刚度矩阵和右端向量
    globalStiffnessMatrix_->Zero();
    solutionVector_->Zero();
    
    // 根据坐标系类型选择组装方法
    bool success = false;
    
    // 检测坐标系类型
    std::string coordinateSystem = detectCoordinateSystem();
    ELMER_INFO("检测到坐标系类型: {}", coordinateSystem);
    
    // 根据坐标系类型选择相应的组装方法
    if (coordinateSystem == "cartesian") {
        success = assembleCartesian();
    } else if (coordinateSystem == "axisymmetric") {
        success = assembleAxisymmetric();
    } else if (coordinateSystem == "general") {
        success = assembleGeneral();
    } else {
        ELMER_ERROR("未知的坐标系类型: {}", coordinateSystem);
        return false;
    }
    
    if (!success) {
        ELMER_ERROR("系统矩阵组装失败");
        return false;
    }
    
    // 组装边界条件
    if (!assembleBoundaryConditions()) {
        ELMER_ERROR("边界条件组装失败");
        return false;
    }
    
    status_ = SolverStatus::ASSEMBLED;
    ELMER_INFO("磁场求解器系统矩阵组装完成");
    return true;
}

bool MagneticSolver::solve() {
    // 基于Fortran版本的MagneticSolve.F90实现非线性迭代求解
    ELMER_INFO("开始求解Maxwell方程...");
    
    if (status_ != SolverStatus::ASSEMBLED) {
        ELMER_ERROR("求解器未组装，无法求解");
        return false;
    }
    
    if (!globalStiffnessMatrix_ || !solutionVector_) {
        ELMER_ERROR("系统矩阵或右端向量未初始化");
        return false;
    }
    
    // 初始化解向量
    if (!solutionVector_) {
        solutionVector_ = Vector::Create(stiffnessMatrix_->GetNumRows());
        solutionVector_->Zero();
    }
    
    // 非线性迭代求解（牛顿迭代法）
    double prevNorm = 0.0;
    double currentNorm = 0.0;
    bool converged = false;
    
    for (int iter = 1; iter <= maxNonlinearIterations_; ++iter) {
        ELMER_DEBUG("非线性迭代: {}/{}", iter, maxNonlinearIterations_);
        
        // 1. 计算残差向量
        auto residual = computeResidualVector();
        if (!residual) {
            ELMER_ERROR("残差向量计算失败");
            return false;
        }
        
        // 2. 计算残差范数
        currentNorm = computeVectorNorm(*residual);
        ELMER_DEBUG("残差范数: {}", currentNorm);
        
        // 3. 检查收敛性
        if (checkConvergence(prevNorm, currentNorm)) {
            ELMER_INFO("非线性迭代收敛于第 {} 次迭代", iter);
            converged = true;
            break;
        }
        
        // 4. 更新雅可比矩阵（简化实现：使用当前刚度矩阵）
        if (!updateJacobianMatrix()) {
            ELMER_WARN("雅可比矩阵更新失败，使用当前刚度矩阵");
        }
        
        // 5. 求解线性系统：J * Δx = -r
        auto deltaX = solveLinearSystem(*residual);
        if (!deltaX) {
            ELMER_ERROR("线性系统求解失败");
            return false;
        }
        
        // 6. 更新解向量：x = x + Δx
        updateSolutionVector(*deltaX);
        
        prevNorm = currentNorm;
        
        if (iter == maxNonlinearIterations_) {
            ELMER_WARN("达到最大非线性迭代次数");
        }
    }
    
    if (!converged) {
        ELMER_WARN("非线性迭代未收敛");
    }
    
    // 更新磁场变量
    updateMagneticFieldFromSolution();
    
    // 计算相关物理量
    if (!computeLorentzForce()) {
        ELMER_WARN("洛伦兹力计算失败");
    }
    
    if (!computeElectricField()) {
        ELMER_WARN("电场计算失败");
    }
    
    status_ = SolverStatus::SOLVED;
    ELMER_INFO("Maxwell方程求解完成");
    return true;
}

std::vector<double> MagneticSolver::getSolution() const {
    // 返回解向量的数据
    if (solutionVector_ && solutionVector_->Size() > 0) {
        std::vector<double> result(solutionVector_->Size());
        for (int i = 0; i < solutionVector_->Size(); ++i) {
            result[i] = (*solutionVector_)[i];
        }
        return result;
    }
    return {};
}

// ===== 磁场求解器特定功能 =====

void MagneticSolver::setTransientSimulation(bool transient) {
    transientSimulation_ = transient;
    std::cout << "设置瞬态仿真标志: " << (transient ? "是" : "否") << std::endl;
}

void MagneticSolver::setTimeStep(double dt) {
    // TODO: 实现时间步长设置
    std::cout << "设置时间步长: " << dt << std::endl;
}

bool MagneticSolver::computeLorentzForce() {
    // TODO: 实现洛伦兹力计算
    std::cout << "开始计算洛伦兹力..." << std::endl;
    
    // 简化实现：使用J × B公式
    size_t numNodes = mesh_->numberOfNodes();
    if (electricCurrent_.size() >= 3 * numNodes && magneticField_.size() >= 3 * numNodes) {
        for (int i = 0; i < numNodes; ++i) {
            // 计算J × B
            // TODO: 实现完整的向量叉乘计算
            lorentzForce_[3 * i] = 0.0;     // Fx
            lorentzForce_[3 * i + 1] = 0.0; // Fy
            lorentzForce_[3 * i + 2] = 0.0; // Fz
        }
    }
    
    std::cout << "洛伦兹力计算完成" << std::endl;
    return true;
}

bool MagneticSolver::computeElectricField() {
    // TODO: 实现电场计算
    std::cout << "开始计算电场..." << std::endl;
    
    // 简化实现：使用E = J / σ公式
    int numNodes = mesh_->numberOfNodes();
    if (electricCurrent_.size() >= 3 * numNodes && conductivity_.size() >= numNodes) {
        for (int i = 0; i < numNodes; ++i) {
            double sigma = conductivity_[i];
            if (sigma > 0.0) {
                electricField_[3 * i] = electricCurrent_[3 * i] / sigma;     // Ex
                electricField_[3 * i + 1] = electricCurrent_[3 * i + 1] / sigma; // Ey
                electricField_[3 * i + 2] = electricCurrent_[3 * i + 2] / sigma; // Ez
            }
        }
    }
    
    std::cout << "电场计算完成" << std::endl;
    return true;
}

std::vector<double> MagneticSolver::getMagneticField() const {
    return magneticField_;
}

std::vector<double> MagneticSolver::getElectricCurrent() const {
    return electricCurrent_;
}

std::vector<double> MagneticSolver::getLorentzForce() const {
    return lorentzForce_;
}

std::vector<double> MagneticSolver::getElectricField() const {
    return electricField_;
}

// ===== 私有辅助函数 =====

bool MagneticSolver::allocateMemory() {
    // TODO: 实现完整的内存分配
    std::cout << "开始分配内存..." << std::endl;
    
    if (!mesh_) {
        std::cerr << "错误: 网格未设置，无法分配内存" << std::endl;
        return false;
    }
    
    int numNodes = mesh_->numberOfNodes();
    
    // 分配场变量内存
    magneticField_.resize(3 * numNodes, 0.0);
    electricCurrent_.resize(3 * numNodes, 0.0);
    lorentzForce_.resize(3 * numNodes, 0.0);
    electricField_.resize(3 * numNodes, 0.0);
    
    // 分配速度场内存
    U_.resize(numNodes, 0.0);
    V_.resize(numNodes, 0.0);
    W_.resize(numNodes, 0.0);
    
    // 分配材料参数内存
    conductivity_.resize(numNodes, 0.0);
    permeability_.resize(numNodes, 0.0);
    
    // 分配施加的磁场内存
    appliedMagneticFieldX_.resize(numNodes, 0.0);
    appliedMagneticFieldY_.resize(numNodes, 0.0);
    appliedMagneticFieldZ_.resize(numNodes, 0.0);
    
    allocationsDone_ = true;
    std::cout << "内存分配完成，节点数: " << numNodes << std::endl;
    return true;
}

bool MagneticSolver::getMaterialParameters() {
    // TODO: 实现从材料数据库获取参数
    std::cout << "开始获取材料参数..." << std::endl;
    
    if (!materialDB_) {
        std::cerr << "错误: 材料数据库未设置" << std::endl;
        return false;
    }
    
    int numNodes = mesh_->numberOfNodes();
    
    // 简化实现：设置默认材料参数
    for (int i = 0; i < numNodes; ++i) {
        conductivity_[i] = 1.0e6;  // 默认电导率：铜的电导率
        permeability_[i] = 1.0;    // 默认磁导率：真空磁导率
    }
    
    std::cout << "材料参数获取完成" << std::endl;
    return true;
}

bool MagneticSolver::getVelocityField() {
    // TODO: 实现速度场获取
    std::cout << "开始获取速度场..." << std::endl;
    
    // 简化实现：设置默认速度场（静止）
    int numNodes = mesh_->numberOfNodes();
    for (int i = 0; i < numNodes; ++i) {
        U_[i] = 0.0;
        V_[i] = 0.0;
        W_[i] = 0.0;
    }
    
    std::cout << "速度场获取完成" << std::endl;
    return true;
}

bool MagneticSolver::assembleElementMatrix(int elementId) {
    // 基于Fortran版本的MagneticSolve.F90中的矩阵组装逻辑实现
    // 实现质量矩阵和刚度矩阵的组装
    
    if (!mesh_ || elementId < 0 || elementId >= mesh_->getBulkElements().size()) {
        ELMER_ERROR("无效的单元索引: {}", elementId);
        return false;
    }
    
    auto& bulkElements = mesh_->getBulkElements();
    auto& element = bulkElements[elementId];
    auto elementNodes = element.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 获取单元节点坐标
    std::vector<double> nodeCoordsX(numElementNodes);
    std::vector<double> nodeCoordsY(numElementNodes);
    std::vector<double> nodeCoordsZ(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        auto node = mesh_->getNodes()[nodeId];
        nodeCoordsX[i] = node.x;
        nodeCoordsY[i] = node.y;
        nodeCoordsZ[i] = node.z;
    }
    
    // 获取单元材料参数
    std::vector<double> elemConductivity(numElementNodes);
    std::vector<double> elemPermeability(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        if (nodeId < conductivity_.size() && nodeId < permeability_.size()) {
            elemConductivity[i] = conductivity_[nodeId];
            elemPermeability[i] = permeability_[nodeId];
        }
    }
    
    // 初始化单元刚度矩阵和质量矩阵
    std::vector<std::vector<double>> elementStiffness(numElementNodes * 3, 
                                                     std::vector<double>(numElementNodes * 3, 0.0));
    std::vector<std::vector<double>> elementMass(numElementNodes * 3, 
                                                std::vector<double>(numElementNodes * 3, 0.0));
    
    // 高斯积分点（根据单元类型选择）
    int numGaussPoints = 8; // 对于3D单元使用8个高斯点
    
    for (int gaussPoint = 0; gaussPoint < numGaussPoints; ++gaussPoint) {
        // 计算高斯点坐标和权重
        double xi, eta, zeta, weight;
        getGaussPoint3D(gaussPoint, xi, eta, zeta, weight);
        
        // 计算形状函数和导数
        std::vector<double> shapeFunctions(numElementNodes);
        std::vector<double> dShapeDXi(numElementNodes);
        std::vector<double> dShapeDEta(numElementNodes);
        std::vector<double> dShapeDZeta(numElementNodes);
        
        computeShapeFunctions3D(xi, eta, zeta, shapeFunctions, dShapeDXi, dShapeDEta, dShapeDZeta);
        
        // 计算雅可比矩阵
        double jacobian[3][3] = {{0}};
        for (int i = 0; i < numElementNodes; ++i) {
            jacobian[0][0] += dShapeDXi[i] * nodeCoordsX[i];
            jacobian[0][1] += dShapeDXi[i] * nodeCoordsY[i];
            jacobian[0][2] += dShapeDXi[i] * nodeCoordsZ[i];
            jacobian[1][0] += dShapeDEta[i] * nodeCoordsX[i];
            jacobian[1][1] += dShapeDEta[i] * nodeCoordsY[i];
            jacobian[1][2] += dShapeDEta[i] * nodeCoordsZ[i];
            jacobian[2][0] += dShapeDZeta[i] * nodeCoordsX[i];
            jacobian[2][1] += dShapeDZeta[i] * nodeCoordsY[i];
            jacobian[2][2] += dShapeDZeta[i] * nodeCoordsZ[i];
        }
        
        // 计算雅可比行列式
        double detJ = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
                    - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
                    + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);
        
        if (detJ <= 0.0) {
            ELMER_ERROR("单元 {} 的雅可比行列式为负或零", elementId);
            return false;
        }
        
        // 计算形状函数在全局坐标系中的导数
        std::vector<double> dShapeDX(numElementNodes);
        std::vector<double> dShapeDY(numElementNodes);
        std::vector<double> dShapeDZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            // 使用雅可比矩阵的逆计算全局导数
            dShapeDX[i] = (1.0/detJ) * (
                (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) * dShapeDXi[i] +
                (jacobian[0][2] * jacobian[2][1] - jacobian[0][1] * jacobian[2][2]) * dShapeDEta[i] +
                (jacobian[0][1] * jacobian[1][2] - jacobian[0][2] * jacobian[1][1]) * dShapeDZeta[i]
            );
            
            dShapeDY[i] = (1.0/detJ) * (
                (jacobian[1][2] * jacobian[2][0] - jacobian[1][0] * jacobian[2][2]) * dShapeDXi[i] +
                (jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0]) * dShapeDEta[i] +
                (jacobian[0][2] * jacobian[1][0] - jacobian[0][0] * jacobian[1][2]) * dShapeDZeta[i]
            );
            
            dShapeDZ[i] = (1.0/detJ) * (
                (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]) * dShapeDXi[i] +
                (jacobian[0][1] * jacobian[2][0] - jacobian[0][0] * jacobian[2][1]) * dShapeDEta[i] +
                (jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0]) * dShapeDZeta[i]
            );
        }
        
        // 基于Fortran版本的矩阵组装逻辑
        // 质量矩阵组装: M(i,i) = M(i,i) + s * Basis(q) * Basis(p)
        // 刚度矩阵组装: A(i,i) = A(i,i) + s * dBasisdx(q,j)*dBasisdx(p,j)/Conductivity
        
        for (int p = 0; p < numElementNodes; ++p) {
            for (int q = 0; q < numElementNodes; ++q) {
                // 计算材料参数在高斯点的平均值
                double avgConductivity = 0.5 * (elemConductivity[p] + elemConductivity[q]);
                double avgPermeability = 0.5 * (elemPermeability[p] + elemPermeability[q]);
                
                // 积分权重
                double s = weight * detJ;
                
                // 组装质量矩阵（对角线项）
                for (int i = 0; i < 3; ++i) {
                    int row = p * 3 + i;
                    int col = q * 3 + i;
                    elementMass[row][col] += s * shapeFunctions[p] * shapeFunctions[q];
                }
                
                // 组装刚度矩阵（扩散项）
                for (int i = 0; i < 3; ++i) {
                    int row = p * 3 + i;
                    int col = q * 3 + i;
                    
                    // 扩散项：A(i,i) = A(i,i) + s * dBasisdx(q,j)*dBasisdx(p,j)/Conductivity
                    double diffusiveTerm = s * (
                        dShapeDX[p] * dShapeDX[q] + 
                        dShapeDY[p] * dShapeDY[q] + 
                        dShapeDZ[p] * dShapeDZ[q]
                    ) / avgConductivity;
                    
                    elementStiffness[row][col] += diffusiveTerm;
                }
                
                // 组装刚度矩阵（涡流项）
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        int row = p * 3 + i;
                        int col = q * 3 + j;
                        
                        // 涡流项：A(i,j) = A(i,j) - s * Basis(q) * dVelodx(i,j) * Basis(p)
                        // 这里简化处理，使用单位速度场
                        double convectionTerm = -s * shapeFunctions[p] * 1.0 * shapeFunctions[q];
                        
                        elementStiffness[row][col] += convectionTerm;
                    }
                }
            }
        }
    }
    
    // 将组装好的矩阵存储到全局矩阵中
    if (!globalStiffnessMatrix_) {
        ELMER_ERROR("全局刚度矩阵未初始化");
        return false;
    }
    
    if (!globalMassMatrix_) {
        ELMER_ERROR("全局质量矩阵未初始化");
        return false;
    }
    
    // 将单元矩阵组装到全局矩阵中
    for (int p = 0; p < numElementNodes; ++p) {
        int globalNodeP = static_cast<int>(elementNodes[p]);
        
        for (int q = 0; q < numElementNodes; ++q) {
            int globalNodeQ = static_cast<int>(elementNodes[q]);
            
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    int localRow = p * 3 + i;
                    int localCol = q * 3 + j;
                    int globalRow = globalNodeP * 3 + i;
                    int globalCol = globalNodeQ * 3 + j;
                    
                    // 添加到全局矩阵
                    globalStiffnessMatrix_->AddToElement(globalRow, globalCol, elementStiffness[localRow][localCol]);
                    globalMassMatrix_->AddToElement(globalRow, globalCol, elementMass[localRow][localCol]);
                }
            }
        }
    }
    
    ELMER_DEBUG("单元 {} 矩阵组装完成，节点数: {}", elementId, numElementNodes);
    return true;
}

bool MagneticSolver::assembleBoundaryConditions() {
    // 实现边界条件组装
    // 基于Fortran版本的MagneticSolve.F90中的边界条件处理
    ELMER_INFO("开始组装边界条件...");
    
    if (!bc_ || !mesh_) {
        ELMER_WARN("边界条件管理器或网格未设置，跳过边界条件组装");
        return true;
    }
    
    // 获取边界单元数量
    int numBoundaryElements = mesh_->getBoundaryElements().size();
    
    if (numBoundaryElements == 0) {
        ELMER_INFO("没有边界单元，跳过边界条件组装");
        return true;
    }
    
    // 遍历所有边界单元
    int processedBC = 0;
    auto& boundaryElements = mesh_->getBoundaryElements();
    
    for (int bcId = 0; bcId < numBoundaryElements; ++bcId) {
        if (bcId >= boundaryElements.size()) {
            ELMER_ERROR("边界单元索引超出范围: {} >= {}", bcId, boundaryElements.size());
            return false;
        }
        
        auto& boundaryElement = boundaryElements[bcId];
        auto elementNodes = boundaryElement.getNodeIndices();
        size_t numElementNodes = elementNodes.size();
        
        // 检测边界条件类型
        std::string bcType = detectBoundaryConditionType(bcId, boundaryElement);
        ELMER_DEBUG("边界条件 {} 类型: {}", bcId, bcType);
        
        bool success = false;
        
        // 根据边界条件类型应用相应的处理
        if (bcType == "dirichlet") {
            success = applyDirichletBoundaryCondition(bcId, boundaryElement);
        } else if (bcType == "neumann") {
            success = applyNeumannBoundaryCondition(bcId, boundaryElement);
        } else if (bcType == "magnetic_force") {
            success = applyMagneticForceBoundaryCondition(bcId, boundaryElement);
        } else if (bcType == "periodic") {
            success = applyPeriodicBoundaryCondition(bcId, boundaryElement);
        } else {
            ELMER_WARN("未知边界条件类型: {}，跳过处理", bcType);
            success = true; // 跳过未知类型
        }
        
        if (!success) {
            ELMER_ERROR("边界条件 {} 应用失败", bcId);
            return false;
        }
        
        processedBC++;
    }
    
    ELMER_INFO("边界条件组装完成，共处理 {} 个边界条件", processedBC);
    return true;
}

std::string MagneticSolver::detectBoundaryConditionType(int bcId, const Element& boundaryElement) {
    // 基于Fortran版本的边界条件检测逻辑
    // 检测边界条件的类型（狄利克雷、诺伊曼、磁力、周期性等）
    
    if (!bc_ || !mesh_) {
        ELMER_WARN("边界条件管理器或网格未设置，使用默认狄利克雷边界条件");
        return "dirichlet";
    }
    
    // 获取边界条件参数
    auto elementNodes = boundaryElement.getNodeIndices();
    
    // 检查是否存在狄利克雷边界条件（固定值）
    bool hasDirichlet = false;
    bool hasNeumann = false;
    bool hasMagneticForce = false;
    bool hasPeriodic = false;
    
    // 简化实现：基于边界条件管理器的检测
    // 实际应该根据边界条件的具体参数进行检测
    
    // 临时简化实现：基于边界条件ID的简单检测
    // TODO: 实现完整的边界条件类型检测逻辑
    if (bcId % 4 == 0) {
        hasDirichlet = true;
    } else if (bcId % 4 == 1) {
        hasNeumann = true;
    } else if (bcId % 4 == 2) {
        hasMagneticForce = true;
    } else if (bcId % 4 == 3) {
        hasPeriodic = true;
    }
    
    // 优先级：狄利克雷 > 诺伊曼 > 磁力 > 周期性 > 默认狄利克雷
    if (hasDirichlet) {
        return "dirichlet";
    } else if (hasNeumann) {
        return "neumann";
    } else if (hasMagneticForce) {
        return "magnetic_force";
    } else if (hasPeriodic) {
        return "periodic";
    }
    
    // 默认返回狄利克雷边界条件
    return "dirichlet";
}

bool MagneticSolver::applyPeriodicBoundaryCondition(int bcId, const Element& boundaryElement) {
    // 实现周期性边界条件
    // 基于Fortran版本的周期性边界条件处理
    
    if (!bc_ || !globalStiffnessMatrix_ || !solutionVector_) {
        ELMER_ERROR("必要的组件未初始化");
        return false;
    }
    
    auto elementNodes = boundaryElement.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 获取主节点和从节点（周期性边界条件对）
    std::vector<int> masterNodes;
    std::vector<int> slaveNodes;
    
    // 简化实现：假设前一半节点是主节点，后一半是从节点
    for (int i = 0; i < numElementNodes / 2; ++i) {
        masterNodes.push_back(static_cast<int>(elementNodes[i]));
    }
    for (int i = numElementNodes / 2; i < numElementNodes; ++i) {
        slaveNodes.push_back(static_cast<int>(elementNodes[i]));
    }
    
    // 应用周期性边界条件：从节点 = 主节点
    for (size_t i = 0; i < masterNodes.size() && i < slaveNodes.size(); ++i) {
        int masterNode = masterNodes[i];
        int slaveNode = slaveNodes[i];
        
        // 对每个自由度应用周期性条件
        for (int dof = 0; dof < 3; ++dof) {
            int masterDof = masterNode * 3 + dof;
            int slaveDof = slaveNode * 3 + dof;
            
            // 设置从节点的解等于主节点的解
            // 这需要在求解过程中特殊处理
            
            // 简化实现：在刚度矩阵中设置约束
            if (masterDof < stiffnessMatrix_->GetNumRows() && 
                slaveDof < stiffnessMatrix_->GetNumRows()) {
                
                // 将主节点的贡献加到从节点
                for (int col = 0; col < stiffnessMatrix_->GetNumCols(); ++col) {
                    double masterValue = stiffnessMatrix_->GetElement(masterDof, col);
                    double slaveValue = stiffnessMatrix_->GetElement(slaveDof, col);
                    
                    // 合并贡献
                    stiffnessMatrix_->SetElement(slaveDof, col, masterValue + slaveValue);
                }
                
                // 清零从节点的行（除了对角线）
                for (int col = 0; col < stiffnessMatrix_->GetNumCols(); ++col) {
                    if (col != slaveDof) {
                        stiffnessMatrix_->SetElement(slaveDof, col, 0.0);
                    }
                }
                
                // 设置从节点对角线为1
                stiffnessMatrix_->SetElement(slaveDof, slaveDof, 1.0);
                
                // 设置右端向量
                (*rhsVector_)[slaveDof] = (*rhsVector_)[masterDof];
            }
        }
    }
    return true;
}

bool MagneticSolver::checkConvergence(double prevNorm, double currentNorm) const {
    // 基于Fortran版本的收敛性检查实现
    // 实现多种收敛准则：相对残差、绝对残差、增量范数
    
    if (prevNorm == 0.0) {
        ELMER_DEBUG("第一次迭代，不检查收敛");
        return false; // 第一次迭代，不检查收敛
    }
    
    // 计算相对残差变化
    double relativeChange = 2.0 * std::abs(prevNorm - currentNorm) / (prevNorm + currentNorm);
    
    // 计算绝对残差
    double absoluteResidual = currentNorm;
    
    // 计算相对残差（相对于初始残差）
    static double initialNorm = prevNorm;
    double relativeResidual = currentNorm / initialNorm;
    
    ELMER_DEBUG("收敛检查: 相对变化={}, 绝对残差={}, 相对残差={}, 容差={}", 
                relativeChange, absoluteResidual, relativeResidual, nonlinearTolerance_);
    
    // 多种收敛准则（满足任一条件即可）
    bool converged = false;
    
    // 1. 相对残差变化准则
    if (relativeChange < nonlinearTolerance_) {
        ELMER_DEBUG("满足相对变化收敛准则");
        converged = true;
    }
    
    // 2. 绝对残差准则
    if (absoluteResidual < 1e-12) {
        ELMER_DEBUG("满足绝对残差收敛准则");
        converged = true;
    }
    
    // 3. 相对残差准则
    if (relativeResidual < nonlinearTolerance_ * 0.1) {
        ELMER_DEBUG("满足相对残差收敛准则");
        converged = true;
    }
    
    // 4. 检查残差是否发散
    if (currentNorm > 1e6 * prevNorm && prevNorm > 1e-12) {
        ELMER_WARN("残差发散: currentNorm={}, prevNorm={}", currentNorm, prevNorm);
        // 不设置收敛，但也不报错，让迭代继续
    }
    
    return converged;
}

bool MagneticSolver::assembleCartesian() {
    // 实现笛卡尔坐标系下的Maxwell方程组装
    // 基于Fortran版本的MagneticSolve.F90中的MaxwellCompose函数
    ELMER_INFO("使用笛卡尔坐标系组装Maxwell方程...");
    
    if (!mesh_ || !stiffnessMatrix_ || !rhsVector_) {
        ELMER_ERROR("必要的组件未初始化");
        return false;
    }
    
    int numElements = mesh_->getBulkElements().size();
    int numNodes = mesh_->numberOfNodes();
    
    // 检查材料参数是否已获取
    if (conductivity_.empty() || permeability_.empty()) {
        ELMER_ERROR("材料参数未获取");
        return false;
    }
    
    // 遍历所有单元进行组装
    for (int elemId = 0; elemId < numElements; ++elemId) {
        // 获取单元信息
        auto& bulkElements = mesh_->getBulkElements();
        if (elemId >= bulkElements.size()) {
            ELMER_ERROR("单元索引超出范围");
            return false;
        }
        
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        int numElementNodes = elementNodes.size();
        
        // 计算单元矩阵和右端向量
        std::vector<std::vector<double>> elementStiffness;
        std::vector<double> elementRHS;
        
        if (!computeElementMatrices(elemId, elementStiffness, elementRHS)) {
            ELMER_ERROR("单元 {} 矩阵计算失败", elemId);
            return false;
        }
        
        // 组装到全局系统
        if (!assembleToGlobalSystem(elemId, elementStiffness, elementRHS)) {
            ELMER_ERROR("单元 {} 组装到全局系统失败", elemId);
            return false;
        }
    }
    
    ELMER_INFO("笛卡尔坐标系组装完成，共处理 {} 个单元", numElements);
    return true;
}

bool MagneticSolver::assembleAxisymmetric() {
    // 基于Fortran版本的MagneticSolve.F90实现柱对称坐标系组装
    ELMER_INFO("使用柱对称坐标系组装Maxwell方程...");
    
    if (!mesh_ || !stiffnessMatrix_ || !rhsVector_) {
        ELMER_ERROR("网格、刚度矩阵或右端向量未初始化");
        return false;
    }
    
    int numElements = mesh_->getBulkElements().size();
    int numNodes = mesh_->numberOfNodes();
    
    // 检查材料参数是否已获取
    if (conductivity_.empty() || permeability_.empty()) {
        ELMER_ERROR("材料参数未获取");
        return false;
    }
    
    // 清零刚度矩阵和右端向量
    stiffnessMatrix_->Zero();
    rhsVector_->Zero();
    
    // 遍历所有单元进行组装
    for (int elemId = 0; elemId < numElements; ++elemId) {
        // 获取单元信息
        auto& bulkElements = mesh_->getBulkElements();
        if (elemId >= bulkElements.size()) {
            ELMER_ERROR("单元索引超出范围");
            return false;
        }
        
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        int numElementNodes = elementNodes.size();
        
        // 计算单元矩阵和右端向量（柱对称坐标系）
        std::vector<std::vector<double>> elementStiffness;
        std::vector<double> elementRHS;
        
        if (!computeAxisymmetricElementMatrices(elemId, elementStiffness, elementRHS)) {
            ELMER_ERROR("单元 {} 柱对称矩阵计算失败", elemId);
            return false;
        }
        
        // 组装到全局系统
        if (!assembleToGlobalSystem(elemId, elementStiffness, elementRHS)) {
            ELMER_ERROR("单元 {} 组装到全局系统失败", elemId);
            return false;
        }
    }
    
    // 组装边界条件
    if (!assembleBoundaryConditions()) {
        ELMER_ERROR("边界条件组装失败");
        return false;
    }
    
    ELMER_INFO("柱对称坐标系组装完成，共处理 {} 个单元", numElements);
    return true;
}

bool MagneticSolver::assembleGeneral() {
    // 基于Fortran版本的MagneticSolve.F90实现一般坐标系组装
    ELMER_INFO("使用一般坐标系组装Maxwell方程...");
    
    if (!mesh_ || !stiffnessMatrix_ || !rhsVector_) {
        ELMER_ERROR("网格、刚度矩阵或右端向量未初始化");
        return false;
    }
    
    int numElements = mesh_->getBulkElements().size();
    int numNodes = mesh_->numberOfNodes();
    
    // 检查材料参数是否已获取
    if (conductivity_.empty() || permeability_.empty()) {
        ELMER_ERROR("材料参数未获取");
        return false;
    }
    
    // 清零刚度矩阵和右端向量
    stiffnessMatrix_->Zero();
    rhsVector_->Zero();
    
    // 遍历所有单元进行组装
    for (int elemId = 0; elemId < numElements; ++elemId) {
        // 获取单元信息
        auto& bulkElements = mesh_->getBulkElements();
        if (elemId >= bulkElements.size()) {
            ELMER_ERROR("单元索引超出范围");
            return false;
        }
        
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        int numElementNodes = elementNodes.size();
        
        // 计算单元矩阵和右端向量（一般坐标系）
        std::vector<std::vector<double>> elementStiffness;
        std::vector<double> elementRHS;
        
        if (!computeGeneralElementMatrices(elemId, elementStiffness, elementRHS)) {
            ELMER_ERROR("单元 {} 一般坐标系矩阵计算失败", elemId);
            return false;
        }
        
        // 组装到全局系统
        if (!assembleToGlobalSystem(elemId, elementStiffness, elementRHS)) {
            ELMER_ERROR("单元 {} 组装到全局系统失败", elemId);
            return false;
        }
    }
    
    // 组装边界条件
    if (!assembleBoundaryConditions()) {
        ELMER_ERROR("边界条件组装失败");
        return false;
    }
    
    ELMER_INFO("一般坐标系组装完成，共处理 {} 个单元", numElements);
    return true;
}

bool MagneticSolver::computeAxisymmetricElementMatrices(int elementId, 
                                                          std::vector<std::vector<double>>& elementStiffness,
                                                          std::vector<double>& elementRHS) {
    // 基于Fortran版本的MagneticSolve.F90实现柱对称坐标系单元矩阵计算
    
    if (!mesh_ || elementId < 0 || elementId >= mesh_->getBulkElements().size()) {
        ELMER_ERROR("无效的单元索引: {}", elementId);
        return false;
    }
    
    auto& bulkElements = mesh_->getBulkElements();
    auto& element = bulkElements[elementId];
    auto elementNodes = element.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 初始化单元矩阵和右端向量
    elementStiffness.resize(numElementNodes * 3, std::vector<double>(numElementNodes * 3, 0.0));
    elementRHS.resize(numElementNodes * 3, 0.0);
    
    // 获取单元材料参数
    std::vector<double> elemConductivity(numElementNodes);
    std::vector<double> elemPermeability(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        if (nodeId < conductivity_.size() && nodeId < permeability_.size()) {
            elemConductivity[i] = conductivity_[nodeId];
            elemPermeability[i] = permeability_[nodeId];
        }
    }
    
    // 获取单元节点坐标
    std::vector<double> nodeCoordsX(numElementNodes);
    std::vector<double> nodeCoordsY(numElementNodes);
    std::vector<double> nodeCoordsZ(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        auto node = mesh_->getNodes()[nodeId];
        nodeCoordsX[i] = node.x;
        nodeCoordsY[i] = node.y;
        nodeCoordsZ[i] = node.z;
    }
    
    // 柱对称坐标系下的特殊处理：径向坐标r = x，轴向坐标z = y
    // 使用2D高斯积分，考虑轴对称特性
    int numGaussPoints = 4; // 对于2D单元使用4个高斯点
    
    for (int gaussPoint = 0; gaussPoint < numGaussPoints; ++gaussPoint) {
        // 计算高斯点坐标和权重
        double xi, eta, weight;
        getGaussPoint2D(gaussPoint, xi, eta, weight);
        
        // 计算形状函数和导数
        std::vector<double> shapeFunctions(numElementNodes);
        std::vector<double> dShapeDXi(numElementNodes);
        std::vector<double> dShapeDEta(numElementNodes);
        
        computeShapeFunctions2D(xi, eta, shapeFunctions, dShapeDXi, dShapeDEta);
        
        // 计算雅可比矩阵（柱对称坐标系）
        double jacobian[2][2] = {{0}};
        double r = 0.0; // 径向坐标
        
        for (int i = 0; i < numElementNodes; ++i) {
            jacobian[0][0] += dShapeDXi[i] * nodeCoordsX[i]; // dr/dξ
            jacobian[0][1] += dShapeDXi[i] * nodeCoordsY[i]; // dz/dξ
            jacobian[1][0] += dShapeDEta[i] * nodeCoordsX[i]; // dr/dη
            jacobian[1][1] += dShapeDEta[i] * nodeCoordsY[i]; // dz/dη
            r += shapeFunctions[i] * nodeCoordsX[i]; // 径向坐标
        }
        
        // 计算雅可比行列式
        double detJ = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
        
        if (detJ <= 0.0 || r <= 0.0) {
            ELMER_ERROR("单元 {} 的雅可比行列式为负或零，或径向坐标无效", elementId);
            return false;
        }
        
        // 计算形状函数在全局坐标系中的导数
        std::vector<double> dShapeDR(numElementNodes);
        std::vector<double> dShapeDZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            dShapeDR[i] = (1.0/detJ) * (jacobian[1][1] * dShapeDXi[i] - jacobian[0][1] * dShapeDEta[i]);
            dShapeDZ[i] = (1.0/detJ) * (-jacobian[1][0] * dShapeDXi[i] + jacobian[0][0] * dShapeDEta[i]);
        }
        
        // 柱对称坐标系下的积分权重因子
        double integrationFactor = weight * detJ * 2.0 * 3.14159265358979323846 * r;
        
        // 计算单元刚度矩阵（考虑轴对称特性）
        for (int i = 0; i < numElementNodes; ++i) {
            for (int j = 0; j < numElementNodes; ++j) {
                // 计算材料参数在高斯点的平均值
                double avgConductivity = 0.5 * (elemConductivity[i] + elemConductivity[j]);
                double avgPermeability = 0.5 * (elemPermeability[i] + elemPermeability[j]);
                
                // 柱对称坐标系下的刚度矩阵项
                double stiffnessTerm = integrationFactor * (
                    avgConductivity * (dShapeDR[i] * dShapeDR[j] + dShapeDZ[i] * dShapeDZ[j]) +
                    (1.0 / avgPermeability) * (dShapeDR[i] * dShapeDR[j] + dShapeDZ[i] * dShapeDZ[j])
                );
                
                // 组装到3x3子矩阵（柱对称坐标系只有2个自由度）
                for (int dofI = 0; dofI < 2; ++dofI) { // r和z方向
                    for (int dofJ = 0; dofJ < 2; ++dofJ) {
                        int row = i * 3 + dofI;
                        int col = j * 3 + dofJ;
                        
                        // 对角线项
                        if (dofI == dofJ) {
                            elementStiffness[row][col] += stiffnessTerm;
                        }
                    }
                }
                
                // 周向分量的特殊处理（轴对称）
                double phiStiffnessTerm = integrationFactor * (
                    avgConductivity * (shapeFunctions[i] * shapeFunctions[j]) / (r * r) +
                    (1.0 / avgPermeability) * (shapeFunctions[i] * shapeFunctions[j]) / (r * r)
                );
                
                int phiRow = i * 3 + 2; // 周向分量索引
                int phiCol = j * 3 + 2;
                elementStiffness[phiRow][phiCol] += phiStiffnessTerm;
            }
        }
        
        // 计算右端向量
        for (int i = 0; i < numElementNodes; ++i) {
            double shapeFunction = shapeFunctions[i];
            
            // 考虑瞬态项
            if (transientSimulation_ && timeStep_ > 0.0) {
                for (int j = 0; j < numElementNodes; ++j) {
                    double massTerm = integrationFactor * shapeFunctions[i] * shapeFunctions[j] * 
                                    elemConductivity[j] / timeStep_;
                    
                    for (int dof = 0; dof < 3; ++dof) {
                        int row = i * 3 + dof;
                        int col = j * 3 + dof;
                        elementStiffness[row][col] += massTerm;
                    }
                }
            }
            
            // 添加源项
            for (int dof = 0; dof < 3; ++dof) {
                int index = i * 3 + dof;
                elementRHS[index] += integrationFactor * shapeFunction * 0.0; // 零源项
            }
        }
    }
    return true;
}

bool MagneticSolver::computeGeneralElementMatrices(int elementId, 
                                                   std::vector<std::vector<double>>& elementStiffness,
                                                   std::vector<double>& elementRHS) {
    // 基于Fortran版本的MagneticSolve.F90实现一般坐标系单元矩阵计算
    
    if (!mesh_ || elementId < 0 || elementId >= mesh_->getBulkElements().size()) {
        ELMER_ERROR("无效的单元索引: {}", elementId);
        return false;
    }
    
    auto& bulkElements = mesh_->getBulkElements();
    auto& element = bulkElements[elementId];
    auto elementNodes = element.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 初始化单元矩阵和右端向量
    elementStiffness.resize(numElementNodes * 3, std::vector<double>(numElementNodes * 3, 0.0));
    elementRHS.resize(numElementNodes * 3, 0.0);
    
    // 获取单元材料参数
    std::vector<double> elemConductivity(numElementNodes);
    std::vector<double> elemPermeability(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        if (nodeId < conductivity_.size() && nodeId < permeability_.size()) {
            elemConductivity[i] = conductivity_[nodeId];
            elemPermeability[i] = permeability_[nodeId];
        }
    }
    
    // 获取单元节点坐标
    std::vector<double> nodeCoordsX(numElementNodes);
    std::vector<double> nodeCoordsY(numElementNodes);
    std::vector<double> nodeCoordsZ(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        auto node = mesh_->getNodes()[nodeId];
        nodeCoordsX[i] = node.x;
        nodeCoordsY[i] = node.y;
        nodeCoordsZ[i] = node.z;
    }
    
    // 一般坐标系下的3D高斯积分
    int numGaussPoints = 8; // 对于3D单元使用8个高斯点
    
    for (int gaussPoint = 0; gaussPoint < numGaussPoints; ++gaussPoint) {
        // 计算高斯点坐标和权重
        double xi, eta, zeta, weight;
        getGaussPoint3D(gaussPoint, xi, eta, zeta, weight);
        
        // 计算形状函数和导数
        std::vector<double> shapeFunctions(numElementNodes);
        std::vector<double> dShapeDXi(numElementNodes);
        std::vector<double> dShapeDEta(numElementNodes);
        std::vector<double> dShapeDZeta(numElementNodes);
        
        computeShapeFunctions3D(xi, eta, zeta, shapeFunctions, dShapeDXi, dShapeDEta, dShapeDZeta);
        
        // 计算雅可比矩阵
        double jacobian[3][3] = {{0}};
        for (int i = 0; i < numElementNodes; ++i) {
            jacobian[0][0] += dShapeDXi[i] * nodeCoordsX[i];
            jacobian[0][1] += dShapeDXi[i] * nodeCoordsY[i];
            jacobian[0][2] += dShapeDXi[i] * nodeCoordsZ[i];
            jacobian[1][0] += dShapeDEta[i] * nodeCoordsX[i];
            jacobian[1][1] += dShapeDEta[i] * nodeCoordsY[i];
            jacobian[1][2] += dShapeDEta[i] * nodeCoordsZ[i];
            jacobian[2][0] += dShapeDZeta[i] * nodeCoordsX[i];
            jacobian[2][1] += dShapeDZeta[i] * nodeCoordsY[i];
            jacobian[2][2] += dShapeDZeta[i] * nodeCoordsZ[i];
        }
        
        // 计算雅可比行列式
        double detJ = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
                    - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
                    + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);
        
        if (detJ <= 0.0) {
            ELMER_ERROR("单元 {} 的雅可比行列式为负或零", elementId);
            return false;
        }
        
        // 计算形状函数在全局坐标系中的导数
        std::vector<double> dShapeDX(numElementNodes);
        std::vector<double> dShapeDY(numElementNodes);
        std::vector<double> dShapeDZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            // 使用雅可比矩阵的逆计算全局导数
            dShapeDX[i] = (1.0/detJ) * (
                (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) * dShapeDXi[i] +
                (jacobian[0][2] * jacobian[2][1] - jacobian[0][1] * jacobian[2][2]) * dShapeDEta[i] +
                (jacobian[0][1] * jacobian[1][2] - jacobian[0][2] * jacobian[1][1]) * dShapeDZeta[i]
            );
            
            dShapeDY[i] = (1.0/detJ) * (
                (jacobian[1][2] * jacobian[2][0] - jacobian[1][0] * jacobian[2][2]) * dShapeDXi[i] +
                (jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0]) * dShapeDEta[i] +
                (jacobian[0][2] * jacobian[1][0] - jacobian[0][0] * jacobian[1][2]) * dShapeDZeta[i]
            );
            
            dShapeDZ[i] = (1.0/detJ) * (
                (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]) * dShapeDXi[i] +
                (jacobian[0][1] * jacobian[2][0] - jacobian[0][0] * jacobian[2][1]) * dShapeDEta[i] +
                (jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0]) * dShapeDZeta[i]
            );
        }
        
        // 积分权重因子
        double integrationFactor = weight * detJ;
        
        // 计算单元刚度矩阵（一般坐标系）
        for (int i = 0; i < numElementNodes; ++i) {
            for (int j = 0; j < numElementNodes; ++j) {
                // 计算材料参数在高斯点的平均值
                double avgConductivity = 0.5 * (elemConductivity[i] + elemConductivity[j]);
                double avgPermeability = 0.5 * (elemPermeability[i] + elemPermeability[j]);
                
                // 一般坐标系下的刚度矩阵项
                double stiffnessTerm = integrationFactor * (
                    avgConductivity * (dShapeDX[i] * dShapeDX[j] + dShapeDY[i] * dShapeDY[j] + dShapeDZ[i] * dShapeDZ[j]) +
                    (1.0 / avgPermeability) * (dShapeDX[i] * dShapeDX[j] + dShapeDY[i] * dShapeDY[j] + dShapeDZ[i] * dShapeDZ[j])
                );
                
                // 组装到3x3子矩阵
                for (int dofI = 0; dofI < 3; ++dofI) {
                    for (int dofJ = 0; dofJ < 3; ++dofJ) {
                        int row = i * 3 + dofI;
                        int col = j * 3 + dofJ;
                        
                        // 对角线项
                        if (dofI == dofJ) {
                            elementStiffness[row][col] += stiffnessTerm;
                        }
                    }
                }
            }
        }
        
        // 计算右端向量
        for (int i = 0; i < numElementNodes; ++i) {
            double shapeFunction = shapeFunctions[i];
            
            // 考虑瞬态项
            if (transientSimulation_ && timeStep_ > 0.0) {
                for (int j = 0; j < numElementNodes; ++j) {
                    double massTerm = integrationFactor * shapeFunctions[i] * shapeFunctions[j] * 
                                    elemConductivity[j] / timeStep_;
                    
                    for (int dof = 0; dof < 3; ++dof) {
                        int row = i * 3 + dof;
                        int col = j * 3 + dof;
                        elementStiffness[row][col] += massTerm;
                    }
                }
            }
            
            // 添加源项
            for (int dof = 0; dof < 3; ++dof) {
                int index = i * 3 + dof;
                elementRHS[index] += integrationFactor * shapeFunction * 0.0; // 零源项
            }
        }
    }
    
    return true;
}

std::string MagneticSolver::detectCoordinateSystem() {
    // 基于Fortran版本的MagneticSolve.F90实现坐标系检测
    // 简化实现：根据网格几何特征检测坐标系类型
    
    if (!mesh_) {
        std::cerr << "错误: 网格未初始化，无法检测坐标系" << std::endl;
        return "cartesian"; // 默认返回笛卡尔坐标系
    }
    
    // 简化实现：直接返回笛卡尔坐标系
    // TODO: 修复网格节点分析代码的编译问题
    return "cartesian";
}

bool MagneticSolver::computeElementMatrices(int elementId, 
                                std::vector<std::vector<double>>& elementStiffness,
                                std::vector<double>& elementRHS) {
    // 基于Fortran版本的MagneticSolve.F90中的MaxwellCompose函数实现
    // 计算单元刚度矩阵和右端向量
    
    if (!mesh_ || elementId < 0 || elementId >= mesh_->getBulkElements().size()) {
        ELMER_ERROR("无效的单元索引: {}", elementId);
        return false;
    }
    
    auto& bulkElements = mesh_->getBulkElements();
    auto& element = bulkElements[elementId];
    auto elementNodes = element.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 初始化单元矩阵和右端向量
    elementStiffness.resize(numElementNodes * 3, std::vector<double>(numElementNodes * 3, 0.0));
    elementRHS.resize(numElementNodes * 3, 0.0);
    
    // 获取单元材料参数
    std::vector<double> elemConductivity(numElementNodes);
    std::vector<double> elemPermeability(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        if (nodeId < conductivity_.size() && nodeId < permeability_.size()) {
            elemConductivity[i] = conductivity_[nodeId];
            elemPermeability[i] = permeability_[nodeId];
        }
    }
    
    // 获取单元节点坐标
    std::vector<double> nodeCoordsX(numElementNodes);
    std::vector<double> nodeCoordsY(numElementNodes);
    std::vector<double> nodeCoordsZ(numElementNodes);
    
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        auto node = mesh_->getNodes()[nodeId];
        nodeCoordsX[i] = node.x;
        nodeCoordsY[i] = node.y;
        nodeCoordsZ[i] = node.z;
    }
    
    // 计算单元雅可比矩阵和形状函数导数
    // 基于Fortran的MaxwellCompose函数实现
    
    // 高斯积分点（根据单元类型选择）
    int numGaussPoints = 8; // 对于3D单元使用8个高斯点
    
    for (int gaussPoint = 0; gaussPoint < numGaussPoints; ++gaussPoint) {
        // 计算高斯点坐标和权重
        double xi, eta, zeta, weight;
        getGaussPoint3D(gaussPoint, xi, eta, zeta, weight);
        
        // 计算形状函数和导数
        std::vector<double> shapeFunctions(numElementNodes);
        std::vector<double> dShapeDXi(numElementNodes);
        std::vector<double> dShapeDEta(numElementNodes);
        std::vector<double> dShapeDZeta(numElementNodes);
        
        computeShapeFunctions3D(xi, eta, zeta, shapeFunctions, dShapeDXi, dShapeDEta, dShapeDZeta);
        
        // 计算雅可比矩阵
        double jacobian[3][3] = {{0}};
        for (int i = 0; i < numElementNodes; ++i) {
            jacobian[0][0] += dShapeDXi[i] * nodeCoordsX[i];
            jacobian[0][1] += dShapeDXi[i] * nodeCoordsY[i];
            jacobian[0][2] += dShapeDXi[i] * nodeCoordsZ[i];
            jacobian[1][0] += dShapeDEta[i] * nodeCoordsX[i];
            jacobian[1][1] += dShapeDEta[i] * nodeCoordsY[i];
            jacobian[1][2] += dShapeDEta[i] * nodeCoordsZ[i];
            jacobian[2][0] += dShapeDZeta[i] * nodeCoordsX[i];
            jacobian[2][1] += dShapeDZeta[i] * nodeCoordsY[i];
            jacobian[2][2] += dShapeDZeta[i] * nodeCoordsZ[i];
        }
        
        // 计算雅可比行列式
        double detJ = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
                    - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
                    + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);
        
        if (detJ <= 0.0) {
            ELMER_ERROR("单元 {} 的雅可比行列式为负或零", elementId);
            return false;
        }
        
        // 计算形状函数在全局坐标系中的导数
        std::vector<double> dShapeDX(numElementNodes);
        std::vector<double> dShapeDY(numElementNodes);
        std::vector<double> dShapeDZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            // 使用雅可比矩阵的逆计算全局导数
            dShapeDX[i] = (1.0/detJ) * (
                (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) * dShapeDXi[i] +
                (jacobian[0][2] * jacobian[2][1] - jacobian[0][1] * jacobian[2][2]) * dShapeDEta[i] +
                (jacobian[0][1] * jacobian[1][2] - jacobian[0][2] * jacobian[1][1]) * dShapeDZeta[i]
            );
            
            dShapeDY[i] = (1.0/detJ) * (
                (jacobian[1][2] * jacobian[2][0] - jacobian[1][0] * jacobian[2][2]) * dShapeDXi[i] +
                (jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0]) * dShapeDEta[i] +
                (jacobian[0][2] * jacobian[1][0] - jacobian[0][0] * jacobian[1][2]) * dShapeDZeta[i]
            );
            
            dShapeDZ[i] = (1.0/detJ) * (
                (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]) * dShapeDXi[i] +
                (jacobian[0][1] * jacobian[2][0] - jacobian[0][0] * jacobian[2][1]) * dShapeDEta[i] +
                (jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0]) * dShapeDZeta[i]
            );
        }
        
        // 计算单元刚度矩阵（基于Maxwell方程）
        for (int i = 0; i < numElementNodes; ++i) {
            for (int j = 0; j < numElementNodes; ++j) {
                // 计算材料参数在高斯点的平均值
                double avgConductivity = 0.5 * (elemConductivity[i] + elemConductivity[j]);
                double avgPermeability = 0.5 * (elemPermeability[i] + elemPermeability[j]);
                
                // 计算刚度矩阵项（涡流项 + 扩散项）
                double stiffnessTerm = weight * detJ * (
                    avgConductivity * (dShapeDX[i] * dShapeDX[j] + dShapeDY[i] * dShapeDY[j] + dShapeDZ[i] * dShapeDZ[j]) +
                    (1.0 / avgPermeability) * (dShapeDX[i] * dShapeDX[j] + dShapeDY[i] * dShapeDY[j] + dShapeDZ[i] * dShapeDZ[j])
                );
                
                // 组装到3x3子矩阵
                for (int dofI = 0; dofI < 3; ++dofI) {
                    for (int dofJ = 0; dofJ < 3; ++dofJ) {
                        int row = i * 3 + dofI;
                        int col = j * 3 + dofJ;
                        
                        // 对角线项
                        if (dofI == dofJ) {
                            elementStiffness[row][col] += stiffnessTerm;
                        }
                    }
                }
            }
        }
        
        // 计算右端向量（源项和边界条件）
        for (int i = 0; i < numElementNodes; ++i) {
            double shapeFunction = shapeFunctions[i];
            
            // 考虑瞬态项（如果启用）
            if (transientSimulation_ && timeStep_ > 0.0) {
                // 添加质量矩阵贡献
                for (int j = 0; j < numElementNodes; ++j) {
                    double massTerm = weight * detJ * shapeFunctions[i] * shapeFunctions[j] * elemConductivity[j] / timeStep_;
                    
                    for (int dof = 0; dof < 3; ++dof) {
                        int row = i * 3 + dof;
                        int col = j * 3 + dof;
                        elementStiffness[row][col] += massTerm;
                    }
                }
            }
            
            // 添加源项（简化实现）
            for (int dof = 0; dof < 3; ++dof) {
                int index = i * 3 + dof;
                elementRHS[index] += weight * detJ * shapeFunction * 0.0; // 零源项
            }
        }
    }
    
    return true;
}

bool MagneticSolver::assembleToGlobalSystem(int elementId,
                                const std::vector<std::vector<double>>& elementStiffness,
                                const std::vector<double>& elementRHS) {
    // 将单元矩阵组装到全局系统
    
    if (!stiffnessMatrix_ || !rhsVector_) {
        std::cerr << "错误: 全局系统未初始化" << std::endl;
        return false;
    }
    
    auto& bulkElements = mesh_->getBulkElements();
    if (elementId < 0 || elementId >= bulkElements.size()) {
        std::cerr << "错误: 无效的单元索引" << std::endl;
        return false;
    }
    
    auto& element = bulkElements[elementId];
    auto elementNodes = element.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 检查矩阵尺寸是否匹配
    if (elementStiffness.size() != numElementNodes * 3 || 
        elementRHS.size() != numElementNodes * 3) {
        std::cerr << "错误: 单元矩阵尺寸不匹配" << std::endl;
        return false;
    }
    
    // 组装到全局刚度矩阵
    for (int i = 0; i < numElementNodes; ++i) {
        int globalNodeI = static_cast<int>(elementNodes[i]);
        
        for (int j = 0; j < numElementNodes; ++j) {
            int globalNodeJ = static_cast<int>(elementNodes[j]);
            
            // 组装3x3子矩阵（对应每个节点的3个自由度）
            for (int dofI = 0; dofI < 3; ++dofI) {
                for (int dofJ = 0; dofJ < 3; ++dofJ) {
                    int localRow = i * 3 + dofI;
                    int localCol = j * 3 + dofJ;
                    
                    double value = elementStiffness[localRow][localCol];
                    
                    // 添加到全局矩阵
                    stiffnessMatrix_->AddToElement(globalNodeI * 3 + dofI, 
                                                   globalNodeJ * 3 + dofJ, 
                                                   value);
                }
            }
        }
    }
    
    // 组装到全局右端向量
    for (int i = 0; i < numElementNodes; ++i) {
        int globalNodeI = static_cast<int>(elementNodes[i]);
        
        for (int dofI = 0; dofI < 3; ++dofI) {
            int localIndex = i * 3 + dofI;
            double value = elementRHS[localIndex];
            
            // 添加到全局右端向量
            (*rhsVector_)[globalNodeI * 3 + dofI] += value;
        }
    }
    
    return true;
}

bool MagneticSolver::applyMagneticForceBoundaryCondition(int bcId, const Element& boundaryElement) {
    // 应用磁力边界条件
    // 基于Fortran版本的MagneticSolve.F90中的磁力边界条件处理
    
    if (!bc_ || !stiffnessMatrix_ || !rhsVector_) {
        std::cerr << "错误: 必要的组件未初始化" << std::endl;
        return false;
    }
    
    auto elementNodes = boundaryElement.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 获取边界条件值
    std::vector<double> forceX(numElementNodes, 0.0);
    std::vector<double> forceY(numElementNodes, 0.0);
    std::vector<double> forceZ(numElementNodes, 0.0);
    
    // 简化实现：设置默认值
    for (int i = 0; i < numElementNodes; ++i) {
        forceX[i] = 0.0;
        forceY[i] = 0.0;
        forceZ[i] = 0.0;
    }
    
    // 应用到右端向量
    for (int i = 0; i < numElementNodes; ++i) {
        int globalNode = static_cast<int>(elementNodes[i]);
        
        (*rhsVector_)[globalNode * 3] += forceX[i];     // Fx
        (*rhsVector_)[globalNode * 3 + 1] += forceY[i]; // Fy
        (*rhsVector_)[globalNode * 3 + 2] += forceZ[i]; // Fz
    }
    
    return true;
}

bool MagneticSolver::applyDirichletBoundaryCondition(int bcId, const Element& boundaryElement) {
    // 应用狄利克雷边界条件（固定值边界条件）
    
    if (!bc_ || !stiffnessMatrix_ || !rhsVector_) {
        std::cerr << "错误: 必要的组件未初始化" << std::endl;
        return false;
    }
    
    auto elementNodes = boundaryElement.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 获取边界条件值
    std::vector<double> fixedValueX(numElementNodes, 0.0);
    std::vector<double> fixedValueY(numElementNodes, 0.0);
    std::vector<double> fixedValueZ(numElementNodes, 0.0);
    
    // 简化实现：设置默认值
    for (int i = 0; i < numElementNodes; ++i) {
        fixedValueX[i] = 0.0;
        fixedValueY[i] = 0.0;
        fixedValueZ[i] = 0.0;
    }
    
    // 应用狄利克雷边界条件
        for (int i = 0; i < numElementNodes; ++i) {
            int globalNode = static_cast<int>(elementNodes[i]);
            
            // 设置刚度矩阵对角线为1，右端向量为固定值
            for (int dof = 0; dof < 3; ++dof) {
                int globalDof = globalNode * 3 + dof;
                
                // 手动清零该行：遍历所有列，设置非对角线元素为0
                for (int colDof = 0; colDof < stiffnessMatrix_->GetNumCols(); ++colDof) {
                    if (colDof != globalDof) {
                        stiffnessMatrix_->SetElement(globalDof, colDof, 0.0);
                    }
                }
                
                // 设置对角线为1
                stiffnessMatrix_->SetElement(globalDof, globalDof, 1.0);
                
                // 设置右端向量
                if (dof == 0) {
                    (*rhsVector_)[globalDof] = fixedValueX[i];
                } else if (dof == 1) {
                    (*rhsVector_)[globalDof] = fixedValueY[i];
                } else {
                    (*rhsVector_)[globalDof] = fixedValueZ[i];
                }
            }
        }
    
    return true;
}

bool MagneticSolver::applyNeumannBoundaryCondition(int bcId, const Element& boundaryElement) {
    // 应用诺伊曼边界条件（通量边界条件）
    
    if (!bc_ || !rhsVector_) {
        std::cerr << "错误: 必要的组件未初始化" << std::endl;
        return false;
    }
    
    auto elementNodes = boundaryElement.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 获取边界条件值
    std::vector<double> fluxX(numElementNodes, 0.0);
    std::vector<double> fluxY(numElementNodes, 0.0);
    std::vector<double> fluxZ(numElementNodes, 0.0);
    
    // 简化实现：设置默认值
    for (int i = 0; i < numElementNodes; ++i) {
        fluxX[i] = 0.0;
        fluxY[i] = 0.0;
        fluxZ[i] = 0.0;
    }
    
    // 应用到右端向量
    for (int i = 0; i < numElementNodes; ++i) {
        int globalNode = static_cast<int>(elementNodes[i]);
        
        (*rhsVector_)[globalNode * 3] += fluxX[i];     // Fx
        (*rhsVector_)[globalNode * 3 + 1] += fluxY[i]; // Fy
        (*rhsVector_)[globalNode * 3 + 2] += fluxZ[i]; // Fz
    }
    
    return true;
}

bool MagneticSolver::computeCurl(const std::vector<double>& fieldX,
                                 const std::vector<double>& fieldY, 
                                 const std::vector<double>& fieldZ,
                                 std::vector<double>& curlX,
                                 std::vector<double>& curlY,
                                 std::vector<double>& curlZ) {
    // 基于Fortran版本的旋度计算实现
    // 使用有限元方法计算向量场的旋度
    
    ELMER_INFO("开始计算旋度...");
    
    if (!mesh_ || fieldX.empty() || fieldY.empty() || fieldZ.empty()) {
        ELMER_ERROR("网格或场量数据未初始化");
        return false;
    }
    
    size_t numNodes = mesh_->numberOfNodes();
    
    // 检查输入场量尺寸
    if (fieldX.size() != numNodes || fieldY.size() != numNodes || fieldZ.size() != numNodes) {
        ELMER_ERROR("场量尺寸与节点数不匹配: fieldX={}, fieldY={}, fieldZ={}, nodes={}", 
                   fieldX.size(), fieldY.size(), fieldZ.size(), numNodes);
        return false;
    }
    
    // 初始化输出向量
    curlX.resize(numNodes, 0.0);
    curlY.resize(numNodes, 0.0);
    curlZ.resize(numNodes, 0.0);
    
    // 获取网格节点坐标
    auto& nodes = mesh_->getNodes();
    
    // 使用有限元方法计算旋度
    // 旋度公式：∇ × F = (∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y)
    
    // 遍历所有单元计算旋度
    auto& bulkElements = mesh_->getBulkElements();
    int numElements = bulkElements.size();
    
    for (int elemId = 0; elemId < numElements; ++elemId) {
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        int numElementNodes = elementNodes.size();
        
        // 获取单元节点坐标
        std::vector<double> nodeCoordsX(numElementNodes);
        std::vector<double> nodeCoordsY(numElementNodes);
        std::vector<double> nodeCoordsZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            int nodeId = static_cast<int>(elementNodes[i]);
            auto node = mesh_->getNodes()[nodeId];
            nodeCoordsX[i] = node.x;
            nodeCoordsY[i] = node.y;
            nodeCoordsZ[i] = node.z;
        }
        
        // 获取单元场量值
        std::vector<double> elemFieldX(numElementNodes);
        std::vector<double> elemFieldY(numElementNodes);
        std::vector<double> elemFieldZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            int nodeId = static_cast<int>(elementNodes[i]);
            elemFieldX[i] = fieldX[nodeId];
            elemFieldY[i] = fieldY[nodeId];
            elemFieldZ[i] = fieldZ[nodeId];
        }
        
        // 计算单元旋度
        std::vector<double> elemCurlX(numElementNodes, 0.0);
        std::vector<double> elemCurlY(numElementNodes, 0.0);
        std::vector<double> elemCurlZ(numElementNodes, 0.0);
        
        if (!computeElementCurl(elemId, nodeCoordsX, nodeCoordsY, nodeCoordsZ,
                               elemFieldX, elemFieldY, elemFieldZ,
                               elemCurlX, elemCurlY, elemCurlZ)) {
            ELMER_ERROR("单元 {} 旋度计算失败", elemId);
            return false;
        }
        
        // 将单元旋度组装到节点旋度
        for (int i = 0; i < numElementNodes; ++i) {
            int nodeId = static_cast<int>(elementNodes[i]);
            
            // 累加旋度值（后续需要平均）
            curlX[nodeId] += elemCurlX[i];
            curlY[nodeId] += elemCurlY[i];
            curlZ[nodeId] += elemCurlZ[i];
        }
    }
    
    // 计算节点平均旋度
    std::vector<int> nodeElementCount(numNodes, 0);
    
    // 统计每个节点所属的单元数
    for (int elemId = 0; elemId < numElements; ++elemId) {
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        
        for (int i = 0; i < elementNodes.size(); ++i) {
            int nodeId = static_cast<int>(elementNodes[i]);
            if (nodeId < numNodes) {
                nodeElementCount[nodeId]++;
            }
        }
    }
    
    // 平均旋度值
    for (size_t i = 0; i < numNodes; ++i) {
        if (nodeElementCount[i] > 0) {
            curlX[i] /= nodeElementCount[i];
            curlY[i] /= nodeElementCount[i];
            curlZ[i] /= nodeElementCount[i];
        }
    }
    
    ELMER_INFO("旋度计算完成，节点数: {}", numNodes);
    return true;
}

bool MagneticSolver::computeElementCurl(int elementId,
                                       const std::vector<double>& nodeCoordsX,
                                       const std::vector<double>& nodeCoordsY,
                                       const std::vector<double>& nodeCoordsZ,
                                       const std::vector<double>& fieldX,
                                       const std::vector<double>& fieldY,
                                       const std::vector<double>& fieldZ,
                                       std::vector<double>& curlX,
                                       std::vector<double>& curlY,
                                       std::vector<double>& curlZ) {
    // 计算单元旋度
    // 使用高斯积分计算旋度在单元内的分布
    
    int numElementNodes = nodeCoordsX.size();
    
    // 初始化输出
    curlX.assign(numElementNodes, 0.0);
    curlY.assign(numElementNodes, 0.0);
    curlZ.assign(numElementNodes, 0.0);
    
    // 高斯积分点
    int numGaussPoints = 8;
    
    for (int gaussPoint = 0; gaussPoint < numGaussPoints; ++gaussPoint) {
        double xi, eta, zeta, weight;
        getGaussPoint3D(gaussPoint, xi, eta, zeta, weight);
        
        // 计算形状函数和导数
        std::vector<double> shapeFunctions;
        std::vector<double> dShapeDXi, dShapeDEta, dShapeDZeta;
        computeShapeFunctions3D(xi, eta, zeta, shapeFunctions, dShapeDXi, dShapeDEta, dShapeDZeta);
        
        // 计算雅可比矩阵
        double jacobian[3][3] = {{0}};
        for (int i = 0; i < numElementNodes; ++i) {
            jacobian[0][0] += dShapeDXi[i] * nodeCoordsX[i];
            jacobian[0][1] += dShapeDXi[i] * nodeCoordsY[i];
            jacobian[0][2] += dShapeDXi[i] * nodeCoordsZ[i];
            jacobian[1][0] += dShapeDEta[i] * nodeCoordsX[i];
            jacobian[1][1] += dShapeDEta[i] * nodeCoordsY[i];
            jacobian[1][2] += dShapeDEta[i] * nodeCoordsZ[i];
            jacobian[2][0] += dShapeDZeta[i] * nodeCoordsX[i];
            jacobian[2][1] += dShapeDZeta[i] * nodeCoordsY[i];
            jacobian[2][2] += dShapeDZeta[i] * nodeCoordsZ[i];
        }
        
        // 计算雅可比行列式
        double detJ = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
                    - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
                    + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);
        
        if (detJ <= 0.0) {
            ELMER_ERROR("单元 {} 的雅可比行列式为负或零", elementId);
            return false;
        }
        
        // 计算形状函数在全局坐标系中的导数
        std::vector<double> dShapeDX(numElementNodes);
        std::vector<double> dShapeDY(numElementNodes);
        std::vector<double> dShapeDZ(numElementNodes);
        
        for (int i = 0; i < numElementNodes; ++i) {
            // 使用雅可比矩阵的逆计算全局导数
            dShapeDX[i] = (1.0/detJ) * (
                (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) * dShapeDXi[i] +
                (jacobian[0][2] * jacobian[2][1] - jacobian[0][1] * jacobian[2][2]) * dShapeDEta[i] +
                (jacobian[0][1] * jacobian[1][2] - jacobian[0][2] * jacobian[1][1]) * dShapeDZeta[i]
            );
            
            dShapeDY[i] = (1.0/detJ) * (
                (jacobian[1][2] * jacobian[2][0] - jacobian[1][0] * jacobian[2][2]) * dShapeDXi[i] +
                (jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0]) * dShapeDEta[i] +
                (jacobian[0][2] * jacobian[1][0] - jacobian[0][0] * jacobian[1][2]) * dShapeDZeta[i]
            );
            
            dShapeDZ[i] = (1.0/detJ) * (
                (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]) * dShapeDXi[i] +
                (jacobian[0][1] * jacobian[2][0] - jacobian[0][0] * jacobian[2][1]) * dShapeDEta[i] +
                (jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0]) * dShapeDZeta[i]
            );
        }
        
        // 计算高斯点的场量值
        double Fx = 0.0, Fy = 0.0, Fz = 0.0;
        for (int i = 0; i < numElementNodes; ++i) {
            Fx += shapeFunctions[i] * fieldX[i];
            Fy += shapeFunctions[i] * fieldY[i];
            Fz += shapeFunctions[i] * fieldZ[i];
        }
        
        // 计算场量的导数（用于旋度计算）
        double dFx_dx = 0.0, dFx_dy = 0.0, dFx_dz = 0.0;
        double dFy_dx = 0.0, dFy_dy = 0.0, dFy_dz = 0.0;
        double dFz_dx = 0.0, dFz_dy = 0.0, dFz_dz = 0.0;
        
        for (int i = 0; i < numElementNodes; ++i) {
            dFx_dx += dShapeDX[i] * fieldX[i];
            dFx_dy += dShapeDY[i] * fieldX[i];
            dFx_dz += dShapeDZ[i] * fieldX[i];
            
            dFy_dx += dShapeDX[i] * fieldY[i];
            dFy_dy += dShapeDY[i] * fieldY[i];
            dFy_dz += dShapeDZ[i] * fieldY[i];
            
            dFz_dx += dShapeDX[i] * fieldZ[i];
            dFz_dy += dShapeDY[i] * fieldZ[i];
            dFz_dz += dShapeDZ[i] * fieldZ[i];
        }
        
        // 计算旋度：∇ × F = (∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y)
        double curlX_gauss = dFz_dy - dFy_dz;
        double curlY_gauss = dFx_dz - dFz_dx;
        double curlZ_gauss = dFy_dx - dFx_dy;
        
        // 将高斯点的旋度分配到节点
        for (int i = 0; i < numElementNodes; ++i) {
            curlX[i] += weight * detJ * shapeFunctions[i] * curlX_gauss;
            curlY[i] += weight * detJ * shapeFunctions[i] * curlY_gauss;
            curlZ[i] += weight * detJ * shapeFunctions[i] * curlZ_gauss;
        }
    }
    
    return true;
}

bool MagneticSolver::computeNodalField() {
    // TODO: 实现节点场计算
    std::cout << "开始计算节点场..." << std::endl;
    
    // 简化实现：输出信息
    std::cout << "节点场计算完成" << std::endl;
    return true;
}

// ===== 非线性迭代求解辅助函数实现 =====

std::shared_ptr<Vector> MagneticSolver::computeResidualVector() {
    // 计算残差向量 r = f(x) - Kx
    if (!stiffnessMatrix_ || !solutionVector_ || !rhsVector_) {
        std::cerr << "错误: 系统矩阵、解向量或右端向量未初始化" << std::endl;
        return nullptr;
    }
    
    auto residual = Vector::Create(solutionVector_->Size());
    
    // 计算 Kx
    auto Kx = Vector::Create(solutionVector_->Size());
    stiffnessMatrix_->Multiply(*solutionVector_, *Kx);
    
    // 计算残差 r = f - Kx
    for (int i = 0; i < residual->Size(); ++i) {
        (*residual)[i] = (*rhsVector_)[i] - (*Kx)[i];
    }
    
    return residual;
}

double MagneticSolver::computeVectorNorm(const Vector& vec) {
    // 计算向量的L2范数
    double norm = 0.0;
    for (int i = 0; i < vec.Size(); ++i) {
        norm += vec[i] * vec[i];
    }
    return std::sqrt(norm);
}

bool MagneticSolver::updateJacobianMatrix() {
    // 基于Fortran版本的雅可比矩阵更新实现
    // 对于非线性问题，需要根据当前解重新计算雅可比矩阵
    
    if (!mesh_ || !stiffnessMatrix_ || !solutionVector_) {
        ELMER_ERROR("必要的组件未初始化");
        return false;
    }
    
    // 检查是否为非线性问题
    bool isNonlinearProblem = checkNonlinearity();
    
    if (!isNonlinearProblem) {
        // 线性问题：雅可比矩阵就是刚度矩阵
        ELMER_DEBUG("线性问题，使用当前刚度矩阵作为雅可比矩阵");
        return true;
    }
    
    // 非线性问题：需要重新计算雅可比矩阵
    ELMER_DEBUG("非线性问题，重新计算雅可比矩阵");
    
    // 保存当前刚度矩阵（作为雅可比矩阵的基础）
    auto originalStiffness = stiffnessMatrix_->Clone();
    
    // 清零当前刚度矩阵（准备重新组装）
    stiffnessMatrix_->Zero();
    
    // 基于当前解重新组装系统矩阵
    int numElements = mesh_->getBulkElements().size();
    
    for (int elemId = 0; elemId < numElements; ++elemId) {
        // 获取单元信息
        auto& bulkElements = mesh_->getBulkElements();
        if (elemId >= bulkElements.size()) {
            ELMER_ERROR("单元索引超出范围");
            return false;
        }
        
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        
        // 计算单元雅可比矩阵（考虑非线性材料）
        std::vector<std::vector<double>> elementJacobian;
        std::vector<double> elementRHS;
        
        if (!computeNonlinearElementJacobian(elemId, elementJacobian, elementRHS)) {
            ELMER_ERROR("单元 {} 非线性雅可比矩阵计算失败", elemId);
            return false;
        }
        
        // 组装到全局雅可比矩阵
        if (!assembleToGlobalSystem(elemId, elementJacobian, elementRHS)) {
            ELMER_ERROR("单元 {} 雅可比矩阵组装失败", elemId);
            return false;
        }
    }
    
    // 重新组装边界条件
    if (!assembleBoundaryConditions()) {
        ELMER_ERROR("边界条件重新组装失败");
        return false;
    }
    
    ELMER_DEBUG("雅可比矩阵更新完成");
    return true;
}

bool MagneticSolver::checkNonlinearity() {
    // 检查问题是否为非线性
    // 基于材料参数和场量的非线性特性
    
    if (!materialDB_ || !mesh_) {
        ELMER_DEBUG("材料数据库或网格未设置，假设为线性问题");
        return false;
    }
    
    // 检查材料非线性（如非线性磁导率）
    int numNodes = mesh_->numberOfNodes();
    
    for (int i = 0; i < numNodes; ++i) {
        // 检查磁导率是否与磁场强度相关（非线性磁性材料）
        double B = 0.0;
        if (magneticField_.size() > i) {
            B = magneticField_[i];
        }
        
        // 简化实现：检查磁导率是否随磁场变化
        // 实际应该根据材料模型进行判断
        if (B > 1.0) { // 磁场强度较大时可能进入非线性区域
            ELMER_DEBUG("检测到非线性材料特性，磁场强度={}", B);
            return true;
        }
    }
    
    // 检查速度场的非线性效应
    for (int i = 0; i < numNodes; ++i) {
        double velocity = 0.0;
        if (U_.size() > i) {
            velocity = std::sqrt(U_[i] * U_[i] + 
                               V_[i] * V_[i] + 
                               W_[i] * W_[i]);
        }
        
        if (velocity > 0.1) { // 速度较大时可能产生非线性对流项
            ELMER_DEBUG("检测到非线性对流效应，速度={}", velocity);
            return true;
        }
    }
    
    ELMER_DEBUG("问题为线性");
    return false;
}

bool MagneticSolver::computeNonlinearElementJacobian(int elementId, 
                                                    std::vector<std::vector<double>>& elementJacobian,
                                                    std::vector<double>& elementRHS) {
    // 计算非线性单元的雅可比矩阵
    // 基于当前解考虑材料非线性
    
    if (!mesh_ || elementId < 0 || elementId >= mesh_->getBulkElements().size()) {
        ELMER_ERROR("无效的单元索引: {}", elementId);
        return false;
    }
    
    auto& bulkElements = mesh_->getBulkElements();
    auto& element = bulkElements[elementId];
    auto elementNodes = element.getNodeIndices();
    int numElementNodes = elementNodes.size();
    
    // 初始化雅可比矩阵和右端向量
    elementJacobian.resize(numElementNodes * 3, std::vector<double>(numElementNodes * 3, 0.0));
    elementRHS.resize(numElementNodes * 3, 0.0);
    
    // 获取单元当前解
    std::vector<double> elementSolution(numElementNodes * 3, 0.0);
    for (int i = 0; i < numElementNodes; ++i) {
        int nodeId = static_cast<int>(elementNodes[i]);
        for (int dof = 0; dof < 3; ++dof) {
            int globalDof = nodeId * 3 + dof;
            if (solutionVector_ && globalDof < solutionVector_->Size()) {
                elementSolution[i * 3 + dof] = (*solutionVector_)[globalDof];
            }
        }
    }
    
    // 计算非线性雅可比矩阵（考虑材料参数随场量的变化）
    // 这里实现简化的非线性模型
    
    for (int gaussPoint = 0; gaussPoint < 8; ++gaussPoint) {
        double xi, eta, zeta, weight;
        getGaussPoint3D(gaussPoint, xi, eta, zeta, weight);
        
        // 计算形状函数和导数
        std::vector<double> shapeFunctions;
        std::vector<double> dShapeDXi, dShapeDEta, dShapeDZeta;
        computeShapeFunctions3D(xi, eta, zeta, shapeFunctions, dShapeDXi, dShapeDEta, dShapeDZeta);
        
        // 计算高斯点的磁场强度
        double Bx = 0.0, By = 0.0, Bz = 0.0;
        for (int i = 0; i < numElementNodes; ++i) {
            Bx += shapeFunctions[i] * elementSolution[i * 3];
            By += shapeFunctions[i] * elementSolution[i * 3 + 1];
            Bz += shapeFunctions[i] * elementSolution[i * 3 + 2];
        }
        
        double B = std::sqrt(Bx * Bx + By * By + Bz * Bz);
        
        // 计算非线性磁导率（简化模型）
        double nonlinearPermeability = 1.0 / (1.0 + 0.1 * B * B); // 饱和磁性材料模型
        
        // 计算雅可比矩阵项（考虑非线性）
        for (int i = 0; i < numElementNodes; ++i) {
            for (int j = 0; j < numElementNodes; ++j) {
                // 考虑非线性磁导率的贡献
                double jacobianTerm = weight * nonlinearPermeability * 
                                    (dShapeDXi[i] * dShapeDXi[j] + 
                                     dShapeDEta[i] * dShapeDEta[j] + 
                                     dShapeDZeta[i] * dShapeDZeta[j]);
                
                for (int dofI = 0; dofI < 3; ++dofI) {
                    for (int dofJ = 0; dofJ < 3; ++dofJ) {
                        int row = i * 3 + dofI;
                        int col = j * 3 + dofJ;
                        
                        if (dofI == dofJ) {
                            elementJacobian[row][col] += jacobianTerm;
                        }
                    }
                }
            }
        }
    }
    
    return true;
}

std::shared_ptr<Vector> MagneticSolver::solveLinearSystem(const Vector& residual) {
    // 求解线性系统 J * Δx = -r
    if (!stiffnessMatrix_) {
        std::cerr << "错误: 雅可比矩阵未初始化" << std::endl;
        return nullptr;
    }
    
    // 创建右端向量 -r
    auto negResidual = Vector::Create(residual.Size());
    for (int i = 0; i < residual.Size(); ++i) {
        (*negResidual)[i] = -residual[i];
    }
    
    // 使用迭代求解器求解线性系统
    // TODO: 实现更高效的线性求解器
    
    // 简化实现：直接返回一个小的增量向量
    auto deltaX = Vector::Create(residual.Size());
    for (int i = 0; i < deltaX->Size(); ++i) {
        (*deltaX)[i] = 0.01 * (*negResidual)[i]; // 简化步长控制
    }
    
    return deltaX;
}

void MagneticSolver::updateSolutionVector(const Vector& deltaX) {
    // 更新解向量 x = x + Δx
    if (!solutionVector_ || deltaX.Size() != solutionVector_->Size()) {
        std::cerr << "错误: 解向量或增量向量尺寸不匹配" << std::endl;
        return;
    }
    
    for (int i = 0; i < solutionVector_->Size(); ++i) {
        (*solutionVector_)[i] += deltaX[i];
    }
}

void MagneticSolver::updateMagneticFieldFromSolution() {
    // 从解向量更新磁场变量
    if (!solutionVector_ || !mesh_) {
        std::cerr << "错误: 解向量或网格未初始化" << std::endl;
        return;
    }
    
    int numNodes = mesh_->numberOfNodes();
    
    // 确保磁场向量大小正确
    if (magneticField_.size() < numNodes) {
        magneticField_.resize(numNodes, 0.0);
    }
    
    // 从解向量提取磁场分量
        // 假设解向量按 [Bx1, By1, Bz1, Bx2, By2, Bz2, ...] 排列
        for (int i = 0; i < numNodes; ++i) {
            if (solutionVector_->Size() >= 3 * (i + 1)) {
                // 计算磁场强度（简化实现：使用X分量）
                magneticField_[i] = std::sqrt(
                    (*solutionVector_)[3 * i] * (*solutionVector_)[3 * i] +
                    (*solutionVector_)[3 * i + 1] * (*solutionVector_)[3 * i + 1] +
                    (*solutionVector_)[3 * i + 2] * (*solutionVector_)[3 * i + 2]
                );
            }
        }
    }

// ===== 高斯积分和形状函数计算 =====

void MagneticSolver::getGaussPoint3D(int gaussPoint, double& xi, double& eta, double& zeta, double& weight) {
    // 3D高斯积分点坐标和权重（8点积分）
    // 基于标准高斯积分表
    
    static const double gaussPoints[8][4] = {
        {-0.577350269189626, -0.577350269189626, -0.577350269189626, 1.0},
        { 0.577350269189626, -0.577350269189626, -0.577350269189626, 1.0},
        {-0.577350269189626,  0.577350269189626, -0.577350269189626, 1.0},
        { 0.577350269189626,  0.577350269189626, -0.577350269189626, 1.0},
        {-0.577350269189626, -0.577350269189626,  0.577350269189626, 1.0},
        { 0.577350269189626, -0.577350269189626,  0.577350269189626, 1.0},
        {-0.577350269189626,  0.577350269189626,  0.577350269189626, 1.0},
        { 0.577350269189626,  0.577350269189626,  0.577350269189626, 1.0}
    };
    
    if (gaussPoint >= 0 && gaussPoint < 8) {
        xi = gaussPoints[gaussPoint][0];
        eta = gaussPoints[gaussPoint][1];
        zeta = gaussPoints[gaussPoint][2];
        weight = gaussPoints[gaussPoint][3];
    } else {
        // 默认值
        xi = 0.0;
        eta = 0.0;
        zeta = 0.0;
        weight = 1.0;
    }
}

void MagneticSolver::computeShapeFunctions3D(double xi, double eta, double zeta,
                                            std::vector<double>& shapeFunctions,
                                            std::vector<double>& dShapeDXi,
                                            std::vector<double>& dShapeDEta,
                                            std::vector<double>& dShapeDZeta) {
    // 计算3D线性六面体单元的形函数和导数
    // 假设单元有8个节点
    
    int numNodes = 8; // 六面体单元有8个节点
    shapeFunctions.resize(numNodes);
    dShapeDXi.resize(numNodes);
    dShapeDEta.resize(numNodes);
    dShapeDZeta.resize(numNodes);
    
    // 标准六面体单元形函数
    shapeFunctions[0] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
    shapeFunctions[1] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta);
    shapeFunctions[2] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta);
    shapeFunctions[3] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta);
    shapeFunctions[4] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta);
    shapeFunctions[5] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta);
    shapeFunctions[6] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta);
    shapeFunctions[7] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta);
    
    // 形函数对ξ的导数
    dShapeDXi[0] = -0.125 * (1.0 - eta) * (1.0 - zeta);
    dShapeDXi[1] =  0.125 * (1.0 - eta) * (1.0 - zeta);
    dShapeDXi[2] =  0.125 * (1.0 + eta) * (1.0 - zeta);
    dShapeDXi[3] = -0.125 * (1.0 + eta) * (1.0 - zeta);
    dShapeDXi[4] = -0.125 * (1.0 - eta) * (1.0 + zeta);
    dShapeDXi[5] =  0.125 * (1.0 - eta) * (1.0 + zeta);
    dShapeDXi[6] =  0.125 * (1.0 + eta) * (1.0 + zeta);
    dShapeDXi[7] = -0.125 * (1.0 + eta) * (1.0 + zeta);
    
    // 形函数对η的导数
    dShapeDEta[0] = -0.125 * (1.0 - xi) * (1.0 - zeta);
    dShapeDEta[1] = -0.125 * (1.0 + xi) * (1.0 - zeta);
    dShapeDEta[2] =  0.125 * (1.0 + xi) * (1.0 - zeta);
    dShapeDEta[3] =  0.125 * (1.0 - xi) * (1.0 - zeta);
    dShapeDEta[4] = -0.125 * (1.0 - xi) * (1.0 + zeta);
    dShapeDEta[5] = -0.125 * (1.0 + xi) * (1.0 + zeta);
    dShapeDEta[6] =  0.125 * (1.0 + xi) * (1.0 + zeta);
    dShapeDEta[7] =  0.125 * (1.0 - xi) * (1.0 + zeta);
    
    // 形函数对ζ的导数
    dShapeDZeta[0] = -0.125 * (1.0 - xi) * (1.0 - eta);
    dShapeDZeta[1] = -0.125 * (1.0 + xi) * (1.0 - eta);
    dShapeDZeta[2] = -0.125 * (1.0 + xi) * (1.0 + eta);
    dShapeDZeta[3] = -0.125 * (1.0 - xi) * (1.0 + eta);
    dShapeDZeta[4] =  0.125 * (1.0 - xi) * (1.0 - eta);
    dShapeDZeta[5] =  0.125 * (1.0 + xi) * (1.0 - eta);
    dShapeDZeta[6] =  0.125 * (1.0 + xi) * (1.0 + eta);
    dShapeDZeta[7] =  0.125 * (1.0 - xi) * (1.0 + eta);
}

// ===== 2D高斯积分和形状函数计算 =====

void MagneticSolver::getGaussPoint2D(int gaussPoint, double& xi, double& eta, double& weight) {
    // 2D高斯积分点坐标和权重（4点积分）
    // 基于标准高斯积分表
    
    static const double gaussPoints[4][3] = {
        {-0.577350269189626, -0.577350269189626, 1.0},
        { 0.577350269189626, -0.577350269189626, 1.0},
        {-0.577350269189626,  0.577350269189626, 1.0},
        { 0.577350269189626,  0.577350269189626, 1.0}
    };
    
    if (gaussPoint >= 0 && gaussPoint < 4) {
        xi = gaussPoints[gaussPoint][0];
        eta = gaussPoints[gaussPoint][1];
        weight = gaussPoints[gaussPoint][2];
    } else {
        // 默认值
        xi = 0.0;
        eta = 0.0;
        weight = 1.0;
    }
}

void MagneticSolver::computeShapeFunctions2D(double xi, double eta,
                                            std::vector<double>& shapeFunctions,
                                            std::vector<double>& dShapeDXi,
                                            std::vector<double>& dShapeDEta) {
    // 计算2D线性四边形单元的形函数和导数
    // 假设单元有4个节点
    
    int numNodes = 4; // 四边形单元有4个节点
    shapeFunctions.resize(numNodes);
    dShapeDXi.resize(numNodes);
    dShapeDEta.resize(numNodes);
    
    // 标准四边形单元形函数
    shapeFunctions[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    shapeFunctions[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    shapeFunctions[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    shapeFunctions[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
    
    // 形函数对ξ的导数
    dShapeDXi[0] = -0.25 * (1.0 - eta);
    dShapeDXi[1] =  0.25 * (1.0 - eta);
    dShapeDXi[2] =  0.25 * (1.0 + eta);
    dShapeDXi[3] = -0.25 * (1.0 + eta);
    
    // 形函数对η的导数
    dShapeDEta[0] = -0.25 * (1.0 - xi);
    dShapeDEta[1] = -0.25 * (1.0 + xi);
    dShapeDEta[2] =  0.25 * (1.0 + xi);
    dShapeDEta[3] =  0.25 * (1.0 - xi);
}

void MagneticSolver::deallocateMemory() {
    // 释放所有分配的内存
    magneticField_.clear();
    electricCurrent_.clear();
    lorentzForce_.clear();
    electricField_.clear();
    
    U_.clear();
    V_.clear();
    W_.clear();
    
    conductivity_.clear();
    permeability_.clear();
    
    appliedMagneticFieldX_.clear();
    appliedMagneticFieldY_.clear();
    appliedMagneticFieldZ_.clear();
    
    allocationsDone_ = false;
}

} // namespace elmer