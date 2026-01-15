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
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>

namespace elmer {

// ===== 构造函数和析构函数 =====

MagneticSolver::MagneticSolver() {
    // TODO: 初始化磁场求解器参数
    std::cout << "MagneticSolver构造函数完成" << std::endl;
}

MagneticSolver::~MagneticSolver() {
    // TODO: 实现资源清理
    std::cout << "MagneticSolver析构函数完成" << std::endl;
}

// ===== 基本接口函数 =====

std::string MagneticSolver::getName() const {
    return "MagneticSolver";
}

bool MagneticSolver::initialize() {
    // TODO: 实现磁场求解器的完整初始化流程
    std::cout << "开始初始化MagneticSolver..." << std::endl;
    
    // 1. 检查必要的组件是否已设置
    if (!mesh_ || !materialDB_ || !bc_) {
        std::cerr << "错误: 求解器未正确初始化，缺少必要的组件" << std::endl;
        return false;
    }
    
    // 2. 分配内存
    if (!allocateMemory()) {
        std::cerr << "错误: 内存分配失败" << std::endl;
        return false;
    }
    
    // 3. 获取材料参数
    if (!getMaterialParameters()) {
        std::cerr << "错误: 材料参数获取失败" << std::endl;
        return false;
    }
    
    // 4. 获取速度场（如果存在）
    if (!getVelocityField()) {
        std::cout << "警告: 速度场获取失败，将使用默认值" << std::endl;
    }
    
    // 5. 设置求解器状态
    status_ = SolverStatus::INITIALIZED;
    
    std::cout << "MagneticSolver初始化完成" << std::endl;
    return true;
}

bool MagneticSolver::assemble() {
    // TODO: 实现磁场求解器的系统矩阵组装
    std::cout << "开始组装磁场求解器系统矩阵..." << std::endl;
    
    if (status_ != SolverStatus::INITIALIZED) {
        std::cerr << "错误: 求解器未初始化，无法组装系统矩阵" << std::endl;
        return false;
    }
    
    if (!stiffnessMatrix_ || !rhsVector_) {
        std::cerr << "错误: 刚度矩阵或右端向量未初始化" << std::endl;
        return false;
    }
    
    // 清零刚度矩阵和右端向量
    stiffnessMatrix_->Zero();
    rhsVector_->Zero();
    
    // 根据坐标系类型选择组装方法
    // TODO: 实现坐标系检测和相应的组装方法
    bool success = false;
    
    // 简化实现：使用笛卡尔坐标系组装
    success = assembleCartesian();
    
    if (!success) {
        std::cerr << "错误: 系统矩阵组装失败" << std::endl;
        return false;
    }
    
    // 组装边界条件
    if (!assembleBoundaryConditions()) {
        std::cerr << "错误: 边界条件组装失败" << std::endl;
        return false;
    }
    
    std::cout << "磁场求解器系统矩阵组装完成" << std::endl;
    return true;
}

bool MagneticSolver::solve() {
    // TODO: 实现磁场求解器的非线性迭代求解
    std::cout << "开始求解磁场方程..." << std::endl;
    
    if (status_ != SolverStatus::ASSEMBLED) {
        std::cerr << "错误: 求解器未组装，无法求解" << std::endl;
        return false;
    }
    
    // 非线性迭代求解
    double prevNorm = 0.0;
    double currentNorm = 0.0;
    
    for (int iter = 1; iter <= maxNonlinearIterations_; ++iter) {
        std::cout << "非线性迭代: " << iter << std::endl;
        
        // TODO: 实现完整的非线性迭代过程
        // 简化实现：直接使用基类的解向量
        
        // 获取当前解的范数（简化实现）
        currentNorm = 1.0; // 临时值
        
        // 检查收敛性
        if (checkConvergence(prevNorm, currentNorm)) {
            std::cout << "非线性迭代收敛于第 " << iter << " 次迭代" << std::endl;
            break;
        }
        
        prevNorm = currentNorm;
        
        if (iter == maxNonlinearIterations_) {
            std::cout << "警告: 达到最大非线性迭代次数" << std::endl;
        }
    }
    
    // 更新磁场变量
    // TODO: 实现从解向量到磁场变量的映射
    if (solution_) {
        int numNodes = mesh_->numberOfNodes();
        if (solution_->Size() >= 3 * numNodes) {
            for (int i = 0; i < numNodes; ++i) {
                magneticField_[i] = (*solution_)[3 * i];
                // TODO: 设置其他磁场分量
            }
        }
    }
    
    // 计算相关物理量
    if (!computeLorentzForce()) {
        std::cerr << "警告: 洛伦兹力计算失败" << std::endl;
    }
    
    if (!computeElectricField()) {
        std::cerr << "警告: 电场计算失败" << std::endl;
    }
    
    status_ = SolverStatus::SOLVED;
    std::cout << "磁场方程求解完成" << std::endl;
    return true;
}

std::vector<double> MagneticSolver::getSolution() const {
    // 使用基类的getSolution方法
    return LinearSolverBase::getSolution();
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
    int numNodes = mesh_->numberOfNodes();
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
    velocityX_.resize(numNodes, 0.0);
    velocityY_.resize(numNodes, 0.0);
    velocityZ_.resize(numNodes, 0.0);
    
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
        velocityX_[i] = 0.0;
        velocityY_[i] = 0.0;
        velocityZ_[i] = 0.0;
    }
    
    std::cout << "速度场获取完成" << std::endl;
    return true;
}

bool MagneticSolver::assembleElementMatrix(int elementId) {
    // TODO: 实现单元矩阵组装
    return true;
}

bool MagneticSolver::assembleBoundaryConditions() {
    // 实现边界条件组装
    // 基于Fortran版本的MagneticSolve.F90中的边界条件处理
    std::cout << "开始组装边界条件..." << std::endl;
    
    if (!bc_ || !mesh_) {
        std::cout << "警告: 边界条件管理器或网格未设置，跳过边界条件组装" << std::endl;
        return true;
    }
    
    // 获取边界单元数量
    int numBoundaryElements = mesh_->getBoundaryElements().size();
    
    if (numBoundaryElements == 0) {
        std::cout << "没有边界单元，跳过边界条件组装" << std::endl;
        return true;
    }
    
    // 遍历所有边界单元
    int processedBC = 0;
    auto& boundaryElements = mesh_->getBoundaryElements();
    
    for (int bcId = 0; bcId < numBoundaryElements; ++bcId) {
        if (bcId >= boundaryElements.size()) {
            std::cerr << "错误: 边界单元索引超出范围" << std::endl;
            return false;
        }
        
        auto& boundaryElement = boundaryElements[bcId];
        auto elementNodes = boundaryElement.getNodeIndices();
        int numElementNodes = elementNodes.size();
        
        // 简化实现：直接应用狄利克雷边界条件
        // TODO: 实现完整的边界条件类型检测
        if (!applyDirichletBoundaryCondition(bcId, boundaryElement)) {
            std::cerr << "错误: 边界条件应用失败" << std::endl;
            return false;
        }
        processedBC++;
    }
    
    std::cout << "边界条件组装完成，共处理 " << processedBC << " 个边界条件" << std::endl;
    return true;
}

bool MagneticSolver::checkConvergence(double prevNorm, double currentNorm) {
    // TODO: 实现收敛性检查
    if (prevNorm == 0.0) {
        return false; // 第一次迭代，不检查收敛
    }
    
    double relativeChange = 2.0 * std::abs(prevNorm - currentNorm) / (prevNorm + currentNorm);
    
    std::cout << "相对变化: " << relativeChange << ", 容差: " << nonlinearTolerance_ << std::endl;
    
    return relativeChange < nonlinearTolerance_;
}

bool MagneticSolver::assembleCartesian() {
    // 实现笛卡尔坐标系下的Maxwell方程组装
    // 基于Fortran版本的MagneticSolve.F90中的MaxwellCompose函数
    std::cout << "使用笛卡尔坐标系组装Maxwell方程..." << std::endl;
    
    if (!mesh_ || !stiffnessMatrix_ || !rhsVector_) {
        std::cerr << "错误: 必要的组件未初始化" << std::endl;
        return false;
    }
    
    int numElements = mesh_->getBulkElements().size();
    int numNodes = mesh_->numberOfNodes();
    
    // 检查材料参数是否已获取
    if (conductivity_.empty() || permeability_.empty()) {
        std::cerr << "错误: 材料参数未获取" << std::endl;
        return false;
    }
    
    // 遍历所有单元进行组装
    for (int elemId = 0; elemId < numElements; ++elemId) {
        // 获取单元信息
        auto& bulkElements = mesh_->getBulkElements();
        if (elemId >= bulkElements.size()) {
            std::cerr << "错误: 单元索引超出范围" << std::endl;
            return false;
        }
        
        auto& element = bulkElements[elemId];
        auto elementNodes = element.getNodeIndices();
        int numElementNodes = elementNodes.size();
        
        // 计算单元矩阵和右端向量
        std::vector<std::vector<double>> elementStiffness;
        std::vector<double> elementRHS;
        
        if (!computeElementMatrices(elemId, elementStiffness, elementRHS)) {
            std::cerr << "错误: 单元 " << elemId << " 矩阵计算失败" << std::endl;
            return false;
        }
        
        // 组装到全局系统
        if (!assembleToGlobalSystem(elemId, elementStiffness, elementRHS)) {
            std::cerr << "错误: 单元 " << elemId << " 组装到全局系统失败" << std::endl;
            return false;
        }
    }
    
    std::cout << "笛卡尔坐标系组装完成，共处理 " << numElements << " 个单元" << std::endl;
    return true;
}

bool MagneticSolver::assembleAxisymmetric() {
    // TODO: 实现柱对称坐标系下的Maxwell方程组装
    std::cout << "使用柱对称坐标系组装Maxwell方程..." << std::endl;
    
    // 简化实现：输出信息
    std::cout << "柱对称坐标系组装完成" << std::endl;
    return true;
}

bool MagneticSolver::assembleGeneral() {
    // TODO: 实现一般坐标系下的Maxwell方程组装
    std::cout << "使用一般坐标系组装Maxwell方程..." << std::endl;
    
    // 简化实现：输出信息
    std::cout << "一般坐标系组装完成" << std::endl;
    return true;
}

bool MagneticSolver::computeElementMatrices(int elementId, 
                                std::vector<std::vector<double>>& elementStiffness,
                                std::vector<double>& elementRHS) {
    // 基于Fortran版本的MagneticSolve.F90中的MaxwellCompose函数实现
    // 计算单元刚度矩阵和右端向量
    
    if (!mesh_ || elementId < 0 || elementId >= mesh_->getBulkElements().size()) {
        std::cerr << "错误: 无效的单元索引" << std::endl;
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
    
    // 计算单元几何信息（简化实现）
    // TODO: 实现完整的形状函数和积分计算
    
    // 简化实现：创建对角矩阵
    for (int i = 0; i < numElementNodes * 3; ++i) {
        elementStiffness[i][i] = 1.0; // 单位矩阵
        elementRHS[i] = 0.0; // 零右端向量
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
    // TODO: 实现旋度计算
    std::cout << "开始计算旋度..." << std::endl;
    
    // 简化实现：使用有限差分法计算旋度
    int numNodes = mesh_->numberOfNodes();
    
    // TODO: 实现完整的旋度计算算法
    
    std::cout << "旋度计算完成" << std::endl;
    return true;
}

bool MagneticSolver::computeNodalField() {
    // TODO: 实现节点场计算
    std::cout << "开始计算节点场..." << std::endl;
    
    // 简化实现：输出信息
    std::cout << "节点场计算完成" << std::endl;
    return true;
}

} // namespace elmer