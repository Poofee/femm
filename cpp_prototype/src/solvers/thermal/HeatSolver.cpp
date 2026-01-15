/**
 * @file HeatSolver.cpp
 * @brief Elmer FEM热传导求解器实现
 * 
 * 实现热传导方程的有限元求解，支持稳态和瞬态热分析
 */

#include "HeatSolver.h"
#include "Mesh.h"
#include "MaterialDatabase.h"
#include "BoundaryConditions.h"
#include "CRSMatrix.h"
#include "LinearAlgebra.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

namespace elmer {

// ===== 构造函数和析构函数 =====

HeatSolver::HeatSolver() 
    : LinearSolverBase("HeatSolver") {
    // 设置默认参数
    parameters_.solverName = "HeatSolver";
    parameters_.tolerance = 1.0e-6;
    parameters_.maxIterations = 1000;
    parameters_.linear = true;
    parameters_.transient = false;
    
    // 设置默认热传导参数
    heatParams_.thermalConductivity = 1.0;
    heatParams_.density = 1.0;
    heatParams_.specificHeat = 1.0;
    heatParams_.heatSource = 0.0;
    heatParams_.initialTemperature = 293.15;
    heatParams_.ambientTemperature = 293.15;
    heatParams_.heatTransferCoefficient = 0.0;
}

// ===== 参数设置函数 =====

void HeatSolver::setHeatParameters(const HeatSolverParameters& params) {
    heatParams_ = params;
}

HeatSolverParameters HeatSolver::getHeatParameters() const {
    return heatParams_;
}

void HeatSolver::setDirichletBoundary(const std::vector<int>& nodes, const std::vector<double>& values) {
    if (nodes.size() != values.size()) {
        throw std::invalid_argument("狄利克雷边界条件节点数和值数不匹配");
    }
    dirichletNodes_ = nodes;
    dirichletValues_ = values;
}

void HeatSolver::setNeumannBoundary(const std::vector<int>& edges, const std::vector<double>& values) {
    if (edges.size() != values.size()) {
        throw std::invalid_argument("诺伊曼边界条件边数和值数不匹配");
    }
    neumannEdges_ = edges;
    neumannValues_ = values;
}

void HeatSolver::setRobinBoundary(const std::vector<int>& edges, const std::vector<double>& coefficients, 
                                 const std::vector<double>& ambientTemps) {
    if (edges.size() != coefficients.size() || edges.size() != ambientTemps.size()) {
        throw std::invalid_argument("罗宾边界条件参数数量不匹配");
    }
    robinEdges_ = edges;
    robinCoefficients_ = coefficients;
    robinAmbientTemps_ = ambientTemps;
}

// ===== 求解器核心函数 =====

bool HeatSolver::initialize() {
    if (!mesh_) {
        std::cerr << "错误: 未设置网格数据" << std::endl;
        return false;
    }
    
    // 初始化温度场
    size_t numNodes = mesh_->numberOfNodes();
    temperatureField_.resize(numNodes, heatParams_.initialTemperature);
    heatFluxField_.resize(numNodes * 3, 0.0); // 每个节点3个分量
    
    // 初始化上一时间步温度场（用于瞬态分析）
    if (parameters_.transient) {
        prevTemperature_ = temperatureField_;
    }
    
    // 初始化矩阵和向量
    stiffnessMatrix_ = std::make_shared<CRSMatrix>(numNodes, numNodes);
    rhsVector_ = Vector::Create(numNodes);
    solution_ = Vector::Create(numNodes);
    
    // 设置初始解
    for (int i = 0; i < numNodes; ++i) {
        (*solution_)[i] = heatParams_.initialTemperature;
    }
    
    status_ = SolverStatus::INITIALIZED;
    std::cout << "HeatSolver初始化完成，节点数: " << numNodes << std::endl;
    
    return true;
}

bool HeatSolver::assemble() {
    if (status_ != SolverStatus::INITIALIZED) {
        std::cerr << "错误: 求解器未初始化" << std::endl;
        return false;
    }
    
    std::cout << "开始组装热传导系统..." << std::endl;
    
    // 组装刚度矩阵
    if (!assembleStiffnessMatrix()) {
        std::cerr << "错误: 刚度矩阵组装失败" << std::endl;
        return false;
    }
    
    // 如果是瞬态分析，组装质量矩阵
    if (parameters_.transient) {
        if (!assembleMassMatrix()) {
            std::cerr << "错误: 质量矩阵组装失败" << std::endl;
            return false;
        }
    }
    
    // 组装右端向量
    if (!assembleRhsVector()) {
        std::cerr << "错误: 右端向量组装失败" << std::endl;
        return false;
    }
    
    // 应用边界条件
    if (!applyBoundaryConditions()) {
        std::cerr << "错误: 边界条件应用失败" << std::endl;
        return false;
    }
    
    status_ = SolverStatus::ASSEMBLED;
    std::cout << "热传导系统组装完成" << std::endl;
    
    return true;
}

bool HeatSolver::solve() {
    if (status_ != SolverStatus::ASSEMBLED) {
        std::cerr << "错误: 系统未组装" << std::endl;
        return false;
    }
    
    std::cout << "开始求解热传导系统..." << std::endl;
    
    // 简化求解：直接使用矩阵求解（实际应该使用线性求解器）
    // TODO: 实现真正的线性求解器
    
    // 模拟求解过程
    size_t numNodes = mesh_->numberOfNodes();
    for (size_t i = 0; i < numNodes; ++i) {
        // 简化处理：直接使用右端向量作为解
        (*solution_)[i] = (*rhsVector_)[i];
    }
    
    // 更新温度场
    for (int i = 0; i < numNodes; ++i) {
        temperatureField_[i] = (*solution_)[i];
    }
    
    // 计算热通量
    computeHeatFlux();
    
    status_ = SolverStatus::SOLVED;
    std::cout << "热传导系统求解完成" << std::endl;
    
    return true;
}

std::vector<double> HeatSolver::getSolution() const {
    return temperatureField_;
}

std::vector<double> HeatSolver::getHeatFlux() const {
    return heatFluxField_;
}

// ===== 辅助函数 =====

double HeatSolver::getMaxTemperature() const {
    if (temperatureField_.empty()) {
        return 0.0;
    }
    return *std::max_element(temperatureField_.begin(), temperatureField_.end());
}

double HeatSolver::getMinTemperature() const {
    if (temperatureField_.empty()) {
        return 0.0;
    }
    return *std::min_element(temperatureField_.begin(), temperatureField_.end());
}

double HeatSolver::getAverageTemperature() const {
    if (temperatureField_.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double temp : temperatureField_) {
        sum += temp;
    }
    return sum / temperatureField_.size();
}

bool HeatSolver::checkConvergence() const {
    // 简化收敛检查
    double residual = getResidual();
    return residual < parameters_.tolerance;
}

double HeatSolver::getResidual() const {
    // 简化残差计算
    if (temperatureField_.empty()) {
        return 1.0;
    }
    
    // 计算温度场的标准差作为残差
    double mean = getAverageTemperature();
    double sumSq = 0.0;
    for (double temp : temperatureField_) {
        double diff = temp - mean;
        sumSq += diff * diff;
    }
    return std::sqrt(sumSq / temperatureField_.size());
}

bool HeatSolver::executeTimeStep(int timeStepIndex, double currentTime) {
    if (!parameters_.transient) {
        std::cerr << "错误: 求解器不支持瞬态分析" << std::endl;
        return false;
    }
    
    std::cout << "执行热传导时间步: " << timeStepIndex << ", 时间: " << currentTime << std::endl;
    
    // 保存上一时间步温度场
    prevTemperature_ = temperatureField_;
    
    // 更新时间积分因子
    timeIntegrationFactor_ = 1.0 / parameters_.timeStep;
    
    // 重新组装和求解系统
    if (!assemble()) {
        return false;
    }
    
    if (!solve()) {
        return false;
    }
    
    return true;
}

// ===== 私有辅助函数 =====

bool HeatSolver::assembleStiffnessMatrix() {
    if (!mesh_) {
        return false;
    }
    
    size_t numElements = mesh_->getBulkElements().size();
    size_t numNodes = mesh_->numberOfNodes();
    
    // 清零刚度矩阵
    stiffnessMatrix_->Zero();
    
    // 遍历所有单元
    for (size_t elemId = 0; elemId < numElements; ++elemId) {
        // 计算单元刚度矩阵
        std::vector<std::vector<double>> elementMatrix;
        computeElementMatrix(elemId, elementMatrix);
        
        // 获取单元节点编号
        auto& bulkElements = mesh_->getBulkElements();
        if (elemId < bulkElements.size()) {
            auto elementNodes = bulkElements[elemId].getNodeIndices();
            size_t numElementNodes = elementNodes.size();
            
            // 组装到全局刚度矩阵
            for (size_t i = 0; i < numElementNodes; ++i) {
                for (size_t j = 0; j < numElementNodes; ++j) {
                    int globalI = static_cast<int>(elementNodes[i]);
                    int globalJ = static_cast<int>(elementNodes[j]);
                    // 使用Matrix类的正确接口
                    stiffnessMatrix_->AddToElement(globalI, globalJ, elementMatrix[i][j]);
                }
            }
        }
    }
    
    return true;
}

bool HeatSolver::assembleMassMatrix() {
    // 简化实现：质量矩阵为单位矩阵乘以系数
    if (!mesh_) {
        return false;
    }
    
    size_t numNodes = mesh_->numberOfNodes();
    double massCoeff = heatParams_.density * heatParams_.specificHeat * timeIntegrationFactor_;
    
    for (size_t i = 0; i < numNodes; ++i) {
        stiffnessMatrix_->AddToElement(i, i, massCoeff);
    }
    
    return true;
}

bool HeatSolver::assembleRhsVector() {
    if (!mesh_) {
        return false;
    }
    
    size_t numNodes = mesh_->numberOfNodes();
    rhsVector_->Zero();
    
    // 添加热源项
    for (size_t i = 0; i < numNodes; ++i) {
        (*rhsVector_)[i] += heatParams_.heatSource;
    }
    
    // 如果是瞬态分析，添加上一时间步的贡献
    if (parameters_.transient && !prevTemperature_.empty()) {
        double massCoeff = heatParams_.density * heatParams_.specificHeat * timeIntegrationFactor_;
        for (int i = 0; i < numNodes; ++i) {
            (*rhsVector_)[i] += massCoeff * prevTemperature_[i];
        }
    }
    
    return true;
}

bool HeatSolver::applyBoundaryConditions() {
    if (!mesh_) {
        return false;
    }
    
    size_t numNodes = mesh_->numberOfNodes();
    
    // 应用狄利克雷边界条件
    for (size_t i = 0; i < dirichletNodes_.size(); ++i) {
        int nodeId = dirichletNodes_[i];
        if (nodeId >= 0 && nodeId < numNodes) {
            // 设置对角线元素为1，右端向量为边界值
            // 简化实现：只设置对角线元素
            for (int j = 0; j < numNodes; ++j) {
                if (j == nodeId) {
                    stiffnessMatrix_->SetElement(nodeId, j, 1.0);
                } else {
                    stiffnessMatrix_->SetElement(nodeId, j, 0.0);
                }
            }
            (*rhsVector_)[nodeId] = dirichletValues_[i];
        }
    }
    
    return true;
}

void HeatSolver::computeHeatFlux() {
    // 简化热通量计算
    if (!mesh_) {
        return;
    }
    
    size_t numNodes = mesh_->numberOfNodes();
    
    // 计算温度梯度（简化实现）
    for (size_t i = 0; i < numNodes; ++i) {
        // 简化处理：假设均匀温度梯度
        double gradX = 0.0;
        double gradY = 0.0;
        double gradZ = 0.0;
        
        // 根据傅里叶定律计算热通量
        heatFluxField_[i * 3] = -heatParams_.thermalConductivity * gradX;
        heatFluxField_[i * 3 + 1] = -heatParams_.thermalConductivity * gradY;
        heatFluxField_[i * 3 + 2] = -heatParams_.thermalConductivity * gradZ;
    }
}

void HeatSolver::computeElementMatrix(int elementId, std::vector<std::vector<double>>& elementMatrix) const {
    // 简化单元矩阵计算
    auto& bulkElements = mesh_->getBulkElements();
    if (elementId < bulkElements.size()) {
        auto elementNodes = bulkElements[elementId].getNodeIndices();
        size_t numElementNodes = elementNodes.size();
        
        elementMatrix.resize(numElementNodes, std::vector<double>(numElementNodes, 0.0));
        
        // 简化处理：使用单位矩阵乘以系数
        double coeff = heatParams_.thermalConductivity;
        for (size_t i = 0; i < numElementNodes; ++i) {
            elementMatrix[i][i] = coeff;
        }
    }
}

void HeatSolver::computeElementMassMatrix(int elementId, std::vector<std::vector<double>>& elementMatrix) const {
    // 简化质量矩阵计算
    auto& bulkElements = mesh_->getBulkElements();
    if (elementId < bulkElements.size()) {
        auto elementNodes = bulkElements[elementId].getNodeIndices();
        size_t numElementNodes = elementNodes.size();
        
        elementMatrix.resize(numElementNodes, std::vector<double>(numElementNodes, 0.0));
        
        // 简化处理：使用单位矩阵乘以系数
        double coeff = heatParams_.density * heatParams_.specificHeat;
        for (int i = 0; i < numElementNodes; ++i) {
            elementMatrix[i][i] = coeff;
        }
    }
}

void HeatSolver::computeElementRhsVector(int elementId, std::vector<double>& elementVector) const {
    // 简化右端向量计算
    auto& bulkElements = mesh_->getBulkElements();
    if (elementId < bulkElements.size()) {
        auto elementNodes = bulkElements[elementId].getNodeIndices();
        size_t numElementNodes = elementNodes.size();
        
        elementVector.resize(numElementNodes, 0.0);
        
        // 简化处理：均匀分布热源
        for (size_t i = 0; i < numElementNodes; ++i) {
            elementVector[i] = heatParams_.heatSource / numElementNodes;
        }
    }
}

void HeatSolver::computeShapeFunctions(double xi, double eta, std::vector<double>& N, 
                                      std::vector<std::vector<double>>& dN) const {
    // 简化形函数计算（线性三角形单元）
    N.resize(3);
    dN.resize(3, std::vector<double>(2));
    
    N[0] = 1.0 - xi - eta;
    N[1] = xi;
    N[2] = eta;
    
    dN[0][0] = -1.0; dN[0][1] = -1.0;
    dN[1][0] = 1.0;  dN[1][1] = 0.0;
    dN[2][0] = 0.0;  dN[2][1] = 1.0;
}

// ===== 状态保存和加载 =====

bool HeatSolver::saveState(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    file << "HeatSolver State" << std::endl;
    file << "TemperatureField:" << std::endl;
    for (double temp : temperatureField_) {
        file << temp << " ";
    }
    file << std::endl;
    
    file.close();
    return true;
}

bool HeatSolver::loadState(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    std::getline(file, line); // 跳过标题行
    
    std::getline(file, line); // 读取温度场标题
    if (line != "TemperatureField:") {
        return false;
    }
    
    std::getline(file, line);
    std::istringstream iss(line);
    temperatureField_.clear();
    
    double temp;
    while (iss >> temp) {
        temperatureField_.push_back(temp);
    }
    
    file.close();
    return true;
}

} // namespace elmer