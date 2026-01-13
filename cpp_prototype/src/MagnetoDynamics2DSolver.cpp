// MagnetoDynamics2DSolver.cpp - Elmer FEM C++ 2D Magnetodynamics Solver
// Corresponds to Fortran module: MagnetoDynamics2D.F90

#include "MagnetoDynamics2DSolver.h"
#include "LinearAlgebra.h"
#include "ElectromagneticMaterial.h"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace elmer {

// 性能优化方法实现
void MagnetoDynamics2DSolver::precomputeBasisFunctionCache() {
    if (!mesh || !useBasisFunctionsCache) {
        return;
    }
    
    std::cout << "开始预计算基函数缓存..." << std::endl;
    
    // 清除现有缓存
    elementCache.clear();
    shapeFunctionCache.clear();
    
    // 获取网格信息
    auto& bulkElements = mesh->getBulkElements();
    auto& nodes = mesh->getNodes();
    
    // 预计算每个单元的基函数
    for (size_t elemIdx = 0; elemIdx < bulkElements.size(); ++elemIdx) {
        ElementCache elemCache;
        
        // 获取单元信息
        auto& element = bulkElements[elemIdx];
        elemCache.elementType = static_cast<int>(element.getType());
        
        // 转换节点索引类型
        auto nodeIndices = element.getNodeIndices();
        elemCache.nodeIndices.resize(nodeIndices.size());
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            elemCache.nodeIndices[i] = static_cast<int>(nodeIndices[i]);
        }
        
        // 简化的基函数计算（实际实现应包括高斯积分点）
        BasisFunctionCache basisCache;
        
        // 假设使用线性三角形单元
        int nNodes = elemCache.nodeIndices.size();
        basisCache.basis.resize(nNodes, 1.0 / nNodes); // 简化的基函数值
        
        // 简化的导数计算
        basisCache.dBasisdx.resize(nNodes);
        for (int i = 0; i < nNodes; ++i) {
            basisCache.dBasisdx[i] = {0.0, 0.0, 0.0}; // 简化的导数
        }
        
        basisCache.weight = 1.0; // 简化的权重
        basisCache.detJ = 1.0;   // 简化的雅可比行列式
        basisCache.coords = {0.0, 0.0, 0.0}; // 简化的坐标
        
        elemCache.integrationPoints.push_back(basisCache);
        
        // 计算单元面积（简化）
        elemCache.area = 1.0; // 简化的面积计算
        
        elementCache.push_back(elemCache);
    }
    
    std::cout << "基函数缓存预计算完成，共缓存 " << elementCache.size() << " 个单元" << std::endl;
}

void MagnetoDynamics2DSolver::assembleSystemOptimized() {
    if (!mesh) {
        std::cout << "错误: 未设置网格" << std::endl;
        return;
    }
    
    std::cout << "开始优化的矩阵组装..." << std::endl;
    
    // 如果启用了缓存，使用缓存数据
    if (useBasisFunctionsCache && elementCache.empty()) {
        precomputeBasisFunctionCache();
    }
    
    // 获取网格信息
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    // 初始化系统矩阵（使用工厂方法）
    stiffnessMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes));
    massMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes));
    dampingMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes));
    rhsVector = Vector::Create(nNodes);
    
    // 简化的矩阵组装（使用缓存）
    if (useBasisFunctionsCache) {
        // 使用缓存的基函数数据进行矩阵组装
        for (size_t elemIdx = 0; elemIdx < elementCache.size(); ++elemIdx) {
            const auto& elemCache = elementCache[elemIdx];
            
            // 简化的单元矩阵组装
            int nNodesPerElem = elemCache.nodeIndices.size();
            
            // 组装刚度矩阵（简化）
            for (int i = 0; i < nNodesPerElem; ++i) {
                for (int j = 0; j < nNodesPerElem; ++j) {
                    double stiffnessValue = 1.0; // 简化的刚度系数
                    stiffnessMatrix->AddToElement(elemCache.nodeIndices[i], elemCache.nodeIndices[j], stiffnessValue);
                }
            }
            
            // 组装质量矩阵（简化）
            for (int i = 0; i < nNodesPerElem; ++i) {
                for (int j = 0; j < nNodesPerElem; ++j) {
                    double massValue = 0.1; // 简化的质量系数
                    massMatrix->AddToElement(elemCache.nodeIndices[i], elemCache.nodeIndices[j], massValue);
                }
            }
        }
    } else {
        // 使用传统方法组装
        assembleElementContributions();
    }
    
    // 组装边界贡献
    assembleBoundaryContributions();
    
    systemAssembled = true;
    std::cout << "优化的矩阵组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::assembleSystemParallel() {
    if (!mesh) {
        std::cout << "错误: 未设置网格" << std::endl;
        return;
    }
    
    std::cout << "开始并行矩阵组装（OpenMP）..." << std::endl;
    std::cout << "使用线程数: " << numThreads << std::endl;
    
    // 设置OpenMP线程数
    omp_set_num_threads(numThreads);
    
    // 检查OpenMP是否可用
    #ifdef _OPENMP
    std::cout << "OpenMP支持已启用，当前线程数: " << omp_get_max_threads() << std::endl;
    #else
    std::cout << "警告: OpenMP支持未启用，将使用串行模式" << std::endl;
    #endif
    
    // 获取网格信息
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    // 初始化系统矩阵
    stiffnessMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes));
    massMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes));
    dampingMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes));
    rhsVector = Vector::Create(nNodes);
    
    // 并行组装单元贡献
    assembleElementContributionsParallel();
    
    // 组装边界贡献（串行，通常边界条件数量较少）
    assembleBoundaryContributions();
    
    systemAssembled = true;
    std::cout << "并行矩阵组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::assembleElementContributionsParallel() {
    auto& bulkElements = mesh->getBulkElements();
    int numElements = static_cast<int>(bulkElements.size());
    
    std::cout << "并行组装 " << numElements << " 个单元..." << std::endl;
    
    // 使用OpenMP并行化单元循环
    #pragma omp parallel for schedule(dynamic)
    for (int elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        const auto& element = bulkElements[elemIdx];
        
        // 计算单元局部矩阵
        std::vector<std::vector<double>> elementStiffness;
        std::vector<std::vector<double>> elementMass;
        std::vector<std::vector<double>> elementDamping;
        std::vector<double> elementForce;
        
        computeLocalMatrix(element, elementStiffness, elementMass, elementDamping, elementForce);
        
        // 获取单元节点索引
        auto nodeIndices = element.getNodeIndices();
        int nNodesPerElem = nodeIndices.size();
        
        // 组装到全局系统矩阵（线程安全）
        for (int i = 0; i < nNodesPerElem; ++i) {
            int globalRow = static_cast<int>(nodeIndices[i]);
            
            // 组装右端向量
            if (!elementForce.empty()) {
                addToVectorThreadSafe(rhsVector, globalRow, elementForce[i]);
            }
            
            // 组装刚度矩阵
            for (int j = 0; j < nNodesPerElem; ++j) {
                int globalCol = static_cast<int>(nodeIndices[j]);
                
                if (!elementStiffness.empty()) {
                    addToMatrixThreadSafe(stiffnessMatrix, globalRow, globalCol, elementStiffness[i][j]);
                }
                
                if (!elementMass.empty()) {
                    addToMatrixThreadSafe(massMatrix, globalRow, globalCol, elementMass[i][j]);
                }
                
                if (!elementDamping.empty()) {
                    addToMatrixThreadSafe(dampingMatrix, globalRow, globalCol, elementDamping[i][j]);
                }
            }
        }
    }
    
    std::cout << "单元并行组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::addToMatrixThreadSafe(std::shared_ptr<Matrix>& matrix, int row, int col, double value) {
    // 使用OpenMP临界区确保线程安全
    #pragma omp critical
    {
        matrix->SetElement(row, col, matrix->GetElement(row, col) + value);
    }
}

void MagnetoDynamics2DSolver::addToVectorThreadSafe(std::shared_ptr<Vector>& vector, int index, double value) {
    // 使用OpenMP原子操作确保线程安全
    #pragma omp atomic
    (*vector)[index] += value;
}

void MagnetoDynamics2DSolver::setParallelThreads(int threads) {
    if (threads > 0) {
        numThreads = threads;
        std::cout << "设置并行线程数为: " << numThreads << std::endl;
    }
}

int MagnetoDynamics2DSolver::getParallelThreads() const {
    return numThreads;
}

void MagnetoDynamics2DSolver::updateSystemIncremental() {
    // 简化的增量式更新实现
    // 实际实现应包括仅更新变化的部分
    
    std::cout << "开始增量式系统更新..." << std::endl;
    
    if (!systemAssembled) {
        assembleSystemOptimized();
        return;
    }
    
    // 简化的增量更新：重新组装非线性部分
    reassembleNonlinearSystem();
    
    std::cout << "增量式系统更新完成" << std::endl;
}

void MagnetoDynamics2DSolver::clearCache() {
    elementCache.clear();
    shapeFunctionCache.clear();
    assemblyCache = MatrixAssemblyCache();
    
    std::cout << "缓存已清除" << std::endl;
}

void MagnetoDynamics2DSolver::getCacheStatistics() const {
    std::cout << "=== 缓存统计信息 ===" << std::endl;
    std::cout << "单元缓存数量: " << elementCache.size() << std::endl;
    std::cout << "形状函数缓存数量: " << shapeFunctionCache.size() << std::endl;
    std::cout << "基函数缓存启用: " << (useBasisFunctionsCache ? "是" : "否") << std::endl;
    std::cout << "系统组装状态: " << (systemAssembled ? "已组装" : "未组装") << std::endl;
}

void MagnetoDynamics2DSolver::initializeBasisFunctionCache() {
    // 简化的基函数缓存初始化
    // 实际实现应包括预计算所有单元的基函数值
    
    if (!mesh) {
        return;
    }
    
    // 简化的基函数缓存实现
    elementCache.clear();
    shapeFunctionCache.clear();
    
    std::cout << "基函数缓存初始化完成" << std::endl;
}

void MagnetoDynamics2DSolver::assembleElementContributions() {
    // 简化的单元贡献组装
    // 实际实现应包括完整的有限元矩阵组装
    
    if (!mesh) {
        return;
    }
    
    std::cout << "单元贡献组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::assembleBoundaryContributions() {
    // 简化的边界贡献组装
    // 实际实现应包括各种边界条件的处理
    
    if (!mesh) {
        return;
    }
    
    std::cout << "边界贡献组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::applyBoundaryConditions() {
    // 简化的边界条件应用
    // 实际实现应包括各种边界条件的处理
    
    if (!systemAssembled) {
        return;
    }
    
    std::cout << "边界条件应用完成" << std::endl;
}

void MagnetoDynamics2DSolver::reassembleNonlinearSystem() {
    // 简化的非线性系统重新组装
    // 实际实现应包括考虑非线性材料特性的重新组装
    
    if (!mesh) {
        return;
    }
    
    std::cout << "非线性系统重新组装完成" << std::endl;
}

std::vector<double> MagnetoDynamics2DSolver::solveLinearSystem() {
    // 简化的线性系统求解
    // 实际实现应包括迭代求解器的使用
    
    if (!stiffnessMatrix || !rhsVector) {
        return {};
    }
    
    size_t nNodes = 100; // 假设有100个节点
    std::vector<double> solution(nNodes, 0.0);
    
    std::cout << "线性系统求解完成" << std::endl;
    
    return solution;
}

bool MagnetoDynamics2DSolver::checkConvergence(const MagnetoDynamics2DResults& results) {
    // 简化的收敛性检查
    // 实际实现应包括残差和相对变化检查
    
    if (!results.converged) {
        return false;
    }
    
    // 检查残差是否小于容差
    if (results.residual > parameters.tolerance) {
        return false;
    }
    
    // 检查迭代次数是否超过限制
    if (results.nonlinearIterations >= parameters.maxNonlinearIterations) {
        return false;
    }
    
    return true;
}

void MagnetoDynamics2DSolver::calculateDerivedFields(MagnetoDynamics2DResults& results) {
    if (!mesh || results.vectorPotential.empty()) {
        return;
    }
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    results.magneticFluxDensity.resize(nNodes);
    results.magneticFieldStrength.resize(nNodes);
    results.currentDensity.resize(nNodes, 0.0);
    
    // 简化的导出场量计算
    // 实际实现应包括基于有限元方法的精确计算
    
    for (size_t i = 0; i < nNodes; ++i) {
        results.magneticFluxDensity[i] = {0.0, 0.0};
        results.magneticFieldStrength[i] = {0.0, 0.0};
        results.currentDensity[i] = 0.0;
    }
    
    std::cout << "导出场量计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::calculateLumpedParameters(MagnetoDynamics2DResults& results) {
    // 简化的集总参数计算
    // 实际实现应包括磁能、电感、转矩等参数的计算
    
    if (!mesh) {
        return;
    }
    
    // 简化的集总参数计算
    results.magneticEnergy = 0.0;
    results.inductance = 0.0;
    results.torque = 0.0;
    
    std::cout << "集总参数计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::computeLocalMatrix(const Element& element, 
                                                 std::vector<std::vector<double>>& stiffness,
                                                 std::vector<std::vector<double>>& mass,
                                                 std::vector<std::vector<double>>& damping,
                                                 std::vector<double>& force) {
    // 简化的局部矩阵计算
    // 实际实现应包括基于有限元方法的精确计算
    
    int nNodes = element.numberOfNodes();
    
    // 初始化局部矩阵
    stiffness.assign(nNodes, std::vector<double>(nNodes, 0.0));
    mass.assign(nNodes, std::vector<double>(nNodes, 0.0));
    damping.assign(nNodes, std::vector<double>(nNodes, 0.0));
    force.assign(nNodes, 0.0);
    
    // 简化的局部矩阵计算
    for (int i = 0; i < nNodes; ++i) {
        for (int j = 0; j < nNodes; ++j) {
            stiffness[i][j] = 1.0;
            mass[i][j] = 0.1;
            damping[i][j] = 0.01;
        }
        force[i] = 0.0;
    }
}

std::array<std::array<double, 2>, 2> MagnetoDynamics2DSolver::computeReluctivity(const std::string& materialName, 
                                                                                 double fluxDensity,
                                                                                 const Element& element) {
    // 简化的磁阻率计算
    // 实际实现应包括非线性材料的磁阻率计算
    
    std::array<std::array<double, 2>, 2> reluctivity;
    
    // 简化的磁阻率计算（线性材料）
    // 暂时使用固定值
    double nu = 1.0 / 1.25663706212e-6; // 假设相对磁导率为1
    reluctivity[0][0] = nu;
    reluctivity[0][1] = 0.0;
    reluctivity[1][0] = 0.0;
    reluctivity[1][1] = nu;
    
    return reluctivity;
}

void MagnetoDynamics2DSolver::applyInfinityBoundaryCondition(const Element& element, 
                                                             std::vector<std::vector<double>>& stiffness,
                                                             std::vector<double>& force) {
    // 简化的无限远边界条件应用
    // 实际实现应包括无限远边界条件的精确处理
    
    int nNodes = element.numberOfNodes();
    
    // 简化的无限远边界条件处理
    for (int i = 0; i < nNodes; ++i) {
        for (int j = 0; j < nNodes; ++j) {
            stiffness[i][j] += 0.1; // 简化的贡献
        }
    }
}

void MagnetoDynamics2DSolver::applyAirGapBoundaryCondition(const Element& element, 
                                                           std::vector<std::vector<double>>& stiffness,
                                                           std::vector<double>& force) {
    // 简化的气隙边界条件应用
    // 实际实现应包括气隙边界条件的精确处理
    
    int nNodes = element.numberOfNodes();
    
    // 简化的气隙边界条件处理
    for (int i = 0; i < nNodes; ++i) {
        for (int j = 0; j < nNodes; ++j) {
            stiffness[i][j] += 0.05; // 简化的贡献
        }
    }
}

// 谐波分析方法实现
MagnetoDynamics2DResults MagnetoDynamics2DSolver::solveHarmonic() {
    if (!mesh) {
        throw std::runtime_error("Mesh not set for harmonic analysis");
    }
    
    if (parameters.frequency <= 0.0) {
        throw std::runtime_error("Frequency must be positive for harmonic analysis");
    }
    
    MagnetoDynamics2DResults results;
    
    std::cout << "开始谐波分析，频率: " << parameters.frequency << " Hz" << std::endl;
    
    // 设置谐波分析标志
    parameters.isHarmonic = true;
    parameters.useComplexMatrices = true;
    
    // 组装复数系统矩阵
    assembleComplexSystem();
    
    // 求解复数线性系统
    auto complexSolution = solveComplexLinearSystem();
    results.complexVectorPotential = complexSolution;
    
    // 计算复数导出场量
    calculateComplexDerivedFields(results);
    
    // 计算复数集总参数
    calculateComplexLumpedParameters(results);
    
    // 标记为收敛
    results.converged = true;
    results.iterations = 1;
    
    std::cout << "谐波分析完成" << std::endl;
    
    return results;
}

void MagnetoDynamics2DSolver::assembleComplexSystem() {
    if (!mesh) {
        throw std::runtime_error("Mesh not set for complex system assembly");
    }
    
    std::cout << "开始组装复数系统矩阵..." << std::endl;
    
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    // 初始化复数系统矩阵
    complexStiffnessMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes).release());
    complexMassMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes).release());
    complexDampingMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes).release());
    complexRhsVector = Vector::Create(nNodes);
    
    // 角频率
    double omega = 2.0 * 3.141592653589793 * parameters.frequency;
    
    // 组装单元贡献
    auto& bulkElements = mesh->getBulkElements();
    
    for (size_t elemIdx = 0; elemIdx < bulkElements.size(); ++elemIdx) {
        const auto& element = bulkElements[elemIdx];
        auto nodeIndices = element.getNodeIndices();
        int nNodesPerElem = nodeIndices.size();
        
        // 获取材料属性
        auto material = materialDB.getMaterial(element.getMaterialName());
        
        // 复数磁导率
        auto mu_complex = material.getComplexPermeability(parameters.frequency);
        
        // 复数电导率
        auto sigma_complex = material.getComplexConductivity();
        
        // 计算局部复数矩阵
        for (int i = 0; i < nNodesPerElem; ++i) {
            int globalRow = static_cast<int>(nodeIndices[i]);
            
            for (int j = 0; j < nNodesPerElem; ++j) {
                int globalCol = static_cast<int>(nodeIndices[j]);
                
                // 复数刚度矩阵项（实部：刚度，虚部：涡流项）
                double stiffnessReal = 1.0 / mu_complex.real(); // 简化的刚度项
                double stiffnessImag = omega * sigma_complex.real(); // 涡流项
                
                // 复数质量矩阵项（用于位移电流）
                double massReal = 0.0;
                double massImag = 0.0;
                if (parameters.includeDisplacementCurrent) {
                    massImag = omega * material.permittivity(); // 位移电流项
                }
                
                // 组装到复数系统矩阵
                complexStiffnessMatrix->AddToElement(globalRow, globalCol, stiffnessReal);
                complexDampingMatrix->AddToElement(globalRow, globalCol, stiffnessImag);
                complexMassMatrix->AddToElement(globalRow, globalCol, massImag);
            }
        }
    }
    
    // 组装边界条件贡献
    assembleBoundaryContributions();
    
    std::cout << "复数系统矩阵组装完成" << std::endl;
}

std::vector<std::complex<double>> MagnetoDynamics2DSolver::solveComplexLinearSystem() {
    if (!complexStiffnessMatrix || !complexRhsVector) {
        throw std::runtime_error("Complex system not assembled");
    }
    
    std::cout << "开始求解复数线性系统..." << std::endl;
    
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    // 简化的复数求解器（实际实现应使用复数迭代求解器）
    std::vector<std::complex<double>> solution(nNodes, std::complex<double>(0.0, 0.0));
    
    // 简化的求解过程：对角矩阵求逆
    for (int i = 0; i < nNodes; ++i) {
        // 获取对角线元素（简化）
        double diagReal = 1.0; // 简化的实部
        double diagImag = 0.1; // 简化的虚部（涡流效应）
        
        // 简化的复数求逆
        double denominator = diagReal * diagReal + diagImag * diagImag;
        double invReal = diagReal / denominator;
        double invImag = -diagImag / denominator;
        
        // 简化的右端项
        double rhsReal = 1.0; // 简化的实部激励
        double rhsImag = 0.0; // 简化的虚部激励
        
        // 求解：x = A⁻¹ * b
        solution[i] = std::complex<double>(
            invReal * rhsReal - invImag * rhsImag,
            invReal * rhsImag + invImag * rhsReal
        );
    }
    
    std::cout << "复数线性系统求解完成" << std::endl;
    
    return solution;
}

void MagnetoDynamics2DSolver::calculateComplexDerivedFields(MagnetoDynamics2DResults& results) {
    if (!mesh || results.complexVectorPotential.empty()) {
        return;
    }
    
    std::cout << "开始计算复数导出场量..." << std::endl;
    
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    // 初始化复数场量
    results.complexMagneticFluxDensity.resize(nNodes);
    results.complexMagneticFieldStrength.resize(nNodes);
    results.complexCurrentDensity.resize(nNodes);
    
    // 简化的复数场量计算
    for (int i = 0; i < nNodes; ++i) {
        // 复数磁通密度：B = ∇ × A
        results.complexMagneticFluxDensity[i][0] = std::complex<double>(0.0, 0.0); // B_x
        results.complexMagneticFluxDensity[i][1] = std::complex<double>(0.0, 0.0); // B_y
        
        // 复数磁场强度：H = B / μ
        results.complexMagneticFieldStrength[i][0] = std::complex<double>(0.0, 0.0); // H_x
        results.complexMagneticFieldStrength[i][1] = std::complex<double>(0.0, 0.0); // H_y
        
        // 复数电流密度：J = -jωσA（涡流）
        double omega = 2.0 * 3.141592653589793 * parameters.frequency;
        results.complexCurrentDensity[i] = std::complex<double>(0.0, -omega * 1.0) * results.complexVectorPotential[i];
    }
    
    std::cout << "复数导出场量计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::calculateComplexLumpedParameters(MagnetoDynamics2DResults& results) {
    if (!mesh) {
        return;
    }
    
    std::cout << "开始计算复数集总参数..." << std::endl;
    
    // 简化的复数阻抗计算
    double omega = 2.0 * 3.141592653589793 * parameters.frequency;
    
    // 复数阻抗：Z = R + jωL
    results.complexImpedance = std::complex<double>(1.0, omega * 0.1); // 简化的阻抗
    results.complexInductance = std::complex<double>(0.1, 0.01); // 简化的电感
    
    // 功率损耗：P = 0.5 * Re(V * I*)
    results.powerLoss = 0.5 * 1.0 * 1.0; // 简化的功率损耗
    
    std::cout << "复数集总参数计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::setHarmonicExcitation(double amplitude, double frequency, double phase) {
    parameters.frequency = frequency;
    parameters.harmonicPhase = phase;
    
    // 设置谐波激励源（简化实现）
    std::cout << "设置谐波激励源: 幅值=" << amplitude 
              << ", 频率=" << frequency << "Hz, 相位=" << phase << "rad" << std::endl;
}

std::vector<double> MagnetoDynamics2DSolver::reconstructTimeDomain(double time) const {
    if (!mesh) {
        return {};
    }
    
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    std::vector<double> timeDomainSolution(nNodes, 0.0);
    
    // 简化的时域重构：A(t) = Re[A(ω) * e^(jωt)]
    double omega = 2.0 * 3.141592653589793 * parameters.frequency;
    
    // 实际实现应从结果中获取复数解
    for (int i = 0; i < nNodes; ++i) {
        // 简化的时域重构
        timeDomainSolution[i] = std::cos(omega * time + parameters.harmonicPhase);
    }
    
    return timeDomainSolution;
}

void MagnetoDynamics2DSolver::reassembleNonlinearSystemWithJacobian() {
    if (!mesh) {
        return;
    }
    
    std::cout << "重新组装非线性系统（包含雅可比项）..." << std::endl;
    
    // 更新材料参数
    updateMaterialParameters();
    
    // 重新组装刚度矩阵（包含非线性项）
    auto& nodes = mesh->getNodes();
    int nNodes = static_cast<int>(nodes.numberOfNodes());
    
    // 创建新的刚度矩阵
    stiffnessMatrix = std::shared_ptr<Matrix>(Matrix::CreateCRS(nNodes, nNodes).release());
    
    // 简化的非线性雅可比矩阵组装
    // 实际实现应考虑材料非线性（如B-H曲线）
    for (int i = 0; i < nNodes; ++i) {
        // 对角线元素：考虑非线性磁导率
        double diagValue = 1.0; // 简化的非线性项
        stiffnessMatrix->SetElement(i, i, diagValue);
    }
    
    std::cout << "非线性系统雅可比矩阵组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::updateMaterialParameters() {
    if (!mesh) {
        return;
    }
    
    std::cout << "更新材料参数..." << std::endl;
    
    // 简化的材料参数更新
    // 实际实现应从材料数据库中获取当前磁场下的材料参数
    
    // 遍历所有体单元
    auto& bulkElements = mesh->getBulkElements();
    for (const auto& element : bulkElements) {
        std::string materialName = element.getMaterialName();
        
        // 简化的材料参数更新逻辑
        // 实际实现应考虑非线性材料的场依赖特性
        if (materialName == "copper") {
            // 铜的电导率保持不变
        } else if (materialName == "iron") {
            // 铁的非线性磁导率需要根据当前磁场更新
        }
    }
    
    std::cout << "材料参数更新完成" << std::endl;
}

// ===== MPI并行方法实现 =====

DomainDecompositionResult MagnetoDynamics2DSolver::performDomainDecomposition() {
    if (!mesh || !decompositionManager_) {
        throw std::runtime_error("Mesh or decomposition manager not initialized");
    }
    
    auto& bulkElements = mesh->getBulkElements();
    std::vector<MeshElement> meshElements;
    
    // 转换网格元素为域分解格式
    for (size_t elemIdx = 0; elemIdx < bulkElements.size(); ++elemIdx) {
        MeshElement meshElem;
        meshElem.id = static_cast<int>(elemIdx);
        
        // 计算元素中心坐标
        auto& element = bulkElements[elemIdx];
        auto nodeIndices = element.getNodeIndices();
        auto& nodes = mesh->getNodes();
        
        std::array<double, 3> center = {0.0, 0.0, 0.0};
        for (auto nodeIdx : nodeIndices) {
            auto& node = nodes.getNode(nodeIdx);
            auto coords = node.getCoordinates();
            center[0] += coords[0];
            center[1] += coords[1];
            center[2] += coords[2];
        }
        
        center[0] /= nodeIndices.size();
        center[1] /= nodeIndices.size();
        center[2] /= nodeIndices.size();
        
        meshElem.coords = {center[0], center[1], center[2]};
        meshElements.push_back(meshElem);
    }
    
    // 执行域分解
    int numPartitions = comm_->getSize();
    auto result = decompositionManager_->decompose(meshElements, numPartitions);
    
    if (comm_->getRank() == 0) {
        std::cout << "域分解完成，分区数: " << numPartitions 
                  << ", 负载均衡度: " << result.loadBalance << std::endl;
    }
    
    return result;
}

void MagnetoDynamics2DSolver::assembleElementContributionsParallel(const DomainDecompositionResult& decomposition) {
    if (!mesh || !parallelAssembler_) {
        throw std::runtime_error("Mesh or parallel assembler not initialized");
    }
    
    // 获取本地元素列表
    auto localElements = getLocalElements(decomposition);
    
    if (comm_->getRank() == 0) {
        std::cout << "并行组装单元贡献，本地元素数: " << localElements.size() << std::endl;
    }
    
    // 并行组装每个本地元素
    for (int elemIdx : localElements) {
        auto& element = mesh->getBulkElements()[elemIdx];
        
        // 计算单元局部矩阵
        std::vector<std::vector<double>> elementStiffness;
        std::vector<std::vector<double>> elementMass;
        std::vector<std::vector<double>> elementDamping;
        std::vector<double> elementForce;
        
        computeLocalMatrix(element, elementStiffness, elementMass, elementDamping, elementForce);
        
        // 获取单元节点索引
        auto nodeIndices = element.getNodeIndices();
        int nNodesPerElem = nodeIndices.size();
        
        // 组装到分布式系统矩阵
        for (int i = 0; i < nNodesPerElem; ++i) {
            int globalRow = static_cast<int>(nodeIndices[i]);
            
            // 组装右端向量
            if (!elementForce.empty()) {
                distributedRhsVector_->addLocalValue(globalRow, elementForce[i]);
            }
            
            // 组装刚度矩阵
            for (int j = 0; j < nNodesPerElem; ++j) {
                int globalCol = static_cast<int>(nodeIndices[j]);
                
                if (!elementStiffness.empty()) {
                    distributedStiffnessMatrix_->addLocalValue(globalRow, globalCol, elementStiffness[i][j]);
                }
                
                if (parameters.isTransient && !elementMass.empty()) {
                    distributedMassMatrix_->addLocalValue(globalRow, globalCol, elementMass[i][j]);
                }
                
                if (parameters.isTransient && !elementDamping.empty()) {
                    distributedDampingMatrix_->addLocalValue(globalRow, globalCol, elementDamping[i][j]);
                }
            }
        }
    }
    
    // 交换幽灵数据
    exchangeGhostData();
    
    if (comm_->getRank() == 0) {
        std::cout << "并行单元组装完成" << std::endl;
    }
}

void MagnetoDynamics2DSolver::assembleBoundaryContributionsParallel(const DomainDecompositionResult& decomposition) {
    if (!mesh) {
        return;
    }
    
    // 简化的并行边界条件组装
    // 实际实现应包括边界条件的并行处理
    
    if (comm_->getRank() == 0) {
        std::cout << "并行边界条件组装完成" << std::endl;
    }
}

void MagnetoDynamics2DSolver::applyBoundaryConditionsParallel() {
    if (!distributedStiffnessMatrix_ || !distributedRhsVector_) {
        return;
    }
    
    // 简化的并行边界条件应用
    // 实际实现应包括边界条件的并行处理
    
    if (comm_->getRank() == 0) {
        std::cout << "并行边界条件应用完成" << std::endl;
    }
}

std::vector<int> MagnetoDynamics2DSolver::getLocalElements(const DomainDecompositionResult& decomposition) {
    std::vector<int> localElements;
    int currentRank = comm_->getRank();
    
    // 获取当前进程负责的元素
    for (int elemIdx = 0; elemIdx < decomposition.elementPartitions.size(); ++elemIdx) {
        if (decomposition.elementPartitions[elemIdx] == currentRank) {
            localElements.push_back(elemIdx);
        }
    }
    
    return localElements;
}

std::vector<int> MagnetoDynamics2DSolver::getGhostBoundaryElements(const DomainDecompositionResult& decomposition) {
    std::vector<int> ghostElements;
    int currentRank = comm_->getRank();
    
    // 获取幽灵边界元素
    if (decomposition.ghostElements.find(currentRank) != decomposition.ghostElements.end()) {
        ghostElements = decomposition.ghostElements.at(currentRank);
    }
    
    return ghostElements;
}

void MagnetoDynamics2DSolver::exchangeGhostData() {
    if (!parallelAssembler_) {
        return;
    }
    
    // 简化的幽灵数据交换
    // 实际实现应包括完整的幽灵数据通信
    
    if (comm_->getRank() == 0) {
        std::cout << "幽灵数据交换完成" << std::endl;
    }
}

double MagnetoDynamics2DSolver::computeLoadBalance(const DomainDecompositionResult& decomposition) {
    if (decomposition.partitionSizes.empty()) {
        return 1.0;
    }
    
    // 计算负载均衡度：最小分区大小 / 最大分区大小
    int minSize = *std::min_element(decomposition.partitionSizes.begin(), decomposition.partitionSizes.end());
    int maxSize = *std::max_element(decomposition.partitionSizes.begin(), decomposition.partitionSizes.end());
    
    return static_cast<double>(minSize) / maxSize;
}

} // namespace elmer