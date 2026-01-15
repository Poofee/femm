// ParallelMatrixAssembly.cpp - MPI并行矩阵组装器实现
// 对应Fortran模块: ParallelAssembly.F90

#include "ParallelMatrixAssembly.h"
#include "MPIConfig.h"
#include "MatrixAssembly.h"
#include "BoundaryConditions.h"
#include <chrono>
#include <algorithm>
#include <numeric>

namespace elmer {

// ===== ParallelMatrixAssembler实现 =====

ParallelMatrixAssembler::ParallelMatrixAssembler(
    std::shared_ptr<MPICommunicator> comm,
    std::shared_ptr<DomainDecompositionManager> decompositionManager)
    : comm_(comm ? comm : MPIUtils::getDefaultComm())
    , decompositionManager_(decompositionManager ? decompositionManager : 
        std::make_shared<DomainDecompositionManager>(comm_)) {
}

void ParallelMatrixAssembler::initialize(int globalSize, 
                                        const DomainDecompositionResult& decompositionResult) {
    // 创建分布式矩阵和向量
    localStiffnessMatrix_ = std::make_shared<DistributedMatrix>(globalSize, comm_);
    localMassMatrix_ = std::make_shared<DistributedMatrix>(globalSize, comm_);
    localDampingMatrix_ = std::make_shared<DistributedMatrix>(globalSize, comm_);
    localRhsVector_ = std::make_shared<DistributedVector>(globalSize, comm_);
    
    // 初始化矩阵和向量
    localStiffnessMatrix_->initialize();
    localMassMatrix_->initialize();
    localDampingMatrix_->initialize();
    localRhsVector_->initialize();
    
    // 重置幽灵数据
    ghostStiffnessData_.clear();
    ghostMassData_.clear();
    ghostDampingData_.clear();
    ghostRhsData_.clear();
    
    isAssembled_ = false;
}

void ParallelMatrixAssembler::assembleElement(
    int elementId,
    const std::vector<std::vector<double>>& elementStiffness,
    const std::vector<std::vector<double>>& elementMass,
    const std::vector<std::vector<double>>& elementDamping,
    const std::vector<double>& elementRhs,
    const std::vector<int>& nodeIndices) {
    
    // 组装刚度矩阵
    if (!elementStiffness.empty()) {
        assembleElementToLocalMatrix(localStiffnessMatrix_, elementStiffness, nodeIndices);
    }
    
    // 组装质量矩阵
    if (!elementMass.empty()) {
        assembleElementToLocalMatrix(localMassMatrix_, elementMass, nodeIndices);
    }
    
    // 组装阻尼矩阵
    if (!elementDamping.empty()) {
        assembleElementToLocalMatrix(localDampingMatrix_, elementDamping, nodeIndices);
    }
    
    // 组装右端向量
    if (!elementRhs.empty()) {
        assembleElementToLocalVector(localRhsVector_, elementRhs, nodeIndices);
    }
}

void ParallelMatrixAssembler::exchangeGhostData() {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // 准备幽灵数据交换
    prepareGhostDataExchange();
    
    // 执行非阻塞通信
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> recvRequests;
    
    // 发送本地幽灵数据到其他进程
    for (const auto& ghostData : ghostStiffnessData_) {
        int destRank = ghostData.first;
        const auto& data = ghostData.second;
        
        MPI_Request request;
        comm_->isend(destRank, 0, data.data(), data.size(), MPI_DOUBLE, &request);
        sendRequests.push_back(request);
    }
    
    // 接收其他进程的幽灵数据
    for (const auto& ghostData : ghostStiffnessData_) {
        int srcRank = ghostData.first;
        auto& recvBuffer = ghostStiffnessData_[srcRank];
        
        MPI_Request request;
        comm_->irecv(srcRank, 0, recvBuffer.data(), recvBuffer.size(), MPI_DOUBLE, &request);
        recvRequests.push_back(request);
    }
    
    // 等待所有通信完成
    if (!sendRequests.empty()) {
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
    }
    
    if (!recvRequests.empty()) {
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);
    }
    
    // 处理接收到的幽灵数据
    processReceivedGhostData();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    
    // 记录通信时间（简化实现）
    // 在实际应用中应该使用更精确的计时
}

void ParallelMatrixAssembler::finalizeAssembly() {
    // 交换幽灵数据
    exchangeGhostData();
    
    // 完成矩阵组装
    localStiffnessMatrix_->finalizeAssembly();
    localMassMatrix_->finalizeAssembly();
    localDampingMatrix_->finalizeAssembly();
    localRhsVector_->finalizeAssembly();
    
    isAssembled_ = true;
}

void ParallelMatrixAssembler::applyBoundaryConditions(
    const std::vector<BoundaryCondition>& boundaryConditions) {
    
    // 应用边界条件到本地系统
    for (const auto& bc : boundaryConditions) {
        if (bc.type == BoundaryConditionType::DIRICHLET) {
            // 检查边界条件是否在本地分区内
            if (localRhsVector_->isLocalIndex(bc.dofIndex)) {
                // 应用Dirichlet边界条件
                localStiffnessMatrix_->applyDirichletBC(bc.dofIndex, bc.value);
                localMassMatrix_->applyDirichletBC(bc.dofIndex, bc.value);
                localDampingMatrix_->applyDirichletBC(bc.dofIndex, bc.value);
                localRhsVector_->applyDirichletBC(bc.dofIndex, bc.value);
            }
        }
    }
}

void ParallelMatrixAssembler::reset() {
    localStiffnessMatrix_.reset();
    localMassMatrix_.reset();
    localDampingMatrix_.reset();
    localRhsVector_.reset();
    
    ghostStiffnessData_.clear();
    ghostMassData_.clear();
    ghostDampingData_.clear();
    ghostRhsData_.clear();
    
    isAssembled_ = false;
}

ParallelMatrixAssembler::AssemblyStatistics ParallelMatrixAssembler::getStatistics() const {
    AssemblyStatistics stats;
    
    if (decompositionManager_) {
        // 获取本地元素统计
        auto localElements = decompositionManager_->getLocalElements(
            DomainDecompositionResult()); // 需要传入实际的分解结果
        stats.localElements = localElements.size();
        
        auto ghostElements = decompositionManager_->getGhostElements(
            DomainDecompositionResult());
        stats.ghostElements = ghostElements.size();
        
        auto boundaryElements = decompositionManager_->getBoundaryElements(
            DomainDecompositionResult());
        stats.boundaryElements = boundaryElements.size();
    }
    
    // 简化实现：在实际应用中应该使用更精确的计时
    stats.assemblyTime = 0.0;
    stats.communicationTime = 0.0;
    stats.loadBalance = 1.0;
    
    return stats;
}

// ===== 私有方法实现 =====

void ParallelMatrixAssembler::assembleElementToLocalMatrix(
    std::shared_ptr<DistributedMatrix>& localMatrix,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<int>& nodeIndices) {
    
    int numNodes = nodeIndices.size();
    
    for (int i = 0; i < numNodes; ++i) {
        int globalRow = nodeIndices[i];
        
        for (int j = 0; j < numNodes; ++j) {
            int globalCol = nodeIndices[j];
            double value = elementMatrix[i][j];
            
            if (value != 0.0) {
                localMatrix->addToElement(globalRow, globalCol, value);
            }
        }
    }
}

void ParallelMatrixAssembler::assembleElementToLocalVector(
    std::shared_ptr<DistributedVector>& localVector,
    const std::vector<double>& elementVector,
    const std::vector<int>& nodeIndices) {
    
    int numNodes = nodeIndices.size();
    
    for (int i = 0; i < numNodes; ++i) {
        int globalIndex = nodeIndices[i];
        double value = elementVector[i];
        
        if (value != 0.0) {
            localVector->addToElement(globalIndex, value);
        }
    }
}

void ParallelMatrixAssembler::prepareGhostDataExchange() {
    // 简化实现：在实际应用中需要更复杂的幽灵数据管理
    
    // 获取幽灵元素列表
    auto ghostElements = decompositionManager_->getGhostElements(
        DomainDecompositionResult()); // 需要传入实际的分解结果
    
    // 为每个幽灵元素准备数据
    for (int ghostElemId : ghostElements) {
        // 在实际应用中，这里需要获取幽灵元素的矩阵贡献
        // 简化实现：创建空的幽灵数据缓冲区
        ghostStiffnessData_[ghostElemId] = std::vector<double>();
        ghostMassData_[ghostElemId] = std::vector<double>();
        ghostDampingData_[ghostElemId] = std::vector<double>();
        ghostRhsData_[ghostElemId] = std::vector<double>();
    }
}

void ParallelMatrixAssembler::processReceivedGhostData() {
    // 简化实现：在实际应用中需要处理接收到的幽灵数据
    
    // 更新本地矩阵的幽灵数据
    updateLocalMatrixWithGhostData(localStiffnessMatrix_, ghostStiffnessData_);
    updateLocalMatrixWithGhostData(localMassMatrix_, ghostMassData_);
    updateLocalMatrixWithGhostData(localDampingMatrix_, ghostDampingData_);
    
    // 更新本地向量的幽灵数据
    updateLocalVectorWithGhostData(localRhsVector_, ghostRhsData_);
}

void ParallelMatrixAssembler::updateLocalMatrixWithGhostData(
    std::shared_ptr<DistributedMatrix>& localMatrix,
    const std::map<int, std::vector<double>>& ghostData) {
    
    // 简化实现：在实际应用中需要将幽灵数据添加到本地矩阵
    for (const auto& ghostEntry : ghostData) {
        int ghostElemId = ghostEntry.first;
        const auto& data = ghostEntry.second;
        
        // 在实际应用中，这里需要解析幽灵数据并添加到矩阵
        // 简化实现：跳过具体实现
    }
}

void ParallelMatrixAssembler::updateLocalVectorWithGhostData(
    std::shared_ptr<DistributedVector>& localVector,
    const std::map<int, std::vector<double>>& ghostData) {
    
    // 简化实现：在实际应用中需要将幽灵数据添加到本地向量
    for (const auto& ghostEntry : ghostData) {
        int ghostElemId = ghostEntry.first;
        const auto& data = ghostEntry.second;
        
        // 在实际应用中，这里需要解析幽灵数据并添加到向量
        // 简化实现：跳过具体实现
    }
}

// ===== ParallelMatrixAssemblyManager实现 =====

ParallelMatrixAssemblyManager::ParallelMatrixAssemblyManager(
    std::shared_ptr<MPICommunicator> comm,
    std::shared_ptr<DomainDecompositionManager> decompositionManager)
    : comm_(comm ? comm : MPIUtils::getDefaultComm())
    , decompositionManager_(decompositionManager ? decompositionManager : 
        std::make_shared<DomainDecompositionManager>(comm_)) {
}

std::shared_ptr<ParallelMatrixAssembler> ParallelMatrixAssemblyManager::createAssembler(
    int globalSize, const DomainDecompositionResult& decompositionResult) {
    
    auto assembler = std::make_shared<ParallelMatrixAssembler>(comm_, decompositionManager_);
    assembler->initialize(globalSize, decompositionResult);
    return assembler;
}

std::shared_ptr<DistributedLinearSystem> ParallelMatrixAssemblyManager::assembleSystem(
    std::shared_ptr<Mesh> mesh,
    const MaterialDatabase& materialDB,
    const MagnetoDynamics2DParameters& parameters) {
    
    if (!mesh) {
        throw std::runtime_error("Mesh not provided for parallel assembly");
    }
    
    // 获取网格信息
    auto& nodes = mesh->getNodes();
    int globalSize = static_cast<int>(nodes.numberOfNodes());
    
    // 执行域分解
    std::vector<MeshElement> elements;
    // 在实际应用中需要从网格创建MeshElement列表
    
    auto decompositionResult = decompositionManager_->decompose(elements, comm_->getSize());
    
    // 创建并行组装器
    auto assembler = createAssembler(globalSize, decompositionResult);
    
    // 获取本地元素
    auto localElements = decompositionManager_->getLocalElements(decompositionResult);
    
    // 组装本地元素
    for (int elemId : localElements) {
        // 在实际应用中需要计算单元矩阵
        std::vector<std::vector<double>> elementStiffness;
        std::vector<std::vector<double>> elementMass;
        std::vector<std::vector<double>> elementDamping;
        std::vector<double> elementRhs;
        std::vector<int> nodeIndices;
        
        // 简化实现：跳过具体的单元矩阵计算
        assembler->assembleElement(elemId, elementStiffness, elementMass, 
                                 elementDamping, elementRhs, nodeIndices);
    }
    
    // 完成组装
    assembler->finalizeAssembly();
    
    // 创建分布式线性系统
    auto linearSystem = std::make_shared<DistributedLinearSystem>(
        assembler->getStiffnessMatrix(),
        assembler->getMassMatrix(),
        assembler->getDampingMatrix(),
        assembler->getRhsVector(),
        comm_);
    
    return linearSystem;
}

std::shared_ptr<ParallelMatrixAssemblyManager> 
ParallelMatrixAssemblyManager::getDefaultManager(std::shared_ptr<MPICommunicator> comm) {
    return std::make_shared<ParallelMatrixAssemblyManager>(comm);
}

} // namespace elmer