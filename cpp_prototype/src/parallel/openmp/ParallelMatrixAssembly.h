// ParallelMatrixAssembly.h - MPI并行矩阵组装�?// 对应Fortran模块: ParallelAssembly.F90

#pragma once

#include "../../parallel/mpi/DistributedLinearAlgebra.h"
#include "DomainDecomposition.h"
#include "MatrixAssembly.h"
#include <memory>
#include <vector>
#include <map>

namespace elmer {

/**
 * @brief MPI并行矩阵组装�? * 
 * 负责在分布式内存环境中组装全局系统矩阵�? * 支持域分解和幽灵数据交换�? */
class ParallelMatrixAssembler {
private:
    std::shared_ptr<MPICommunicator> comm_;           // MPI通信�?    std::shared_ptr<DomainDecompositionManager> decompositionManager_; // 域分解管理器
    
    // 本地矩阵数据
    std::shared_ptr<DistributedMatrix> localStiffnessMatrix_;   // 本地刚度矩阵
    std::shared_ptr<DistributedMatrix> localMassMatrix_;        // 本地质量矩阵
    std::shared_ptr<DistributedMatrix> localDampingMatrix_;     // 本地阻尼矩阵
    std::shared_ptr<DistributedVector> localRhsVector_;         // 本地右端向量
    
    // 幽灵数据管理
    std::map<int, std::vector<double>> ghostStiffnessData_;     // 幽灵刚度数据
    std::map<int, std::vector<double>> ghostMassData_;          // 幽灵质量数据
    std::map<int, std::vector<double>> ghostDampingData_;       // 幽灵阻尼数据
    std::map<int, std::vector<double>> ghostRhsData_;           // 幽灵右端数据
    
    // 组装状�?    bool isAssembled_ = false;
    
public:
    /**
     * @brief 构造函�?     * @param comm MPI通信�?     * @param decompositionManager 域分解管理器
     */
    ParallelMatrixAssembler(
        std::shared_ptr<MPICommunicator> comm = nullptr,
        std::shared_ptr<DomainDecompositionManager> decompositionManager = nullptr);
    
    virtual ~ParallelMatrixAssembler() = default;
    
    /**
     * @brief 初始化并行组装器
     * @param globalSize 全局系统大小
     * @param decompositionResult 域分解结�?     */
    void initialize(int globalSize, const DomainDecompositionResult& decompositionResult);
    
    /**
     * @brief 组装单元矩阵到本地系�?     * @param elementId 单元ID
     * @param elementStiffness 单元刚度矩阵
     * @param elementMass 单元质量矩阵
     * @param elementDamping 单元阻尼矩阵
     * @param elementRhs 单元右端向量
     * @param nodeIndices 单元节点索引
     */
    void assembleElement(
        int elementId,
        const std::vector<std::vector<double>>& elementStiffness,
        const std::vector<std::vector<double>>& elementMass,
        const std::vector<std::vector<double>>& elementDamping,
        const std::vector<double>& elementRhs,
        const std::vector<int>& nodeIndices);
    
    /**
     * @brief 交换幽灵数据
     */
    void exchangeGhostData();
    
    /**
     * @brief 完成矩阵组装
     */
    void finalizeAssembly();
    
    /**
     * @brief 应用边界条件
     * @param boundaryConditions 边界条件列表
     */
    void applyBoundaryConditions(const std::vector<BoundaryCondition>& boundaryConditions);
    
    // 获取组装结果
    std::shared_ptr<DistributedMatrix> getStiffnessMatrix() const { return localStiffnessMatrix_; }
    std::shared_ptr<DistributedMatrix> getMassMatrix() const { return localMassMatrix_; }
    std::shared_ptr<DistributedMatrix> getDampingMatrix() const { return localDampingMatrix_; }
    std::shared_ptr<DistributedVector> getRhsVector() const { return localRhsVector_; }
    
    /**
     * @brief 检查组装状�?     */
    bool isAssembled() const { return isAssembled_; }
    
    /**
     * @brief 重置组装器状�?     */
    void reset();
    
    /**
     * @brief 获取组装统计信息
     */
    struct AssemblyStatistics {
        int localElements;        // 本地元素数量
        int ghostElements;        // 幽灵元素数量
        int boundaryElements;     // 边界元素数量
        double assemblyTime;      // 组装时间（秒�?        double communicationTime; // 通信时间（秒�?        double loadBalance;       // 负载均衡�?    };
    
    AssemblyStatistics getStatistics() const;
    
private:
    /**
     * @brief 组装单元矩阵到本地矩�?     */
    void assembleElementToLocalMatrix(
        std::shared_ptr<DistributedMatrix>& localMatrix,
        const std::vector<std::vector<double>>& elementMatrix,
        const std::vector<int>& nodeIndices);
    
    /**
     * @brief 组装单元向量到本地向�?     */
    void assembleElementToLocalVector(
        std::shared_ptr<DistributedVector>& localVector,
        const std::vector<double>& elementVector,
        const std::vector<int>& nodeIndices);
    
    /**
     * @brief 准备幽灵数据交换
     */
    void prepareGhostDataExchange();
    
    /**
     * @brief 处理接收到的幽灵数据
     */
    void processReceivedGhostData();
    
    /**
     * @brief 更新本地矩阵的幽灵数�?     */
    void updateLocalMatrixWithGhostData(
        std::shared_ptr<DistributedMatrix>& localMatrix,
        const std::map<int, std::vector<double>>& ghostData);
    
    /**
     * @brief 更新本地向量的幽灵数�?     */
    void updateLocalVectorWithGhostData(
        std::shared_ptr<DistributedVector>& localVector,
        const std::map<int, std::vector<double>>& ghostData);
};

/**
 * @brief 并行矩阵组装管理�? */
class ParallelMatrixAssemblyManager {
private:
    std::shared_ptr<MPICommunicator> comm_;
    std::shared_ptr<DomainDecompositionManager> decompositionManager_;
    std::shared_ptr<ParallelMatrixAssembler> assembler_;
    
public:
    ParallelMatrixAssemblyManager(
        std::shared_ptr<MPICommunicator> comm = nullptr,
        std::shared_ptr<DomainDecompositionManager> decompositionManager = nullptr);
    
    /**
     * @brief 创建并行矩阵组装�?     */
    std::shared_ptr<ParallelMatrixAssembler> createAssembler(int globalSize, 
                                                             const DomainDecompositionResult& decompositionResult);
    
    /**
     * @brief 执行并行矩阵组装
     * @param mesh 网格
     * @param materialDB 材料数据�?     * @param parameters 求解器参�?     * @return 组装好的分布式线性系�?     */
    std::shared_ptr<DistributedLinearSystem> assembleSystem(
        std::shared_ptr<Mesh> mesh,
        const MaterialDatabase& materialDB,
        const MagnetoDynamics2DParameters& parameters);
    
    /**
     * @brief 获取默认管理�?     */
    static std::shared_ptr<ParallelMatrixAssemblyManager> getDefaultManager(
        std::shared_ptr<MPICommunicator> comm = nullptr);
};

} // namespace elmer

