// DomainDecomposition.h - 域分解算法定义
// 对应Fortran模块: DomainDecomposition.F90, Partitioning.F90

#pragma once

#include "MPIConfig.h"
#include <memory>
#include <vector>
#include <map>
#include <set>

namespace elmer {

// ===== 网格元素信息 =====
struct MeshElement {
    int id;                    // 元素ID
    std::vector<double> coords; // 元素坐标（中心点或顶点）
    int partition;             // 分区编号
    std::set<int> neighbors;   // 相邻元素
    
    MeshElement(int elemId = -1) : id(elemId), partition(-1) {}
};

// ===== 域分解结果 =====
struct DomainDecompositionResult {
    int numPartitions;                         // 分区数量
    std::map<int, std::vector<int>> partitionMap; // 分区到元素的映射
    std::map<int, std::vector<int>> ghostElements; // 幽灵元素映射
    std::map<int, std::vector<int>> boundaryElements; // 边界元素映射
    std::vector<int> elementPartitions;        // 元素到分区的映射
    
    // 分区负载信息
    std::vector<int> partitionSizes;           // 每个分区的元素数量
    double loadBalance;                        // 负载均衡度
    
    DomainDecompositionResult(int nParts = 1) : numPartitions(nParts), loadBalance(1.0) {
        partitionSizes.resize(nParts, 0);
    }
};

// ===== 域分解算法基类 =====
class DomainDecomposer {
protected:
    std::shared_ptr<MPICommunicator> comm_;    // MPI通信器
    
public:
    DomainDecomposer(std::shared_ptr<MPICommunicator> comm = nullptr)
        : comm_(comm ? comm : MPIUtils::getDefaultComm()) {}
    
    virtual ~DomainDecomposer() = default;
    
    /**
     * @brief 执行域分解
     * @param elements 网格元素列表
     * @param numPartitions 分区数量
     * @return 域分解结果
     */
    virtual DomainDecompositionResult decompose(
        const std::vector<MeshElement>& elements, 
        int numPartitions) = 0;
    
    /**
     * @brief 评估分解质量
     * @param result 域分解结果
     * @param elements 网格元素列表
     * @return 质量指标（越小越好）
     */
    virtual double evaluateQuality(
        const DomainDecompositionResult& result,
        const std::vector<MeshElement>& elements) = 0;
    
    // 获取通信器
    std::shared_ptr<MPICommunicator> getCommunicator() const { return comm_; }
};

// ===== 坐标分解算法 =====
class CoordinateDecomposer : public DomainDecomposer {
private:
    int coordinateAxis_;  // 坐标轴（0=x, 1=y, 2=z）
    
public:
    CoordinateDecomposer(int axis = 0, std::shared_ptr<MPICommunicator> comm = nullptr)
        : DomainDecomposer(comm), coordinateAxis_(axis) {}
    
    DomainDecompositionResult decompose(
        const std::vector<MeshElement>& elements, 
        int numPartitions) override;
    
    double evaluateQuality(
        const DomainDecompositionResult& result,
        const std::vector<MeshElement>& elements) override;
    
    void setCoordinateAxis(int axis) { coordinateAxis_ = axis; }
};

// ===== 递归坐标分解算法 =====
class RecursiveCoordinateDecomposer : public DomainDecomposer {
private:
    int maxDepth_;  // 最大递归深度
    
public:
    RecursiveCoordinateDecomposer(int maxDepth = 3, std::shared_ptr<MPICommunicator> comm = nullptr)
        : DomainDecomposer(comm), maxDepth_(maxDepth) {}
    
    DomainDecompositionResult decompose(
        const std::vector<MeshElement>& elements, 
        int numPartitions) override;
    
    double evaluateQuality(
        const DomainDecompositionResult& result,
        const std::vector<MeshElement>& elements) override;
    
    void setMaxDepth(int depth) { maxDepth_ = depth; }
};

// ===== 图分解算法（基于METIS） =====
class GraphDecomposer : public DomainDecomposer {
private:
    double imbalanceTolerance_;  // 负载不平衡容忍度
    
public:
    GraphDecomposer(double imbalanceTol = 1.05, std::shared_ptr<MPICommunicator> comm = nullptr)
        : DomainDecomposer(comm), imbalanceTolerance_(imbalanceTol) {}
    
    DomainDecompositionResult decompose(
        const std::vector<MeshElement>& elements, 
        int numPartitions) override;
    
    double evaluateQuality(
        const DomainDecompositionResult& result,
        const std::vector<MeshElement>& elements) override;
    
    void setImbalanceTolerance(double tol) { imbalanceTolerance_ = tol; }
};

// ===== 域分解管理器 =====
class DomainDecompositionManager {
private:
    std::shared_ptr<DomainDecomposer> decomposer_;
    std::shared_ptr<MPICommunicator> comm_;
    
public:
    DomainDecompositionManager(std::shared_ptr<DomainDecomposer> decomposer = nullptr,
                              std::shared_ptr<MPICommunicator> comm = nullptr)
        : decomposer_(decomposer), comm_(comm ? comm : MPIUtils::getDefaultComm()) {}
    
    /**
     * @brief 执行域分解
     * @param elements 网格元素列表
     * @param numPartitions 分区数量（默认使用进程数）
     * @return 域分解结果
     */
    DomainDecompositionResult decompose(
        const std::vector<MeshElement>& elements, 
        int numPartitions = -1);
    
    /**
     * @brief 获取本地分区信息
     * @param result 域分解结果
     * @return 本地元素列表
     */
    std::vector<int> getLocalElements(const DomainDecompositionResult& result) const;
    
    /**
     * @brief 获取幽灵元素列表
     * @param result 域分解结果
     * @return 本地幽灵元素列表
     */
    std::vector<int> getGhostElements(const DomainDecompositionResult& result) const;
    
    /**
     * @brief 获取边界元素列表
     * @param result 域分解结果
     * @return 本地边界元素列表
     */
    std::vector<int> getBoundaryElements(const DomainDecompositionResult& result) const;
    
    /**
     * @brief 设置分解器
     */
    void setDecomposer(std::shared_ptr<DomainDecomposer> decomposer) {
        decomposer_ = decomposer;
    }
    
    /**
     * @brief 获取默认分解器（基于可用性选择）
     */
    static std::shared_ptr<DomainDecomposer> getDefaultDecomposer(
        std::shared_ptr<MPICommunicator> comm = nullptr);
};

// ===== 域分解工具函数 =====
namespace DomainDecompositionUtils {
    
    /**
     * @brief 创建网格邻接关系
     * @param elements 网格元素列表
     * @param tolerance 邻接容差
     */
    void createElementAdjacency(std::vector<MeshElement>& elements, double tolerance = 1e-6);
    
    /**
     * @brief 计算分区负载均衡度
     * @param result 域分解结果
     * @return 负载均衡度（1.0表示完美均衡）
     */
    double calculateLoadBalance(const DomainDecompositionResult& result);
    
    /**
     * @brief 计算分区边界长度
     * @param result 域分解结果
     * @param elements 网格元素列表
     * @return 边界长度指标
     */
    double calculateBoundaryLength(
        const DomainDecompositionResult& result,
        const std::vector<MeshElement>& elements);
    
    /**
     * @brief 优化域分解结果
     * @param result 域分解结果
     * @param elements 网格元素列表
     * @param maxIterations 最大迭代次数
     * @return 优化后的结果
     */
    DomainDecompositionResult optimizeDecomposition(
        const DomainDecompositionResult& result,
        const std::vector<MeshElement>& elements,
        int maxIterations = 100);
    
} // namespace DomainDecompositionUtils

} // namespace elmer