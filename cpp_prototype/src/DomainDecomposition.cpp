// DomainDecomposition.cpp - 域分解算法实现
// 对应Fortran模块: DomainDecomposition.F90, Partitioning.F90

#include "DomainDecomposition.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>

namespace elmer {

// ===== CoordinateDecomposer实现 =====
DomainDecompositionResult CoordinateDecomposer::decompose(
    const std::vector<MeshElement>& elements, 
    int numPartitions) {
    
    if (numPartitions <= 0) {
        numPartitions = comm_->getSize();
    }
    
    DomainDecompositionResult result(numPartitions);
    
    // 计算坐标范围
    double minCoord = std::numeric_limits<double>::max();
    double maxCoord = std::numeric_limits<double>::lowest();
    
    for (const auto& elem : elements) {
        if (elem.coords.size() > coordinateAxis_) {
            double coord = elem.coords[coordinateAxis_];
            minCoord = std::min(minCoord, coord);
            maxCoord = std::max(maxCoord, coord);
        }
    }
    
    // 广播坐标范围（如果使用MPI）
    if (comm_->getSize() > 1) {
        comm_->allreduce(&minCoord, &minCoord, 1, MPI_DOUBLE, MPI_MIN);
        comm_->allreduce(&maxCoord, &maxCoord, 1, MPI_DOUBLE, MPI_MAX);
    }
    
    double coordRange = maxCoord - minCoord;
    if (coordRange < 1e-12) {
        // 所有元素在同一个位置，使用简单轮询分配
        for (size_t i = 0; i < elements.size(); ++i) {
            int partition = i % numPartitions;
            result.partitionMap[partition].push_back(elements[i].id);
            result.elementPartitions.push_back(partition);
        }
    } else {
        // 基于坐标的分解
        double partitionWidth = coordRange / numPartitions;
        
        for (const auto& elem : elements) {
            if (elem.coords.size() > coordinateAxis_) {
                double coord = elem.coords[coordinateAxis_];
                int partition = static_cast<int>((coord - minCoord) / partitionWidth);
                partition = std::min(partition, numPartitions - 1);
                
                result.partitionMap[partition].push_back(elem.id);
                result.elementPartitions.push_back(partition);
            } else {
                // 没有坐标信息，使用默认分配
                int partition = elem.id % numPartitions;
                result.partitionMap[partition].push_back(elem.id);
                result.elementPartitions.push_back(partition);
            }
        }
    }
    
    // 计算分区大小
    for (int i = 0; i < numPartitions; ++i) {
        result.partitionSizes[i] = result.partitionMap[i].size();
    }
    
    // 识别幽灵元素和边界元素
    if (comm_->getSize() > 1) {
        // 简化实现：在实际应用中需要更复杂的邻接分析
        for (const auto& elem : elements) {
            for (int neighborId : elem.neighbors) {
                if (neighborId < result.elementPartitions.size()) {
                    int elemPartition = result.elementPartitions[elem.id];
                    int neighborPartition = result.elementPartitions[neighborId];
                    
                    if (elemPartition != neighborPartition) {
                        // 边界元素
                        result.boundaryElements[elemPartition].push_back(elem.id);
                        
                        // 幽灵元素
                        result.ghostElements[elemPartition].push_back(neighborId);
                    }
                }
            }
        }
    }
    
    // 计算负载均衡度
    result.loadBalance = DomainDecompositionUtils::calculateLoadBalance(result);
    
    return result;
}

double CoordinateDecomposer::evaluateQuality(
    const DomainDecompositionResult& result,
    const std::vector<MeshElement>& elements) {
    
    // 质量指标：负载不均衡度 + 边界长度
    double imbalance = 1.0 - result.loadBalance;
    double boundaryLength = DomainDecompositionUtils::calculateBoundaryLength(result, elements);
    
    // 归一化边界长度（假设最大边界长度为元素数量的平方根）
    double maxBoundary = std::sqrt(static_cast<double>(elements.size()));
    double normalizedBoundary = boundaryLength / maxBoundary;
    
    return imbalance + normalizedBoundary;
}

// ===== RecursiveCoordinateDecomposer实现 =====
DomainDecompositionResult RecursiveCoordinateDecomposer::decompose(
    const std::vector<MeshElement>& elements, 
    int numPartitions) {
    
    if (numPartitions <= 0) {
        numPartitions = comm_->getSize();
    }
    
    // 递归坐标分解算法
    DomainDecompositionResult result(numPartitions);
    
    // 初始分配：所有元素属于分区0
    for (const auto& elem : elements) {
        result.partitionMap[0].push_back(elem.id);
        result.elementPartitions.push_back(0);
    }
    
    // 递归分解
    int currentAxis = 0;
    int currentDepth = 0;
    
    while (result.partitionMap.size() < static_cast<size_t>(numPartitions) && 
           currentDepth < maxDepth_) {
        
        std::map<int, std::vector<int>> newPartitionMap;
        
        for (const auto& partition : result.partitionMap) {
            int partitionId = partition.first;
            const auto& partitionElements = partition.second;
            
            if (partitionElements.size() <= 1) {
                // 分区太小，不进行进一步分解
                newPartitionMap[partitionId] = partitionElements;
                continue;
            }
            
            // 计算当前分区的坐标范围
            double minCoord = std::numeric_limits<double>::max();
            double maxCoord = std::numeric_limits<double>::lowest();
            
            for (int elemId : partitionElements) {
                if (elemId < elements.size() && 
                    elements[elemId].coords.size() > currentAxis) {
                    double coord = elements[elemId].coords[currentAxis];
                    minCoord = std::min(minCoord, coord);
                    maxCoord = std::max(maxCoord, coord);
                }
            }
            
            double coordRange = maxCoord - minCoord;
            if (coordRange < 1e-12) {
                // 坐标范围太小，切换到下一个轴
                currentAxis = (currentAxis + 1) % 3;
                newPartitionMap[partitionId] = partitionElements;
                continue;
            }
            
            // 将当前分区分成两个子分区
            double midCoord = (minCoord + maxCoord) / 2.0;
            
            std::vector<int> leftPartition, rightPartition;
            for (int elemId : partitionElements) {
                if (elemId < elements.size() && 
                    elements[elemId].coords.size() > currentAxis) {
                    double coord = elements[elemId].coords[currentAxis];
                    if (coord <= midCoord) {
                        leftPartition.push_back(elemId);
                    } else {
                        rightPartition.push_back(elemId);
                    }
                } else {
                    // 没有坐标信息，随机分配
                    if (rand() % 2 == 0) {
                        leftPartition.push_back(elemId);
                    } else {
                        rightPartition.push_back(elemId);
                    }
                }
            }
            
            // 分配新的分区ID
            int leftPartitionId = partitionId;
            int rightPartitionId = newPartitionMap.size();
            
            newPartitionMap[leftPartitionId] = leftPartition;
            newPartitionMap[rightPartitionId] = rightPartition;
            
            // 更新元素分区映射
            for (int elemId : leftPartition) {
                result.elementPartitions[elemId] = leftPartitionId;
            }
            for (int elemId : rightPartition) {
                result.elementPartitions[elemId] = rightPartitionId;
            }
        }
        
        result.partitionMap = newPartitionMap;
        currentAxis = (currentAxis + 1) % 3;
        currentDepth++;
    }
    
    // 计算分区大小
    result.partitionSizes.clear();
    result.partitionSizes.resize(numPartitions, 0);
    for (const auto& partition : result.partitionMap) {
        if (partition.first < numPartitions) {
            result.partitionSizes[partition.first] = partition.second.size();
        }
    }
    
    // 识别幽灵元素和边界元素（简化实现）
    DomainDecompositionUtils::createElementAdjacency(const_cast<std::vector<MeshElement>&>(elements));
    
    // 计算负载均衡度
    result.loadBalance = DomainDecompositionUtils::calculateLoadBalance(result);
    
    return result;
}

double RecursiveCoordinateDecomposer::evaluateQuality(
    const DomainDecompositionResult& result,
    const std::vector<MeshElement>& elements) {
    
    // 使用与坐标分解相同的质量评估方法
    CoordinateDecomposer coordDecomposer(0, comm_);
    return coordDecomposer.evaluateQuality(result, elements);
}

// ===== GraphDecomposer实现 =====
DomainDecompositionResult GraphDecomposer::decompose(
    const std::vector<MeshElement>& elements, 
    int numPartitions) {
    
    if (numPartitions <= 0) {
        numPartitions = comm_->getSize();
    }
    
    DomainDecompositionResult result(numPartitions);
    
    // 简化实现：使用基于度的贪心算法
    // 在实际应用中应该使用METIS等专业图分解库
    
    // 创建邻接矩阵
    std::vector<std::vector<int>> adjacency(elements.size());
    for (const auto& elem : elements) {
        for (int neighborId : elem.neighbors) {
            if (neighborId < elements.size()) {
                adjacency[elem.id].push_back(neighborId);
            }
        }
    }
    
    // 计算每个元素的度
    std::vector<int> degrees(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        degrees[i] = adjacency[i].size();
    }
    
    // 贪心分配：优先分配高度元素
    std::vector<bool> assigned(elements.size(), false);
    std::vector<int> partitionLoads(numPartitions, 0);
    
    // 按度降序排序元素
    std::vector<int> elementIds(elements.size());
    std::iota(elementIds.begin(), elementIds.end(), 0);
    std::sort(elementIds.begin(), elementIds.end(), 
              [&](int a, int b) { return degrees[a] > degrees[b]; });
    
    for (int elemId : elementIds) {
        if (assigned[elemId]) continue;
        
        // 找到负载最小的分区
        int minLoadPartition = 0;
        int minLoad = partitionLoads[0];
        for (int i = 1; i < numPartitions; ++i) {
            if (partitionLoads[i] < minLoad) {
                minLoad = partitionLoads[i];
                minLoadPartition = i;
            }
        }
        
        // 分配当前元素
        result.partitionMap[minLoadPartition].push_back(elemId);
        result.elementPartitions[elemId] = minLoadPartition;
        partitionLoads[minLoadPartition]++;
        assigned[elemId] = true;
        
        // 尝试分配相邻元素到同一分区（减少边界）
        for (int neighborId : adjacency[elemId]) {
            if (!assigned[neighborId] && partitionLoads[minLoadPartition] < 
                elements.size() / numPartitions * imbalanceTolerance_) {
                
                result.partitionMap[minLoadPartition].push_back(neighborId);
                result.elementPartitions[neighborId] = minLoadPartition;
                partitionLoads[minLoadPartition]++;
                assigned[neighborId] = true;
            }
        }
    }
    
    // 分配剩余元素
    for (size_t i = 0; i < elements.size(); ++i) {
        if (!assigned[i]) {
            int minLoadPartition = 0;
            int minLoad = partitionLoads[0];
            for (int j = 1; j < numPartitions; ++j) {
                if (partitionLoads[j] < minLoad) {
                    minLoad = partitionLoads[j];
                    minLoadPartition = j;
                }
            }
            
            result.partitionMap[minLoadPartition].push_back(i);
            result.elementPartitions[i] = minLoadPartition;
            partitionLoads[minLoadPartition]++;
        }
    }
    
    // 计算分区大小
    for (int i = 0; i < numPartitions; ++i) {
        result.partitionSizes[i] = result.partitionMap[i].size();
    }
    
    // 识别幽灵元素和边界元素
    for (const auto& elem : elements) {
        for (int neighborId : elem.neighbors) {
            if (neighborId < result.elementPartitions.size()) {
                int elemPartition = result.elementPartitions[elem.id];
                int neighborPartition = result.elementPartitions[neighborId];
                
                if (elemPartition != neighborPartition) {
                    result.boundaryElements[elemPartition].push_back(elem.id);
                    result.ghostElements[elemPartition].push_back(neighborId);
                }
            }
        }
    }
    
    // 计算负载均衡度
    result.loadBalance = DomainDecompositionUtils::calculateLoadBalance(result);
    
    return result;
}

double GraphDecomposer::evaluateQuality(
    const DomainDecompositionResult& result,
    const std::vector<MeshElement>& elements) {
    
    // 图分解质量指标：边界切割数量
    int edgeCut = 0;
    for (const auto& elem : elements) {
        for (int neighborId : elem.neighbors) {
            if (neighborId < result.elementPartitions.size()) {
                if (result.elementPartitions[elem.id] != result.elementPartitions[neighborId]) {
                    edgeCut++;
                }
            }
        }
    }
    
    // 归一化边界切割数量
    double totalEdges = 0;
    for (const auto& elem : elements) {
        totalEdges += elem.neighbors.size();
    }
    
    double normalizedEdgeCut = edgeCut / totalEdges;
    double imbalance = 1.0 - result.loadBalance;
    
    return normalizedEdgeCut + imbalance;
}

// ===== DomainDecompositionManager实现 =====
DomainDecompositionResult DomainDecompositionManager::decompose(
    const std::vector<MeshElement>& elements, 
    int numPartitions) {
    
    if (!decomposer_) {
        decomposer_ = getDefaultDecomposer(comm_);
    }
    
    if (numPartitions <= 0) {
        numPartitions = comm_->getSize();
    }
    
    return decomposer_->decompose(elements, numPartitions);
}

std::vector<int> DomainDecompositionManager::getLocalElements(
    const DomainDecompositionResult& result) const {
    
    int rank = comm_->getRank();
    if (result.partitionMap.find(rank) != result.partitionMap.end()) {
        return result.partitionMap.at(rank);
    }
    return {};
}

std::vector<int> DomainDecompositionManager::getGhostElements(
    const DomainDecompositionResult& result) const {
    
    int rank = comm_->getRank();
    if (result.ghostElements.find(rank) != result.ghostElements.end()) {
        return result.ghostElements.at(rank);
    }
    return {};
}

std::vector<int> DomainDecompositionManager::getBoundaryElements(
    const DomainDecompositionResult& result) const {
    
    int rank = comm_->getRank();
    if (result.boundaryElements.find(rank) != result.boundaryElements.end()) {
        return result.boundaryElements.at(rank);
    }
    return {};
}

std::shared_ptr<DomainDecomposer> DomainDecompositionManager::getDefaultDecomposer(
    std::shared_ptr<MPICommunicator> comm) {
    
    // 根据进程数量选择默认分解器
    int numProcesses = comm ? comm->getSize() : 1;
    
    if (numProcesses <= 8) {
        // 小规模并行，使用坐标分解
        return std::make_shared<CoordinateDecomposer>(0, comm);
    } else if (numProcesses <= 64) {
        // 中等规模并行，使用递归坐标分解
        return std::make_shared<RecursiveCoordinateDecomposer>(3, comm);
    } else {
        // 大规模并行，使用图分解
        return std::make_shared<GraphDecomposer>(1.05, comm);
    }
}

// ===== DomainDecompositionUtils实现 =====
namespace DomainDecompositionUtils {
    
void createElementAdjacency(std::vector<MeshElement>& elements, double tolerance) {
    // 简化实现：基于坐标距离创建邻接关系
    // 在实际应用中应该使用网格拓扑信息
    
    for (size_t i = 0; i < elements.size(); ++i) {
        for (size_t j = i + 1; j < elements.size(); ++j) {
            if (elements[i].coords.size() == elements[j].coords.size()) {
                double distance = 0.0;
                for (size_t k = 0; k < elements[i].coords.size(); ++k) {
                    double diff = elements[i].coords[k] - elements[j].coords[k];
                    distance += diff * diff;
                }
                distance = std::sqrt(distance);
                
                if (distance < tolerance) {
                    elements[i].neighbors.insert(j);
                    elements[j].neighbors.insert(i);
                }
            }
        }
    }
}

double calculateLoadBalance(const DomainDecompositionResult& result) {
    if (result.partitionSizes.empty()) return 1.0;
    
    int maxSize = *std::max_element(result.partitionSizes.begin(), result.partitionSizes.end());
    int minSize = *std::min_element(result.partitionSizes.begin(), result.partitionSizes.end());
    
    if (maxSize == 0) return 0.0;
    
    return static_cast<double>(minSize) / maxSize;
}

double calculateBoundaryLength(
    const DomainDecompositionResult& result,
    const std::vector<MeshElement>& elements) {
    
    double boundaryLength = 0.0;
    
    for (const auto& boundary : result.boundaryElements) {
        boundaryLength += boundary.second.size();
    }
    
    return boundaryLength;
}

DomainDecompositionResult optimizeDecomposition(
    const DomainDecompositionResult& result,
    const std::vector<MeshElement>& elements,
    int maxIterations) {
    
    // 简化实现：使用局部交换优化
    DomainDecompositionResult optimizedResult = result;
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        bool improved = false;
        
        // 尝试交换边界元素
        for (const auto& boundary : result.boundaryElements) {
            int partitionId = boundary.first;
            
            for (int elemId : boundary.second) {
                // 查找相邻分区
                for (int neighborId : elements[elemId].neighbors) {
                    if (neighborId < result.elementPartitions.size()) {
                        int neighborPartition = result.elementPartitions[neighborId];
                        
                        if (neighborPartition != partitionId) {
                            // 尝试交换元素到相邻分区
                            // 简化实现：在实际应用中需要更复杂的成本函数
                            
                            // 检查交换是否改善负载均衡
                            int currentPartitionSize = optimizedResult.partitionSizes[partitionId];
                            int neighborPartitionSize = optimizedResult.partitionSizes[neighborPartition];
                            
                            if (currentPartitionSize > neighborPartitionSize + 1) {
                                // 执行交换
                                auto& currentPartition = optimizedResult.partitionMap[partitionId];
                                auto& neighborPartition = optimizedResult.partitionMap[neighborPartition];
                                
                                // 从当前分区移除元素
                                currentPartition.erase(
                                    std::remove(currentPartition.begin(), currentPartition.end(), elemId),
                                    currentPartition.end()
                                );
                                
                                // 添加到相邻分区
                                neighborPartition.push_back(elemId);
                                
                                // 更新分区大小
                                optimizedResult.partitionSizes[partitionId]--;
                                optimizedResult.partitionSizes[neighborPartition]++;
                                
                                // 更新元素分区映射
                                optimizedResult.elementPartitions[elemId] = neighborPartition;
                                
                                improved = true;
                                break;
                            }
                        }
                    }
                }
                
                if (improved) break;
            }
            
            if (improved) break;
        }
        
        if (!improved) break;
    }
    
    // 重新计算负载均衡度
    optimizedResult.loadBalance = calculateLoadBalance(optimizedResult);
    
    return optimizedResult;
}

} // namespace DomainDecompositionUtils

} // namespace elmer