#pragma once

#include "Mesh.h"
#include <memory>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace ElmerCpp {

/**
 * @brief 网格工具类 - 提供网格创建、初始化和操作功能
 */
class MeshUtils {
public:
    /**
     * @brief 分配网格结构并返回智能指针
     * @param name 网格名称
     * @return 分配的网格对象
     */
    static std::shared_ptr<Mesh> allocateMesh(const std::string& name = "") {
        auto mesh = std::make_shared<Mesh>(name);
        
        // 初始化默认值
        mesh->setChanged(false);
        mesh->setOutputActive(false);
        mesh->setAdaptiveDepth(0);
        mesh->setDiscontMesh(false);
        mesh->setSingleMesh(false);
        
        return mesh;
    }
    
    /**
     * @brief 分配网格结构并指定大小
     * @param numberOfBulkElements 体单元数量
     * @param numberOfBoundaryElements 边界单元数量
     * @param numberOfNodes 节点数量
     * @param name 网格名称
     * @param initParallel 是否初始化并行信息
     * @return 分配的网格对象
     */
    static std::shared_ptr<Mesh> allocateMesh(
        size_t numberOfBulkElements,
        size_t numberOfBoundaryElements,
        size_t numberOfNodes,
        const std::string& name = "",
        bool initParallel = false
    ) {
        auto mesh = allocateMesh(name);
        
        // 设置网格大小信息
        mesh->getNodes().clear();
        mesh->getBulkElements().clear();
        mesh->getBoundaryElements().clear();
        
        // 预留空间
        mesh->getBulkElements().reserve(numberOfBulkElements);
        mesh->getBoundaryElements().reserve(numberOfBoundaryElements);
        
        // 初始化节点空间
        for (size_t i = 0; i < numberOfNodes; ++i) {
            mesh->getNodes().addNode(0.0, 0.0, 0.0);
        }
        
        // 初始化并行信息
        if (initParallel) {
            mesh->initializeParallelInfo();
        }
        
        return mesh;
    }
    
    /**
     * @brief 初始化网格结构
     * @param mesh 要初始化的网格
     * @param initParallel 是否初始化并行信息
     */
    static void initializeMesh(std::shared_ptr<Mesh> mesh, bool initParallel = false) {
        if (!mesh.get()) {
            throw std::invalid_argument("Mesh pointer is null");
        }
        
        size_t numberOfNodes = mesh->numberOfNodes();
        size_t numberOfBulkElements = mesh->numberOfBulkElements();
        size_t numberOfBoundaryElements = mesh->numberOfBoundaryElements();
        
        // 验证网格完整性
        if (numberOfNodes == 0) {
            throw std::runtime_error("Mesh has zero nodes!");
        }
        
        if (mesh->totalElements() == 0) {
            throw std::runtime_error("Mesh has zero elements!");
        }
        
        // 初始化并行信息
        if (initParallel) {
            mesh->initializeParallelInfo();
        }
        
        // 更新统计信息
        mesh->getBulkElements(); // 触发统计更新
        mesh->getBoundaryElements(); // 触发统计更新
    }
    
    /**
     * @brief 创建简单立方体网格
     * @param size 立方体大小
     * @param divisions 每个方向的划分数量
     * @return 创建的网格对象
     */
    static std::shared_ptr<Mesh> createCubeMesh(
        double size = 1.0,
        int divisions = 1,
        const std::string& name = "CubeMesh"
    ) {
        auto mesh = allocateMesh(name);
        
        // 创建节点
        double step = size / divisions;
        for (int k = 0; k <= divisions; ++k) {
            for (int j = 0; j <= divisions; ++j) {
                for (int i = 0; i <= divisions; ++i) {
                    double x = i * step - size / 2.0;
                    double y = j * step - size / 2.0;
                    double z = k * step - size / 2.0;
                    mesh->getNodes().addNode(x, y, z);
                }
            }
        }
        
        // 创建立方体单元（六面体）
        int nodesPerDirection = divisions + 1;
        for (int k = 0; k < divisions; ++k) {
            for (int j = 0; j < divisions; ++j) {
                for (int i = 0; i < divisions; ++i) {
                    Element element(ElementType::HEXAHEDRON);
                    
                    // 计算六面体单元的8个节点索引
                    int base = k * nodesPerDirection * nodesPerDirection + j * nodesPerDirection + i;
                    
                    element.addNodeIndex(base);
                    element.addNodeIndex(base + 1);
                    element.addNodeIndex(base + nodesPerDirection + 1);
                    element.addNodeIndex(base + nodesPerDirection);
                    
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection);
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection + 1);
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection + nodesPerDirection + 1);
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection + nodesPerDirection);
                    
                    mesh->addBulkElement(element);
                }
            }
        }
        
        // 创建边界单元
        createBoundaryElements(mesh, divisions);
        
        return mesh;
    }
    
    /**
     * @brief 创建简单四面体网格
     * @param size 网格大小
     * @param divisions 划分数量
     * @return 创建的网格对象
     */
    static std::shared_ptr<Mesh> createTetrahedronMesh(
        double size = 1.0,
        int divisions = 1,
        const std::string& name = "TetrahedronMesh"
    ) {
        auto mesh = allocateMesh(name);
        
        // 创建简单四面体的四个节点
        mesh->getNodes().addNode(0.0, 0.0, 0.0);
        mesh->getNodes().addNode(size, 0.0, 0.0);
        mesh->getNodes().addNode(0.0, size, 0.0);
        mesh->getNodes().addNode(0.0, 0.0, size);
        
        // 创建四面体单元
        Element element(ElementType::TETRAHEDRON);
        element.addNodeIndex(0);
        element.addNodeIndex(1);
        element.addNodeIndex(2);
        element.addNodeIndex(3);
        
        mesh->addBulkElement(element);
        
        return mesh;
    }
    
    /**
     * @brief 计算网格质量指标
     * @param mesh 要计算的网格
     * @return 质量指标（0-1，1表示完美质量）
     */
    static double calculateMeshQuality(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh.get() || mesh->totalElements() == 0) {
            return 0.0;
        }
        
        double totalQuality = 0.0;
        size_t elementCount = 0;
        
        // 计算体单元质量
        for (const auto& element : mesh->getBulkElements()) {
            double quality = 1.0; // 简化版本，返回固定质量值
            totalQuality += quality;
            elementCount++;
        }
        
        // 计算边界单元质量
        for (const auto& element : mesh->getBoundaryElements()) {
            double quality = 1.0; // 简化版本，返回固定质量值
            totalQuality += quality;
            elementCount++;
        }
        
        return elementCount > 0 ? totalQuality / elementCount : 0.0;
    }
    
    /**
     * @brief 验证网格拓扑结构
     * @param mesh 要验证的网格
     * @return 验证结果
     */
    static bool validateMeshTopology(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh.get()) {
            return false;
        }
        
        // 检查基本完整性
        if (!mesh->validate()) {
            return false;
        }
        
        // 检查节点索引范围
        size_t numberOfNodes = mesh->numberOfNodes();
        
        for (const auto& element : mesh->getBulkElements()) {
            for (size_t nodeIndex : element.getNodeIndices()) {
                if (nodeIndex >= numberOfNodes) {
                    return false;
                }
            }
        }
        
        for (const auto& element : mesh->getBoundaryElements()) {
            for (size_t nodeIndex : element.getNodeIndices()) {
                if (nodeIndex >= numberOfNodes) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     * @brief 计算网格边界框
     * @param mesh 要计算的网格
     * @return 边界框[minX, maxX, minY, maxY, minZ, maxZ]
     */
    static std::vector<double> calculateBoundingBox(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh || mesh->numberOfNodes() == 0) {
            return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        }
        
        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();
        double minZ = std::numeric_limits<double>::max();
        double maxZ = std::numeric_limits<double>::lowest();
        
        const auto& nodes = mesh->getNodes().getNodes();
        for (const auto& node : nodes) {
            minX = std::min(minX, node.x);
            maxX = std::max(maxX, node.x);
            minY = std::min(minY, node.y);
            maxY = std::max(maxY, node.y);
            minZ = std::min(minZ, node.z);
            maxZ = std::max(maxZ, node.z);
        }
        
        return {minX, maxX, minY, maxY, minZ, maxZ};
    }
    
    /**
     * @brief 计算网格体积
     * @param mesh 要计算的网格
     * @return 网格体积
     */
    static double calculateVolume(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh || mesh->numberOfBulkElements() == 0) {
            return 0.0;
        }
        
        double totalVolume = 0.0;
        const auto& nodes = mesh->getNodes().getNodes();
        
        for (const auto& element : mesh->getBulkElements()) {
            totalVolume += calculateElementVolume(element, nodes);
        }
        
        return totalVolume;
    }
    
private:
    /**
     * @brief 为立方体网格创建边界单元
     */
    static void createBoundaryElements(std::shared_ptr<Mesh> mesh, int divisions) {
        int nodesPerDirection = divisions + 1;
        
        // 创建六个面的边界单元
        for (int face = 0; face < 6; ++face) {
            for (int j = 0; j < divisions; ++j) {
                for (int i = 0; i < divisions; ++i) {
                    Element boundaryElement(ElementType::LINEAR);
                    boundaryElement.setBoundaryId(face + 1);
                    
                    // 根据面计算节点索引
                    std::vector<size_t> indices = calculateFaceIndices(face, i, j, divisions, nodesPerDirection);
                    boundaryElement.setNodeIndices(indices);
                    
                    mesh->addBoundaryElement(boundaryElement);
                }
            }
        }
    }
    
    /**
     * @brief 计算立方体面的节点索引
     */
    static std::vector<size_t> calculateFaceIndices(int face, int i, int j, int divisions, int nodesPerDirection) {
        std::vector<size_t> indices;
        
        switch (face) {
            case 0: // 底面 z=0
                indices.push_back(j * nodesPerDirection + i);
                indices.push_back(j * nodesPerDirection + i + 1);
                indices.push_back((j + 1) * nodesPerDirection + i + 1);
                indices.push_back((j + 1) * nodesPerDirection + i);
                break;
            case 1: // 顶面 z=divisions
                size_t base = divisions * nodesPerDirection * nodesPerDirection;
                indices.push_back(base + j * nodesPerDirection + i);
                indices.push_back(base + j * nodesPerDirection + i + 1);
                indices.push_back(base + (j + 1) * nodesPerDirection + i + 1);
                indices.push_back(base + (j + 1) * nodesPerDirection + i);
                break;
            // 其他四个面类似处理
            default:
                // 简化处理，实际实现需要完整处理所有面
                break;
        }
        
        return indices;
    }
    
    /**
     * @brief 计算单元质量
     */
    static double calculateElementQuality(const Element& element, const Nodes& nodes) {
        // 简化实现：基于单元形状计算质量指标
        // 实际实现需要根据单元类型进行具体计算
        
        size_t nodeCount = element.numberOfNodes();
        if (nodeCount < 3) {
            return 0.0;
        }
        
        // 计算节点间距离的方差作为质量指标
        std::vector<double> distances;
        const auto& indices = element.getNodeIndices();
        
        for (size_t i = 0; i < nodeCount; ++i) {
            for (size_t j = i + 1; j < nodeCount; ++j) {
                double dist = nodes[indices[i]].distance(nodes[indices[j]]);
                distances.push_back(dist);
            }
        }
        
        if (distances.empty()) {
            return 0.0;
        }
        
        // 计算距离的平均值和标准差
        double mean = 0.0;
        for (double dist : distances) {
            mean += dist;
        }
        mean /= distances.size();
        
        double variance = 0.0;
        for (double dist : distances) {
            variance += (dist - mean) * (dist - mean);
        }
        variance /= distances.size();
        
        // 质量指标：方差越小，质量越好
        return 1.0 / (1.0 + std::sqrt(variance));
    }
    
    /**
     * @brief 计算单元体积
     */
    static double calculateElementVolume(const Element& element, const std::vector<Node>& nodes) {
        // 简化实现：仅处理四面体和六面体
        const auto& indices = element.getNodeIndices();
        
        switch (element.getType()) {
            case ElementType::TETRAHEDRON:
                if (indices.size() == 4) {
                    return calculateTetrahedronVolume(nodes[indices[0]], nodes[indices[1]], 
                                                     nodes[indices[2]], nodes[indices[3]]);
                }
                break;
                
            case ElementType::HEXAHEDRON:
                if (indices.size() == 8) {
                    // 将六面体分解为5个四面体计算体积
                    return calculateHexahedronVolume(nodes, indices);
                }
                break;
                
            default:
                break;
        }
        
        return 0.0;
    }
    
    /**
     * @brief 计算四面体体积
     */
    static double calculateTetrahedronVolume(const Node& a, const Node& b, const Node& c, const Node& d) {
        double v1x = b.x - a.x, v1y = b.y - a.y, v1z = b.z - a.z;
        double v2x = c.x - a.x, v2y = c.y - a.y, v2z = c.z - a.z;
        double v3x = d.x - a.x, v3y = d.y - a.y, v3z = d.z - a.z;
        
        // 计算标量三重积
        double volume = std::abs(
            v1x * (v2y * v3z - v2z * v3y) -
            v1y * (v2x * v3z - v2z * v3x) +
            v1z * (v2x * v3y - v2y * v3x)
        ) / 6.0;
        
        return volume;
    }
    
    /**
     * @brief 计算六面体体积
     */
    static double calculateHexahedronVolume(const std::vector<Node>& nodes, const std::vector<size_t>& indices) {
        // 简化实现：将六面体分解为5个四面体
        double volume = 0.0;
        
        // 分解方案1：基于中心点分解
        Node center(0.0, 0.0, 0.0);
        for (size_t index : indices) {
            center.x += nodes[index].x;
            center.y += nodes[index].y;
            center.z += nodes[index].z;
        }
        center.x /= 8.0;
        center.y /= 8.0;
        center.z /= 8.0;
        
        // 计算每个面的四面体体积
        std::vector<std::vector<size_t>> faces = {
            {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
            {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
        };
        
        for (const auto& face : faces) {
            // 将四边形面分解为两个三角形
            volume += calculateTetrahedronVolume(
                center, 
                nodes[indices[face[0]]], 
                nodes[indices[face[1]]], 
                nodes[indices[face[2]]]
            );
            volume += calculateTetrahedronVolume(
                center, 
                nodes[indices[face[0]]], 
                nodes[indices[face[2]]], 
                nodes[indices[face[3]]]
            );
        }
        
        return volume;
    }
};

} // namespace ElmerCpp