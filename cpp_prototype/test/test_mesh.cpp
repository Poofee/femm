/**
 * @file test_mesh.cpp
 * @brief Mesh模块测试文件
 * 
 * 测试Mesh模块的核心功能，包括节点管理、元素操作、网格验证等。
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "../src/Mesh.h"

using namespace elmer;

/**
 * @brief 自定义断言宏，用于测试
 */
#define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "断言失败: " << message << std::endl; \
            std::cerr << "位置: " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

/**
 * @brief 自定义近似相等断言宏
 */
#define ASSERT_NEAR(a, b, tolerance) \
    do { \
        if (std::abs((a) - (b)) > (tolerance)) { \
            std::cerr << "近似相等断言失败: " << (a) << " != " << (b) << std::endl; \
            std::cerr << "容差: " << (tolerance) << std::endl; \
            std::cerr << "位置: " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

/**
 * @brief 测试节点基本功能
 */
bool TestNodeBasic() {
    std::cout << "测试节点基本功能..." << std::endl;
    
    // 测试默认构造函数
    Node node1;
    ASSERT(node1.x == 0.0 && node1.y == 0.0 && node1.z == 0.0, "默认构造函数应初始化坐标为零");
    
    // 测试参数化构造函数
    Node node2(1.0, 2.0, 3.0);
    ASSERT(node2.x == 1.0 && node2.y == 2.0 && node2.z == 3.0, "参数化构造函数应正确设置坐标");
    
    // 测试距离计算
    Node node3(4.0, 5.0, 6.0);
    double distance = node2.distance(node3);
    double expected = std::sqrt(9.0 + 9.0 + 9.0); // sqrt(3^2 + 3^2 + 3^2)
    ASSERT_NEAR(distance, expected, 1e-12);
    
    std::cout << "节点基本功能测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试节点集合管理
 */
bool TestNodesCollection() {
    std::cout << "测试节点集合管理..." << std::endl;
    
    Nodes nodes;
    
    // 测试添加节点
    nodes.addNode(Node(1.0, 1.0, 1.0));
    nodes.addNode(2.0, 2.0, 2.0);
    nodes.addNode(3.0, 3.0, 3.0);
    
    ASSERT(nodes.numberOfNodes() == 3, "节点数量应为3");
    
    // 测试索引访问
    Node& node = nodes[1];
    ASSERT(node.x == 2.0 && node.y == 2.0 && node.z == 2.0, "索引访问应返回正确的节点");
    
    // 测试常量访问
    const Nodes& constNodes = nodes;
    const Node& constNode = constNodes[0];
    ASSERT(constNode.x == 1.0 && constNode.y == 1.0 && constNode.z == 1.0, "常量索引访问应返回正确的节点");
    
    // 测试清空
    nodes.clear();
    ASSERT(nodes.numberOfNodes() == 0, "清空后节点数量应为0");
    
    std::cout << "节点集合管理测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试元素基本功能
 */
bool TestElementBasic() {
    std::cout << "测试元素基本功能..." << std::endl;
    
    // 创建节点集合
    Nodes nodes;
    nodes.addNode(0.0, 0.0, 0.0);
    nodes.addNode(1.0, 0.0, 0.0);
    nodes.addNode(1.0, 1.0, 0.0);
    nodes.addNode(0.0, 1.0, 0.0);
    
    // 创建四边形元素
    Element element(ElementType::LINEAR, 1);
    element.addNodeIndex(0);
    element.addNodeIndex(1);
    element.addNodeIndex(2);
    element.addNodeIndex(3);
    
    ASSERT(element.numberOfNodes() == 4, "元素节点数应为4");
    ASSERT(element.getType() == ElementType::LINEAR, "元素类型应为LINEAR");
    ASSERT(element.getBodyId() == 1, "体ID应为1");
    ASSERT(element.isBulkElement(), "应为体元素");
    
    // 测试质心计算
    Node centroid = element.calculateCentroid(nodes);
    ASSERT_NEAR(centroid.x, 0.5, 1e-12);
    ASSERT_NEAR(centroid.y, 0.5, 1e-12);
    ASSERT_NEAR(centroid.z, 0.0, 1e-12);
    
    // 测试边界框计算
    Node minCorner, maxCorner;
    element.getBoundingBox(nodes, minCorner, maxCorner);
    ASSERT_NEAR(minCorner.x, 0.0, 1e-12);
    ASSERT_NEAR(minCorner.y, 0.0, 1e-12);
    ASSERT_NEAR(maxCorner.x, 1.0, 1e-12);
    ASSERT_NEAR(maxCorner.y, 1.0, 1e-12);
    
    std::cout << "元素基本功能测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试网格基本功能
 */
bool TestMeshBasic() {
    std::cout << "测试网格基本功能..." << std::endl;
    
    Mesh mesh("TestMesh");
    
    // 添加节点
    mesh.getNodes().addNode(0.0, 0.0, 0.0);
    mesh.getNodes().addNode(1.0, 0.0, 0.0);
    mesh.getNodes().addNode(1.0, 1.0, 0.0);
    mesh.getNodes().addNode(0.0, 1.0, 0.0);
    
    ASSERT(mesh.numberOfNodes() == 4, "网格节点数应为4");
    
    // 添加体元素
    Element bulkElement(ElementType::LINEAR, 1);
    bulkElement.addNodeIndex(0);
    bulkElement.addNodeIndex(1);
    bulkElement.addNodeIndex(2);
    bulkElement.addNodeIndex(3);
    mesh.addBulkElement(bulkElement);
    
    ASSERT(mesh.numberOfBulkElements() == 1, "体元素数应为1");
    
    // 添加边界元素
    Element boundaryElement(ElementType::LINEAR, 0);
    boundaryElement.setBoundaryId(1);
    boundaryElement.addNodeIndex(0);
    boundaryElement.addNodeIndex(1);
    mesh.addBoundaryElement(boundaryElement);
    
    ASSERT(mesh.numberOfBoundaryElements() == 1, "边界元素数应为1");
    ASSERT(mesh.totalElements() == 2, "总元素数应为2");
    
    // 测试网格验证
    ASSERT(mesh.validate(), "网格验证应通过");
    
    // 测试元素质心计算
    Node centroid = mesh.calculateElementCentroid(bulkElement);
    ASSERT_NEAR(centroid.x, 0.5, 1e-12);
    ASSERT_NEAR(centroid.y, 0.5, 1e-12);
    ASSERT_NEAR(centroid.z, 0.0, 1e-12);
    
    // 测试边界框计算
    Node minCorner, maxCorner;
    mesh.getBoundingBox(minCorner, maxCorner);
    ASSERT_NEAR(minCorner.x, 0.0, 1e-12);
    ASSERT_NEAR(minCorner.y, 0.0, 1e-12);
    ASSERT_NEAR(maxCorner.x, 1.0, 1e-12);
    ASSERT_NEAR(maxCorner.y, 1.0, 1e-12);
    
    std::cout << "网格基本功能测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试网格统计功能
 */
bool TestMeshStatistics() {
    std::cout << "测试网格统计功能..." << std::endl;
    
    Mesh mesh("StatisticsTest");
    
    // 添加节点
    mesh.getNodes().addNode(0.0, 0.0, 0.0);
    mesh.getNodes().addNode(1.0, 0.0, 0.0);
    mesh.getNodes().addNode(1.0, 1.0, 0.0);
    mesh.getNodes().addNode(0.0, 1.0, 0.0);
    
    // 添加不同节点数的元素
    Element element4(ElementType::LINEAR, 1);
    element4.addNodeIndex(0);
    element4.addNodeIndex(1);
    element4.addNodeIndex(2);
    element4.addNodeIndex(3);
    mesh.addBulkElement(element4);
    
    Element element3(ElementType::LINEAR, 1);
    element3.addNodeIndex(0);
    element3.addNodeIndex(1);
    element3.addNodeIndex(2);
    mesh.addBulkElement(element3);
    
    ASSERT(mesh.getMaxElementNodes() == 4, "最大元素节点数应为4");
    
    // 测试统计信息输出（不会崩溃即可）
    mesh.getMeshStatistics(std::cout);
    
    std::cout << "网格统计功能测试通过" << std::endl;
    return true;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始Mesh模块测试..." << std::endl;
    
    bool allPassed = true;
    
    try {
        allPassed &= TestNodeBasic();
        allPassed &= TestNodesCollection();
        allPassed &= TestElementBasic();
        allPassed &= TestMeshBasic();
        allPassed &= TestMeshStatistics();
        
        if (allPassed) {
            std::cout << "\n✅ 所有Mesh模块测试通过！" << std::endl;
            return 0;
        } else {
            std::cout << "\n❌ 部分Mesh模块测试失败！" << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "测试异常: " << e.what() << std::endl;
        return 1;
    }
}