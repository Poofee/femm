#include "../src/Mesh.h"
#include "../src/MeshUtils.h"
#include "../src/MeshIO.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace elmer;

/**
 * @brief 测试节点基本功能
 */
void testNode() {
    std::cout << "测试节点功能..." << std::endl;
    
    // 测试默认构造函数
    Node node1;
    assert(node1.x == 0.0);
    assert(node1.y == 0.0);
    assert(node1.z == 0.0);
    
    // 测试参数化构造函数
    Node node2(1.0, 2.0, 3.0);
    assert(node2.x == 1.0);
    assert(node2.y == 2.0);
    assert(node2.z == 3.0);
    
    // 测试距离计算
    Node node3(4.0, 5.0, 6.0);
    double distance = node2.distance(node3);
    double expected = std::sqrt(9.0 + 9.0 + 9.0); // sqrt(3^2 + 3^2 + 3^2)
    assert(std::abs(distance - expected) < 1e-10);
    
    std::cout << "节点测试通过!" << std::endl;
}

/**
 * @brief 测试节点集合功能
 */
void testNodes() {
    std::cout << "测试节点集合功能..." << std::endl;
    
    Nodes nodes;
    
    // 测试添加节点
    nodes.addNode(0.0, 0.0, 0.0);
    nodes.addNode(1.0, 0.0, 0.0);
    nodes.addNode(0.0, 1.0, 0.0);
    
    assert(nodes.numberOfNodes() == 3);
    
    // 测试索引访问
    assert(nodes[0].x == 0.0);
    assert(nodes[1].y == 0.0);
    assert(nodes[2].z == 0.0);
    
    // 测试修改节点
    nodes[0].x = 5.0;
    assert(nodes[0].x == 5.0);
    
    // 测试清空
    nodes.clear();
    assert(nodes.numberOfNodes() == 0);
    
    std::cout << "节点集合测试通过!" << std::endl;
}

/**
 * @brief 测试单元功能
 */
void testElement() {
    std::cout << "测试单元功能..." << std::endl;
    
    // 测试默认构造函数
    Element element1;
    assert(element1.getType() == ElementType::LINEAR);
    assert(element1.getBodyId() == 0);
    assert(element1.getBoundaryId() == 0);
    assert(element1.numberOfNodes() == 0);
    assert(element1.isBulkElement());
    assert(!element1.isBoundaryElement());
    
    // 测试参数化构造函数
    Element element2(ElementType::TETRAHEDRON, 1);
    assert(element2.getType() == ElementType::TETRAHEDRON);
    assert(element2.getBodyId() == 1);
    assert(element2.isBulkElement());
    
    // 测试边界单元
    Element boundary(ElementType::LINEAR);
    boundary.setBoundaryId(2);
    assert(boundary.isBoundaryElement());
    assert(!boundary.isBulkElement());
    
    // 测试节点索引操作
    element2.addNodeIndex(0);
    element2.addNodeIndex(1);
    element2.addNodeIndex(2);
    element2.addNodeIndex(3);
    
    assert(element2.numberOfNodes() == 4);
    const auto& indices = element2.getNodeIndices();
    assert(indices.size() == 4);
    assert(indices[0] == 0);
    assert(indices[3] == 3);
    
    std::cout << "单元测试通过!" << std::endl;
}

/**
 * @brief 测试网格基本功能
 */
void testMesh() {
    std::cout << "测试网格功能..." << std::endl;
    
    Mesh mesh("TestMesh");
    
    // 测试名称
    assert(mesh.getName() == "TestMesh");
    mesh.setName("NewMesh");
    assert(mesh.getName() == "NewMesh");
    
    // 测试节点操作
    mesh.getNodes().addNode(0.0, 0.0, 0.0);
    mesh.getNodes().addNode(1.0, 0.0, 0.0);
    mesh.getNodes().addNode(0.0, 1.0, 0.0);
    
    assert(mesh.numberOfNodes() == 3);
    
    // 测试体单元操作
    Element bulkElement(ElementType::LINEAR);
    bulkElement.addNodeIndex(0);
    bulkElement.addNodeIndex(1);
    bulkElement.setBodyId(1);
    
    mesh.addBulkElement(bulkElement);
    assert(mesh.numberOfBulkElements() == 1);
    
    // 测试边界单元操作
    Element boundaryElement(ElementType::LINEAR);
    boundaryElement.addNodeIndex(0);
    boundaryElement.addNodeIndex(2);
    boundaryElement.setBoundaryId(1);
    
    mesh.addBoundaryElement(boundaryElement);
    assert(mesh.numberOfBoundaryElements() == 1);
    
    // 测试总单元数
    assert(mesh.totalElements() == 2);
    
    // 测试网格验证
    assert(mesh.validate());
    
    // 测试清空
    mesh.clear();
    assert(mesh.numberOfNodes() == 0);
    assert(mesh.numberOfBulkElements() == 0);
    assert(mesh.numberOfBoundaryElements() == 0);
    
    std::cout << "网格测试通过!" << std::endl;
}

/**
 * @brief 测试网格工具功能
 */
void testMeshUtils() {
    std::cout << "测试网格工具功能..." << std::endl;
    
    // 测试网格分配
    auto mesh1 = MeshUtils::allocateMesh("TestMesh1");
    assert(mesh1->getName() == "TestMesh1");
    assert(!mesh1->isChanged());
    assert(!mesh1->isOutputActive());
    
    // 测试带大小参数的网格分配
    auto mesh2 = MeshUtils::allocateMesh(10, 5, 20, "TestMesh2");
    assert(mesh2->numberOfNodes() == 20);
    assert(mesh2->getBulkElements().capacity() >= 10);
    assert(mesh2->getBoundaryElements().capacity() >= 5);
    
    // 测试立方体网格创建
    auto cubeMesh = MeshUtils::createCubeMesh(1.0, 2, "CubeMesh");
    assert(cubeMesh->numberOfNodes() == 27); // (2+1)^3 = 27
    assert(cubeMesh->numberOfBulkElements() == 8); // 2^3 = 8
    assert(cubeMesh->numberOfBoundaryElements() > 0);
    
    // 测试四面体网格创建
    auto tetraMesh = MeshUtils::createTetrahedronMesh(1.0, 1, "TetraMesh");
    assert(tetraMesh->numberOfNodes() == 4);
    assert(tetraMesh->numberOfBulkElements() == 1);
    
    // 测试网格质量计算
    double quality = MeshUtils::calculateMeshQuality(tetraMesh);
    assert(quality > 0.0 && quality <= 1.0);
    
    // 测试网格拓扑验证
    assert(MeshUtils::validateMeshTopology(tetraMesh));
    
    // 测试边界框计算
    auto bbox = MeshUtils::calculateBoundingBox(tetraMesh);
    assert(bbox.size() == 6);
    assert(bbox[0] <= bbox[1]); // minX <= maxX
    assert(bbox[2] <= bbox[3]); // minY <= maxY
    assert(bbox[4] <= bbox[5]); // minZ <= maxZ
    
    // 测试体积计算
    double volume = MeshUtils::calculateVolume(tetraMesh);
    assert(volume > 0.0);
    
    std::cout << "网格工具测试通过!" << std::endl;
}

/**
 * @brief 测试网格IO功能
 */
void testMeshIO() {
    std::cout << "测试网格IO功能..." << std::endl;
    
    // 创建测试网格
    auto testMesh = MeshIO::createTestMesh();
    assert(testMesh->numberOfNodes() == 4);
    assert(testMesh->numberOfBulkElements() == 1);
    assert(testMesh->numberOfBoundaryElements() == 4);
    
    // 测试网格保存和加载
    const std::string testFile = "test_mesh.mesh";
    
    try {
        // 保存网格
        MeshIO::saveMesh(testMesh, testFile);
        
        // 加载网格
        auto loadedMesh = MeshIO::loadMesh(testFile);
        
        // 验证加载的网格
        assert(loadedMesh->numberOfNodes() == testMesh->numberOfNodes());
        assert(loadedMesh->numberOfBulkElements() == testMesh->numberOfBulkElements());
        assert(loadedMesh->numberOfBoundaryElements() == testMesh->numberOfBoundaryElements());
        
        // 验证节点坐标
        const auto& originalNodes = testMesh->getNodes().getNodes();
        const auto& loadedNodes = loadedMesh->getNodes().getNodes();
        
        for (size_t i = 0; i < originalNodes.size(); ++i) {
            assert(std::abs(originalNodes[i].x - loadedNodes[i].x) < 1e-10);
            assert(std::abs(originalNodes[i].y - loadedNodes[i].y) < 1e-10);
            assert(std::abs(originalNodes[i].z - loadedNodes[i].z) < 1e-10);
        }
        
        // 测试VTK导出
        MeshIO::exportToVTK(testMesh, "test_mesh.vtk");
        
        std::cout << "网格IO测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "网格IO测试失败: " << e.what() << std::endl;
        throw;
    }
}

/**
 * @brief 测试并行信息功能
 */
void testParallelInfo() {
    std::cout << "测试并行信息功能..." << std::endl;
    
    ParallelInfo parallelInfo;
    
    // 测试初始化
    parallelInfo.initialize(10);
    assert(parallelInfo.globalDOFs.size() == 10);
    assert(parallelInfo.gInterface.size() == 10);
    assert(parallelInfo.neighbourList.size() == 10);
    
    // 测试清空
    parallelInfo.clear();
    assert(parallelInfo.globalDOFs.empty());
    assert(parallelInfo.gInterface.empty());
    assert(parallelInfo.neighbourList.empty());
    assert(parallelInfo.numberOfIfDOFs == 0);
    
    std::cout << "并行信息测试通过!" << std::endl;
}

/**
 * @brief 综合测试：创建复杂网格并验证
 */
void testComplexMesh() {
    std::cout << "测试复杂网格功能..." << std::endl;
    
    // 创建复杂立方体网格
    auto complexMesh = MeshUtils::createCubeMesh(2.0, 3, "ComplexCube");
    
    // 验证网格统计信息
    assert(complexMesh->numberOfNodes() == 64); // (3+1)^3 = 64
    assert(complexMesh->numberOfBulkElements() == 27); // 3^3 = 27
    assert(complexMesh->numberOfBoundaryElements() == 54); // 6个面 * 3^2 = 54
    
    // 验证网格完整性
    assert(complexMesh->validate());
    
    // 验证网格质量
    double quality = MeshUtils::calculateMeshQuality(complexMesh);
    assert(quality > 0.0 && quality <= 1.0);
    
    // 验证边界框
    auto bbox = MeshUtils::calculateBoundingBox(complexMesh);
    assert(std::abs(bbox[0] - (-1.0)) < 1e-10); // minX = -1.0
    assert(std::abs(bbox[1] - 1.0) < 1e-10);    // maxX = 1.0
    
    // 验证体积
    double volume = MeshUtils::calculateVolume(complexMesh);
    assert(std::abs(volume - 8.0) < 0.1); // 2^3 = 8
    
    std::cout << "复杂网格测试通过!" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始ElmerCpp网格模块测试..." << std::endl;
    std::cout << "============================" << std::endl;
    
    try {
        testNode();
        testNodes();
        testElement();
        testMesh();
        testParallelInfo();
        testMeshUtils();
        testMeshIO();
        testComplexMesh();
        
        std::cout << "============================" << std::endl;
        std::cout << "所有网格模块测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "测试失败: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}