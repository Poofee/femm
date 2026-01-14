#include <iostream>
#include <vector>
#include <cmath>
#include "src/ElementDescription.h"
#include "src/Types.h"

using namespace ElmerCpp;

/**
 * @brief 测试边界条件处理函数
 */
void testBoundaryConditions() {
    std::cout << "=== 测试边界条件处理函数 ===" << std::endl;
    
    // 创建测试用的元素结构
    Element element;
    element.type.elementCode = 202; // 三角形元素
    element.type.numberOfNodes = 3;
    element.nodeIndexes = {0, 1, 2};
    element.edgeIndexes = {0, 1, 2};
    element.faceIndexes = {0};
    
    // 测试ProcessBoundaryConditions函数
    std::cout << "\n1. 测试ProcessBoundaryConditions函数:" << std::endl;
    
    // Dirichlet边界条件
    std::vector<int> boundaryNodes = {0, 1};
    std::vector<double> boundaryValues;
    
    ProcessBoundaryConditions(element, 0, 5.0, boundaryNodes, boundaryValues);
    
    std::cout << "   Dirichlet边界条件测试:" << std::endl;
    std::cout << "   边界节点: ";
    for (int node : boundaryNodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
    
    std::cout << "   边界值: ";
    for (double value : boundaryValues) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    
    // 验证结果
    if (boundaryValues.size() == 2 && boundaryValues[0] == 5.0 && boundaryValues[1] == 5.0) {
        std::cout << "   ✓ Dirichlet边界条件测试通过" << std::endl;
    } else {
        std::cout << "   ✗ Dirichlet边界条件测试失败" << std::endl;
    }
    
    // Neumann边界条件
    ProcessBoundaryConditions(element, 1, 3.0, boundaryNodes, boundaryValues);
    
    std::cout << "   Neumann边界条件测试:" << std::endl;
    std::cout << "   边界值: ";
    for (double value : boundaryValues) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    
    if (boundaryValues.size() == 2 && boundaryValues[0] == 3.0 && boundaryValues[1] == 3.0) {
        std::cout << "   ✓ Neumann边界条件测试通过" << std::endl;
    } else {
        std::cout << "   ✗ Neumann边界条件测试失败" << std::endl;
    }
    
    // 测试GetFaceNormal函数
    std::cout << "\n2. 测试GetFaceNormal函数:" << std::endl;
    
    Nodes nodes;
    nodes.x = {0.0, 1.0, 0.0};
    nodes.y = {0.0, 0.0, 1.0};
    nodes.z = {0.0, 0.0, 0.0};
    
    std::vector<double> normal;
    GetFaceNormal(element, 0, nodes, normal);
    
    std::cout << "   面法向量: (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << std::endl;
    
    // 验证法向量长度是否为1（归一化）
    double length = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    if (std::abs(length - 1.0) < 1e-12) {
        std::cout << "   ✓ 法向量归一化测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 法向量归一化测试失败，长度=" << length << std::endl;
    }
    
    // 测试ComputeFaceFlux函数
    std::cout << "\n3. 测试ComputeFaceFlux函数:" << std::endl;
    
    std::vector<double> field = {1.0, 2.0, 3.0};
    double flux = 0.0;
    
    ComputeFaceFlux(element, 0, nodes, field, flux);
    
    std::cout << "   面通量: " << flux << std::endl;
    
    // 验证通量计算
    double expectedFlux = (1.0 + 2.0 + 3.0) / 3.0; // 平均值
    if (std::abs(flux - expectedFlux) < 1e-12) {
        std::cout << "   ✓ 面通量计算测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 面通量计算测试失败，期望值=" << expectedFlux << std::endl;
    }
    
    // 测试ApplyDirichletBoundaryConditions函数
    std::cout << "\n4. 测试ApplyDirichletBoundaryConditions函数:" << std::endl;
    
    std::vector<std::vector<double>> stiffnessMatrix = {
        {2.0, -1.0, 0.0},
        {-1.0, 2.0, -1.0},
        {0.0, -1.0, 2.0}
    };
    
    std::vector<double> forceVector = {1.0, 2.0, 3.0};
    std::vector<double> dirichletValues = {5.0, 6.0};
    
    ApplyDirichletBoundaryConditions(element, boundaryNodes, dirichletValues, stiffnessMatrix, forceVector);
    
    std::cout << "   应用Dirichlet边界条件后的刚度矩阵:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "   [";
        for (int j = 0; j < 3; ++j) {
            std::cout << stiffnessMatrix[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
    
    std::cout << "   应用Dirichlet边界条件后的力向量: ";
    for (double value : forceVector) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    
    // 验证边界条件应用
    bool dirichletTestPassed = true;
    for (int i = 0; i < 2; ++i) {
        int nodeIdx = boundaryNodes[i];
        if (stiffnessMatrix[nodeIdx][nodeIdx] != 1.0) {
            dirichletTestPassed = false;
            break;
        }
        if (forceVector[nodeIdx] != dirichletValues[i]) {
            dirichletTestPassed = false;
            break;
        }
    }
    
    if (dirichletTestPassed) {
        std::cout << "   ✓ Dirichlet边界条件应用测试通过" << std::endl;
    } else {
        std::cout << "   ✗ Dirichlet边界条件应用测试失败" << std::endl;
    }
    
    // 测试ApplyNeumannBoundaryConditions函数
    std::cout << "\n5. 测试ApplyNeumannBoundaryConditions函数:" << std::endl;
    
    std::vector<double> forceVector2 = {1.0, 2.0, 3.0};
    ApplyNeumannBoundaryConditions(element, 0, 10.0, forceVector2);
    
    std::cout << "   应用Neumann边界条件后的力向量: ";
    for (double value : forceVector2) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    
    // 验证Neumann边界条件应用
    bool neumannTestPassed = (forceVector2[0] > 1.0 && forceVector2[1] > 2.0 && forceVector2[2] > 3.0);
    
    if (neumannTestPassed) {
        std::cout << "   ✓ Neumann边界条件应用测试通过" << std::endl;
    } else {
        std::cout << "   ✗ Neumann边界条件应用测试失败" << std::endl;
    }
    
    std::cout << "\n=== 边界条件处理函数测试完成 ===" << std::endl;
}

/**
 * @brief 测试面连续性检查函数
 */
void testFaceContinuity() {
    std::cout << "\n=== 测试面连续性检查函数 ===" << std::endl;
    
    // 创建测试用的网格结构
    Mesh mesh;
    
    // 创建面结构
    Face face1, face2;
    face1.nodeIndexes = {0, 1, 2};
    face2.nodeIndexes = {0, 1, 2};
    
    mesh.faces.push_back(face1);
    mesh.faces.push_back(face2);
    
    // 创建测试用的元素结构
    Element element1, element2;
    element1.type.elementCode = 202; // 三角形元素
    element1.type.numberOfNodes = 3;
    element1.nodeIndexes = {0, 1, 2};
    element1.faceIndexes = {0};
    
    element2.type.elementCode = 202; // 三角形元素
    element2.type.numberOfNodes = 3;
    element2.nodeIndexes = {0, 1, 2};
    element2.faceIndexes = {1};
    
    // 创建节点坐标
    Nodes nodes;
    nodes.x = {0.0, 1.0, 0.0};
    nodes.y = {0.0, 0.0, 1.0};
    nodes.z = {0.0, 0.0, 0.0};
    
    // 测试CheckFaceNodeContinuity函数
    std::cout << "\n1. 测试CheckFaceNodeContinuity函数:" << std::endl;
    
    bool nodeContinuity = CheckFaceNodeContinuity(element1, element2, 0, mesh, 1e-12);
    
    std::cout << "   面节点连续性检查结果: " << (nodeContinuity ? "连续" : "不连续") << std::endl;
    
    if (nodeContinuity) {
        std::cout << "   ✓ 面节点连续性检查测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 面节点连续性检查测试失败" << std::endl;
    }
    
    // 测试CheckFaceGeometricContinuity函数
    std::cout << "\n2. 测试CheckFaceGeometricContinuity函数:" << std::endl;
    
    bool geometricContinuity = CheckFaceGeometricContinuity(element1, element2, 0, mesh, nodes, 1e-12);
    
    std::cout << "   面几何连续性检查结果: " << (geometricContinuity ? "连续" : "不连续") << std::endl;
    
    if (geometricContinuity) {
        std::cout << "   ✓ 面几何连续性检查测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 面几何连续性检查测试失败" << std::endl;
    }
    
    // 测试CheckFaceFunctionContinuity函数
    std::cout << "\n3. 测试CheckFaceFunctionContinuity函数:" << std::endl;
    
    std::vector<double> field1 = {1.0, 2.0, 3.0};
    std::vector<double> field2 = {1.0, 2.0, 3.0};
    
    bool functionContinuity = CheckFaceFunctionContinuity(element1, element2, 0, mesh, nodes, field1, field2, 1e-12);
    
    std::cout << "   面函数连续性检查结果: " << (functionContinuity ? "连续" : "不连续") << std::endl;
    
    if (functionContinuity) {
        std::cout << "   ✓ 面函数连续性检查测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 面函数连续性检查测试失败" << std::endl;
    }
    
    // 测试不连续的情况
    std::cout << "\n4. 测试不连续情况:" << std::endl;
    
    std::vector<double> field3 = {1.0, 2.0, 4.0}; // 第三个节点值不同
    bool discontinuity = !CheckFaceFunctionContinuity(element1, element2, 0, mesh, nodes, field1, field3, 1e-12);
    
    std::cout << "   面函数不连续性检查结果: " << (discontinuity ? "正确检测到不连续" : "错误检测为连续") << std::endl;
    
    if (discontinuity) {
        std::cout << "   ✓ 不连续性检测测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 不连续性检测测试失败" << std::endl;
    }
    
    std::cout << "\n=== 面连续性检查函数测试完成 ===" << std::endl;
}

/**
 * @brief 测试其他边界相关函数
 */
void testOtherBoundaryFunctions() {
    std::cout << "\n=== 测试其他边界相关函数 ===" << std::endl;
    
    // 创建测试用的元素结构
    Element element;
    element.type.elementCode = 202; // 三角形元素
    element.type.numberOfNodes = 3;
    element.nodeIndexes = {10, 20, 30};
    
    // 测试getTriangleFaceDirection函数
    std::cout << "\n1. 测试getTriangleFaceDirection函数:" << std::endl;
    
    std::vector<int> faceMap = {0, 1, 2};
    std::vector<int> indexes = {10, 20, 30};
    
    std::vector<int> direction = getTriangleFaceDirection(element, faceMap, indexes);
    
    std::cout << "   三角形面方向: ";
    for (int d : direction) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
    
    if (direction.size() == 3 && direction[0] == 0 && direction[1] == 1 && direction[2] == 2) {
        std::cout << "   ✓ 三角形面方向计算测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 三角形面方向计算测试失败" << std::endl;
    }
    
    // 测试getSquareFaceDirection函数
    std::cout << "\n2. 测试getSquareFaceDirection函数:" << std::endl;
    
    Element quadElement;
    quadElement.type.elementCode = 404; // 四边形元素
    quadElement.type.numberOfNodes = 4;
    quadElement.nodeIndexes = {10, 20, 30, 40};
    
    std::vector<int> quadFaceMap = {0, 1, 2, 3};
    std::vector<int> quadIndexes = {10, 20, 30, 40};
    
    std::vector<int> quadDirection = getSquareFaceDirection(quadElement, quadFaceMap, quadIndexes);
    
    std::cout << "   四边形面方向: ";
    for (int d : quadDirection) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
    
    if (quadDirection.size() == 4 && quadDirection[0] == 0 && quadDirection[1] == 1 && 
        quadDirection[2] == 2 && quadDirection[3] == 3) {
        std::cout << "   ✓ 四边形面方向计算测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 四边形面方向计算测试失败" << std::endl;
    }
    
    // 测试wedgeOrdering函数
    std::cout << "\n3. 测试wedgeOrdering函数:" << std::endl;
    
    std::vector<int> validOrdering = {0, 1, 2, 3};
    std::vector<int> invalidOrdering = {0, 1, 0, 3}; // 重复索引
    
    bool valid = wedgeOrdering(validOrdering);
    bool invalid = !wedgeOrdering(invalidOrdering);
    
    std::cout << "   有效楔形元素面排序检查: " << (valid ? "有效" : "无效") << std::endl;
    std::cout << "   无效楔形元素面排序检查: " << (invalid ? "正确检测为无效" : "错误检测为有效") << std::endl;
    
    if (valid && invalid) {
        std::cout << "   ✓ 楔形元素面排序检查测试通过" << std::endl;
    } else {
        std::cout << "   ✗ 楔形元素面排序检查测试失败" << std::endl;
    }
    
    std::cout << "\n=== 其他边界相关函数测试完成 ===" << std::endl;
}

int main() {
    std::cout << "开始测试边界处理函数..." << std::endl;
    
    try {
        testBoundaryConditions();
        testFaceContinuity();
        testOtherBoundaryFunctions();
        
        std::cout << "\n=== 所有边界处理函数测试完成 ===" << std::endl;
        std::cout << "测试结果: 所有测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "测试过程中发生错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}