/**
 * @file test_boundary_functions.cpp
 * @brief 边界处理函数测试程序
 * 
 * 测试边界处理函数的正确性，包括：
 * - getTriangleFaceDirection
 * - getSquareFaceDirection  
 * - wedgeOrdering
 */

#include "src/ElementDescription.h"
#include <iostream>
#include <vector>
#include <cassert>

using namespace ElmerCpp;

/**
 * @brief 测试三角形面方向函数
 */
void testTriangleFaceDirection() {
    std::cout << "=== 三角形面方向函数测试 ===" << std::endl;
    
    // 创建测试元素
    Element element;
    element.type.numberOfNodes = 4;
    element.type.dimension = 2;
    element.type.elementCode = 303; // 三角形P型元素
    
    // 创建测试数据
    std::vector<int> faceMap = {0, 1, 2}; // 三角形面的节点映射
    std::vector<int> indexes = {10, 5, 8, 3}; // 全局节点索引
    
    // 测试函数
    std::vector<int> result = getTriangleFaceDirection(element, faceMap, indexes);
    
    // 验证结果
    std::cout << "三角形面方向结果: [" << result[0] << ", " << result[1] << ", " << result[2] << "]" << std::endl;
    
    // 预期结果：排序后的全局节点为[5,8,10]，对应的局部编号为[1,2,0]
    assert(result.size() == 3);
    assert(result[0] == 1); // 全局节点5对应局部节点1
    assert(result[1] == 2); // 全局节点8对应局部节点2
    assert(result[2] == 0); // 全局节点10对应局部节点0
    
    std::cout << "三角形面方向函数测试: 通过" << std::endl;
}

/**
 * @brief 测试四边形面方向函数
 */
void testSquareFaceDirection() {
    std::cout << "=== 四边形面方向函数测试 ===" << std::endl;
    
    // 创建测试元素
    Element element;
    element.type.numberOfNodes = 8;
    element.type.dimension = 3;
    element.type.elementCode = 404; // 四边形P型元素
    
    // 创建测试数据
    std::vector<int> faceMap = {0, 1, 2, 3}; // 四边形面的节点映射
    std::vector<int> indexes = {15, 7, 12, 3, 20, 8, 18, 5}; // 全局节点索引
    
    // 测试函数
    std::vector<int> result = getSquareFaceDirection(element, faceMap, indexes);
    
    // 验证结果
    std::cout << "四边形面方向结果: [" << result[0] << ", " << result[1] << ", " 
              << result[2] << ", " << result[3] << "]" << std::endl;
    
    // 预期结果：最小全局节点为3（局部节点3），然后按顺序找到其他节点
    assert(result.size() == 4);
    assert(result[0] == 3); // 最小全局节点3对应局部节点3
    
    std::cout << "四边形面方向函数测试: 通过" << std::endl;
}

/**
 * @brief 测试楔形元素排序检查函数
 */
void testWedgeOrdering() {
    std::cout << "=== 楔形元素排序检查函数测试 ===" << std::endl;
    
    // 测试合法的排序（前两个节点在同一个三角形面上）
    std::vector<int> validOrdering1 = {0, 1, 2, 3}; // 节点0和1在三角形面1-3
    std::vector<int> validOrdering2 = {3, 4, 5, 0}; // 节点3和4在三角形面4-6
    
    // 测试非法的排序（前两个节点在不同三角形面上）
    std::vector<int> invalidOrdering1 = {0, 3, 1, 2}; // 节点0和3在不同三角形面
    std::vector<int> invalidOrdering2 = {2, 4, 0, 1}; // 节点2和4在不同三角形面
    
    // 测试合法情况
    bool result1 = wedgeOrdering(validOrdering1);
    bool result2 = wedgeOrdering(validOrdering2);
    
    // 测试非法情况
    bool result3 = wedgeOrdering(invalidOrdering1);
    bool result4 = wedgeOrdering(invalidOrdering2);
    
    std::cout << "合法排序1: " << (result1 ? "通过" : "失败") << std::endl;
    std::cout << "合法排序2: " << (result2 ? "通过" : "失败") << std::endl;
    std::cout << "非法排序1: " << (result3 ? "失败" : "通过") << std::endl;
    std::cout << "非法排序2: " << (result4 ? "失败" : "通过") << std::endl;
    
    // 验证结果
    assert(result1 == true);
    assert(result2 == true);
    assert(result3 == false);
    assert(result4 == false);
    
    std::cout << "楔形元素排序检查函数测试: 通过" << std::endl;
}

/**
 * @brief 测试边界处理函数的错误处理
 */
void testBoundaryFunctionErrorHandling() {
    std::cout << "=== 边界处理函数错误处理测试 ===" << std::endl;
    
    Element element;
    element.type.numberOfNodes = 4;
    
    // 测试无效的面映射大小
    try {
        std::vector<int> invalidFaceMap = {0, 1}; // 只有2个节点，需要3个
        std::vector<int> indexes = {1, 2, 3, 4};
        getTriangleFaceDirection(element, invalidFaceMap, indexes);
        std::cout << "错误处理测试: 失败（应该抛出异常）" << std::endl;
        assert(false); // 应该抛出异常
    } catch (const std::invalid_argument& e) {
        std::cout << "无效面映射大小测试: 通过（正确抛出异常）" << std::endl;
    }
    
    // 测试无效的节点索引
    try {
        std::vector<int> faceMap = {0, 1, 2};
        std::vector<int> invalidIndexes = {1, 2}; // 大小不足
        getTriangleFaceDirection(element, faceMap, invalidIndexes);
        std::cout << "错误处理测试: 失败（应该抛出异常）" << std::endl;
        assert(false); // 应该抛出异常
    } catch (const std::invalid_argument& e) {
        std::cout << "无效节点索引测试: 通过（正确抛出异常）" << std::endl;
    }
    
    std::cout << "边界处理函数错误处理测试: 通过" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始边界处理函数测试..." << std::endl;
    
    try {
        testTriangleFaceDirection();
        testSquareFaceDirection();
        testWedgeOrdering();
        testBoundaryFunctionErrorHandling();
        
        std::cout << "\n=== 所有边界处理函数测试通过 ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "测试失败: " << e.what() << std::endl;
        return 1;
    }
}