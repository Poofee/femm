/**
 * @file test_boundary_simple.cpp
 * @brief 边界处理函数简单验证程序
 * 
 * 用于验证边界处理函数的正确性，不依赖复杂的构建系统。
 */

#include <iostream>
#include <vector>
#include <algorithm>

// 简化版本的Element结构
struct ElementType {
    int numberOfNodes;
    int dimension;
    int elementCode;
    
    ElementType() : numberOfNodes(0), dimension(0), elementCode(0) {}
    ElementType(int nodes, int dim, int code) : 
        numberOfNodes(nodes), dimension(dim), elementCode(code) {}
};

struct Element {
    ElementType type;
    
    Element() {}
    Element(const ElementType& t) : type(t) {}
};

/**
 * @brief 辅助排序函数
 */
void sort(std::vector<int>& arr) {
    std::sort(arr.begin(), arr.end());
}

/**
 * @brief 获取三角形面的全局方向
 */
std::vector<int> getTriangleFaceDirection(const Element& element, 
                                        const std::vector<int>& faceMap,
                                        const std::vector<int>& indexes) {
    // 检查输入参数
    if (faceMap.size() != 3) {
        throw std::invalid_argument("getTriangleFaceDirection: 面映射必须包含3个节点");
    }
    
    if (indexes.size() < static_cast<size_t>(element.type.numberOfNodes)) {
        throw std::invalid_argument("getTriangleFaceDirection: 索引数组大小不足");
    }
    
    // 将面的全局节点放入排序后的数组中
    std::vector<int> nodes(3);
    for (int i = 0; i < 3; ++i) {
        int nodeIndex = faceMap[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(indexes.size())) {
            throw std::invalid_argument("getTriangleFaceDirection: 无效的面映射节点索引");
        }
        nodes[i] = indexes[nodeIndex];
    }
    
    // 对节点进行排序
    sort(nodes);
    
    // 查找排序后节点的局部编号
    std::vector<int> globalDir(3, 0);
    for (int i = 0; i < element.type.numberOfNodes; ++i) {
        if (nodes[0] == indexes[i]) {
            globalDir[0] = i;
        } else if (nodes[1] == indexes[i]) {
            globalDir[1] = i;
        } else if (nodes[2] == indexes[i]) {
            globalDir[2] = i;
        }
    }
    
    return globalDir;
}

/**
 * @brief 检查楔形元素面的局部编号是否合法
 */
bool wedgeOrdering(const std::vector<int>& ordering) {
    if (ordering.size() != 4) {
        return false;
    }
    
    // 检查是否包含重复节点
    std::vector<int> sorted = ordering;
    sort(sorted);
    for (size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i] == sorted[i-1]) {
            return false;
        }
    }
    
    // 检查节点编号范围
    for (int node : ordering) {
        if (node < 0 || node >= 6) {
            return false;
        }
    }
    
    return true;
}

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
    if (result.size() == 3 && result[0] == 1 && result[1] == 2 && result[2] == 0) {
        std::cout << "✓ 三角形面方向函数测试: 通过" << std::endl;
    } else {
        std::cout << "✗ 三角形面方向函数测试: 失败" << std::endl;
    }
}

/**
 * @brief 测试楔形元素排序函数
 */
void testWedgeOrdering() {
    std::cout << "\n=== 楔形元素排序函数测试 ===" << std::endl;
    
    // 测试合法排序
    std::vector<int> validOrdering = {0, 1, 2, 3};
    bool validResult = wedgeOrdering(validOrdering);
    std::cout << "合法排序测试: " << (validResult ? "通过" : "失败") << std::endl;
    
    // 测试非法排序（重复节点）
    std::vector<int> invalidOrdering1 = {0, 1, 1, 3};
    bool invalidResult1 = wedgeOrdering(invalidOrdering1);
    std::cout << "重复节点测试: " << (!invalidResult1 ? "通过" : "失败") << std::endl;
    
    // 测试非法排序（超出范围）
    std::vector<int> invalidOrdering2 = {0, 1, 2, 10};
    bool invalidResult2 = wedgeOrdering(invalidOrdering2);
    std::cout << "超出范围测试: " << (!invalidResult2 ? "通过" : "失败") << std::endl;
    
    // 测试非法排序（大小错误）
    std::vector<int> invalidOrdering3 = {0, 1, 2};
    bool invalidResult3 = wedgeOrdering(invalidOrdering3);
    std::cout << "大小错误测试: " << (!invalidResult3 ? "通过" : "失败") << std::endl;
}

/**
 * @brief 测试错误处理
 */
void testErrorHandling() {
    std::cout << "\n=== 错误处理测试 ===" << std::endl;
    
    Element element;
    element.type.numberOfNodes = 4;
    
    try {
        // 测试面映射大小错误
        std::vector<int> invalidFaceMap = {0, 1}; // 只有2个节点
        std::vector<int> indexes = {1, 2, 3, 4};
        getTriangleFaceDirection(element, invalidFaceMap, indexes);
        std::cout << "✗ 面映射大小错误测试: 失败（未抛出异常）" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "✓ 面映射大小错误测试: 通过" << std::endl;
    }
    
    try {
        // 测试索引数组大小不足
        std::vector<int> validFaceMap = {0, 1, 2};
        std::vector<int> smallIndexes = {1, 2}; // 只有2个索引
        getTriangleFaceDirection(element, validFaceMap, smallIndexes);
        std::cout << "✗ 索引数组大小不足测试: 失败（未抛出异常）" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "✓ 索引数组大小不足测试: 通过" << std::endl;
    }
}

int main() {
    std::cout << "边界处理函数验证程序" << std::endl;
    std::cout << "========================" << std::endl;
    
    testTriangleFaceDirection();
    testWedgeOrdering();
    testErrorHandling();
    
    std::cout << "\n========================" << std::endl;
    std::cout << "验证程序执行完成" << std::endl;
    
    return 0;
}