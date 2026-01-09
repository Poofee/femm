/**
 * @file test_element_utils.cpp
 * @brief ElementUtils模块测试套件
 * 
 * 测试ElementUtils模块的各种功能，包括矩阵创建、初始化、释放等操作。
 */

#include <iostream>
#include <vector>
#include <memory>
#include "../src/ElementUtils.h"

using namespace elmer;

// 工具函数：检查两个数值是否接近
bool isClose(double a, double b, double tolerance = 1e-12) {
    return std::abs(a - b) < tolerance;
}

// 测试矩阵格式验证
void testMatrixFormatValidation() {
    std::cout << "=== 测试矩阵格式验证 ===" << std::endl;
    
    // 测试有效格式
    if (!ElementUtils::validateMatrixFormat(ElementUtils::MatrixFormat::CRS)) {
        std::cout << "FAILED: CRS格式验证失败" << std::endl;
        return;
    }
    
    if (!ElementUtils::validateMatrixFormat(ElementUtils::MatrixFormat::BAND)) {
        std::cout << "FAILED: BAND格式验证失败" << std::endl;
        return;
    }
    
    if (!ElementUtils::validateMatrixFormat(ElementUtils::MatrixFormat::SBAND)) {
        std::cout << "FAILED: SBAND格式验证失败" << std::endl;
        return;
    }
    
    std::cout << "PASSED: 矩阵格式验证测试" << std::endl;
}

// 测试矩阵创建功能
void testMatrixCreation() {
    std::cout << "=== 测试矩阵创建功能 ===" << std::endl;
    
    try {
        // 测试CRS矩阵创建
        auto crsMatrix = ElementUtils::createMatrix(
            nullptr, nullptr, nullptr,
            std::vector<int>{}, // 空排列
            1,                  // 1个自由度
            ElementUtils::MatrixFormat::CRS,
            false               // 不重新排序
        );
        
        if (!crsMatrix) {
            std::cout << "FAILED: CRS矩阵创建失败" << std::endl;
            return;
        }
        
        std::cout << "PASSED: CRS矩阵创建测试" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 矩阵创建异常: " << e.what() << std::endl;
        return;
    }
}

// 测试错误处理
void testErrorHandling() {
    std::cout << "=== 测试错误处理 ===" << std::endl;
    
    // 测试无效的自由度数量
    try {
        auto matrix = ElementUtils::createMatrix(
            nullptr, nullptr, nullptr,
            std::vector<int>{},
            0,  // 无效的自由度数量
            ElementUtils::MatrixFormat::CRS,
            false
        );
        
        std::cout << "FAILED: 应该对无效自由度数量抛出异常" << std::endl;
        return;
        
    } catch (const std::exception& e) {
        std::cout << "PASSED: 正确捕获无效自由度数量异常" << std::endl;
    }
    
    // 测试空矩阵指针检查
    try {
        ElementUtils::freeMatrix(nullptr);
        std::cout << "PASSED: 空矩阵指针释放处理正常" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 空矩阵指针释放异常: " << e.what() << std::endl;
        return;
    }
}

// 测试DG辐射索引计算
void testDgRadiationIndexes() {
    std::cout << "=== 测试DG辐射索引计算 ===" << std::endl;
    
    try {
        std::vector<int> elemInds;
        
        // 测试正常情况
        ElementUtils::dgRadiationIndexes(nullptr, 5, elemInds, false);
        
        if (elemInds.size() != 5) {
            std::cout << "FAILED: 索引数组大小不正确" << std::endl;
            return;
        }
        
        // 验证索引值
        for (int i = 0; i < 5; ++i) {
            if (elemInds[i] != i) {
                std::cout << "FAILED: 索引值不正确" << std::endl;
                return;
            }
        }
        
        std::cout << "PASSED: DG辐射索引计算测试" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: DG辐射索引计算异常: " << e.what() << std::endl;
        return;
    }
}

// 测试打包排序添加操作
void testPackSortAdd() {
    std::cout << "=== 测试打包排序添加操作 ===" << std::endl;
    
    try {
        std::vector<int> ind, perm;
        
        // 测试正常情况
        ElementUtils::packSortAdd(3, ind, perm);
        
        if (ind.size() != 3 || perm.size() != 3) {
            std::cout << "FAILED: 数组大小不正确" << std::endl;
            return;
        }
        
        std::cout << "PASSED: 打包排序添加操作测试" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 打包排序添加操作异常: " << e.what() << std::endl;
        return;
    }
}

// 测试列表矩阵创建
void testListMatrixCreation() {
    std::cout << "=== 测试列表矩阵创建 ===" << std::endl;
    
    try {
        // 测试列表矩阵创建
        auto listMatrix = ElementUtils::makeListMatrix(
            nullptr, nullptr, nullptr, nullptr,
            false,  // 不重新排序
            2,      // 2个自由度
            ElementUtils::MatrixFormat::CRS
        );
        
        if (!listMatrix) {
            std::cout << "FAILED: 列表矩阵创建失败" << std::endl;
            return;
        }
        
        std::cout << "PASSED: 列表矩阵创建测试" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 列表矩阵创建异常: " << e.what() << std::endl;
        return;
    }
}

// 测试列表矩阵数组创建
void testListMatrixArrayCreation() {
    std::cout << "=== 测试列表矩阵数组创建 ===" << std::endl;
    
    try {
        // 测试列表矩阵数组创建
        auto matrixArray = ElementUtils::makeListMatrixArray(
            nullptr, nullptr, nullptr, nullptr,
            false,  // 不重新排序
            1,      // 1个自由度
            ElementUtils::MatrixFormat::CRS
        );
        
        if (matrixArray.empty()) {
            std::cout << "FAILED: 列表矩阵数组创建失败" << std::endl;
            return;
        }
        
        std::cout << "PASSED: 列表矩阵数组创建测试" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 列表矩阵数组创建异常: " << e.what() << std::endl;
        return;
    }
}

int main() {
    std::cout << "开始ElementUtils模块测试..." << std::endl;
    std::cout << "==================================" << std::endl;
    
    testMatrixFormatValidation();
    testMatrixCreation();
    testErrorHandling();
    testDgRadiationIndexes();
    testPackSortAdd();
    testListMatrixCreation();
    testListMatrixArrayCreation();
    
    std::cout << "==================================" << std::endl;
    std::cout << "所有测试完成！" << std::endl;
    
    return 0;
}