// test_crs_matrix.cpp - CRSMatrix模块测试
// 验证CRSMatrix类的核心功能

#include "../src/Types.h"
#include "../src/CRSMatrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <memory>

// 测试辅助函数
bool isClose(double a, double b, double tolerance = 1e-12) {
    return std::abs(a - b) < tolerance;
}

void printTestResult(const std::string& testName, bool passed) {
    std::cout << (passed ? "PASSED" : "FAILED") << ": " << testName << std::endl;
}

// 测试矩阵创建和基本属性
void testMatrixCreation() {
    std::cout << "=== 测试矩阵创建和基本属性 ===" << std::endl;
    
    try {
        auto matrix = elmer::CreateCRSMatrix(5, 5);
        
        if (matrix->GetNumRows() != 5 || matrix->GetNumCols() != 5) {
            printTestResult("矩阵维度检查", false);
            return;
        }
        
        printTestResult("矩阵维度检查", true);
        
        // 测试零矩阵
        matrix->Zero();
        
        // 测试获取元素（应该返回0）
        bool allZero = true;
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (!isClose(matrix->GetElement(i, j), 0.0)) {
                    allZero = false;
                    break;
                }
            }
        }
        
        printTestResult("零矩阵初始化", allZero);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 矩阵创建异常 - " << e.what() << std::endl;
    }
}

// 测试设置和获取矩阵元素
void testSetGetElements() {
    std::cout << "=== 测试设置和获取矩阵元素 ===" << std::endl;
    
    try {
        auto matrix = elmer::CreateCRSMatrix(3, 3);
        
        // 设置对角线元素
        matrix->SetElement(0, 0, 1.0);
        matrix->SetElement(1, 1, 2.0);
        matrix->SetElement(2, 2, 3.0);
        
        // 设置非对角线元素
        matrix->SetElement(0, 1, 0.5);
        matrix->SetElement(1, 0, 0.5);
        
        // 验证设置的元素
        bool elementsCorrect = true;
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(0, 0), 1.0);
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(1, 1), 2.0);
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(2, 2), 3.0);
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(0, 1), 0.5);
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(1, 0), 0.5);
        
        // 验证未设置的元素为0
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(0, 2), 0.0);
        elementsCorrect = elementsCorrect && isClose(matrix->GetElement(2, 0), 0.0);
        
        printTestResult("设置和获取矩阵元素", elementsCorrect);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 元素设置异常 - " << e.what() << std::endl;
    }
}

// 测试矩阵加法操作
void testAddToElement() {
    std::cout << "=== 测试矩阵加法操作 ===" << std::endl;
    
    try {
        auto matrix = elmer::CreateCRSMatrix(2, 2);
        
        // 初始设置
        matrix->SetElement(0, 0, 1.0);
        matrix->SetElement(0, 1, 0.5);
        
        // 加法操作
        matrix->AddToElement(0, 0, 0.5);
        matrix->AddToElement(0, 1, 0.5);
        matrix->AddToElement(1, 1, 1.0);  // 新元素
        
        // 验证结果
        bool addCorrect = true;
        addCorrect = addCorrect && isClose(matrix->GetElement(0, 0), 1.5);
        addCorrect = addCorrect && isClose(matrix->GetElement(0, 1), 1.0);
        addCorrect = addCorrect && isClose(matrix->GetElement(1, 1), 1.0);
        addCorrect = addCorrect && isClose(matrix->GetElement(1, 0), 0.0);
        
        printTestResult("矩阵元素加法操作", addCorrect);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 加法操作异常 - " << e.what() << std::endl;
    }
}

// 测试矩阵-向量乘法
void testMatrixVectorMultiplication() {
    std::cout << "=== 测试矩阵-向量乘法 ===" << std::endl;
    
    try {
        auto matrix = elmer::CreateCRSMatrix(3, 3);
        auto x = elmer::Vector::Create(3);
        auto y = elmer::Vector::Create(3);
        
        // 创建对称矩阵
        matrix->SetElement(0, 0, 2.0);
        matrix->SetElement(1, 1, 3.0);
        matrix->SetElement(2, 2, 4.0);
        matrix->SetElement(0, 1, 1.0);
        matrix->SetElement(1, 0, 1.0);
        
        // 设置向量
        (*x)[0] = 1.0;
        (*x)[1] = 2.0;
        (*x)[2] = 3.0;
        
        // 执行乘法
        matrix->Multiply(*x, *y);
        
        // 验证结果: y = A * x
        // y[0] = 2*1 + 1*2 + 0*3 = 4
        // y[1] = 1*1 + 3*2 + 0*3 = 7  
        // y[2] = 0*1 + 0*2 + 4*3 = 12
        bool multCorrect = true;
        multCorrect = multCorrect && isClose((*y)[0], 4.0);
        multCorrect = multCorrect && isClose((*y)[1], 7.0);
        multCorrect = multCorrect && isClose((*y)[2], 12.0);
        
        printTestResult("矩阵-向量乘法", multCorrect);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 矩阵乘法异常 - " << e.what() << std::endl;
    }
}

// 测试转置矩阵-向量乘法
void testTransposeMultiplication() {
    std::cout << "=== 测试转置矩阵-向量乘法 ===" << std::endl;
    
    try {
        auto matrix = elmer::CreateCRSMatrix(2, 3);  // 非方阵
        auto x = elmer::Vector::Create(2);
        auto y = elmer::Vector::Create(3);
        
        // 创建矩阵
        matrix->SetElement(0, 0, 1.0);
        matrix->SetElement(0, 1, 2.0);
        matrix->SetElement(1, 2, 3.0);
        
        // 设置向量
        (*x)[0] = 2.0;
        (*x)[1] = 3.0;
        
        // 执行转置乘法: y = A^T * x
         auto crsMatrix = dynamic_cast<elmer::CRSMatrix*>(matrix.get());
         if (crsMatrix) {
             crsMatrix->MultiplyTranspose(*x, *y);
         } else {
             std::cout << "FAILED: 无法转换为CRSMatrix类型" << std::endl;
             return;
         }
        
        // 验证结果
        // y[0] = 1*2 + 0*3 = 2
        // y[1] = 2*2 + 0*3 = 4
        // y[2] = 0*2 + 3*3 = 9
        bool transCorrect = true;
        transCorrect = transCorrect && isClose((*y)[0], 2.0);
        transCorrect = transCorrect && isClose((*y)[1], 4.0);
        transCorrect = transCorrect && isClose((*y)[2], 9.0);
        
        printTestResult("转置矩阵-向量乘法", transCorrect);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 转置乘法异常 - " << e.what() << std::endl;
    }
}

// 测试对称性检查
void testSymmetryCheck() {
    std::cout << "=== 测试对称性检查 ===" << std::endl;
    
    try {
        // 创建对称矩阵
        auto symMatrix = elmer::CreateCRSMatrix(2, 2);
        symMatrix->SetElement(0, 0, 1.0);
        symMatrix->SetElement(1, 1, 2.0);
        symMatrix->SetElement(0, 1, 0.5);
        symMatrix->SetElement(1, 0, 0.5);
        
        // 创建非对称矩阵
        auto asymMatrix = elmer::CreateCRSMatrix(2, 2);
        asymMatrix->SetElement(0, 0, 1.0);
        asymMatrix->SetElement(1, 1, 2.0);
        asymMatrix->SetElement(0, 1, 0.5);
        asymMatrix->SetElement(1, 0, 1.5);  // 不对称
        
        // 检查对称性
        auto symCRS = dynamic_cast<elmer::CRSMatrix*>(symMatrix.get());
        auto asymCRS = dynamic_cast<elmer::CRSMatrix*>(asymMatrix.get());
        
        if (!symCRS || !asymCRS) {
            std::cout << "FAILED: 无法转换为CRSMatrix类型" << std::endl;
            return;
        }
        
        if (!symCRS->IsSymmetric(1e-12)) {
            std::cout << "FAILED: 对称矩阵检查失败" << std::endl;
            return;
        }
        
        if (asymCRS->IsSymmetric(1e-12)) {
            std::cout << "FAILED: 非对称矩阵检查失败" << std::endl;
            return;
        }
        
        printTestResult("矩阵对称性检查", true);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 对称性检查异常 - " << e.what() << std::endl;
    }
}

// 测试错误处理
void testErrorHandling() {
    std::cout << "=== 测试错误处理 ===" << std::endl;
    
    try {
        auto matrix = elmer::CreateCRSMatrix(2, 2);
        
        // 测试越界访问
        bool caughtOutOfRange = false;
        try {
            matrix->GetElement(5, 0);  // 行越界
        } catch (const std::out_of_range&) {
            caughtOutOfRange = true;
        }
        
        // 测试维度不匹配
        bool caughtInvalidArgument = false;
        try {
            auto x = elmer::Vector::Create(3);  // 错误维度
            auto y = elmer::Vector::Create(2);
            matrix->Multiply(*x, *y);
        } catch (const std::invalid_argument&) {
            caughtInvalidArgument = true;
        }
        
        bool errorHandling = caughtOutOfRange && caughtInvalidArgument;
        printTestResult("错误处理", errorHandling);
        
    } catch (const std::exception& e) {
        std::cout << "FAILED: 错误处理测试异常 - " << e.what() << std::endl;
    }
}

// 主测试函数
int main() {
    std::cout << "开始CRSMatrix模块测试..." << std::endl;
    std::cout << "==================================" << std::endl;
    
    testMatrixCreation();
    testSetGetElements();
    testAddToElement();
    testMatrixVectorMultiplication();
    testTransposeMultiplication();
    testSymmetryCheck();
    testErrorHandling();
    
    std::cout << "==================================" << std::endl;
    std::cout << "所有测试完成！" << std::endl;
    
    return 0;
}