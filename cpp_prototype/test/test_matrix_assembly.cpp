// test_matrix_assembly.cpp - 矩阵组装模块测试程序

#include "../src/MatrixAssembly.h"
#include "../src/ElmerCpp.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

using namespace elmer;

// 测试矩阵基本操作
void testBasicMatrixOperations() {
    std::cout << "测试矩阵基本操作..." << std::endl;
    
    // 创建一个10x10的CRS矩阵
    auto matrix = Matrix::CreateCRS(10, 10);
    
    // 测试设置和获取矩阵元素
    MatrixAssembly::SetMatrixElement(matrix, 0, 0, 1.0);
    MatrixAssembly::SetMatrixElement(matrix, 1, 1, 2.0);
    MatrixAssembly::SetMatrixElement(matrix, 2, 2, 3.0);
    
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 0) == 1.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 1, 1) == 2.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 2, 2) == 3.0);
    
    // 测试修改矩阵元素
    Real oldValue = MatrixAssembly::ChangeMatrixElement(matrix, 0, 0, 5.0);
    assert(oldValue == 1.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 0) == 5.0);
    
    // 测试添加到矩阵元素
    MatrixAssembly::AddToMatrixElement(matrix, 0, 0, 1.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 0) == 6.0);
    
    std::cout << "基本矩阵操作测试通过!" << std::endl;
}

// 测试矩阵向量乘法
void testMatrixVectorMultiplication() {
    std::cout << "测试矩阵向量乘法..." << std::endl;
    
    // 创建一个3x3的稠密矩阵
    auto matrix = Matrix::CreateDense(3, 3);
    
    // 设置矩阵元素（单位矩阵）
    MatrixAssembly::SetMatrixElement(matrix, 0, 0, 1.0);
    MatrixAssembly::SetMatrixElement(matrix, 1, 1, 1.0);
    MatrixAssembly::SetMatrixElement(matrix, 2, 2, 1.0);
    
    // 创建测试向量
    std::vector<Real> x = {1.0, 2.0, 3.0};
    std::vector<Real> y(3, 0.0);
    
    // 执行矩阵向量乘法
    MatrixAssembly::MatrixVectorMultiply(matrix, x, y);
    
    // 验证结果
    assert(y[0] == 1.0);
    assert(y[1] == 2.0);
    assert(y[2] == 3.0);
    
    // 测试非单位矩阵
    MatrixAssembly::SetMatrixElement(matrix, 0, 1, 1.0);
    MatrixAssembly::SetMatrixElement(matrix, 1, 2, 1.0);
    
    y = {0.0, 0.0, 0.0};
    MatrixAssembly::MatrixVectorMultiply(matrix, x, y);
    
    assert(y[0] == 1.0 + 2.0);  // 1*1 + 1*2
    assert(y[1] == 2.0 + 3.0);  // 1*2 + 1*3
    assert(y[2] == 3.0);        // 1*3
    
    std::cout << "矩阵向量乘法测试通过!" << std::endl;
}

// 测试单元矩阵组装
void testElementMatrixAssembly() {
    std::cout << "测试单元矩阵组装..." << std::endl;
    
    // 创建一个6x6的全局矩阵
    auto globalMatrix = Matrix::CreateCRS(6, 6);
    
    // 创建一个2x2的单元刚度矩阵
    std::vector<std::vector<Real>> elementMatrix = {
        {2.0, -1.0},
        {-1.0, 2.0}
    };
    
    // 自由度索引（组装到全局矩阵的位置1和3）
    std::vector<Integer> dofIndices = {1, 3};
    
    // 组装单元矩阵
    MatrixAssembly::AssembleElementMatrix(globalMatrix, dofIndices, elementMatrix);
    
    // 验证组装结果
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 1, 1) == 2.0);
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 1, 3) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 3, 1) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 3, 3) == 2.0);
    
    // 测试多次组装（叠加）
    MatrixAssembly::AssembleElementMatrix(globalMatrix, dofIndices, elementMatrix);
    
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 1, 1) == 4.0);
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 1, 3) == -2.0);
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 3, 1) == -2.0);
    assert(MatrixAssembly::GetMatrixElement(globalMatrix, 3, 3) == 4.0);
    
    std::cout << "单元矩阵组装测试通过!" << std::endl;
}

// 测试单元向量组装
void testElementVectorAssembly() {
    std::cout << "测试单元向量组装..." << std::endl;
    
    // 创建一个6维的全局向量
    std::vector<Real> globalVector(6, 0.0);
    
    // 创建一个2维的单元载荷向量
    std::vector<Real> elementVector = {1.5, 2.5};
    
    // 自由度索引
    std::vector<Integer> dofIndices = {2, 4};
    
    // 组装单元向量
    MatrixAssembly::AssembleElementVector(globalVector, dofIndices, elementVector);
    
    // 验证组装结果
    assert(globalVector[2] == 1.5);
    assert(globalVector[4] == 2.5);
    
    // 测试多次组装（叠加）
    MatrixAssembly::AssembleElementVector(globalVector, dofIndices, elementVector);
    
    assert(globalVector[2] == 3.0);
    assert(globalVector[4] == 5.0);
    
    std::cout << "单元向量组装测试通过!" << std::endl;
}

// 测试Dirichlet边界条件应用
void testDirichletBoundaryConditions() {
    std::cout << "测试Dirichlet边界条件..." << std::endl;
    
    // 创建一个3x3的矩阵
    auto matrix = Matrix::CreateDense(3, 3);
    
    // 设置矩阵元素（对称正定矩阵）
    MatrixAssembly::SetMatrixElement(matrix, 0, 0, 4.0);
    MatrixAssembly::SetMatrixElement(matrix, 0, 1, -1.0);
    MatrixAssembly::SetMatrixElement(matrix, 0, 2, -1.0);
    MatrixAssembly::SetMatrixElement(matrix, 1, 0, -1.0);
    MatrixAssembly::SetMatrixElement(matrix, 1, 1, 4.0);
    MatrixAssembly::SetMatrixElement(matrix, 1, 2, -1.0);
    MatrixAssembly::SetMatrixElement(matrix, 2, 0, -1.0);
    MatrixAssembly::SetMatrixElement(matrix, 2, 1, -1.0);
    MatrixAssembly::SetMatrixElement(matrix, 2, 2, 4.0);
    
    // 创建右端项
    std::vector<Real> rhs = {1.0, 2.0, 3.0};
    
    // 应用Dirichlet边界条件（第0个自由度固定为5.0）
    MatrixAssembly::ApplyDirichletBC(matrix, rhs, 0, 5.0);
    
    // 验证边界条件应用
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 0) == 1.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 1) == 0.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 2) == 0.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 1, 0) == 0.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 2, 0) == 0.0);
    assert(rhs[0] == 5.0);
    
    // 验证其他行保持不变
    assert(MatrixAssembly::GetMatrixElement(matrix, 1, 1) == 4.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 1, 2) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 2, 1) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(matrix, 2, 2) == 4.0);
    
    std::cout << "Dirichlet边界条件测试通过!" << std::endl;
}

// 测试矩阵属性检查
void testMatrixProperties() {
    std::cout << "测试矩阵属性检查..." << std::endl;
    
    // 创建一个对称矩阵
    auto symmetricMatrix = Matrix::CreateDense(3, 3);
    MatrixAssembly::SetMatrixElement(symmetricMatrix, 0, 0, 4.0);
    MatrixAssembly::SetMatrixElement(symmetricMatrix, 0, 1, 1.0);
    MatrixAssembly::SetMatrixElement(symmetricMatrix, 1, 0, 1.0);
    MatrixAssembly::SetMatrixElement(symmetricMatrix, 1, 1, 4.0);
    MatrixAssembly::SetMatrixElement(symmetricMatrix, 2, 2, 4.0);
    
    // 测试对称性检查
    assert(MatrixAssembly::IsSymmetric(symmetricMatrix));
    
    // 创建一个非对称矩阵
    auto nonSymmetricMatrix = Matrix::CreateDense(2, 2);
    MatrixAssembly::SetMatrixElement(nonSymmetricMatrix, 0, 1, 1.0);
    MatrixAssembly::SetMatrixElement(nonSymmetricMatrix, 1, 0, 2.0);
    
    assert(!MatrixAssembly::IsSymmetric(nonSymmetricMatrix));
    
    // 测试正定性检查（简化实现）
    auto positiveDefiniteMatrix = Matrix::CreateDense(2, 2);
    MatrixAssembly::SetMatrixElement(positiveDefiniteMatrix, 0, 0, 4.0);
    MatrixAssembly::SetMatrixElement(positiveDefiniteMatrix, 1, 1, 4.0);
    
    assert(MatrixAssembly::IsPositiveDefinite(positiveDefiniteMatrix));
    
    std::cout << "矩阵属性检查测试通过!" << std::endl;
}

// 测试复数矩阵操作
void testComplexMatrixOperations() {
    std::cout << "测试复数矩阵操作..." << std::endl;
    
    // 创建一个4x4的矩阵（用于存储2x2复数矩阵）
    auto matrix = Matrix::CreateDense(4, 4);
    
    // 添加复数矩阵元素（实部=2.0，虚部=1.0）
    MatrixAssembly::AddToComplexMatrixElement(matrix, 0, 0, 2.0, 1.0);
    
    // 验证复数矩阵存储格式
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 0) == 2.0);  // 实部
    assert(MatrixAssembly::GetMatrixElement(matrix, 0, 1) == -1.0); // -虚部
    assert(MatrixAssembly::GetMatrixElement(matrix, 1, 0) == 1.0);  // 虚部
    assert(MatrixAssembly::GetMatrixElement(matrix, 1, 1) == 2.0);  // 实部
    
    std::cout << "复数矩阵操作测试通过!" << std::endl;
}

// 主测试函数
int main() {
    std::cout << "开始矩阵组装模块测试..." << std::endl;
    std::cout << "=======================================" << std::endl;
    
    try {
        testBasicMatrixOperations();
        testMatrixVectorMultiplication();
        testElementMatrixAssembly();
        testElementVectorAssembly();
        testDirichletBoundaryConditions();
        testMatrixProperties();
        testComplexMatrixOperations();
        
        std::cout << "=======================================" << std::endl;
        std::cout << "所有矩阵组装模块测试通过!" << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "测试失败: " << e.what() << std::endl;
        return 1;
    }
}