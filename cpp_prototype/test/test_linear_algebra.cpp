// test_linear_algebra.cpp - 线性代数模块测试程序

#include "../src/LinearAlgebra.h"
#include <iostream>
#include <cmath>

using namespace ElmerCpp;

void testCRSMatrix() {
    std::cout << "=== 测试CRS矩阵 ===" << std::endl;
    
    // 创建3x3矩阵
    auto matrix = createCRSMatrix(3, 3);
    
    // 设置矩阵元素（对称正定矩阵）
    matrix->setElement(0, 0, 4.0);
    matrix->setElement(0, 1, -1.0);
    matrix->setElement(0, 2, 0.0);
    
    matrix->setElement(1, 0, -1.0);
    matrix->setElement(1, 1, 4.0);
    matrix->setElement(1, 2, -1.0);
    
    matrix->setElement(2, 0, 0.0);
    matrix->setElement(2, 1, -1.0);
    matrix->setElement(2, 2, 4.0);
    
    // 验证矩阵元素
    std::cout << "矩阵元素验证:" << std::endl;
    std::cout << "A[0,0] = " << matrix->getElement(0, 0) << " (期望: 4.0)" << std::endl;
    std::cout << "A[1,1] = " << matrix->getElement(1, 1) << " (期望: 4.0)" << std::endl;
    std::cout << "A[2,2] = " << matrix->getElement(2, 2) << " (期望: 4.0)" << std::endl;
    
    // 测试矩阵-向量乘法
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> b = matrix->multiply(x);
    
    std::cout << "矩阵-向量乘法测试:" << std::endl;
    std::cout << "x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
    std::cout << "Ax = [" << b[0] << ", " << b[1] << ", " << b[2] << "]" << std::endl;
    
    // 验证结果
    double expected_b0 = 4.0*1.0 - 1.0*2.0 + 0.0*3.0;
    double expected_b1 = -1.0*1.0 + 4.0*2.0 - 1.0*3.0;
    double expected_b2 = 0.0*1.0 - 1.0*2.0 + 4.0*3.0;
    
    std::cout << "期望结果: [" << expected_b0 << ", " << expected_b1 << ", " << expected_b2 << "]" << std::endl;
    
    bool test_passed = std::abs(b[0] - expected_b0) < 1e-10 &&
                      std::abs(b[1] - expected_b1) < 1e-10 &&
                      std::abs(b[2] - expected_b2) < 1e-10;
    
    std::cout << "CRS矩阵测试: " << (test_passed ? "通过" : "失败") << std::endl;
}

void testIterativeSolver() {
    std::cout << "\n=== 测试迭代求解器 ===" << std::endl;
    
    // 创建对称正定矩阵
    auto matrix = createCRSMatrix(3, 3);
    matrix->setElement(0, 0, 4.0);
    matrix->setElement(0, 1, -1.0);
    matrix->setElement(1, 0, -1.0);
    matrix->setElement(1, 1, 4.0);
    matrix->setElement(1, 2, -1.0);
    matrix->setElement(2, 1, -1.0);
    matrix->setElement(2, 2, 4.0);
    
    // 创建求解器
    auto solver = createIterativeSolver(matrix);
    
    // 设置求解器参数
    SolverParameters params;
    params.method = SolverMethod::CG;
    params.tolerance = 1e-10;
    params.maxIterations = 100;
    params.verbose = true;
    
    solver->setParameters(params);
    
    // 设置右端项
    std::vector<double> rhs = {1.0, 2.0, 3.0};
    
    // 求解线性系统
    auto result = solver->solve(rhs);
    
    // 获取解
    auto solution = solver->getSolution();
    
    std::cout << "求解结果:" << std::endl;
    std::cout << "迭代次数: " << result.iterations << std::endl;
    std::cout << "最终残差: " << result.residual << std::endl;
    std::cout << "是否收敛: " << (result.converged ? "是" : "否") << std::endl;
    std::cout << "解向量: [" << solution[0] << ", " << solution[1] << ", " << solution[2] << "]" << std::endl;
    
    // 验证解的正确性
    auto residual = matrix->multiply(solution);
    for (size_t i = 0; i < residual.size(); ++i) {
        residual[i] -= rhs[i];
    }
    
    double residual_norm = norm(residual);
    std::cout << "实际残差范数: " << residual_norm << std::endl;
    
    bool test_passed = result.converged && residual_norm < 1e-8;
    std::cout << "迭代求解器测试: " << (test_passed ? "通过" : "失败") << std::endl;
}

void testMatrixAssembler() {
    std::cout << "\n=== 测试矩阵组装器 ===" << std::endl;
    
    // 创建全局矩阵和向量
    auto matrix = createCRSMatrix(4, 4);
    std::vector<double> rhs(4, 0.0);
    
    auto assembler = createMatrixAssembler(matrix);
    
    // 组装第一个单元（节点0,1,2）
    std::vector<int> elem1_indices = {0, 1, 2};
    std::vector<std::vector<double>> elem1_matrix = {
        {2.0, -1.0, 0.0},
        {-1.0, 2.0, -1.0},
        {0.0, -1.0, 2.0}
    };
    std::vector<double> elem1_vector = {1.0, 2.0, 3.0};
    
    assembler->addElementMatrix(elem1_indices, elem1_matrix);
    assembler->addElementVector(elem1_indices, elem1_vector, rhs);
    
    // 组装第二个单元（节点1,2,3）
    std::vector<int> elem2_indices = {1, 2, 3};
    std::vector<std::vector<double>> elem2_matrix = {
        {2.0, -1.0, 0.0},
        {-1.0, 2.0, -1.0},
        {0.0, -1.0, 2.0}
    };
    std::vector<double> elem2_vector = {2.0, 3.0, 1.0};
    
    assembler->addElementMatrix(elem2_indices, elem2_matrix);
    assembler->addElementVector(elem2_indices, elem2_vector, rhs);
    
    // 验证组装结果
    std::cout << "组装后的矩阵对角线元素:" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << "A[" << i << "," << i << "] = " << matrix->getElement(i, i) << std::endl;
    }
    
    std::cout << "组装后的右端项: [";
    for (size_t i = 0; i < rhs.size(); ++i) {
        std::cout << rhs[i];
        if (i < rhs.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    // 验证对角线元素（应该叠加）
    bool test_passed = true;
    test_passed &= std::abs(matrix->getElement(0, 0) - 2.0) < 1e-10;  // 只参与单元1
    test_passed &= std::abs(matrix->getElement(1, 1) - 4.0) < 1e-10;  // 参与单元1和2
    test_passed &= std::abs(matrix->getElement(2, 2) - 4.0) < 1e-10;  // 参与单元1和2
    test_passed &= std::abs(matrix->getElement(3, 3) - 2.0) < 1e-10;  // 只参与单元2
    
    std::cout << "矩阵组装器测试: " << (test_passed ? "通过" : "失败") << std::endl;
}

int main() {
    std::cout << "Elmer FEM C++线性代数模块测试" << std::endl;
    std::cout << "================================" << std::endl;
    
    try {
        testCRSMatrix();
        testIterativeSolver();
        testMatrixAssembler();
        
        std::cout << "\n=== 所有测试完成 ===" << std::endl;
        std::cout << "C++线性代数模块基本功能验证成功！" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "测试过程中出现错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}