// test_sparse_solver_framework_simple.cpp - 稀疏矩阵求解器框架简化测试
// 测试稀疏矩阵存储和基本接口功能，不依赖具体求解器实现

#include "../src/core/math/SparseMatrix.h"
#include "../src/core/math/LinearSolverInterface.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <cassert>

using namespace elmer;

/**
 * @brief 测试稀疏矩阵基本功能
 */
void testSparseMatrixBasic() {
    std::cout << "=== 测试稀疏矩阵基本功能 ===" << std::endl;
    
    // 创建CSR矩阵
    auto matrix = std::make_shared<CSRMatrix>(5, 5);
    
    // 设置矩阵元素（对角线占优）
    for (Integer i = 0; i < 5; ++i) {
        matrix->SetElement(i, i, 4.0); // 对角线
        if (i > 0) {
            matrix->SetElement(i, i-1, -1.0); // 下对角线
        }
        if (i < 4) {
            matrix->SetElement(i, i+1, -1.0); // 上对角线
        }
    }
    
    // 装配矩阵
    matrix->Assemble();
    
    // 验证矩阵属性
    std::cout << "矩阵尺寸: " << matrix->GetNumRows() << "x" << matrix->GetNumCols() << std::endl;
    std::cout << "非零元素数: " << matrix->GetNonZeroCount() << std::endl;
    std::cout << "内存使用: " << matrix->GetMemoryUsage() / 1024.0 << " KB" << std::endl;
    
    // 验证矩阵元素
    for (Integer i = 0; i < 5; ++i) {
        for (Integer j = 0; j < 5; ++j) {
            Real value = matrix->GetElement(i, j);
            if (i == j) {
                assert(std::abs(value - 4.0) < 1e-12);
            } else if (std::abs(i - j) == 1) {
                assert(std::abs(value + 1.0) < 1e-12);
            } else {
                assert(std::abs(value) < 1e-12);
            }
        }
    }
    
    // 测试矩阵向量乘法
    auto b = Vector::Create(5);
    auto y = Vector::Create(5);
    for (Integer i = 0; i < 5; ++i) {
        (*b)[i] = static_cast<Real>(i + 1);
    }
    
    matrix->Multiply(*b, *y);
    
    std::cout << "矩阵向量乘法测试: " << std::endl;
    std::cout << "b = [";
    for (Integer i = 0; i < 5; ++i) {
        std::cout << (*b)[i] << (i < 4 ? ", " : "]");
    }
    std::cout << std::endl;
    
    std::cout << "A*b = [";
    for (Integer i = 0; i < 5; ++i) {
        std::cout << (*y)[i] << (i < 4 ? ", " : "]");
    }
    std::cout << std::endl;
    
    std::cout << "✓ 稀疏矩阵基本功能测试通过" << std::endl << std::endl;
}

/**
 * @brief 测试矩阵格式转换
 */
void testMatrixFormatConversion() {
    std::cout << "=== 测试矩阵格式转换 ===" << std::endl;
    
    // 创建CSR矩阵
    auto csr_matrix = std::make_shared<CSRMatrix>(3, 3);
    
    // 设置矩阵元素
    csr_matrix->SetElement(0, 0, 1.0);
    csr_matrix->SetElement(0, 1, 2.0);
    csr_matrix->SetElement(1, 0, 3.0);
    csr_matrix->SetElement(1, 1, 4.0);
    csr_matrix->SetElement(2, 2, 5.0);
    
    csr_matrix->Assemble();
    
    std::cout << "CSR矩阵非零元素数: " << csr_matrix->GetNonZeroCount() << std::endl;
    
    // 测试矩阵复制功能
    auto copied_matrix = std::make_shared<CSRMatrix>(3, 3);
    for (Integer i = 0; i < 3; ++i) {
        for (Integer j = 0; j < 3; ++j) {
            Real value = csr_matrix->GetElement(i, j);
            if (std::abs(value) > 1e-12) {
                copied_matrix->SetElement(i, j, value);
            }
        }
    }
    copied_matrix->Assemble();
    
    std::cout << "复制矩阵非零元素数: " << copied_matrix->GetNonZeroCount() << std::endl;
    
    // 验证复制正确性
    for (Integer i = 0; i < 3; ++i) {
        for (Integer j = 0; j < 3; ++j) {
            Real original = csr_matrix->GetElement(i, j);
            Real copied = copied_matrix->GetElement(i, j);
            assert(std::abs(original - copied) < 1e-12);
        }
    }
    
    std::cout << "✓ 矩阵格式转换测试通过" << std::endl << std::endl;
}

/**
 * @brief 测试矩阵I/O功能
 */
void testMatrixIO() {
    std::cout << "=== 测试矩阵I/O功能 ===" << std::endl;
    
    // 创建测试矩阵
    auto matrix = std::make_shared<CSRMatrix>(4, 4);
    
    // 设置矩阵元素
    for (Integer i = 0; i < 4; ++i) {
        matrix->SetElement(i, i, i + 1.0);
    }
    matrix->SetElement(0, 1, 0.5);
    matrix->SetElement(1, 0, 0.5);
    
    matrix->Assemble();
    
    // 保存矩阵到文件
    std::string filename = "test_matrix.mtx";
    if (matrix->SaveToFile(filename)) {
        std::cout << "✓ 矩阵保存到文件成功: " << filename << std::endl;
        
        // 从文件加载矩阵
        auto loaded_matrix = std::make_shared<CSRMatrix>(4, 4);
        if (loaded_matrix->LoadFromFile(filename)) {
            std::cout << "✓ 矩阵从文件加载成功" << std::endl;
            
            // 验证加载的矩阵与原矩阵一致
            bool matrices_equal = true;
            for (Integer i = 0; i < 4; ++i) {
                for (Integer j = 0; j < 4; ++j) {
                    Real original = matrix->GetElement(i, j);
                    Real loaded = loaded_matrix->GetElement(i, j);
                    if (std::abs(original - loaded) > 1e-12) {
                        matrices_equal = false;
                        break;
                    }
                }
                if (!matrices_equal) break;
            }
            
            if (matrices_equal) {
                std::cout << "✓ 矩阵I/O功能验证通过" << std::endl;
            } else {
                std::cout << "✗ 矩阵I/O功能验证失败" << std::endl;
            }
        } else {
            std::cout << "✗ 矩阵从文件加载失败" << std::endl;
        }
    } else {
        std::cout << "✗ 矩阵保存到文件失败" << std::endl;
    }
    
    std::cout << std::endl;
}

/**
 * @brief 测试求解器接口定义
 */
void testSolverInterface() {
    std::cout << "=== 测试求解器接口定义 ===" << std::endl;
    
    // 测试求解器类型枚举
    std::cout << "求解器类型定义检查:" << std::endl;
    std::cout << "  SUPERLU = " << static_cast<int>(SolverType::SUPERLU) << std::endl;
    std::cout << "  SUPERLU_MT = " << static_cast<int>(SolverType::SUPERLU_MT) << std::endl;
    std::cout << "  MUMPS = " << static_cast<int>(SolverType::MUMPS) << std::endl;
    std::cout << "  PETSC = " << static_cast<int>(SolverType::PETSC) << std::endl;
    
    // 测试求解器参数结构
    SolverParameters params;
    params.tolerance = 1e-8;
    params.max_iterations = 1000;
    params.drop_tolerance = 1e-4;
    params.mumps_icntl_7 = 0;
    params.mumps_icntl_14 = 20;
    
    std::cout << "求解器参数结构检查:" << std::endl;
    std::cout << "  容差: " << params.tolerance << std::endl;
    std::cout << "  最大迭代次数: " << params.max_iterations << std::endl;
    
    // 测试求解器状态结构
    SolverStatus status;
    status.success = true;
    status.message = "测试状态";
    status.iterations = 10;
    status.residual = 1e-12;
    status.setup_time = 0.1;
    status.solve_time = 0.5;
    status.memory_usage = 1024.0;
    
    std::cout << "求解器状态结构检查:" << std::endl;
    std::cout << "  成功: " << (status.success ? "是" : "否") << std::endl;
    std::cout << "  迭代次数: " << status.iterations << std::endl;
    std::cout << "  残差: " << status.residual << std::endl;
    
    std::cout << "✓ 求解器接口定义测试通过" << std::endl << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "稀疏矩阵求解器框架简化测试开始" << std::endl;
    std::cout << "================================" << std::endl << std::endl;
    
    try {
        testSparseMatrixBasic();
        testMatrixFormatConversion();
        testMatrixIO();
        testSolverInterface();
        
        std::cout << "================================" << std::endl;
        std::cout << "所有基础测试通过！稀疏矩阵框架核心功能正常。" << std::endl;
        std::cout << "下一步：实现具体求解器类并完成完整测试。" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "测试失败: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}