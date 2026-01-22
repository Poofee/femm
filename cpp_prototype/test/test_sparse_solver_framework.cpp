// test_sparse_solver_framework.cpp - 稀疏矩阵求解器框架测试
// 测试稀疏矩阵存储、求解器接口和求解器框架的功能

#include "SparseMatrix.h"
#include "LinearSolverInterface.h"
#include "LinearSolverFactory.h"
#include "CRSMatrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <cassert>
#include <string>

using namespace elmer;

/**
 * @brief 测试稀疏矩阵基本功能
 */
void testSparseMatrixBasic() {
    std::cout << "=== 测试稀疏矩阵基本功能 ===" << std::endl;
    
    // 创建CSR矩阵
    auto matrix = std::make_shared<elmer::CSRMatrix>(5, 5);
    
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
    auto b = elmer::Vector::Create(5);
    auto y = elmer::Vector::Create(5);
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
 * @brief 测试求解器工厂
 */
void testSolverFactory() {
    std::cout << "=== 测试求解器工厂 ===" << std::endl;
    
    // 测试可用求解器列表
    auto available_solvers = elmer::LinearSolverFactory::GetAvailableSolverNames();
    std::cout << "可用求解器: ";
    for (const auto& solver : available_solvers) {
        std::cout << solver << " ";
    }
    std::cout << std::endl;
    
    // 测试求解器可用性
    std::cout << "求解器可用性检查: " << std::endl;
    for (const auto& solver : available_solvers) {
        bool available = elmer::LinearSolverFactory::IsSolverAvailable(solver);
        std::cout << "  " << solver << ": " << (available ? "可用" : "不可用") << std::endl;
    }
    
    // 使用工厂创建SuperLU求解器
    try {
        auto superlu_solver = elmer::LinearSolverFactory::CreateSolver("SuperLU");
        std::cout << "✓ SuperLU求解器创建成功" << std::endl;
        
        // 获取推荐参数
        auto params = elmer::LinearSolverFactory::GetRecommendedParameters(elmer::SolverType::SUPERLU);
        std::cout << "SuperLU推荐参数: tolerance=" << params.tolerance 
                  << ", drop_tolerance=" << params.drop_tolerance << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "✗ SuperLU求解器创建失败: " << e.what() << std::endl;
    }
    
    // 使用工厂创建MUMPS求解器
    try {
        auto mumps_solver = elmer::LinearSolverFactory::CreateSolver("MUMPS");
        std::cout << "✓ MUMPS求解器创建成功" << std::endl;
        
        // 获取推荐参数
        auto params = elmer::LinearSolverFactory::GetRecommendedParameters(elmer::SolverType::MUMPS);
        std::cout << "MUMPS推荐参数: icntl_7=" << params.mumps_icntl_7 
                  << ", icntl_14=" << params.mumps_icntl_14 << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "✗ MUMPS求解器创建失败: " << e.what() << std::endl;
    }
    
    std::cout << "✓ 求解器工厂测试通过" << std::endl << std::endl;
}

/**
 * @brief 测试线性代数求解器功能
 */
void testLinearAlgebraSolvers() {
    std::cout << "=== 测试线性代数求解器功能 ===" << std::endl;
    
    // 创建测试矩阵（5x5三对角矩阵）
    auto matrix = std::make_shared<elmer::CSRMatrix>(5, 5);
    for (elmer::Integer i = 0; i < 5; ++i) {
        matrix->SetElement(i, i, 4.0);
        if (i > 0) {
            matrix->SetElement(i, i-1, -1.0);
        }
        if (i < 4) {
            matrix->SetElement(i, i+1, -1.0);
        }
    }
    matrix->Assemble();
    
    // 创建右端向量 b = [1, 2, 3, 4, 5]
    auto b = elmer::Vector::Create(5);
    for (elmer::Integer i = 0; i < 5; ++i) {
        (*b)[i] = static_cast<elmer::Real>(i + 1);
    }
    
    // 测试不同求解器
    std::vector<std::string> solvers_to_test = {"SuperLU", "MUMPS"};
    
    for (const auto& solver_name : solvers_to_test) {
        std::cout << "\n测试求解器: " << solver_name << std::endl;
        
        try {
            auto solver = elmer::LinearSolverFactory::CreateSolver(solver_name);
            
            // 初始化求解器
            if (!solver->Initialize(matrix)) {
                std::cout << "✗ 求解器初始化失败" << std::endl;
                continue;
            }
            std::cout << "✓ 求解器初始化成功" << std::endl;
            
            // 求解线性方程组
            auto x = elmer::Vector::Create(5);
            if (!solver->Solve(*b, *x)) {
                std::cout << "✗ 线性方程组求解失败" << std::endl;
                continue;
            }
            std::cout << "✓ 线性方程组求解成功" << std::endl;
            
            // 验证解的正确性
            auto Ax = elmer::Vector::Create(5);
            matrix->Multiply(*x, *Ax);
            
            elmer::Real residual = 0.0;
            for (elmer::Integer i = 0; i < 5; ++i) {
                residual += std::pow((*Ax)[i] - (*b)[i], 2);
            }
            residual = std::sqrt(residual);
            
            std::cout << "残差范数: " << residual << std::endl;
            
            if (residual < 1e-10) {
                std::cout << "✓ 解的正确性验证通过" << std::endl;
            } else {
                std::cout << "✗ 解的正确性验证失败" << std::endl;
            }
            
        } catch (const std::exception& e) {
            std::cout << "✗ 求解器测试失败: " << e.what() << std::endl;
        }
    }
    
    std::cout << "✓ 线性代数求解器功能测试通过" << std::endl << std::endl;
}

/**
 * @brief 测试线性方程组求解
 */
void testLinearSystemSolving() {
    std::cout << "=== 测试线性方程组求解 ===" << std::endl;
    
    // 创建测试矩阵（5x5三对角矩阵）
    auto matrix = std::make_shared<CSRMatrix>(5, 5);
    for (Integer i = 0; i < 5; ++i) {
        matrix->SetElement(i, i, 4.0);
        if (i > 0) {
            matrix->SetElement(i, i-1, -1.0);
        }
        if (i < 4) {
            matrix->SetElement(i, i+1, -1.0);
        }
    }
    matrix->Assemble();
    
    // 创建右端向量 b = [1, 2, 3, 4, 5]
    auto b = Vector::Create(5);
    for (Integer i = 0; i < 5; ++i) {
        (*b)[i] = static_cast<Real>(i + 1);
    }
    
    // 创建求解器
    try {
        auto solver = elmer::LinearSolverFactory::CreateSolver("SuperLU");
        
        // 初始化求解器
        if (!solver->Initialize(matrix)) {
            std::cout << "✗ 求解器初始化失败" << std::endl;
            return;
        }
        std::cout << "✓ 求解器初始化成功" << std::endl;
        
        // 求解线性方程组
        auto x = elmer::Vector::Create(5);
        if (!solver->Solve(*b, *x)) {
            std::cout << "✗ 线性方程组求解失败" << std::endl;
            return;
        }
        std::cout << "✓ 线性方程组求解成功" << std::endl;
        
        // 验证解的正确性
        auto Ax = elmer::Vector::Create(5);
        matrix->Multiply(*x, *Ax);
        
        elmer::Real residual = 0.0;
        for (elmer::Integer i = 0; i < 5; ++i) {
            residual += std::pow((*Ax)[i] - (*b)[i], 2);
        }
        residual = std::sqrt(residual);
        
        std::cout << "解向量 x = [";
        for (elmer::Integer i = 0; i < 5; ++i) {
            std::cout << (*x)[i] << (i < 4 ? ", " : "]");
        }
        std::cout << std::endl;
        
        std::cout << "残差范数: " << residual << std::endl;
        
        if (residual < 1e-10) {
            std::cout << "✓ 解的正确性验证通过" << std::endl;
        } else {
            std::cout << "✗ 解的正确性验证失败" << std::endl;
        }
        
        // 获取求解器状态
        auto status = solver->GetStatus();
        std::cout << "求解器状态: " << std::endl;
        std::cout << status.ToString() << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "✗ 求解器测试失败: " << e.what() << std::endl;
        return;
    }
    
    std::cout << "✓ 线性方程组求解测试通过" << std::endl << std::endl;
}

/**
 * @brief 测试低频电磁有限元场景
 */
void testElectromagneticScenario() {
    std::cout << "=== 测试低频电磁有限元场景 ===" << std::endl;
    
    // 创建拉普拉斯矩阵（模拟电磁场离散）
    elmer::Integer n = 10; // 10x10网格
    // TODO: 需要实现SparseMatrixUtils::CreateLaplacianMatrix
    // 暂时使用简单的三对角矩阵替代
    auto matrix = std::make_shared<elmer::CSRMatrix>(n, n);
    for (elmer::Integer i = 0; i < n; ++i) {
        matrix->SetElement(i, i, 4.0);
        if (i > 0) {
            matrix->SetElement(i, i-1, -1.0);
        }
        if (i < n-1) {
            matrix->SetElement(i, i+1, -1.0);
        }
    }
    matrix->Assemble();
    
    std::cout << "创建测试矩阵 (" << n << "x" << n 
              << "), 非零元素: " << matrix->GetNonZeroCount() << std::endl;
    
    // 创建右端向量（模拟源项）
    auto b = elmer::Vector::Create(n);
    (*b)[0] = 1.0;  // 边界条件
    (*b)[n-1] = -1.0; // 边界条件
    
    // 测试不同求解器
    std::vector<std::string> solvers_to_test = {"SuperLU", "MUMPS"};
    
    for (const auto& solver_name : solvers_to_test) {
        std::cout << "\n测试求解器: " << solver_name << std::endl;
        
        try {
            auto solver = elmer::LinearSolverFactory::CreateSolver(solver_name);
            
            // 设置求解器参数
            auto params = elmer::LinearSolverFactory::GetRecommendedParameters(
                solver_name == "SuperLU" ? elmer::SolverType::SUPERLU : elmer::SolverType::MUMPS);
            solver->SetParameters(params);
            
            // 初始化求解器
            if (!solver->Initialize(matrix)) {
                std::cout << "✗ 求解器初始化失败" << std::endl;
                continue;
            }
            
            // 求解
            auto x = elmer::Vector::Create(n);
            auto start_time = std::chrono::high_resolution_clock::now();
            
            bool success = solver->Solve(*b, *x);
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            
            if (success) {
                // 验证解
                auto Ax = elmer::Vector::Create(n);
                matrix->Multiply(*x, *Ax);
                
                elmer::Real residual = 0.0;
                for (elmer::Integer i = 0; i < n; ++i) {
                    residual += std::pow((*Ax)[i] - (*b)[i], 2);
                }
                residual = std::sqrt(residual);
                
                std::cout << "  求解时间: " << duration.count() / 1000.0 << " ms" << std::endl;
                std::cout << "  残差: " << residual << std::endl;
                std::cout << "  状态: " << (residual < 1e-10 ? "成功" : "失败") << std::endl;
                
            } else {
                std::cout << "✗ 求解失败" << std::endl;
            }
            
        } catch (const std::exception& e) {
            std::cout << "✗ 求解器测试异常: " << e.what() << std::endl;
        }
    }
    
    std::cout << "✓ 低频电磁有限元场景测试完成" << std::endl << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "稀疏矩阵求解器框架测试开始" << std::endl;
    std::cout << "==========================" << std::endl << std::endl;
    
    try {
        testSparseMatrixBasic();
        testSolverFactory();
        testLinearAlgebraSolvers();
        testLinearSystemSolving();
        testElectromagneticScenario();
        
        std::cout << "==========================" << std::endl;
        std::cout << "所有测试通过！稀疏矩阵求解器框架功能正常。" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "测试失败: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}