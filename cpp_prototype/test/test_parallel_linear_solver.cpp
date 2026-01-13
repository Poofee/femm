// test_parallel_linear_solver.cpp - 并行线性求解器测试

#include "../src/ParallelLinearSolver.h"
#include "../src/DistributedLinearAlgebra.h"
#include "../src/MPIConfig.h"
#include "../src/CRSMatrix.h"
#include <iostream>
#include <cmath>
#include <memory>

using namespace elmer;

/**
 * @brief 创建简单的测试矩阵
 */
std::shared_ptr<DistributedMatrix> createTestMatrix(std::shared_ptr<MPICommunicator> comm) {
    int globalSize = 100;
    auto matrix = std::make_shared<DistributedMatrix>(globalSize, globalSize, comm);
    
    // 创建对角占优矩阵
    for (int i = 0; i < globalSize; ++i) {
        matrix->SetElement(i, i, 2.0); // 对角线元素
        if (i > 0) {
            matrix->SetElement(i, i-1, -1.0); // 下对角线
        }
        if (i < globalSize - 1) {
            matrix->SetElement(i, i+1, -1.0); // 上对角线
        }
    }
    
    return matrix;
}

/**
 * @brief 创建测试右端向量
 */
std::shared_ptr<DistributedVector> createTestRhs(std::shared_ptr<MPICommunicator> comm) {
    int globalSize = 100;
    auto rhs = std::make_shared<DistributedVector>(globalSize, comm);
    
    // 创建简单的右端向量：b_i = i
    for (int i = 0; i < globalSize; ++i) {
        rhs->SetElement(i, static_cast<double>(i));
    }
    
    return rhs;
}

/**
 * @brief 测试并行共轭梯度法
 */
void testParallelCG() {
    auto comm = MPIUtils::getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    
    if (rank == 0) {
        std::cout << "=== 测试并行共轭梯度法 ===" << std::endl;
        std::cout << "进程数: " << size << std::endl;
    }
    
    // 创建测试系统
    auto A = createTestMatrix(comm);
    auto b = createTestRhs(comm);
    auto x = std::make_shared<DistributedVector>(A->getGlobalRows(), comm);
    
    // 创建线性系统
    auto linearSystem = std::make_shared<DistributedLinearSystem>(A, b, x, comm);
    
    // 创建并行CG求解器
    auto solver = std::make_shared<ParallelConjugateGradientSolver>(comm);
    solver->setLinearSystem(linearSystem);
    solver->setParameters(1e-8, 1000, 1); // 容差1e-8，最大迭代1000，输出详细
    
    // 求解
    bool success = solver->solve(x);
    
    if (rank == 0) {
        std::cout << "求解结果: " << (success ? "成功" : "失败") << std::endl;
        std::cout << "迭代次数: " << solver->getIterations() << std::endl;
        std::cout << "最终残差: " << solver->getResidual() << std::endl;
        std::cout << "收敛状态: " << (solver->isConverged() ? "收敛" : "未收敛") << std::endl;
    }
    
    // 验证解的正确性
    auto residual = std::make_shared<DistributedVector>(A->getGlobalRows(), comm);
    A->multiply(*x, *residual);
    residual->axpy(-1.0, *b);
    
    double residual_norm = residual->norm();
    double global_residual = 0.0;
    comm->allReduce(&residual_norm, &global_residual, 1, MPI_DOUBLE, MPI_SUM);
    global_residual = std::sqrt(global_residual);
    
    if (rank == 0) {
        std::cout << "验证残差: " << global_residual << std::endl;
        std::cout << "测试结果: " << (global_residual < 1e-6 ? "通过" : "失败") << std::endl;
    }
}

/**
 * @brief 测试并行GMRES求解器
 */
void testParallelGMRES() {
    auto comm = MPIUtils::getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    
    if (rank == 0) {
        std::cout << "\n=== 测试并行GMRES求解器 ===" << std::endl;
        std::cout << "进程数: " << size << std::endl;
    }
    
    // 创建测试系统
    auto A = createTestMatrix(comm);
    auto b = createTestRhs(comm);
    auto x = std::make_shared<DistributedVector>(A->getGlobalRows(), comm);
    
    // 创建线性系统
    auto linearSystem = std::make_shared<DistributedLinearSystem>(A, b, x, comm);
    
    // 创建并行GMRES求解器
    auto solver = std::make_shared<ParallelGMRESSolver>(comm);
    solver->setLinearSystem(linearSystem);
    solver->setParameters(1e-8, 1000, 1);
    
    // 求解
    bool success = solver->solve(x);
    
    if (rank == 0) {
        std::cout << "求解结果: " << (success ? "成功" : "失败") << std::endl;
        std::cout << "迭代次数: " << solver->getIterations() << std::endl;
        std::cout << "最终残差: " << solver->getResidual() << std::endl;
        std::cout << "收敛状态: " << (solver->isConverged() ? "收敛" : "未收敛") << std::endl;
    }
    
    // 验证解的正确性
    auto residual = std::make_shared<DistributedVector>(A->getGlobalRows(), comm);
    A->multiply(*x, *residual);
    residual->axpy(-1.0, *b);
    
    double residual_norm = residual->norm();
    double global_residual = 0.0;
    comm->allReduce(&residual_norm, &global_residual, 1, MPI_DOUBLE, MPI_SUM);
    global_residual = std::sqrt(global_residual);
    
    if (rank == 0) {
        std::cout << "验证残差: " << global_residual << std::endl;
        std::cout << "测试结果: " << (global_residual < 1e-6 ? "通过" : "失败") << std::endl;
    }
}

/**
 * @brief 测试求解器工厂
 */
void testSolverFactory() {
    auto comm = MPIUtils::getDefaultComm();
    int rank = comm->getRank();
    
    if (rank == 0) {
        std::cout << "\n=== 测试求解器工厂 ===" << std::endl;
    }
    
    // 创建测试系统
    auto A = createTestMatrix(comm);
    auto b = createTestRhs(comm);
    auto x = std::make_shared<DistributedVector>(A->getGlobalRows(), comm);
    
    // 创建线性系统
    auto linearSystem = std::make_shared<DistributedLinearSystem>(A, b, x, comm);
    
    // 测试CG求解器
    auto cgSolver = ParallelSolverFactory::createSolver(
        ParallelSolverFactory::CG, comm);
    cgSolver->setLinearSystem(linearSystem);
    
    bool cgSuccess = cgSolver->solve(x);
    
    if (rank == 0) {
        std::cout << "CG求解器测试: " << (cgSuccess ? "通过" : "失败") << std::endl;
    }
    
    // 测试GMRES求解器
    x->setZero(); // 重置解向量
    auto gmresSolver = ParallelSolverFactory::createSolver(
        ParallelSolverFactory::GMRES, comm);
    gmresSolver->setLinearSystem(linearSystem);
    
    bool gmresSuccess = gmresSolver->solve(x);
    
    if (rank == 0) {
        std::cout << "GMRES求解器测试: " << (gmresSuccess ? "通过" : "失败") << std::endl;
    }
}

int main(int argc, char** argv) {
    // 初始化MPI
    MPI_Init(&argc, &argv);
    
    try {
        // 运行测试
        testParallelCG();
        testParallelGMRES();
        testSolverFactory();
        
        auto comm = MPIUtils::getDefaultComm();
        int rank = comm->getRank();
        if (rank == 0) {
            std::cout << "\n=== 所有测试完成 ===" << std::endl;
        }
    }
    catch (const std::exception& e) {
        auto comm = MPIUtils::getDefaultComm();
        int rank = comm->getRank();
        if (rank == 0) {
            std::cerr << "测试失败: " << e.what() << std::endl;
        }
    }
    
    // 清理MPI
    MPI_Finalize();
    
    return 0;
}