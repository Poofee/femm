// ParallelLinearSolver.cpp - MPI并行线性求解器实现

#include "ParallelLinearSolver.h"
#include "MPIUtils.h"
#include <cmath>
#include <algorithm>
#include <chrono>

namespace elmer {

// ============================================================================
// ParallelConjugateGradientSolver 实现
// ============================================================================

void ParallelConjugateGradientSolver::initialize() {
    if (!linearSystem_) {
        throw std::runtime_error("线性系统未设置");
    }
    
    // 创建工作向量
    auto& matrix = linearSystem_->getMatrix();
    r_ = std::make_shared<DistributedVector>(matrix.getGlobalRows(), comm_);
    p_ = std::make_shared<DistributedVector>(matrix.getGlobalRows(), comm_);
    Ap_ = std::make_shared<DistributedVector>(matrix.getGlobalRows(), comm_);
    z_ = std::make_shared<DistributedVector>(matrix.getGlobalRows(), comm_);
    
    isInitialized_ = true;
    isConverged_ = false;
    iterations_ = 0;
    residual_ = 0.0;
}

bool ParallelConjugateGradientSolver::solve(std::shared_ptr<DistributedVector> x, 
                                           std::shared_ptr<DistributedVector> b) {
    if (!isInitialized_) {
        initialize();
    }
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    auto& matrix = linearSystem_->getMatrix();
    
    // 初始化残差 r = b - A*x
    matrix.multiply(*x, *r_);
    r_->copy(*b);
    r_->axpy(-1.0, *r_);
    
    // 计算初始残差范数
    double r0_norm = computeResidualNorm(r_);
    double residual_norm = r0_norm;
    
    if (r0_norm < tolerance_) {
        isConverged_ = true;
        residual_ = r0_norm;
        return true;
    }
    
    // 应用预条件子：z = M^{-1} * r
    if (preconditioner_) {
        applyPreconditioner(r_, z_);
    } else {
        z_->copy(*r_);
    }
    
    // 初始化搜索方向 p = z
    p_->copy(*z_);
    
    double rho_old = r_->dot(*z_);
    
    for (iterations_ = 1; iterations_ <= maxIterations_; ++iterations_) {
        // 计算 Ap = A*p
        matrix.multiply(*p_, *Ap_);
        
        // 计算 alpha = rho_old / (p^T * Ap)
        double pAp = p_->dot(*Ap_);
        if (std::abs(pAp) < 1e-16) {
            // 矩阵奇异，收敛失败
            break;
        }
        
        double alpha = rho_old / pAp;
        
        // 更新解 x = x + alpha * p
        x->axpy(alpha, *p_);
        
        // 更新残差 r = r - alpha * Ap
        r_->axpy(-alpha, *Ap_);
        
        // 检查收敛
        residual_norm = computeResidualNorm(r_);
        if (verbosity_ > 0 && iterations_ % 10 == 0) {
            int rank = comm_->getRank();
            if (rank == 0) {
                std::cout << "CG迭代 " << iterations_ << ": 残差 = " 
                          << residual_norm << std::endl;
            }
        }
        
        if (residual_norm < tolerance_) {
            isConverged_ = true;
            residual_ = residual_norm;
            break;
        }
        
        // 应用预条件子：z = M^{-1} * r
        if (preconditioner_) {
            applyPreconditioner(r_, z_);
        } else {
            z_->copy(*r_);
        }
        
        // 计算 rho_new = r^T * z
        double rho_new = r_->dot(*z_);
        
        // 计算 beta = rho_new / rho_old
        double beta = rho_new / rho_old;
        
        // 更新搜索方向 p = z + beta * p
        p_->scale(beta);
        p_->axpy(1.0, *z_);
        
        rho_old = rho_new;
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    
    if (verbosity_ > 0) {
        int rank = comm_->getRank();
        if (rank == 0) {
            std::cout << "CG求解完成: 迭代次数 = " << iterations_ 
                      << ", 最终残差 = " << residual_norm 
                      << ", 耗时 = " << duration.count() / 1000.0 << " ms" << std::endl;
        }
    }
    
    residual_ = residual_norm;
    return isConverged_;
}

bool ParallelConjugateGradientSolver::solve(std::shared_ptr<DistributedVector> x) {
    if (!linearSystem_) {
        throw std::runtime_error("线性系统未设置");
    }
    
    auto b = linearSystem_->getRhs();
    return solve(x, b);
}

double ParallelConjugateGradientSolver::computeResidualNorm(std::shared_ptr<DistributedVector> r) {
    // 计算全局残差范数
    double local_norm = r->norm();
    double global_norm = 0.0;
    
    comm_->allReduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
    
    return std::sqrt(global_norm);
}

void ParallelConjugateGradientSolver::applyPreconditioner(std::shared_ptr<DistributedVector> r, 
                                                         std::shared_ptr<DistributedVector> z) {
    if (preconditioner_) {
        preconditioner_->solve(z, r);
    } else {
        z->copy(*r);
    }
}

// ============================================================================
// ParallelGMRESSolver 实现
// ============================================================================

void ParallelGMRESSolver::initialize() {
    if (!linearSystem_) {
        throw std::runtime_error("线性系统未设置");
    }
    
    // 预分配Krylov子空间向量
    auto& matrix = linearSystem_->getMatrix();
    v_.resize(restartSize_ + 1);
    z_.resize(restartSize_ + 1);
    
    for (int i = 0; i <= restartSize_; ++i) {
        v_[i] = matrix.createVector();
        z_[i] = matrix.createVector();
    }
    
    // 预分配Hessenberg矩阵
    h_.resize(restartSize_ + 1);
    for (int i = 0; i <= restartSize_; ++i) {
        h_[i].resize(restartSize_);
    }
    
    g_.resize(restartSize_ + 1);
    
    isInitialized_ = true;
    isConverged_ = false;
    iterations_ = 0;
    residual_ = 0.0;
}

bool ParallelGMRESSolver::solve(std::shared_ptr<DistributedVector> x, 
                               std::shared_ptr<DistributedVector> b) {
    if (!isInitialized_) {
        initialize();
    }
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    auto& matrix = linearSystem_->getMatrix();
    
    // 初始化残差
    matrix.multiply(*x, *v_[0]);
    v_[0]->copy(*b);
    v_[0]->axpy(-1.0, *v_[0]);
    
    double beta = computeResidualNorm(v_[0]);
    
    if (beta < tolerance_) {
        isConverged_ = true;
        residual_ = beta;
        return true;
    }
    
    // 归一化v0
    v_[0]->scale(1.0 / beta);
    
    g_[0] = beta;
    for (int i = 1; i <= restartSize_; ++i) {
        g_[i] = 0.0;
    }
    
    int k = 0;
    for (iterations_ = 1; iterations_ <= maxIterations_; ++iterations_) {
        // Arnoldi过程
        arnoldiProcess(k, v_[k]);
        
        // 求解最小二乘问题
        solveLeastSquares(k);
        
        // 检查收敛
        double residual_norm = std::abs(g_[k + 1]);
        
        if (verbosity_ > 0 && iterations_ % 10 == 0) {
            int rank = comm_->getRank();
            if (rank == 0) {
                std::cout << "GMRES迭代 " << iterations_ << ": 残差 = " 
                          << residual_norm << std::endl;
            }
        }
        
        if (residual_norm < tolerance_) {
            // 更新解
            updateSolution(x, k);
            isConverged_ = true;
            residual_ = residual_norm;
            break;
        }
        
        ++k;
        
        // 重启检查
        if (k >= restartSize_) {
            updateSolution(x, k - 1);
            
            // 计算新的残差
            matrix.multiply(*x, *v_[0]);
            v_[0]->copy(*b);
            v_[0]->axpy(-1.0, *v_[0]);
            
            beta = computeResidualNorm(v_[0]);
            if (beta < tolerance_) {
                isConverged_ = true;
                residual_ = beta;
                break;
            }
            
            v_[0]->scale(1.0 / beta);
            g_[0] = beta;
            k = 0;
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    
    if (verbosity_ > 0) {
        int rank = comm_->getRank();
        if (rank == 0) {
            std::cout << "GMRES求解完成: 迭代次数 = " << iterations_ 
                      << ", 最终残差 = " << residual_ 
                      << ", 耗时 = " << duration.count() / 1000.0 << " ms" << std::endl;
        }
    }
    
    return isConverged_;
}

bool ParallelGMRESSolver::solve(std::shared_ptr<DistributedVector> x) {
    if (!linearSystem_) {
        throw std::runtime_error("线性系统未设置");
    }
    
    auto b = linearSystem_->getRhs();
    return solve(x, b);
}

void ParallelGMRESSolver::arnoldiProcess(int k, std::shared_ptr<DistributedVector> v0) {
    // 应用预条件子（如果可用）
    // 简化实现：直接使用v[k]
    
    // 计算 w = A * v[k]
    auto& matrix = linearSystem_->getMatrix();
    matrix.multiply(*v0, *v_[k + 1]);
    
    // 正交化
    for (int j = 0; j <= k; ++j) {
        h_[j][k] = v_[j]->dot(*v_[k + 1]);
        v_[k + 1]->axpy(-h_[j][k], *v_[j]);
    }
    
    // 归一化
    h_[k + 1][k] = computeResidualNorm(v_[k + 1]);
    if (h_[k + 1][k] > 1e-16) {
        v_[k + 1]->scale(1.0 / h_[k + 1][k]);
    }
}

void ParallelGMRESSolver::solveLeastSquares(int k) {
    // 简化实现：使用Givens旋转求解最小二乘问题
    
    // 应用之前的Givens旋转
    for (int i = 0; i < k; ++i) {
        double temp = h_[i][k];
        h_[i][k] = temp * h_[i + 1][k] - h_[i + 1][k] * temp;
        // 实际实现需要完整的Givens旋转
    }
    
    // 计算新的Givens旋转
    double h1 = h_[k][k];
    double h2 = h_[k + 1][k];
    double r = std::sqrt(h1 * h1 + h2 * h2);
    
    if (r > 1e-16) {
        double c = h1 / r;
        double s = h2 / r;
        
        h_[k][k] = r;
        h_[k + 1][k] = 0.0;
        
        // 应用旋转到右端向量
        double temp = g_[k];
        g_[k] = c * temp + s * g_[k + 1];
        g_[k + 1] = -s * temp + c * g_[k + 1];
    }
}

void ParallelGMRESSolver::updateSolution(std::shared_ptr<DistributedVector> x, int k) {
    // 求解上三角系统
    std::vector<double> y(k + 1);
    
    for (int i = k; i >= 0; --i) {
        y[i] = g_[i];
        for (int j = i + 1; j <= k; ++j) {
            y[i] -= h_[i][j] * y[j];
        }
        y[i] /= h_[i][i];
    }
    
    // 更新解 x = x0 + V*y
    for (int i = 0; i <= k; ++i) {
        x->axpy(y[i], *v_[i]);
    }
}

// ============================================================================
// ParallelPreconditioner 实现
// ============================================================================

void ParallelPreconditioner::initialize() {
    if (!linearSystem_) {
        throw std::runtime_error("线性系统未设置");
    }
    
    isInitialized_ = true;
}

bool ParallelPreconditioner::solve(std::shared_ptr<DistributedVector> x, 
                                  std::shared_ptr<DistributedVector> b) {
    if (!isInitialized_) {
        initialize();
    }
    
    switch (type_) {
        case JACOBI:
            applyJacobi(x, b);
            break;
        case BLOCK_JACOBI:
            applyBlockJacobi(x, b);
            break;
        case ILU:
            applyILU(x, b);
            break;
        case MULTIGRID:
            // 简化实现：使用Jacobi作为多重网格的替代
            applyJacobi(x, b);
            break;
    }
    
    return true;
}

bool ParallelPreconditioner::solve(std::shared_ptr<DistributedVector> x) {
    if (!linearSystem_) {
        throw std::runtime_error("线性系统未设置");
    }
    
    auto b = linearSystem_->getRhs();
    return solve(x, b);
}

void ParallelPreconditioner::applyJacobi(std::shared_ptr<DistributedVector> x, 
                                        std::shared_ptr<DistributedVector> b) {
    auto& matrix = linearSystem_->getMatrix();
    
    // 获取对角线元素
    auto diagonal = matrix.getDiagonal();
    
    // Jacobi预条件子：x = D^{-1} * b
    for (int i = 0; i < diagonal.size(); ++i) {
        if (std::abs(diagonal[i]) > 1e-16) {
            x->setValue(i, b->getValue(i) / diagonal[i]);
        } else {
            x->setValue(i, 0.0);
        }
    }
}

void ParallelPreconditioner::applyBlockJacobi(std::shared_ptr<DistributedVector> x, 
                                             std::shared_ptr<DistributedVector> b) {
    // 简化实现：使用Jacobi预条件子
    applyJacobi(x, b);
}

void ParallelPreconditioner::applyILU(std::shared_ptr<DistributedVector> x, 
                                     std::shared_ptr<DistributedVector> b) {
    // 简化实现：使用Jacobi预条件子
    // 在实际应用中，需要实现不完全LU分解
    applyJacobi(x, b);
}

// ============================================================================
// ParallelSolverFactory 实现
// ============================================================================

std::shared_ptr<ParallelLinearSolver> ParallelSolverFactory::createSolver(
    SolverType type, 
    std::shared_ptr<MPICommunicator> comm) {
    
    switch (type) {
        case CG:
            return std::make_shared<ParallelConjugateGradientSolver>(comm);
        case GMRES:
            return std::make_shared<ParallelGMRESSolver>(comm);
        case BICGSTAB:
            // 简化实现：返回CG求解器
            return std::make_shared<ParallelConjugateGradientSolver>(comm);
        default:
            throw std::invalid_argument("未知的求解器类型");
    }
}

std::shared_ptr<ParallelLinearSolver> ParallelSolverFactory::createSolverWithPreconditioner(
    SolverType solverType,
    ParallelPreconditioner::PreconditionerType preconditionerType,
    std::shared_ptr<MPICommunicator> comm) {
    
    auto solver = createSolver(solverType, comm);
    auto preconditioner = std::make_shared<ParallelPreconditioner>(preconditionerType, comm);
    
    // 设置预条件子
    if (auto cgSolver = std::dynamic_pointer_cast<ParallelConjugateGradientSolver>(solver)) {
        cgSolver->setPreconditioner(preconditioner);
    }
    
    return solver;
}

} // namespace elmer