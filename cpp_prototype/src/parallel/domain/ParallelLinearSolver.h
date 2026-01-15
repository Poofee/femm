// ParallelLinearSolver.h - MPI并行线性求解器
// 对应Fortran模块: ParallelSolver.F90, ParallelIterativeSolver.F90

#pragma once

#include "../../parallel/mpi/DistributedLinearAlgebra.h"
#include "../../core/math/IterativeSolver.h"
#include "LinearAlgebra.h"
#include <memory>
#include <vector>
#include <functional>

namespace elmer {

/**
 * @brief 并行线性求解器基类
 */
class ParallelLinearSolver {
protected:
    std::shared_ptr<MPICommunicator> comm_;           // MPI通信�?    std::shared_ptr<DistributedLinearSystem> linearSystem_; // 分布式线性系�?    
    // 求解器参�?    double tolerance_ = 1.0e-8;        // 收敛容差
    int maxIterations_ = 1000;         // 最大迭代次�?    int verbosity_ = 0;                // 输出详细程度
    
    // 求解状�?    bool isInitialized_ = false;
    bool isConverged_ = false;
    int iterations_ = 0;
    double residual_ = 0.0;
    
public:
    ParallelLinearSolver(std::shared_ptr<MPICommunicator> comm = nullptr)
        : comm_(comm ? comm : MPIUtils::getDefaultComm()) {}
    
    virtual ~ParallelLinearSolver() = default;
    
    /**
     * @brief 设置线性系�?     */
    void setLinearSystem(std::shared_ptr<DistributedLinearSystem> linearSystem) {
        linearSystem_ = linearSystem;
        isInitialized_ = false;
    }
    
    /**
     * @brief 设置求解器参�?     */
    void setParameters(double tolerance, int maxIterations, int verbosity = 0) {
        tolerance_ = tolerance;
        maxIterations_ = maxIterations;
        verbosity_ = verbosity;
    }
    
    /**
     * @brief 初始化求解器
     */
    virtual void initialize() = 0;
    
    /**
     * @brief 求解线性系�?     * @param x 解向量（输出�?     * @param b 右端向量
     * @return 求解是否成功
     */
    virtual bool solve(std::shared_ptr<DistributedVector> x, 
                      std::shared_ptr<DistributedVector> b) = 0;
    
    /**
     * @brief 求解线性系统（使用线性系统内部的右端向量�?     * @param x 解向量（输出�?     * @return 求解是否成功
     */
    virtual bool solve(std::shared_ptr<DistributedVector> x) = 0;
    
    /**
     * @brief 获取求解状�?     */
    bool isConverged() const { return isConverged_; }
    int getIterations() const { return iterations_; }
    double getResidual() const { return residual_; }
    
    /**
     * @brief 重置求解器状�?     */
    virtual void reset() {
        isInitialized_ = false;
        isConverged_ = false;
        iterations_ = 0;
        residual_ = 0.0;
    }
    
    /**
     * @brief 获取求解器统计信�?     */
    struct SolverStatistics {
        int iterations;         // 迭代次数
        double residual;        // 最终残�?        double solveTime;       // 求解时间（秒�?        double communicationTime; // 通信时间（秒�?        bool converged;         // 收敛状�?    };
    
    virtual SolverStatistics getStatistics() const {
        SolverStatistics stats;
        stats.iterations = iterations_;
        stats.residual = residual_;
        stats.solveTime = 0.0;
        stats.communicationTime = 0.0;
        stats.converged = isConverged_;
        return stats;
    }
};

/**
 * @brief 并行共轭梯度法求解器
 */
class ParallelConjugateGradientSolver : public ParallelLinearSolver {
private:
    // 预条件子
    std::shared_ptr<ParallelLinearSolver> preconditioner_;
    
    // 工作向量
    std::shared_ptr<DistributedVector> r_;  // 残差向量
    std::shared_ptr<DistributedVector> p_;  // 搜索方向向量
    std::shared_ptr<DistributedVector> Ap_; // 矩阵向量乘积结果
    std::shared_ptr<DistributedVector> z_;  // 预条件残差向�?    
public:
    ParallelConjugateGradientSolver(std::shared_ptr<MPICommunicator> comm = nullptr)
        : ParallelLinearSolver(comm) {}
    
    /**
     * @brief 设置预条件子
     */
    void setPreconditioner(std::shared_ptr<ParallelLinearSolver> preconditioner) {
        preconditioner_ = preconditioner;
    }
    
    void initialize() override;
    bool solve(std::shared_ptr<DistributedVector> x, 
               std::shared_ptr<DistributedVector> b) override;
    bool solve(std::shared_ptr<DistributedVector> x) override;
    
private:
    /**
     * @brief 计算残差范数
     */
    double computeResidualNorm(std::shared_ptr<DistributedVector> r);
    
    /**
     * @brief 应用预条件子
     */
    void applyPreconditioner(std::shared_ptr<DistributedVector> r, 
                            std::shared_ptr<DistributedVector> z);
};

/**
 * @brief 并行GMRES求解�? */
class ParallelGMRESSolver : public ParallelLinearSolver {
private:
    int restartSize_ = 30;  // 重启大小
    
    // Krylov子空间向�?    std::vector<std::shared_ptr<DistributedVector>> v_; // 正交基向�?    std::vector<std::shared_ptr<DistributedVector>> z_; // 预条件基向量
    
    // Hessenberg矩阵
    std::vector<std::vector<double>> h_; // Hessenberg矩阵
    std::vector<double> g_;              // 右端向量
    
public:
    ParallelGMRESSolver(std::shared_ptr<MPICommunicator> comm = nullptr)
        : ParallelLinearSolver(comm) {}
    
    /**
     * @brief 设置重启大小
     */
    void setRestartSize(int restartSize) {
        restartSize_ = restartSize;
    }
    
    void initialize() override;
    bool solve(std::shared_ptr<DistributedVector> x, 
               std::shared_ptr<DistributedVector> b) override;
    bool solve(std::shared_ptr<DistributedVector> x) override;
    
private:
    /**
     * @brief Arnoldi过程
     */
    void arnoldiProcess(int k, std::shared_ptr<DistributedVector> v0);
    
    /**
     * @brief 求解最小二乘问�?     */
    void solveLeastSquares(int k);
    
    /**
     * @brief 更新解向�?     */
    void updateSolution(std::shared_ptr<DistributedVector> x, int k);
};

/**
 * @brief 并行预条件子
 */
class ParallelPreconditioner : public ParallelLinearSolver {
public:
    enum PreconditionerType {
        JACOBI,         // Jacobi预条件子
        BLOCK_JACOBI,   // 块Jacobi预条件子
        ILU,            // 不完全LU分解
        MULTIGRID       // 多重网格预条件子
    };
    
private:
    PreconditionerType type_;
    
public:
    ParallelPreconditioner(PreconditionerType type, std::shared_ptr<MPICommunicator> comm = nullptr)
        : ParallelLinearSolver(comm), type_(type) {}
    
    void initialize() override;
    bool solve(std::shared_ptr<DistributedVector> x, 
               std::shared_ptr<DistributedVector> b) override;
    bool solve(std::shared_ptr<DistributedVector> x) override;
    
private:
    /**
     * @brief Jacobi预条件子
     */
    void applyJacobi(std::shared_ptr<DistributedVector> x, 
                    std::shared_ptr<DistributedVector> b);
    
    /**
     * @brief 块Jacobi预条件子
     */
    void applyBlockJacobi(std::shared_ptr<DistributedVector> x, 
                         std::shared_ptr<DistributedVector> b);
    
    /**
     * @brief ILU预条件子
     */
    void applyILU(std::shared_ptr<DistributedVector> x, 
                 std::shared_ptr<DistributedVector> b);
};

/**
 * @brief 并行线性求解器工厂
 */
class ParallelSolverFactory {
public:
    enum SolverType {
        CG,     // 共轭梯度�?        GMRES,  // GMRES方法
        BICGSTAB // 双共轭梯度稳定法
    };
    
    /**
     * @brief 创建并行求解�?     */
    static std::shared_ptr<ParallelLinearSolver> createSolver(
        SolverType type, 
        std::shared_ptr<MPICommunicator> comm = nullptr);
    
    /**
     * @brief 创建带预条件子的并行求解�?     */
    static std::shared_ptr<ParallelLinearSolver> createSolverWithPreconditioner(
        SolverType solverType,
        ParallelPreconditioner::PreconditionerType preconditionerType,
        std::shared_ptr<MPICommunicator> comm = nullptr);
};

} // namespace elmer

