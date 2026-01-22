// LinearSolverInterface.h - 线性求解器统一接口
// 提供SuperLU、MUMPS等求解器的统一抽象层

#pragma once

#include "SparseMatrix.h"
#include "Types.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <sstream>

namespace elmer {

/**
 * @brief 求解器类型枚举
 */
enum class SolverType {
    SUPERLU,        ///< SuperLU串行求解器
    SUPERLU_MT,     ///< SuperLU多线程求解器
    MUMPS,          ///< MUMPS并行求解器
    PETSC,          ///< PETSc求解器
    UMFPACK,        ///< UMFPACK求解器
    EIGEN,          ///< Eigen内置求解器
    CUSTOM          ///< 自定义求解器
};

/**
 * @brief 求解器参数配置结构
 */
struct SolverParameters {
    // 通用参数
    Real tolerance = 1e-12;                    ///< 收敛容差
    Integer max_iterations = 1000;             ///< 最大迭代次数
    bool verbose = false;                      ///< 详细输出
    bool check_symmetry = false;               ///< 检查矩阵对称性
    
    // SuperLU特定参数
    Real drop_tolerance = 1e-4;                ///< 丢弃容差
    Integer panel_size = 10;                   ///< 面板大小
    Integer relax = 1;                         ///< 松弛参数
    
    // MUMPS特定参数
    Integer mumps_icntl_7 = 0;                 ///< 排序策略
    Integer mumps_icntl_14 = 20;               ///< 内存百分比
    Integer mumps_icntl_28 = 0;                ///< 并行策略
    Integer mumps_icntl_29 = 0;                ///< 并行计算节点
    
    // 性能调优参数
    bool use_metis_ordering = true;            ///< 使用METIS排序
    bool use_scaling = true;                   ///< 使用缩放
    bool use_iterative_refinement = true;      ///< 使用迭代精化
    
    SolverParameters() = default;
    
    /**
     * @brief 设置参数值
     */
    void SetParameter(const std::string& name, const std::string& value);
    
    /**
     * @brief 获取参数值
     */
    std::string GetParameter(const std::string& name) const;
};

/**
 * @brief 求解器状态信息
 */
struct SolverStatus {
    bool success = false;                      ///< 求解是否成功
    std::string message;                       ///< 状态消息
    Integer iterations = 0;                    ///< 迭代次数
    Real residual = 0.0;                       ///< 最终残差
    Real setup_time = 0.0;                     ///< 设置时间(秒)
    Real solve_time = 0.0;                     ///< 求解时间(秒)
    Real memory_usage = 0.0;                   ///< 内存使用(MB)
    
    /**
     * @brief 重置状态
     */
    void Reset();
    
    /**
     * @brief 转换为字符串
     */
    std::string ToString() const;
};

/**
 * @brief 线性求解器抽象基类
 * 
 * 提供统一的接口，支持不同的稀疏矩阵求解器
 */
class LinearSolver {
protected:
    SolverType solver_type_;                   ///< 求解器类型
    SolverParameters parameters_;              ///< 求解器参数
    SolverStatus status_;                      ///< 求解器状态
    bool is_initialized_ = false;              ///< 是否已初始化
    
public:
    /**
     * @brief 构造函数
     * @param type 求解器类型
     */
    explicit LinearSolver(SolverType type) : solver_type_(type) {}
    
    virtual ~LinearSolver() = default;
    
    // 基本信息获取
    SolverType GetSolverType() const { return solver_type_; }
    const SolverParameters& GetParameters() const { return parameters_; }
    const SolverStatus& GetStatus() const { return status_; }
    bool IsInitialized() const { return is_initialized_; }
    
    // 参数配置
    void SetParameters(const SolverParameters& params) { parameters_ = params; }
    void SetParameter(const std::string& name, const std::string& value);
    
    // 纯虚函数接口
    
    /**
     * @brief 初始化求解器
     * @param matrix 稀疏矩阵
     * @return 是否初始化成功
     */
    virtual bool Initialize(std::shared_ptr<SparseMatrix> matrix) = 0;
    
    /**
     * @brief 求解线性方程组 Ax = b
     * @param b 右端向量
     * @param x 解向量（输出）
     * @return 求解状态
     */
    virtual bool Solve(const Vector& b, Vector& x) = 0;
    
    /**
     * @brief 求解多个右端项
     * @param B 右端矩阵（列向量组成）
     * @param X 解矩阵（输出）
     * @return 求解状态
     */
    virtual bool SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) = 0;
    
    /**
     * @brief 重新分解矩阵（当矩阵结构不变，仅值变化时使用）
     * @param matrix 新的矩阵（可选，如果为nullptr则使用原矩阵）
     * @return 是否成功
     */
    virtual bool Refactorize(std::shared_ptr<SparseMatrix> matrix = nullptr) = 0;
    
    /**
     * @brief 清理求解器资源
     */
    virtual void Cleanup() = 0;
    
    /**
     * @brief 检查求解器是否支持特定功能
     */
    virtual bool SupportsFeature(const std::string& feature) const = 0;
    
    /**
     * @brief 获取求解器信息
     */
    virtual std::string GetSolverInfo() const = 0;
    
    /**
     * @brief 验证求解器配置
     */
    virtual bool ValidateConfiguration() const = 0;
    
protected:
    /**
     * @brief 更新求解器状态
     */
    void UpdateStatus(bool success, const std::string& message = "");
    
    /**
     * @brief 记录时间信息
     */
    void RecordTime(bool is_setup);
};

/**
 * @brief 直接求解器基类
 * 
 * 提供直接求解器的通用功能
 */
class DirectSolver : public LinearSolver {
protected:
    std::shared_ptr<SparseMatrix> matrix_;     ///< 当前矩阵
    bool matrix_factorized_ = false;           ///< 矩阵是否已分解
    
public:
    DirectSolver(SolverType type) : LinearSolver(type) {}
    
    // LinearSolver接口实现
    virtual bool Initialize(std::shared_ptr<SparseMatrix> matrix) = 0;
    virtual bool Solve(const Vector& b, Vector& x) = 0;
    virtual bool SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) = 0;
    virtual bool Refactorize(std::shared_ptr<SparseMatrix> matrix) = 0;
    virtual void Cleanup() = 0;
    
    // 直接求解器特定方法
    
    /**
     * @brief 矩阵分解
     */
    virtual bool Factorize() = 0;
    
    /**
     * @brief 前向/后向代换
     */
    virtual bool BackSubstitute(const Vector& b, Vector& x) = 0;
    
    /**
     * @brief 计算行列式
     */
    virtual Real ComputeDeterminant() const = 0;
    
    /**
     * @brief 计算条件数估计
     */
    virtual Real EstimateConditionNumber() const = 0;
    
    /**
     * @brief 获取矩阵的L和U因子（如果可用）
     */
    virtual void GetLUFactors(std::shared_ptr<SparseMatrix>& L, 
                             std::shared_ptr<SparseMatrix>& U) const = 0;
};

/**
 * @brief 迭代求解器基类
 * 
 * 提供迭代求解器的通用功能
 */
class IterativeSolver : public LinearSolver {
protected:
    std::shared_ptr<SparseMatrix> matrix_;     ///< 当前矩阵
    std::shared_ptr<LinearSolver> preconditioner_; ///< 预条件子
    
public:
    IterativeSolver(SolverType type) : LinearSolver(type) {}
    
    // LinearSolver接口实现
    bool Initialize(std::shared_ptr<SparseMatrix> matrix) override;
    bool Solve(const Vector& b, Vector& x) override;
    bool SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) override;
    bool Refactorize(std::shared_ptr<SparseMatrix> matrix) override;
    void Cleanup() override;
    
    // 迭代求解器特定方法
    
    /**
     * @brief 设置预条件子
     */
    virtual void SetPreconditioner(std::shared_ptr<LinearSolver> preconditioner);
    
    /**
     * @brief 获取当前残差历史
     */
    virtual std::vector<Real> GetResidualHistory() const = 0;
    
    /**
     * @brief 重启求解器
     */
    virtual void Restart() = 0;
    
    /**
     * @brief 设置初始猜测
     */
    virtual void SetInitialGuess(const Vector& x0) = 0;
};

/**
 * @brief 求解器异常类
 */
class SolverException : public std::runtime_error {
private:
    SolverType solver_type_;
    std::string detailed_message_;
    
public:
    SolverException(SolverType type, const std::string& message, 
                   const std::string& details = "")
        : std::runtime_error(message), solver_type_(type), detailed_message_(details) {}
    
    SolverType GetSolverType() const { return solver_type_; }
    const std::string& GetDetailedMessage() const { return detailed_message_; }
    
    /**
     * @brief 获取完整的错误信息
     */
    std::string GetFullMessage() const;
};

/**
 * @brief 求解器工厂类
 */
class SolverFactory {
public:
    /**
     * @brief 创建求解器实例
     * @param type 求解器类型
     * @return 求解器指针
     */
    static std::unique_ptr<LinearSolver> CreateSolver(SolverType type);
    
    /**
     * @brief 从字符串创建求解器
     * @param solver_name 求解器名称（"superlu", "mumps", 等）
     * @return 求解器指针
     */
    static std::unique_ptr<LinearSolver> CreateSolver(const std::string& solver_name);
    
    /**
     * @brief 获取可用的求解器列表
     */
    static std::vector<std::string> GetAvailableSolvers();
    
    /**
     * @brief 检查求解器是否可用
     */
    static bool IsSolverAvailable(SolverType type);
    static bool IsSolverAvailable(const std::string& solver_name);
    
    /**
     * @brief 获取求解器推荐配置
     */
    static SolverParameters GetRecommendedParameters(SolverType type);
};

/**
 * @brief 求解器管理器类
 * 
 * 管理多个求解器实例，提供统一的接口
 */
class SolverManager {
private:
    std::unordered_map<std::string, std::unique_ptr<LinearSolver>> solvers_;
    std::string default_solver_;
    
public:
    SolverManager() = default;
    
    /**
     * @brief 添加求解器
     * @param name 求解器名称
     * @param solver 求解器实例
     */
    void AddSolver(const std::string& name, std::unique_ptr<LinearSolver> solver);
    
    /**
     * @brief 获取求解器
     * @param name 求解器名称
     * @return 求解器指针
     */
    LinearSolver* GetSolver(const std::string& name);
    
    /**
     * @brief 移除求解器
     */
    void RemoveSolver(const std::string& name);
    
    /**
     * @brief 设置默认求解器
     */
    void SetDefaultSolver(const std::string& name);
    
    /**
     * @brief 使用默认求解器求解
     */
    bool Solve(const Vector& b, Vector& x);
    
    /**
     * @brief 使用指定求解器求解
     */
    bool Solve(const std::string& solver_name, const Vector& b, Vector& x);
    
    /**
     * @brief 获取所有求解器名称
     */
    std::vector<std::string> GetSolverNames() const;
    
    /**
     * @brief 清理所有求解器
     */
    void CleanupAll();
    
    /**
     * @brief 获取求解器统计信息
     */
    std::string GetStatistics() const;
};

// ============================================================================
// 内联函数实现
// ============================================================================

inline void SolverStatus::Reset() {
    success = false;
    message.clear();
    iterations = 0;
    residual = 0.0;
    setup_time = 0.0;
    solve_time = 0.0;
    memory_usage = 0.0;
}

inline std::string SolverStatus::ToString() const {
    std::ostringstream oss;
    oss << "状态: " << (success ? "成功" : "失败") << "\n"
        << "消息: " << message << "\n"
        << "迭代次数: " << iterations << "\n"
        << "残差: " << std::scientific << residual << "\n"
        << "设置时间: " << setup_time << "秒\n"
        << "求解时间: " << solve_time << "秒\n"
        << "内存使用: " << memory_usage << "MB";
    return oss.str();
}

inline std::string SolverException::GetFullMessage() const {
    std::ostringstream oss;
    oss << "求解器错误 (" << static_cast<int>(solver_type_) << "): " 
        << what();
    if (!detailed_message_.empty()) {
        oss << "\n详细信息: " << detailed_message_;
    }
    return oss.str();
}

} // namespace elmer