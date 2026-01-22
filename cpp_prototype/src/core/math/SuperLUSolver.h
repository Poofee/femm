// SuperLUSolver.h - SuperLU求解器封装
// 封装SuperLU 5.3.0串行双精度版本

#pragma once

#include "LinearSolverInterface.h"
#include "SparseMatrix.h"
#include <memory>

// SuperLU头文件（实际使用时需要包含正确的路径）
// #include <slu_ddefs.h>

namespace elmer {

/**
 * @brief SuperLU求解器封装类
 * 
 * 封装SuperLU 5.3.0串行双精度版本，提供高效的稀疏矩阵直接求解
 */
class SuperLUSolver : public DirectSolver {
private:
    // SuperLU数据结构（模拟定义，实际使用时需要包含SuperLU头文件）
    struct SuperLUData {
        void* L = nullptr;          ///< L因子
        void* U = nullptr;          ///< U因子
        void* perm_r = nullptr;     ///< 行置换
        void* perm_c = nullptr;     ///< 列置换
        void* etree = nullptr;      ///< 消去树
        void* R = nullptr;          ///< 行缩放因子
        void* C = nullptr;          ///< 列缩放因子
        void* ferr = nullptr;       ///< 前向误差界
        void* berr = nullptr;       ///< 后向误差界
        void* mem_usage = nullptr;  ///< 内存使用统计
        void* stat = nullptr;       ///< 统计信息
        
        Integer info = 0;           ///< 返回信息
        Real rpg = 0.0;             ///< 增长因子
        Real rcond = 0.0;           ///< 条件数倒数
    };
    
    std::unique_ptr<SuperLUData> superlu_data_; ///< SuperLU内部数据
    
    // 矩阵转换相关
    std::vector<Integer> row_pointers_;         ///< 行指针（0-based）
    std::vector<Integer> col_indices_;          ///< 列索引（0-based）
    std::vector<Real> values_;                  ///< 非零值
    
    // 工作数组
    std::vector<Integer> etree_;                ///< 消去树
    std::vector<Integer> perm_r_;               ///< 行置换
    std::vector<Integer> perm_c_;               ///< 列置换
    std::vector<Real> R_;                       ///< 行缩放
    std::vector<Real> C_;                       ///< 列缩放
    
    bool use_metis_ = true;                     ///< 使用METIS排序
    bool use_scaling_ = true;                   ///< 使用缩放
    bool use_iterative_refinement_ = true;      ///< 使用迭代精化
    
public:
    /**
     * @brief 构造函数
     */
    SuperLUSolver();
    
    /**
     * @brief 析构函数
     */
    ~SuperLUSolver() override;
    
    // DirectSolver接口实现
    bool Factorize() override;
    bool BackSubstitute(const Vector& b, Vector& x) override;
    Real ComputeDeterminant() const override;
    Real EstimateConditionNumber() const override;
    void GetLUFactors(std::shared_ptr<SparseMatrix>& L, 
                     std::shared_ptr<SparseMatrix>& U) const override;
    
    // LinearSolver接口实现
    bool Initialize(std::shared_ptr<SparseMatrix> matrix) override;
    bool Solve(const Vector& b, Vector& x) override;
    bool SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) override;
    bool Refactorize(std::shared_ptr<SparseMatrix> matrix) override;
    void Cleanup() override;
    
    bool SupportsFeature(const std::string& feature) const override;
    std::string GetSolverInfo() const override;
    bool ValidateConfiguration() const override;
    
    // SuperLU特定方法
    
    /**
     * @brief 设置排序选项
     * @param use_metis 是否使用METIS排序
     */
    void SetOrdering(bool use_metis);
    
    /**
     * @brief 设置缩放选项
     * @param use_scaling 是否使用缩放
     */
    void SetScaling(bool use_scaling);
    
    /**
     * @brief 设置迭代精化选项
     * @param use_refinement 是否使用迭代精化
     */
    void SetIterativeRefinement(bool use_refinement);
    
    /**
     * @brief 获取SuperLU统计信息
     */
    std::string GetSuperLUStatistics() const;
    
    /**
     * @brief 获取内存使用详情
     */
    std::string GetMemoryUsageDetails() const;
    
    /**
     * @brief 检查矩阵是否适合SuperLU求解
     */
    bool IsMatrixSuitable() const;
    
private:
    /**
     * @brief 转换矩阵到SuperLU格式
     */
    bool ConvertMatrixToSuperLUFormat();
    
    /**
     * @brief 初始化SuperLU选项
     */
    void InitializeSuperLUOptions();
    
    /**
     * @brief 执行矩阵分解
     */
    bool PerformFactorization();
    
    /**
     * @brief 执行前向/后向代换
     */
    bool PerformBackSubstitution(const Vector& b, Vector& x);
    
    /**
     * @brief 检查SuperLU返回代码
     */
    bool CheckSuperLUInfo(Integer info) const;
    
    /**
     * @brief 获取错误信息
     */
    std::string GetErrorMessage(Integer info) const;
    
    /**
     * @brief 计算残差
     */
    Real ComputeResidual(const Vector& b, const Vector& x) const;
    
    /**
     * @brief 验证解的正确性
     */
    bool ValidateSolution(const Vector& b, const Vector& x, Real tolerance = 1e-12) const;
};

} // namespace elmer