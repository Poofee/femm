// MumpsSolver.h - MUMPS求解器封装
// 封装MUMPS 5.6.2 MPI并行版本，支持64位整数

#pragma once

#include "LinearSolverInterface.h"
#include "SparseMatrix.h"
#include <memory>

// MUMPS头文件（实际使用时需要包含正确的路径）
// #include <dmumps_c.h>

namespace elmer {

/**
 * @brief MUMPS求解器封装类
 * 
 * 封装MUMPS 5.6.2 MPI并行版本，支持大规模稀疏矩阵并行求解
 */
class MumpsSolver : public DirectSolver {
private:
    // MUMPS数据结构（模拟定义，实际使用时需要包含MUMPS头文件）
    struct MumpsData {
        void* id = nullptr;           ///< MUMPS主结构
        Integer comm_fortran = 0;     ///< Fortran通信器
        Integer sym = 0;              ///< 对称性标志
        Integer par = 1;              ///< 并行模式
        Integer job = -1;             ///< 作业类型
        
        // 输入数据
        Integer n = 0;                ///< 矩阵阶数
        Integer nz = 0;               ///< 非零元素数
        Integer* irn = nullptr;       ///< 行索引（1-based）
        Integer* jcn = nullptr;       ///< 列索引（1-based）
        Real* a = nullptr;            ///< 非零值
        
        // 输出数据
        Real* rhs = nullptr;          ///< 右端项
        Real* sol = nullptr;          ///< 解向量
        Integer info[40] = {0};       ///< 信息数组
        Integer infog[40] = {0};      ///< 全局信息数组
        Real rinfo[40] = {0.0};       ///< 实数信息数组
        Real rinfog[40] = {0.0};      ///< 全局实数信息数组
        
        ~MumpsData();
    };
    
    std::unique_ptr<MumpsData> mumps_data_;    ///< MUMPS内部数据
    
    // MPI相关
    bool use_mpi_ = false;                    ///< 是否使用MPI
    Integer mpi_rank_ = 0;                    ///< MPI秩
    Integer mpi_size_ = 1;                    ///< MPI大小
    
    // 矩阵数据
    std::vector<Integer> irn_;                ///< 行索引（1-based）
    std::vector<Integer> jcn_;                ///< 列索引（1-based）
    std::vector<Real> a_;                     ///< 非零值
    
    // 并行配置
    Integer icntl_7_ = 0;                     ///< 排序策略
    Integer icntl_14_ = 20;                   ///< 内存百分比
    Integer icntl_28_ = 0;                    ///< 并行策略
    Integer icntl_29_ = 0;                    ///< 并行计算节点
    
    bool centralized_ = true;                 ///< 集中式输入
    bool distributed_rhs_ = false;            ///< 分布式右端项
    
public:
    /**
     * @brief 构造函数（串行模式）
     */
    MumpsSolver();
    
    /**
     * @brief 构造函数（MPI并行模式）
     * @param comm MPI通信器
     */
    explicit MumpsSolver(void* comm);
    
    /**
     * @brief 析构函数
     */
    ~MumpsSolver() override;
    
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
    
    // MUMPS特定方法
    
    /**
     * @brief 设置排序策略
     * @param strategy 排序策略（0=AMD, 1=AMF, 2=SCOTCH, 3=PORD, 4=PT-SCOTCH, 5=ParMETIS）
     */
    void SetOrderingStrategy(Integer strategy);
    
    /**
     * @brief 设置内存百分比
     * @param percent 内存百分比（默认20）
     */
    void SetMemoryPercentage(Integer percent);
    
    /**
     * @brief 设置并行策略
     * @param strategy 并行策略（0=自动, 1=集中式, 2=分布式）
     */
    void SetParallelStrategy(Integer strategy);
    
    /**
     * @brief 设置并行计算节点
     * @param nodes 计算节点数
     */
    void SetParallelNodes(Integer nodes);
    
    /**
     * @brief 设置矩阵对称性
     * @param symmetric 是否对称
     */
    void SetMatrixSymmetric(bool symmetric);
    
    /**
     * @brief 设置输入模式
     * @param centralized 是否集中式输入
     */
    void SetInputMode(bool centralized);
    
    /**
     * @brief 设置右端项分布模式
     * @param distributed 是否分布式右端项
     */
    void SetRhsDistribution(bool distributed);
    
    /**
     * @brief 获取MUMPS统计信息
     */
    std::string GetMumpsStatistics() const;
    
    /**
     * @brief 获取并行效率
     */
    Real GetParallelEfficiency() const;
    
    /**
     * @brief 检查MPI环境
     */
    bool CheckMPIEnvironment() const;
    
    /**
     * @brief 获取MPI秩和大小
     */
    void GetMPIInfo(Integer& rank, Integer& size) const;
    
private:
    /**
     * @brief 初始化MUMPS环境
     */
    bool InitializeMUMPSEnvironment();
    
    /**
     * @brief 转换矩阵到MUMPS格式
     */
    bool ConvertMatrixToMUMPSFormat();
    
    /**
     * @brief 设置MUMPS控制参数
     */
    void SetMUMPSControlParameters();
    
    /**
     * @brief 执行MUMPS作业
     */
    bool PerformMUMPSJob(Integer job);
    
    /**
     * @brief 分析阶段
     */
    bool AnalysisPhase();
    
    /**
     * @brief 分解阶段
     */
    bool FactorizationPhase();
    
    /**
     * @brief 求解阶段
     */
    bool SolutionPhase(const Vector& b, Vector& x);
    
    /**
     * @brief 检查MUMPS返回代码
     */
    bool CheckMUMPSInfo() const;
    
    /**
     * @brief 获取错误信息
     */
    std::string GetErrorMessage() const;
    
    /**
     * @brief 分布式矩阵处理
     */
    bool ProcessDistributedMatrix();
    
    /**
     * @brief 分布式右端项处理
     */
    bool ProcessDistributedRhs(const Vector& b, Vector& x);
    
    /**
     * @brief 收集分布式解
     */
    bool GatherDistributedSolution(Vector& x);
    
    /**
     * @brief 验证并行一致性
     */
    bool ValidateParallelConsistency() const;
};

} // namespace elmer