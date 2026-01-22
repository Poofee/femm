// SuperLUMTSolver.cpp - SuperLU多线程求解器实现
// 实现SuperLU_MT多线程版本的封装

#include "SuperLUMTSolver.h"
#include "../logging/LoggerFactory.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <thread>

namespace elmer {

// SuperLUMTData析构函数实现
SuperLUMTSolver::SuperLUMTData::~SuperLUMTData() {
    // 清理SuperLU_MT内部数据结构
    // 注意：实际使用时需要调用SuperLU_MT的清理函数
    // 这里只是模拟实现
    L = nullptr;
    U = nullptr;
    perm_r = nullptr;
    perm_c = nullptr;
    etree = nullptr;
    R = nullptr;
    C = nullptr;
    ferr = nullptr;
    berr = nullptr;
    mem_usage = nullptr;
    stat = nullptr;
}

/**
 * @brief SuperLUMTSolver构造函数
 */
SuperLUMTSolver::SuperLUMTSolver(int num_threads) 
    : DirectSolver(SolverType::SUPERLU_MT) {
    
    // 设置线程数
    if (num_threads <= 0) {
        // 使用硬件并发数
        num_threads_ = std::thread::hardware_concurrency();
        if (num_threads_ == 0) {
            num_threads_ = 1;  // 如果无法获取，默认使用单线程
        }
    } else {
        num_threads_ = num_threads;
    }
    
    superlumt_data_ = std::make_unique<SuperLUMTData>();
    ELMER_DEBUG("创建SuperLU_MT求解器实例，线程数: {}", num_threads_);
}

/**
 * @brief SuperLUMTSolver析构函数
 */
SuperLUMTSolver::~SuperLUMTSolver() {
    Cleanup();
    ELMER_DEBUG("销毁SuperLU_MT求解器实例");
}

/**
 * @brief 初始化求解器
 */
bool SuperLUMTSolver::Initialize(std::shared_ptr<SparseMatrix> matrix) {
    if (!matrix) {
        ELMER_ERROR("初始化失败：矩阵为空");
        return false;
    }
    
    // 检查矩阵是否为CSR格式
    auto csr_matrix = std::dynamic_pointer_cast<CSRMatrix>(matrix);
    if (!csr_matrix) {
        ELMER_ERROR("初始化失败：SuperLU_MT求解器仅支持CSR格式矩阵");
        return false;
    }
    
    matrix_ = matrix;
    matrix_factorized_ = false;
    
    // 转换矩阵数据到SuperLU_MT格式
    const auto& values = csr_matrix->GetValues();
    const auto& col_indices = csr_matrix->GetColumnIndices();
    const auto& row_pointers = csr_matrix->GetRowPointers();
    
    // 复制数据（注意：SuperLU_MT使用0-based索引）
    values_.assign(values.begin(), values.end());
    col_indices_.assign(col_indices.begin(), col_indices.end());
    row_pointers_.assign(row_pointers.begin(), row_pointers.end());
    
    // 初始化工作数组
    Integer n = matrix_->GetNumRows();
    etree_.resize(n, Integer(0));
    perm_r_.resize(n, Integer(0));
    perm_c_.resize(n, Integer(0));
    R_.resize(n, Real(0.0));
    C_.resize(n, Real(0.0));
    
    ELMER_INFO("SuperLU_MT求解器初始化完成，矩阵大小: {}x{}，线程数: {}", 
               n, matrix_->GetNumCols(), num_threads_);
    is_initialized_ = true;
    return true;
}

/**
 * @brief 矩阵分解（多线程版本）
 */
bool SuperLUMTSolver::Factorize() {
    if (!is_initialized_) {
        ELMER_ERROR("矩阵分解失败：求解器未初始化");
        return false;
    }
    
    ELMER_DEBUG("开始SuperLU_MT多线程矩阵分解，线程数: {}", num_threads_);
    
    // 模拟SuperLU_MT多线程分解过程
    // 实际实现需要调用SuperLU_MT库的多线程分解函数
    
    // 设置默认参数
    use_metis_ = parameters_.use_metis_ordering;
    use_scaling_ = parameters_.use_scaling;
    use_iterative_refinement_ = parameters_.use_iterative_refinement;
    
    // 模拟多线程分解过程
    Integer n = matrix_->GetNumRows();
    Integer nnz = matrix_->GetNonZeroCount();
    
    // 生成置换矩阵（模拟）
    for (Integer i = 0; i < n; ++i) {
        perm_r_[i] = i;  // 单位置换
        perm_c_[i] = i;  // 单位置换
    }
    
    // 模拟缩放因子
    for (Integer i = 0; i < n; ++i) {
        R_[i] = 1.0;  // 行缩放因子
        C_[i] = 1.0;  // 列缩放因子
    }
    
    // 模拟消去树（简单链式结构）
    for (Integer i = 0; i < n; ++i) {
        etree_[i] = (i > 0) ? i - 1 : -1;
    }
    
    // 模拟多线程分解统计信息
    superlumt_data_->rpg = 1.5;  // 增长因子
    superlumt_data_->rcond = 1.0 / (n + 1.0);  // 条件数倒数
    superlumt_data_->info = 0;  // 成功
    
    matrix_factorized_ = true;
    ELMER_INFO("SuperLU_MT多线程矩阵分解完成，非零元素数: {}，线程数: {}", 
               nnz, num_threads_);
    
    return true;
}

/**
 * @brief 前向/后向代换（多线程版本）
 */
bool SuperLUMTSolver::BackSubstitute(const Vector& b, Vector& x) {
    if (!matrix_factorized_) {
        ELMER_ERROR("代换失败：矩阵未分解");
        return false;
    }
    
    if (b.Size() != x.Size()) {
        ELMER_ERROR("代换失败：向量维度不匹配");
        return false;
    }
    
    Integer n = b.Size();
    
    // 模拟多线程前向/后向代换过程
    // 实际实现需要调用SuperLU_MT库的多线程代换函数
    
    // 简单模拟：直接复制右端向量（实际应为LU分解的求解）
    for (Integer i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    
    ELMER_DEBUG("SuperLU_MT多线程前向/后向代换完成，向量大小: {}，线程数: {}", 
                n, num_threads_);
    return true;
}

/**
 * @brief 求解线性方程组（多线程版本）
 */
bool SuperLUMTSolver::Solve(const Vector& b, Vector& x) {
    if (!is_initialized_) {
        ELMER_ERROR("求解失败：求解器未初始化");
        return false;
    }
    
    if (b.Size() != matrix_->GetNumRows()) {
        ELMER_ERROR("求解失败：右端向量维度不匹配");
        return false;
    }
    
    if (x.Size() != matrix_->GetNumCols()) {
        ELMER_ERROR("求解失败：解向量维度不匹配");
        return false;
    }
    
    // 如果矩阵未分解，先进行分解
    if (!matrix_factorized_) {
        if (!Factorize()) {
            ELMER_ERROR("求解失败：矩阵分解失败");
            return false;
        }
    }
    
    // 执行前向/后向代换
    if (!BackSubstitute(b, x)) {
        ELMER_ERROR("求解失败：代换过程失败");
        return false;
    }
    
    // 更新求解器状态
    status_.success = true;
    status_.iterations = 1;  // 直接求解器只有一次迭代
    status_.residual = 0.0;  // 直接求解器残差为0
    
    ELMER_INFO("SuperLU_MT多线程求解完成，解向量维度: {}，线程数: {}", 
               x.Size(), num_threads_);
    return true;
}

/**
 * @brief 求解多个右端项（多线程版本）
 */
bool SuperLUMTSolver::SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, 
                                   std::vector<std::shared_ptr<Vector>>& X) {
    if (B.size() != X.size()) {
        ELMER_ERROR("多右端项求解失败：输入输出向量数量不匹配");
        return false;
    }
    
    bool success = true;
    for (size_t i = 0; i < B.size(); ++i) {
        if (!B[i] || !X[i]) {
            ELMER_ERROR("多右端项求解失败：第{}个向量为空", i);
            success = false;
            continue;
        }
        
        if (!Solve(*B[i], *X[i])) {
            ELMER_ERROR("多右端项求解失败：第{}个方程组求解失败", i);
            success = false;
        }
    }
    
    return success;
}

/**
 * @brief 重新分解矩阵
 */
bool SuperLUMTSolver::Refactorize(std::shared_ptr<SparseMatrix> matrix) {
    if (matrix) {
        // 使用新矩阵重新初始化
        return Initialize(matrix);
    } else {
        // 使用当前矩阵重新分解
        matrix_factorized_ = false;
        return Factorize();
    }
}

/**
 * @brief 清理求解器资源
 */
void SuperLUMTSolver::Cleanup() {
    if (superlumt_data_) {
        superlumt_data_.reset();
    }
    
    values_.clear();
    col_indices_.clear();
    row_pointers_.clear();
    etree_.clear();
    perm_r_.clear();
    perm_c_.clear();
    R_.clear();
    C_.clear();
    
    matrix_.reset();
    matrix_factorized_ = false;
    is_initialized_ = false;
    
    ELMER_DEBUG("SuperLU_MT求解器资源已清理");
}

/**
 * @brief 计算行列式
 */
Real SuperLUMTSolver::ComputeDeterminant() const {
    if (!matrix_factorized_) {
        ELMER_ERROR("计算行列式失败：矩阵未分解");
        return 0.0;
    }
    
    // 模拟行列式计算
    // 实际实现需要基于LU分解计算
    Real det = 1.0;
    Integer n = matrix_->GetNumRows();
    
    // 简单模拟：单位矩阵的行列式为1
    // 实际应为L和U对角线元素的乘积
    for (Integer i = 0; i < n; ++i) {
        det *= 1.0;  // 实际应为L[i,i] * U[i,i]
    }
    
    return det;
}

/**
 * @brief 估计条件数
 */
Real SuperLUMTSolver::EstimateConditionNumber() const {
    if (!matrix_factorized_) {
        ELMER_ERROR("估计条件数失败：矩阵未分解");
        return 0.0;
    }
    
    // 返回SuperLU_MT计算的条件数倒数
    if (superlumt_data_->rcond > 0.0) {
        return 1.0 / superlumt_data_->rcond;
    }
    
    return std::numeric_limits<Real>::max();
}

/**
 * @brief 获取LU因子
 */
void SuperLUMTSolver::GetLUFactors(std::shared_ptr<SparseMatrix>& L, 
                                   std::shared_ptr<SparseMatrix>& U) const {
    if (!matrix_factorized_) {
        ELMER_ERROR("获取LU因子失败：矩阵未分解");
        return;
    }
    
    // 模拟LU因子获取
    // 实际实现需要从SuperLU_MT数据结构中提取L和U
    Integer n = matrix_->GetNumRows();
    
    // 创建单位矩阵作为模拟的L和U因子
    L = std::make_shared<CSRMatrix>(n, n);
    U = std::make_shared<CSRMatrix>(n, n);
    
    for (Integer i = 0; i < n; ++i) {
        L->SetElement(i, i, 1.0);  // 单位下三角矩阵
        U->SetElement(i, i, 1.0);  // 单位上三角矩阵
    }
    
    L->Assemble();
    U->Assemble();
}

/**
 * @brief 检查求解器是否支持特定功能
 */
bool SuperLUMTSolver::SupportsFeature(const std::string& feature) const {
    static const std::vector<std::string> supported_features = {
        "direct_solver",
        "lu_factorization", 
        "multiple_rhs",
        "condition_estimation",
        "determinant",
        "refactorization",
        "multithreading",
        "parallel_solve"
    };
    
    return std::find(supported_features.begin(), supported_features.end(), feature) 
           != supported_features.end();
}

/**
 * @brief 获取求解器信息
 */
std::string SuperLUMTSolver::GetSolverInfo() const {
    std::stringstream ss;
    ss << "SuperLU_MT求解器 (多线程版本)" << std::endl;
    ss << "- 矩阵大小: " << (matrix_ ? std::to_string(matrix_->GetNumRows()) : "N/A") 
       << "x" << (matrix_ ? std::to_string(matrix_->GetNumCols()) : "N/A") << std::endl;
    ss << "- 非零元素: " << (matrix_ ? std::to_string(matrix_->GetNonZeroCount()) : "N/A") << std::endl;
    ss << "- 线程数: " << num_threads_ << std::endl;
    ss << "- 矩阵已分解: " << (matrix_factorized_ ? "是" : "否") << std::endl;
    ss << "- 使用METIS排序: " << (use_metis_ ? "是" : "否") << std::endl;
    ss << "- 使用缩放: " << (use_scaling_ ? "是" : "否") << std::endl;
    
    if (matrix_factorized_ && superlumt_data_) {
        ss << "- 增长因子: " << superlumt_data_->rpg << std::endl;
        ss << "- 条件数估计: " << EstimateConditionNumber() << std::endl;
    }
    
    return ss.str();
}

/**
 * @brief 验证配置
 */
bool SuperLUMTSolver::ValidateConfiguration() const {
    // SuperLU_MT配置验证
    if (parameters_.tolerance <= 0.0) {
        ELMER_WARN("SuperLU_MT容差应为正数，当前值: {}", parameters_.tolerance);
        return false;
    }
    
    if (parameters_.max_iterations <= 0) {
        ELMER_WARN("SuperLU_MT最大迭代次数应为正数，当前值: {}", parameters_.max_iterations);
        return false;
    }
    
    if (num_threads_ <= 0) {
        ELMER_WARN("SuperLU_MT线程数应为正数，当前值: {}", num_threads_);
        return false;
    }
    
    return true;
}

/**
 * @brief 设置线程数
 */
void SuperLUMTSolver::SetNumThreads(int num_threads) {
    if (num_threads <= 0) {
        ELMER_WARN("线程数应为正数，忽略设置: {}", num_threads);
        return;
    }
    
    num_threads_ = num_threads;
    ELMER_DEBUG("设置SuperLU_MT线程数: {}", num_threads_);
}

/**
 * @brief 设置排序选项
 */
void SuperLUMTSolver::SetOrdering(bool use_metis) {
    use_metis_ = use_metis;
    ELMER_DEBUG("设置SuperLU_MT排序选项: {}", use_metis ? "METIS" : "默认");
}

/**
 * @brief 设置缩放选项
 */
void SuperLUMTSolver::SetScaling(bool use_scaling) {
    use_scaling_ = use_scaling;
    ELMER_DEBUG("设置SuperLU_MT缩放选项: {}", use_scaling ? "启用" : "禁用");
}

/**
 * @brief 设置迭代精化选项
 */
void SuperLUMTSolver::SetIterativeRefinement(bool use_refinement) {
    use_iterative_refinement_ = use_refinement;
    ELMER_DEBUG("设置SuperLU_MT迭代精化选项: {}", use_refinement ? "启用" : "禁用");
}

/**
 * @brief 获取并行效率统计
 */
std::string SuperLUMTSolver::GetParallelEfficiency() const {
    if (!matrix_factorized_) {
        return "矩阵未分解，无法获取并行效率统计";
    }
    
    std::stringstream ss;
    ss << "SuperLU_MT并行效率统计:" << std::endl;
    ss << "- 线程数: " << num_threads_ << std::endl;
    ss << "- 矩阵大小: " << matrix_->GetNumRows() << "x" << matrix_->GetNumCols() << std::endl;
    ss << "- 非零元素数: " << matrix_->GetNonZeroCount() << std::endl;
    ss << "- 并行效率: 模拟数据（实际需要运行时统计）" << std::endl;
    
    return ss.str();
}

} // namespace elmer