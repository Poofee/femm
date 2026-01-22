// LinearSolverFactory.cpp - 求解器工厂和管理器实现
// 提供求解器的动态创建和管理功能

#include "LinearSolverInterface.h"
#include "SuperLUSolver.h"
#include "SuperLUMTSolver.h"
#include "MumpsSolver.h"
#include <map>
#include <algorithm>
#include <chrono>
#include <sstream>

namespace elmer {

// ============================================================================
// SolverFactory 实现
// ============================================================================

std::unique_ptr<LinearSolver> SolverFactory::CreateSolver(SolverType type) {
    switch (type) {
        case SolverType::SUPERLU:
            return std::make_unique<SuperLUSolver>();
            
        case SolverType::SUPERLU_MT:
            return std::make_unique<SuperLUMTSolver>();
            
        case SolverType::MUMPS:
            return std::make_unique<MumpsSolver>();
            
        case SolverType::PETSC:
            // TODO: 实现PETSc求解器
            throw SolverException(type, "PETSc求解器尚未实现");
            
        case SolverType::UMFPACK:
            // TODO: 实现UMFPACK求解器
            throw SolverException(type, "UMFPACK求解器尚未实现");
            
        case SolverType::EIGEN:
            // TODO: 实现Eigen求解器
            throw SolverException(type, "Eigen求解器尚未实现");
            
        case SolverType::CUSTOM:
            throw SolverException(type, "自定义求解器需要用户提供实现");
            
        default:
            throw SolverException(type, "未知的求解器类型");
    }
}

std::unique_ptr<LinearSolver> SolverFactory::CreateSolver(const std::string& solver_name) {
    static const std::map<std::string, SolverType> solver_map = {
        {"superlu", SolverType::SUPERLU},
        {"superlu_mt", SolverType::SUPERLU_MT},
        {"mumps", SolverType::MUMPS},
        {"petsc", SolverType::PETSC},
        {"umfpack", SolverType::UMFPACK},
        {"eigen", SolverType::EIGEN},
        {"custom", SolverType::CUSTOM}
    };
    
    std::string lower_name = solver_name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
    
    auto it = solver_map.find(lower_name);
    if (it != solver_map.end()) {
        return CreateSolver(it->second);
    }
    
    throw SolverException(SolverType::CUSTOM, "未知的求解器名称: " + solver_name);
}

std::vector<std::string> SolverFactory::GetAvailableSolvers() {
    return {"superlu", "superlu_mt", "mumps", "petsc", "umfpack", "eigen"};
}

bool SolverFactory::IsSolverAvailable(SolverType type) {
    switch (type) {
        case SolverType::SUPERLU:
        case SolverType::SUPERLU_MT:
        case SolverType::MUMPS:
            // 这些求解器已实现
            return true;
            
        case SolverType::PETSC:
        case SolverType::UMFPACK:
        case SolverType::EIGEN:
            // 这些求解器尚未实现
            return false;
            
        default:
            return false;
    }
}

bool SolverFactory::IsSolverAvailable(const std::string& solver_name) {
    static const std::map<std::string, SolverType> solver_map = {
        {"superlu", SolverType::SUPERLU},
        {"superlu_mt", SolverType::SUPERLU_MT},
        {"mumps", SolverType::MUMPS},
        {"petsc", SolverType::PETSC},
        {"umfpack", SolverType::UMFPACK},
        {"eigen", SolverType::EIGEN}
    };
    
    std::string lower_name = solver_name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
    
    auto it = solver_map.find(lower_name);
    if (it != solver_map.end()) {
        return IsSolverAvailable(it->second);
    }
    
    return false;
}

SolverParameters SolverFactory::GetRecommendedParameters(SolverType type) {
    SolverParameters params;
    
    switch (type) {
        case SolverType::SUPERLU:
            params.tolerance = 1e-12;
            params.drop_tolerance = 1e-4;
            params.panel_size = 10;
            params.relax = 1;
            params.use_metis_ordering = true;
            params.use_scaling = true;
            params.use_iterative_refinement = true;
            break;
            
        case SolverType::SUPERLU_MT:
            params.tolerance = 1e-12;
            params.drop_tolerance = 1e-4;
            params.panel_size = 16; // 多线程优化
            params.relax = 2;
            params.use_metis_ordering = true;
            params.use_scaling = true;
            params.use_iterative_refinement = true;
            break;
            
        case SolverType::MUMPS:
            params.tolerance = 1e-12;
            params.mumps_icntl_7 = 0; // AMD排序
            params.mumps_icntl_14 = 20; // 20%额外内存
            params.mumps_icntl_28 = 0; // 自动并行策略
            params.mumps_icntl_29 = 0; // 自动计算节点
            params.use_metis_ordering = true;
            params.use_scaling = true;
            params.use_iterative_refinement = true;
            break;
            
        default:
            // 使用默认参数
            break;
    }
    
    return params;
}

// ============================================================================
// SolverManager 实现
// ============================================================================

void SolverManager::AddSolver(const std::string& name, std::unique_ptr<LinearSolver> solver) {
    if (solvers_.find(name) != solvers_.end()) {
        throw std::runtime_error("求解器名称已存在: " + name);
    }
    
    solvers_[name] = std::move(solver);
    
    // 如果没有默认求解器，设置第一个添加的求解器为默认
    if (default_solver_.empty()) {
        default_solver_ = name;
    }
}

LinearSolver* SolverManager::GetSolver(const std::string& name) {
    auto it = solvers_.find(name);
    if (it == solvers_.end()) {
        throw std::runtime_error("求解器不存在: " + name);
    }
    
    return it->second.get();
}

void SolverManager::RemoveSolver(const std::string& name) {
    auto it = solvers_.find(name);
    if (it != solvers_.end()) {
        // 如果移除的是默认求解器，需要重新设置默认求解器
        if (default_solver_ == name) {
            default_solver_.clear();
            if (!solvers_.empty()) {
                default_solver_ = solvers_.begin()->first;
            }
        }
        
        solvers_.erase(it);
    }
}

void SolverManager::SetDefaultSolver(const std::string& name) {
    if (solvers_.find(name) == solvers_.end()) {
        throw std::runtime_error("求解器不存在: " + name);
    }
    
    default_solver_ = name;
}

bool SolverManager::Solve(const Vector& b, Vector& x) {
    if (default_solver_.empty()) {
        throw std::runtime_error("没有设置默认求解器");
    }
    
    return Solve(default_solver_, b, x);
}

bool SolverManager::Solve(const std::string& solver_name, const Vector& b, Vector& x) {
    LinearSolver* solver = GetSolver(solver_name);
    
    if (!solver->IsInitialized()) {
        throw std::runtime_error("求解器未初始化: " + solver_name);
    }
    
    return solver->Solve(b, x);
}

std::vector<std::string> SolverManager::GetSolverNames() const {
    std::vector<std::string> names;
    for (const auto& pair : solvers_) {
        names.push_back(pair.first);
    }
    return names;
}

void SolverManager::CleanupAll() {
    for (auto& pair : solvers_) {
        pair.second->Cleanup();
    }
}

std::string SolverManager::GetStatistics() const {
    std::ostringstream oss;
    
    oss << "求解器管理器统计信息:\n";
    oss << "求解器数量: " << solvers_.size() << "\n";
    oss << "默认求解器: " << (default_solver_.empty() ? "无" : default_solver_) << "\n\n";
    
    for (const auto& pair : solvers_) {
        const auto& solver = pair.second;
        const auto& status = solver->GetStatus();
        
        oss << "求解器: " << pair.first << "\n";
        oss << "类型: " << static_cast<int>(solver->GetSolverType()) << "\n";
        oss << "初始化: " << (solver->IsInitialized() ? "是" : "否") << "\n";
        oss << "成功求解次数: " << (status.success ? 1 : 0) << "\n";
        oss << "总求解时间: " << status.solve_time << "秒\n";
        oss << "内存使用: " << status.memory_usage << "MB\n";
        oss << "---\n";
    }
    
    return oss.str();
}

// ============================================================================
// SolverParameters 实现
// ============================================================================

void SolverParameters::SetParameter(const std::string& name, const std::string& value) {
    std::string lower_name = name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
    
    if (lower_name == "tolerance") {
        tolerance = std::stod(value);
    } else if (lower_name == "max_iterations") {
        max_iterations = std::stoi(value);
    } else if (lower_name == "verbose") {
        verbose = (value == "true" || value == "1");
    } else if (lower_name == "check_symmetry") {
        check_symmetry = (value == "true" || value == "1");
    } else if (lower_name == "drop_tolerance") {
        drop_tolerance = std::stod(value);
    } else if (lower_name == "panel_size") {
        panel_size = std::stoi(value);
    } else if (lower_name == "relax") {
        relax = std::stoi(value);
    } else if (lower_name == "mumps_icntl_7") {
        mumps_icntl_7 = std::stoi(value);
    } else if (lower_name == "mumps_icntl_14") {
        mumps_icntl_14 = std::stoi(value);
    } else if (lower_name == "mumps_icntl_28") {
        mumps_icntl_28 = std::stoi(value);
    } else if (lower_name == "mumps_icntl_29") {
        mumps_icntl_29 = std::stoi(value);
    } else if (lower_name == "use_metis_ordering") {
        use_metis_ordering = (value == "true" || value == "1");
    } else if (lower_name == "use_scaling") {
        use_scaling = (value == "true" || value == "1");
    } else if (lower_name == "use_iterative_refinement") {
        use_iterative_refinement = (value == "true" || value == "1");
    } else {
        throw std::invalid_argument("未知的参数名称: " + name);
    }
}

std::string SolverParameters::GetParameter(const std::string& name) const {
    std::string lower_name = name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
    
    if (lower_name == "tolerance") {
        return std::to_string(tolerance);
    } else if (lower_name == "max_iterations") {
        return std::to_string(max_iterations);
    } else if (lower_name == "verbose") {
        return verbose ? "true" : "false";
    } else if (lower_name == "check_symmetry") {
        return check_symmetry ? "true" : "false";
    } else if (lower_name == "drop_tolerance") {
        return std::to_string(drop_tolerance);
    } else if (lower_name == "panel_size") {
        return std::to_string(panel_size);
    } else if (lower_name == "relax") {
        return std::to_string(relax);
    } else if (lower_name == "mumps_icntl_7") {
        return std::to_string(mumps_icntl_7);
    } else if (lower_name == "mumps_icntl_14") {
        return std::to_string(mumps_icntl_14);
    } else if (lower_name == "mumps_icntl_28") {
        return std::to_string(mumps_icntl_28);
    } else if (lower_name == "mumps_icntl_29") {
        return std::to_string(mumps_icntl_29);
    } else if (lower_name == "use_metis_ordering") {
        return use_metis_ordering ? "true" : "false";
    } else if (lower_name == "use_scaling") {
        return use_scaling ? "true" : "false";
    } else if (lower_name == "use_iterative_refinement") {
        return use_iterative_refinement ? "true" : "false";
    } else {
        throw std::invalid_argument("未知的参数名称: " + name);
    }
}

// ============================================================================
// LinearSolver 基类实现
// ============================================================================

void LinearSolver::SetParameter(const std::string& name, const std::string& value) {
    parameters_.SetParameter(name, value);
}

void LinearSolver::UpdateStatus(bool success, const std::string& message) {
    status_.success = success;
    status_.message = message;
}

void LinearSolver::RecordTime(bool is_setup) {
    static auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double seconds = duration.count() / 1000000.0;
    
    if (is_setup) {
        status_.setup_time += seconds;
    } else {
        status_.solve_time += seconds;
    }
    
    start_time = end_time;
}

// ============================================================================
// DirectSolver 基类实现
// ============================================================================

// DirectSolver 是一个抽象基类，具体实现由 SuperLUSolver 和 MumpsSolver 提供
// 这些函数应该由具体的求解器类实现

// ============================================================================
// IterativeSolver 基类实现
// ============================================================================

bool IterativeSolver::Initialize(std::shared_ptr<SparseMatrix> matrix) {
    if (!matrix) {
        UpdateStatus(false, "矩阵指针为空");
        return false;
    }
    
    matrix_ = matrix;
    is_initialized_ = true;
    
    UpdateStatus(true, "迭代求解器初始化成功");
    return true;
}

bool IterativeSolver::Solve(const Vector& b, Vector& x) {
    // 迭代求解器的具体实现在子类中完成
    UpdateStatus(false, "迭代求解器基类不能直接求解");
    return false;
}

bool IterativeSolver::SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) {
    // 多右端项求解需要特殊处理
    // 迭代求解器通常不支持多右端项求解，需要由具体实现类处理
    throw SolverException(GetSolverType(), "迭代求解器多右端项求解尚未实现");
}

bool IterativeSolver::Refactorize(std::shared_ptr<SparseMatrix> matrix) {
    if (matrix) {
        matrix_ = matrix;
    }
    
    UpdateStatus(true, "迭代求解器无需重新分解");
    return true;
}

void IterativeSolver::Cleanup() {
    matrix_.reset();
    preconditioner_.reset();
    is_initialized_ = false;
    status_.Reset();
}

void IterativeSolver::SetPreconditioner(std::shared_ptr<LinearSolver> preconditioner) {
    preconditioner_ = preconditioner;
}

// ============================================================================
// 具体求解器类的SolveMultiple实现
// ============================================================================

bool SuperLUSolver::SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) {
    if (!is_initialized_) {
        UpdateStatus(false, "求解器未初始化");
        return false;
    }
    
    if (B.empty()) {
        UpdateStatus(false, "右端项为空");
        return false;
    }
    
    if (!matrix_factorized_) {
        if (!Factorize()) {
            UpdateStatus(false, "矩阵分解失败");
            return false;
        }
    }
    
    X.resize(B.size());
    bool all_success = true;
    
    RecordTime(false); // 开始求解计时
    
    for (size_t i = 0; i < B.size(); ++i) {
        // 检查向量大小是否匹配
        if (!X[i] || X[i]->Size() != B[i]->Size()) {
            // 如果大小不匹配或向量未初始化，需要重新创建向量
            X[i] = Vector::Create(B[i]->Size());
        }
        if (!BackSubstitute(*B[i], *X[i])) {
            all_success = false;
            break;
        }
    }
    
    RecordTime(false); // 结束求解计时
    
    if (all_success) {
        UpdateStatus(true, "SuperLU多右端项求解成功");
    } else {
        UpdateStatus(false, "SuperLU多右端项求解失败");
    }
    
    return all_success;
}

bool SuperLUMTSolver::SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) {
    if (!is_initialized_) {
        UpdateStatus(false, "求解器未初始化");
        return false;
    }
    
    if (B.empty()) {
        UpdateStatus(false, "右端项为空");
        return false;
    }
    
    if (!matrix_factorized_) {
        if (!Factorize()) {
            UpdateStatus(false, "矩阵分解失败");
            return false;
        }
    }
    
    X.resize(B.size());
    bool all_success = true;
    
    RecordTime(false); // 开始求解计时
    
    for (size_t i = 0; i < B.size(); ++i) {
        // 检查向量大小是否匹配
        if (!X[i] || X[i]->Size() != B[i]->Size()) {
            // 如果大小不匹配或向量未初始化，需要重新创建向量
            X[i] = Vector::Create(B[i]->Size());
        }
        if (!BackSubstitute(*B[i], *X[i])) {
            all_success = false;
            break;
        }
    }
    
    RecordTime(false); // 结束求解计时
    
    if (all_success) {
        UpdateStatus(true, "SuperLU_MT多右端项求解成功");
    } else {
        UpdateStatus(false, "SuperLU_MT多右端项求解失败");
    }
    
    return all_success;
}

bool MumpsSolver::SolveMultiple(const std::vector<std::shared_ptr<Vector>>& B, std::vector<std::shared_ptr<Vector>>& X) {
    if (!is_initialized_) {
        UpdateStatus(false, "求解器未初始化");
        return false;
    }
    
    if (B.empty()) {
        UpdateStatus(false, "右端项为空");
        return false;
    }
    
    if (!matrix_factorized_) {
        if (!Factorize()) {
            UpdateStatus(false, "矩阵分解失败");
            return false;
        }
    }
    
    X.resize(B.size());
    bool all_success = true;
    
    RecordTime(false); // 开始求解计时
    
    for (size_t i = 0; i < B.size(); ++i) {
        // 检查向量大小是否匹配
        if (!X[i] || X[i]->Size() != B[i]->Size()) {
            // 如果大小不匹配或向量未初始化，需要重新创建向量
            X[i] = Vector::Create(B[i]->Size());
        }
        if (!BackSubstitute(*B[i], *X[i])) {
            all_success = false;
            break;
        }
    }
    
    RecordTime(false); // 结束求解计时
    
    if (all_success) {
        UpdateStatus(true, "MUMPS多右端项求解成功");
    } else {
        UpdateStatus(false, "MUMPS多右端项求解失败");
    }
    
    return all_success;
}

} // namespace elmer