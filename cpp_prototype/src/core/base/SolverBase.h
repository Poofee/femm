/**
 * @file SolverBase.h
 * @brief Elmer FEM求解器基类接口定义
 * 
 * 提供统一的求解器接口，支持多物理场耦合和并行计算
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "Mesh.h"
#include "MaterialDatabase.h"
#include "BoundaryConditions.h"

namespace elmer {

/**
 * @brief 求解器状态枚举
 */
enum class SolverStatus {
    INITIALIZED,    ///< 已初始化
    ASSEMBLED,      ///< 矩阵已组装
    SOLVED,         ///< 已求解
    FAILED          ///< 求解失败
};

/**
 * @brief 求解器参数结构体
 */
struct SolverBaseParameters {
    std::string solverName;           ///< 求解器名称
    double tolerance = 1.0e-6;        ///< 收敛容差
    int maxIterations = 1000;         ///< 最大迭代次数
    bool linear = true;               ///< 是否为线性问题
    bool transient = false;           ///< 是否为瞬态问题
    double timeStep = 0.0;            ///< 时间步长
    int numTimeSteps = 1;             ///< 时间步数
    
    // 并行参数（MPI支持已移除）
    
    // 输出控制
    bool verbose = true;              ///< 详细输出
    int outputInterval = 1;           ///< 输出间隔
    
    SolverBaseParameters() = default;
    SolverBaseParameters(const std::string& name) : solverName(name) {}
};

/**
 * @brief 求解器基类
 * 
 * 所有Elmer FEM求解器的基类，提供统一的接口和功能
 */
class SolverBase {
protected:
    std::string name_;                    ///< 求解器名称
    SolverBaseParameters parameters_;         ///< 求解器参数
    SolverStatus status_;                 ///< 求解器状态
    
    std::shared_ptr<Mesh> mesh_;          ///< 网格数据
    std::shared_ptr<MaterialDatabase> materialDB_; ///< 材料数据库
    std::shared_ptr<BoundaryConditionManager> bc_; ///< 边界条件管理器
#ifdef USE_MPI
    std::shared_ptr<MPICommunicator> comm_; ///< MPI通信器
#endif
    
    // 求解器变量
    std::map<std::string, std::vector<double>> variables_; ///< 求解器变量
    
public:
    /**
     * @brief 构造函数
     * @param name 求解器名称
     */
    explicit SolverBase(const std::string& name = "UnnamedSolver")
        : name_(name), status_(SolverStatus::INITIALIZED) {}
    
    virtual ~SolverBase() = default;
    
    /**
     * @brief 获取求解器名称
     */
    virtual std::string getName() const { return name_; }
    
    /**
     * @brief 设置求解器参数
     */
    virtual void setParameters(const SolverBaseParameters& params) {
        parameters_ = params;
    }
    
    /**
     * @brief 获取求解器参数
     */
    virtual SolverBaseParameters getParameters() const {
        return parameters_;
    }
    
    /**
     * @brief 设置网格数据
     */
    virtual void setMesh(std::shared_ptr<Mesh> mesh) {
        mesh_ = mesh;
    }
    
    /**
     * @brief 设置材料数据库
     */
    virtual void setMaterialDatabase(std::shared_ptr<MaterialDatabase> materialDB) {
        materialDB_ = materialDB;
    }
    
    /**
     * @brief 设置边界条件
     */
    virtual void setBoundaryConditions(std::shared_ptr<BoundaryConditionManager> bc) {
        bc_ = bc;
    }
    
    /**
     * @brief 设置MPI通信器
     */
#ifdef USE_MPI
    virtual void setMPICommunicator(std::shared_ptr<MPICommunicator> comm) {
        comm_ = comm;
    }
#else
    virtual void setMPICommunicator(void* comm) {
        (void)comm; // 避免未使用参数警告
    }
#endif
    
    /**
     * @brief 初始化求解器
     */
    virtual bool initialize() {
        status_ = SolverStatus::INITIALIZED;
        return true;
    }
    
    /**
     * @brief 组装系统矩阵
     */
    virtual bool assemble() = 0;
    
    /**
     * @brief 求解系统
     */
    virtual bool solve() = 0;
    
    /**
     * @brief 获取求解结果
     */
    virtual std::vector<double> getSolution() const = 0;
    
    /**
     * @brief 获取求解器状态
     */
    virtual SolverStatus getStatus() const {
        return status_;
    }
    
    /**
     * @brief 设置变量值
     */
    virtual void setVariable(const std::string& name, const std::vector<double>& values) {
        variables_[name] = values;
    }
    
    /**
     * @brief 获取变量值
     */
    virtual std::vector<double> getVariable(const std::string& name) const {
        auto it = variables_.find(name);
        if (it != variables_.end()) {
            return it->second;
        }
        return {};
    }
    
    /**
     * @brief 检查是否支持并行计算
     */
    virtual bool supportsParallel() const {
        return false; // 默认不支持，子类可重写
    }
    
    /**
     * @brief 检查是否支持瞬态计算
     */
    virtual bool supportsTransient() const {
        return false; // 默认不支持，子类可重写
    }
    
    /**
     * @brief 获取求解器信息
     */
    virtual std::string getInfo() const {
        return "Solver: " + name_ + " (Status: " + 
               std::to_string(static_cast<int>(status_)) + ")";
    }
};

/**
 * @brief 线性求解器基类
 */
class LinearSolverBase : public SolverBase {
protected:
    std::shared_ptr<Matrix> stiffnessMatrix_;   ///< 刚度矩阵
    std::shared_ptr<Vector> rhsVector_;         ///< 右端向量
    std::shared_ptr<Vector> solution_;          ///< 解向量
    
public:
    LinearSolverBase(const std::string& name = "LinearSolver")
        : SolverBase(name) {}
    
    virtual ~LinearSolverBase() = default;
    
    /**
     * @brief 获取刚度矩阵
     */
    virtual std::shared_ptr<Matrix> getStiffnessMatrix() const {
        return stiffnessMatrix_;
    }
    
    /**
     * @brief 获取右端向量
     */
    virtual std::shared_ptr<Vector> getRhsVector() const {
        return rhsVector_;
    }
    
    /**
     * @brief 获取解向量
     */
    virtual std::vector<double> getSolution() const override {
        if (solution_) {
            std::vector<double> result(solution_->Size());
            for (int i = 0; i < solution_->Size(); ++i) {
                result[i] = (*solution_)[i];
            }
            return result;
        }
        return {};
    }
};

/**
 * @brief 非线性求解器基类
 */
class NonlinearSolverBase : public SolverBase {
protected:
    int nonlinearIterations_ = 0;               ///< 非线性迭代次数
    double nonlinearTolerance_ = 1.0e-6;       ///< 非线性收敛容差
    
public:
    NonlinearSolverBase(const std::string& name = "NonlinearSolver")
        : SolverBase(name) {}
    
    virtual ~NonlinearSolverBase() = default;
    
    /**
     * @brief 设置非线性求解参数
     */
    virtual void setNonlinearParameters(int maxIterations, double tolerance) {
        nonlinearIterations_ = maxIterations;
        nonlinearTolerance_ = tolerance;
    }
    
    /**
     * @brief 检查非线性收敛
     */
    virtual bool checkNonlinearConvergence(const std::vector<double>& residual) {
        double norm = 0.0;
        for (const auto& r : residual) {
            norm += r * r;
        }
        norm = std::sqrt(norm);
        return norm < nonlinearTolerance_;
    }
};

/**
 * @brief 瞬态求解器基类
 */
class TransientSolverBase : public SolverBase {
protected:
    double currentTime_ = 0.0;                  ///< 当前时间
    double timeStep_ = 0.0;                     ///< 时间步长
    int timeStepIndex_ = 0;                     ///< 时间步索引
    
public:
    TransientSolverBase(const std::string& name = "TransientSolver")
        : SolverBase(name) {}
    
    virtual ~TransientSolverBase() = default;
    
    /**
     * @brief 设置时间步参数
     */
    virtual void setTimeParameters(double startTime, double timeStep) {
        currentTime_ = startTime;
        timeStep_ = timeStep;
        timeStepIndex_ = 0;
    }
    
    /**
     * @brief 推进时间步
     */
    virtual void advanceTimeStep() {
        currentTime_ += timeStep_;
        timeStepIndex_++;
    }
    
    /**
     * @brief 获取当前时间
     */
    virtual double getCurrentTime() const {
        return currentTime_;
    }
    
    /**
     * @brief 获取时间步索引
     */
    virtual int getTimeStepIndex() const {
        return timeStepIndex_;
    }
    
    bool supportsTransient() const override {
        return true;
    }
};

} // namespace elmer