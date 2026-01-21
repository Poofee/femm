/**
 * @file SolverRegistry.h
 * @brief Elmer FEM求解器注册机制
 * 
 * 提供动态求解器注册、发现和实例化功能
 */

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include <vector>
#include "SolverBase.h"

namespace elmer {

/**
 * @brief 求解器工厂函数类型
 */
using SolverFactoryFunction = std::function<std::shared_ptr<SolverBase>()>;

/**
 * @brief 求解器注册信息
 */
struct SolverRegistration {
    std::string name;                    ///< 求解器名称
    std::string description;             ///< 求解器描述
    SolverFactoryFunction factory;       ///< 工厂函数
    std::vector<std::string> dependencies; ///< 依赖的求解器
    bool isParallel = false;             ///< 是否支持并行
    bool isTransient = false;            ///< 是否支持瞬态
    
    SolverRegistration() = default;
    
    SolverRegistration(const std::string& n, const std::string& desc, 
                      SolverFactoryFunction fac)
        : name(n), description(desc), factory(fac) {}
};

/**
 * @brief 求解器注册器单例类
 */
class SolverRegistry {
private:
    std::unordered_map<std::string, SolverRegistration> solvers_;
    
    // 单例模式
    SolverRegistry() = default;
    ~SolverRegistry() = default;
    
    SolverRegistry(const SolverRegistry&) = delete;
    SolverRegistry& operator=(const SolverRegistry&) = delete;
    
public:
    /**
     * @brief 获取单例实例
     */
    static SolverRegistry& getInstance() {
        static SolverRegistry instance;
        return instance;
    }
    
    /**
     * @brief 注册求解器
     * @param registration 求解器注册信息
     * @return 是否注册成功
     */
    bool registerSolver(const SolverRegistration& registration) {
        auto it = solvers_.find(registration.name);
        if (it != solvers_.end()) {
            // 已存在同名求解器
            return false;
        }
        
        solvers_[registration.name] = registration;
        return true;
    }
    
    /**
     * @brief 注销求解器
     * @param name 求解器名称
     * @return 是否注销成功
     */
    bool unregisterSolver(const std::string& name) {
        return solvers_.erase(name) > 0;
    }
    
    /**
     * @brief 创建求解器实例
     * @param name 求解器名称
     * @return 求解器实例指针，失败返回nullptr
     */
    std::shared_ptr<SolverBase> createSolver(const std::string& name) {
        auto it = solvers_.find(name);
        if (it != solvers_.end()) {
            return it->second.factory();
        }
        return nullptr;
    }
    
    /**
     * @brief 获取所有已注册的求解器名称
     */
    std::vector<std::string> getRegisteredSolvers() const {
        std::vector<std::string> names;
        for (const auto& pair : solvers_) {
            names.push_back(pair.first);
        }
        return names;
    }
    
    /**
     * @brief 检查求解器是否已注册
     */
    bool isSolverRegistered(const std::string& name) const {
        return solvers_.find(name) != solvers_.end();
    }
    
    /**
     * @brief 获取求解器信息
     */
    SolverRegistration getSolverInfo(const std::string& name) const {
        auto it = solvers_.find(name);
        if (it != solvers_.end()) {
            return it->second;
        }
        throw std::runtime_error("Solver not found: " + name);
    }
    
    /**
     * @brief 清空所有注册的求解器
     */
    void clear() {
        solvers_.clear();
    }
    
    /**
     * @brief 获取支持并行的求解器列表
     */
    std::vector<std::string> getParallelSolvers() const {
        std::vector<std::string> parallelSolvers;
        for (const auto& pair : solvers_) {
            if (pair.second.isParallel) {
                parallelSolvers.push_back(pair.first);
            }
        }
        return parallelSolvers;
    }
    
    /**
     * @brief 获取支持瞬态的求解器列表
     */
    std::vector<std::string> getTransientSolvers() const {
        std::vector<std::string> transientSolvers;
        for (const auto& pair : solvers_) {
            if (pair.second.isTransient) {
                transientSolvers.push_back(pair.first);
            }
        }
        return transientSolvers;
    }
};

/**
 * @brief 求解器注册宏
 * 
 * 用于简化求解器注册过程
 */
#define REGISTER_SOLVER(ClassName, SolverName, Description) \
    class ClassName##RegistryHelper { \
    public: \
        ClassName##RegistryHelper() { \
            SolverRegistry::getInstance().registerSolver( \
                SolverRegistration(SolverName, Description, []() { \
                    return std::make_shared<ClassName>(); \
                }) \
            ); \
        } \
    }; \
    static ClassName##RegistryHelper ClassName##RegistryInstance;

/**
 * @brief 求解器注册宏（带参数）
 */
#define REGISTER_SOLVER_WITH_PARAMS(ClassName, SolverName, Description, ...) \
    class ClassName##RegistryHelper { \
    public: \
        ClassName##RegistryHelper() { \
            SolverRegistry::getInstance().registerSolver( \
                SolverRegistration(SolverName, Description, []() { \
                    return std::make_shared<ClassName>(__VA_ARGS__); \
                }) \
            ); \
        } \
    }; \
    static ClassName##RegistryHelper ClassName##RegistryInstance;

/**
 * @brief 求解器管理器类
 * 
 * 负责管理求解器实例的生命周期和执行顺序
 */
class SolverManager {
private:
    std::vector<std::shared_ptr<SolverBase>> solvers_;
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<MaterialDatabase> materialDB_;
    std::shared_ptr<BoundaryConditionManager> bc_;
    
public:
    SolverManager() = default;
    
    /**
     * @brief 设置网格数据
     */
    void setMesh(std::shared_ptr<Mesh> mesh) {
        mesh_ = mesh;
    }
    
    /**
     * @brief 设置材料数据库
     */
    void setMaterialDatabase(std::shared_ptr<MaterialDatabase> materialDB) {
        materialDB_ = materialDB;
    }
    
    /**
     * @brief 设置边界条件管理器
     */
    void setBoundaryConditions(std::shared_ptr<BoundaryConditionManager> bc) {
        bc_ = bc;
    }
    
    // MPI支持已移除
    
    /**
     * @brief 添加求解器
     */
    void addSolver(const std::string& solverName) {
        auto solver = SolverRegistry::getInstance().createSolver(solverName);
        if (solver) {
            // 设置公共资源
            if (mesh_) solver->setMesh(mesh_);
            if (materialDB_) solver->setMaterialDatabase(materialDB_);
            if (bc_) solver->setBoundaryConditions(bc_);
            
            solvers_.push_back(solver);
        } else {
            throw std::runtime_error("Failed to create solver: " + solverName);
        }
    }
    
    /**
     * @brief 添加求解器实例
     */
    void addSolver(std::shared_ptr<SolverBase> solver) {
        solvers_.push_back(solver);
    }
    
    /**
     * @brief 移除求解器
     */
    void removeSolver(const std::string& solverName) {
        solvers_.erase(
            std::remove_if(solvers_.begin(), solvers_.end(),
                [&](const std::shared_ptr<SolverBase>& solver) {
                    return solver->getName() == solverName;
                }),
            solvers_.end()
        );
    }
    
    /**
     * @brief 获取求解器列表
     */
    std::vector<std::shared_ptr<SolverBase>> getSolvers() const {
        return solvers_;
    }
    
    /**
     * @brief 获取指定名称的求解器
     */
    std::shared_ptr<SolverBase> getSolver(const std::string& solverName) const {
        for (const auto& solver : solvers_) {
            if (solver->getName() == solverName) {
                return solver;
            }
        }
        return nullptr;
    }
    
    /**
     * @brief 初始化所有求解器
     */
    bool initializeAll() {
        bool success = true;
        for (auto& solver : solvers_) {
            if (!solver->initialize()) {
                success = false;
            }
        }
        return success;
    }
    
    /**
     * @brief 按顺序执行所有求解器
     */
    bool executeAll() {
        bool success = true;
        for (auto& solver : solvers_) {
            if (!solver->assemble()) {
                success = false;
                continue;
            }
            if (!solver->solve()) {
                success = false;
                continue;
            }
        }
        return success;
    }
    
    /**
     * @brief 执行稳态仿真
     */
    bool executeSteadyState() {
        return executeAll();
    }
    
    /**
     * @brief 执行瞬态仿真
     */
    bool executeTransient(double startTime, double endTime, double timeStep) {
        bool success = true;
        
        // 设置瞬态参数
        for (auto& solver : solvers_) {
            auto params = solver->getParameters();
            params.transient = true;
            params.timeStep = timeStep;
            solver->setParameters(params);
        }
        
        // 执行时间步进
        double currentTime = startTime;
        int timeStepIndex = 0;
        
        while (currentTime < endTime) {
            std::cout << "时间步 " << timeStepIndex << ", 时间 = " << currentTime << std::endl;
            
            // 执行当前时间步
            if (!executeAll()) {
                success = false;
                break;
            }
            
            currentTime += timeStep;
            timeStepIndex++;
            
            // 检查是否到达结束时间
            if (currentTime >= endTime) {
                break;
            }
        }
        
        return success;
    }
    
    /**
     * @brief 初始化求解器管理器
     */
    bool initialize() {
        return initializeAll();
    }
    
    /**
     * @brief 获取所有求解器的解
     */
    std::map<std::string, std::vector<double>> getSolutions() const {
        std::map<std::string, std::vector<double>> solutions;
        for (const auto& solver : solvers_) {
            solutions[solver->getName()] = solver->getSolution();
        }
        return solutions;
    }
    
    /**
     * @brief 清理资源
     */
    void cleanup() {
        for (auto& solver : solvers_) {
            // 求解器应该有清理方法
            // 这里可以添加清理逻辑
        }
        solvers_.clear();
    }
    
    /**
     * @brief 清空所有求解器
     */
    void clear() {
        solvers_.clear();
    }
    
    /**
     * @brief 获取求解器数量
     */
    size_t size() const {
        return solvers_.size();
    }
    
    /**
     * @brief 检查是否为空
     */
    bool empty() const {
        return solvers_.empty();
    }
};

} // namespace elmer