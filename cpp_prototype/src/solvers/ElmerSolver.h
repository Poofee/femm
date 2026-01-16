#pragma once

#include "../core/math/LinearAlgebra.h"
#include "Mesh.h"
#include "../core/base/Types.h"
#include "../boundary/BoundaryConditions.h"
#include <algorithm>
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace elmer {

/**
 * @brief Elmer求解器基类
 * 
 * 提供求解器的基本功能，包括网格管理、边界条件处理和求解过程
 */
class ElmerSolver {
public:
    ElmerSolver() = default;
    virtual ~ElmerSolver() = default;
    
    /**
     * @brief 初始化求解器
     */
    virtual bool initialize() = 0;
    
    /**
     * @brief 加载网格数据
     */
    virtual bool loadMesh(std::shared_ptr<Mesh> mesh) = 0;
    
    /**
     * @brief 设置求解参数
     */
    virtual void setParameters(const std::unordered_map<std::string, double>& params) = 0;
    
    /**
     * @brief 添加边界条件
     */
    virtual void addBoundaryCondition(std::shared_ptr<BoundaryCondition> bc) = 0;
    
    /**
     * @brief 组装系统矩阵
     */
    virtual bool assembleSystem() = 0;
    
    /**
     * @brief 求解系统
     */
    virtual bool solve() = 0;
    
    /**
     * @brief 获取求解结果
     */
    virtual std::vector<double> getSolution() const = 0;
    
    /**
     * @brief 保存结果
     */
    virtual bool saveResults(const std::string& filename) = 0;
    
    /**
     * @brief 获取求解器名称
     */
    virtual std::string getName() const = 0;
    
    /**
     * @brief 获取求解器类型
     */
    virtual std::string getType() const = 0;
    
    /**
     * @brief 检查求解器是否已初始化
     */
    virtual bool isInitialized() const = 0;
    
    /**
     * @brief 检查网格是否已加载
     */
    virtual bool isMeshLoaded() const = 0;
    
    /**
     * @brief 获取求解器参数
     */
    virtual std::unordered_map<std::string, double> getParameters() const = 0;
    
    /**
     * @brief 获取边界条件管理器
     */
    virtual std::shared_ptr<BoundaryConditionManager> getBoundaryConditionManager() const = 0;
    
    /**
     * @brief 获取网格对象
     */
    virtual std::shared_ptr<Mesh> getMesh() const = 0;
    
    /**
     * @brief 获取系统矩阵
     */
    virtual std::shared_ptr<Matrix> getSystemMatrix() const = 0;
    
    /**
     * @brief 获取右端向量
     */
    virtual std::vector<double> getRHSVector() const = 0;
    
    /**
     * @brief 获取求解器状态信息
     */
    virtual std::string getStatus() const = 0;
    
    /**
     * @brief 获取求解器统计信息
     */
    virtual std::unordered_map<std::string, double> getStatistics() const = 0;
    
    /**
     * @brief 重置求解器状态
     */
    virtual void reset() = 0;
    
    /**
     * @brief 清理求解器资源
     */
    virtual void cleanup() = 0;

protected:
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<BoundaryConditionManager> bcManager_;
    std::unordered_map<std::string, double> parameters_;
    bool initialized_ = false;
    bool meshLoaded_ = false;
    int timeStepIndex_ = 0;
    double currentTime_ = 0.0;
};

/**
 * @brief 求解器管理器
 * 
 * 管理多个求解器的注册、创建和执行
 */
class SolverManager {
public:
    SolverManager() = default;
    
    /**
     * @brief 注册求解器工厂
     */
    void registerSolverFactory(const std::string& solverType, 
                              std::function<std::shared_ptr<ElmerSolver>()> factory);
    
    /**
     * @brief 创建求解器实例
     */
    std::shared_ptr<ElmerSolver> createSolver(const std::string& solverType);
    
    /**
     * @brief 获取已注册的求解器类型
     */
    std::vector<std::string> getRegisteredSolverTypes() const;
    
    /**
     * @brief 检查求解器类型是否已注册
     */
    bool isSolverTypeRegistered(const std::string& solverType) const;
    
    /**
     * @brief 执行多物理场耦合求解
     */
    bool solveCoupledSystem(const std::vector<std::shared_ptr<ElmerSolver>>& solvers);
    
    /**
     * @brief 执行顺序求解
     */
    bool solveSequentialSystem(const std::vector<std::shared_ptr<ElmerSolver>>& solvers);

private:
    std::unordered_map<std::string, std::function<std::shared_ptr<ElmerSolver>()>> solverFactories_;
};

} // namespace elmer