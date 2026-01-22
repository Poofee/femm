/**
 * @file LinearSolverFactory.h
 * @brief 线性代数求解器工厂类
 * 
 * 提供稀疏矩阵线性代数求解器的创建和管理功能。
 */

#pragma once

#include "LinearSolverInterface.h"
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <map>

namespace elmer {

/**
 * @brief 线性代数求解器工厂类
 * 
 * 专门用于创建和管理稀疏矩阵线性代数求解器（SuperLU、MUMPS等）
 */
class LinearSolverFactory {
public:
    /**
     * @brief 创建求解器实例
     * @param type 求解器类型
     * @return 求解器实例
     */
    static std::unique_ptr<LinearSolver> CreateSolver(SolverType type);
    
    /**
     * @brief 创建求解器实例（带参数）
     * @param type 求解器类型
     * @param params 求解器参数
     * @return 求解器实例
     */
    static std::unique_ptr<LinearSolver> CreateSolver(SolverType type, 
                                                     const SolverParameters& params);
    
    /**
     * @brief 创建求解器实例（通过名称）
     * @param solver_name 求解器名称
     * @return 求解器实例
     */
    static std::unique_ptr<LinearSolver> CreateSolver(const std::string& solver_name);
    
    /**
     * @brief 检查求解器是否可用
     * @param type 求解器类型
     * @return 是否可用
     */
    static bool IsSolverAvailable(SolverType type);
    
    /**
     * @brief 检查求解器是否可用（通过名称）
     * @param solver_name 求解器名称
     * @return 是否可用
     */
    static bool IsSolverAvailable(const std::string& solver_name);
    
    /**
     * @brief 获取所有可用求解器类型
     * @return 求解器类型列表
     */
    static std::vector<SolverType> GetAvailableSolverTypes();
    
    /**
     * @brief 获取所有可用求解器名称
     * @return 求解器名称列表
     */
    static std::vector<std::string> GetAvailableSolverNames();
    
    /**
     * @brief 获取求解器描述信息
     * @param type 求解器类型
     * @return 描述信息
     */
    static std::string GetSolverDescription(SolverType type);
    
    /**
     * @brief 获取默认求解器参数
     * @param type 求解器类型
     * @return 默认参数
     */
    static SolverParameters GetDefaultParameters(SolverType type);
    
    /**
     * @brief 获取推荐求解器参数
     * @param type 求解器类型
     * @return 推荐参数
     */
    static SolverParameters GetRecommendedParameters(SolverType type);
};

} // namespace elmer