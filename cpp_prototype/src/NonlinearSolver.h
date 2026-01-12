#pragma once

#include <memory>
#include <vector>
#include <cmath>
#include <iostream>
#include <functional>
#include <chrono>
#include "LinearAlgebra.h"
#include "ElectromagneticMaterial.h"
#include "Mesh.h"

namespace elmer {

/**
 * @brief 非线性求解器配置参数
 */
struct NonlinearSolverParameters {
    double tolerance = 1.0e-8;        ///< 收敛容差
    int maxIterations = 50;           ///< 最大迭代次数
    bool useNewtonRaphson = true;     ///< 使用牛顿-拉夫逊方法
    bool useLineSearch = true;        ///< 使用线搜索
    double relaxationFactor = 1.0;    ///< 松弛因子
    
    // 收敛标准
    double residualTolerance = 1.0e-8;    ///< 残差容差
    double solutionTolerance = 1.0e-8;    ///< 解变化容差
    
    // 输出控制
    bool verbose = true;              ///< 详细输出
    int printFrequency = 5;           ///< 输出频率
    
    NonlinearSolverParameters() = default;
};

/**
 * @brief 非线性求解器结果
 */
struct NonlinearSolverResults {
    bool converged = false;           ///< 是否收敛
    int iterations = 0;               ///< 迭代次数
    double residualNorm = 0.0;        ///< 最终残差范数
    double solutionNorm = 0.0;        ///< 解范数
    std::vector<double> residuals;    ///< 迭代残差历史
    std::vector<double> corrections;  ///< 修正量历史
    
    // 性能统计
    double assemblyTime = 0.0;        ///< 矩阵组装时间(ms)
    double solveTime = 0.0;           ///< 线性求解时间(ms)
    double totalTime = 0.0;           ///< 总时间(ms)
    
    NonlinearSolverResults() = default;
};

/**
 * @brief 非线性求解器基类
 */
class NonlinearSolver {
protected:
    NonlinearSolverParameters parameters;
    
public:
    virtual ~NonlinearSolver() = default;
    
    // 设置参数
    void setParameters(const NonlinearSolverParameters& params) {
        parameters = params;
    }
    
    // 获取参数
    const NonlinearSolverParameters& getParameters() const {
        return parameters;
    }
    
    /**
     * @brief 求解非线性系统
     * 
     * @param initialGuess 初始猜测
     * @param residualFunction 残差函数
     * @param jacobianFunction 雅可比函数
     * @return NonlinearSolverResults 求解结果
     */
    virtual NonlinearSolverResults solve(
        const std::shared_ptr<Vector>& initialGuess,
        std::function<std::shared_ptr<Vector>(const std::shared_ptr<Vector>&)> residualFunction,
        std::function<std::shared_ptr<Matrix>(const std::shared_ptr<Vector>&)> jacobianFunction) = 0;
    
    /**
     * @brief 检查收敛性
     */
    virtual bool checkConvergence(double residualNorm, double correctionNorm, 
                                 double initialResidual, int iteration) const {
        // 残差收敛标准
        if (residualNorm < parameters.residualTolerance) {
            if (parameters.verbose) {
                std::cout << "收敛: 残差范数 " << residualNorm 
                         << " < " << parameters.residualTolerance << std::endl;
            }
            return true;
        }
        
        // 解变化收敛标准
        if (correctionNorm < parameters.solutionTolerance) {
            if (parameters.verbose) {
                std::cout << "收敛: 解变化范数 " << correctionNorm 
                         << " < " << parameters.solutionTolerance << std::endl;
            }
            return true;
        }
        
        // 最大迭代次数
        if (iteration >= parameters.maxIterations) {
            if (parameters.verbose) {
                std::cout << "未收敛: 达到最大迭代次数 " << parameters.maxIterations << std::endl;
            }
            return false;
        }
        
        return false;
    }
};

/**
 * @brief 牛顿-拉夫逊非线性求解器
 */
class NewtonRaphsonSolver : public NonlinearSolver {
private:
    // 计算向量范数的辅助函数
    double vectorNorm(const std::shared_ptr<Vector>& v) {
        double norm = 0.0;
        for (Integer i = 0; i < v->Size(); ++i) {
            norm += (*v)[i] * (*v)[i];
        }
        return std::sqrt(norm);
    }
    
    // 线搜索方法
    double lineSearch(const std::shared_ptr<Vector>& x, 
                     const std::shared_ptr<Vector>& direction,
                     std::function<std::shared_ptr<Vector>(const std::shared_ptr<Vector>&)> residualFunction) {
        if (!parameters.useLineSearch) {
            return parameters.relaxationFactor;
        }
        
        double alpha = parameters.relaxationFactor;
        double initialResidualNorm = vectorNorm(residualFunction(x));
        
        // 简单回溯线搜索
        for (int i = 0; i < 10; ++i) {
            auto trialX = std::shared_ptr<Vector>(Vector::Create(x->Size()));
            for (int j = 0; j < x->Size(); ++j) {
                (*trialX)[j] = (*x)[j] + alpha * (*direction)[j];
            }
            
            double trialResidualNorm = vectorNorm(residualFunction(trialX));
            
            // 如果残差减小，接受该步长
            if (trialResidualNorm < initialResidualNorm) {
                return alpha;
            }
            
            // 否则减小步长
            alpha *= 0.5;
        }
        
        return parameters.relaxationFactor; // 返回默认值
    }
    
public:
    NonlinearSolverResults solve(
        const std::shared_ptr<Vector>& initialGuess,
        std::function<std::shared_ptr<Vector>(const std::shared_ptr<Vector>&)> residualFunction,
        std::function<std::shared_ptr<Matrix>(const std::shared_ptr<Vector>&)> jacobianFunction) override {
        
        NonlinearSolverResults results;
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // 初始化解向量
        auto x = std::shared_ptr<Vector>(Vector::Create(initialGuess->Size()));
        for (int i = 0; i < initialGuess->Size(); ++i) {
            (*x)[i] = (*initialGuess)[i];
        }
        
        // 计算初始残差
        auto residual = residualFunction(x);
        double residualNorm = vectorNorm(residual);
        double initialResidualNorm = residualNorm;
        
        if (parameters.verbose) {
            std::cout << "牛顿-拉夫逊求解器开始" << std::endl;
            std::cout << "初始残差范数: " << residualNorm << std::endl;
        }
        
        results.residuals.push_back(residualNorm);
        
        // 牛顿迭代
        for (int iteration = 0; iteration < parameters.maxIterations; ++iteration) {
            auto iterationStart = std::chrono::high_resolution_clock::now();
            
            // 计算雅可比矩阵
            auto jacobian = jacobianFunction(x);
            auto assemblyTime = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - iterationStart).count() / 1000.0;
            
            // 求解线性系统: J * dx = -r
            auto rhs = std::shared_ptr<Vector>(Vector::Create(residual->Size()));
            for (size_t i = 0; i < residual->Size(); ++i) {
                (*rhs)[i] = -(*residual)[i];
            }
            
            auto solveStart = std::chrono::high_resolution_clock::now();
            auto dx = std::shared_ptr<Vector>(Vector::Create(rhs->Size()));
            // 这里需要实现线性求解器，暂时使用简单的对角求解
            for (int i = 0; i < rhs->Size(); ++i) {
                (*dx)[i] = (*rhs)[i] / jacobian->GetElement(i, i);
            }
            auto solveTime = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - solveStart).count() / 1000.0;
            
            // 线搜索
            double alpha = lineSearch(x, dx, residualFunction);
            
            // 更新解
            double correctionNorm = 0.0;
            for (int i = 0; i < x->Size(); ++i) {
                double correction = alpha * (*dx)[i];
                (*x)[i] += correction;
                correctionNorm += correction * correction;
            }
            correctionNorm = std::sqrt(correctionNorm);
            
            // 计算新残差
            residual = residualFunction(x);
            residualNorm = vectorNorm(residual);
            
            results.residuals.push_back(residualNorm);
            results.corrections.push_back(correctionNorm);
            results.assemblyTime += assemblyTime;
            results.solveTime += solveTime;
            
            // 输出迭代信息
            if (parameters.verbose && (iteration % parameters.printFrequency == 0 || iteration == 0)) {
                std::cout << "迭代 " << iteration << ": 残差=" << residualNorm 
                         << ", 修正量=" << correctionNorm 
                         << ", 步长=" << alpha << std::endl;
            }
            
            // 检查收敛性
            if (checkConvergence(residualNorm, correctionNorm, initialResidualNorm, iteration)) {
                results.converged = true;
                results.iterations = iteration + 1;
                results.residualNorm = residualNorm;
                results.solutionNorm = vectorNorm(x);
                break;
            }
        }
        
        auto totalTime = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;
        
        results.totalTime = totalTime;
        
        if (parameters.verbose) {
            std::cout << "求解完成: " << (results.converged ? "收敛" : "未收敛") 
                     << ", 迭代次数=" << results.iterations 
                     << ", 总时间=" << totalTime << " ms" << std::endl;
        }
        
        return results;
    }
};

/**
 * @brief 非线性材料更新器
 */
class NonlinearMaterialUpdater {
private:
    NonlinearMaterialDatabase& materialDB;
    
public:
    NonlinearMaterialUpdater(NonlinearMaterialDatabase& db) : materialDB(db) {}
    
    /**
     * @brief 更新材料参数基于当前磁场强度
     */
    void updateMaterialParameters(const std::vector<double>& H_values, 
                                 const std::vector<std::string>& materialNames) {
        for (size_t i = 0; i < materialNames.size(); ++i) {
            if (materialDB.hasMaterial(materialNames[i])) {
                auto& material = materialDB.getMaterial(materialNames[i]);
                if (material.isNonlinear) {
                    // 基于当前H值更新材料参数
                    // 这里可以添加更复杂的更新逻辑
                    double H_effective = std::abs(H_values[i]);
                    // 材料参数会在求解过程中动态更新
                }
            }
        }
    }
    
    /**
     * @brief 计算磁场强度从向量势
     */
    std::vector<double> computeHFromPotential(const std::shared_ptr<Vector>& potential, 
                                             const std::shared_ptr<Mesh>& mesh) {
        std::vector<double> H_values;
        // 实现磁场强度计算逻辑
        // 这里需要根据具体的有限元公式实现
        return H_values;
    }
};

} // namespace elmer