/**
 * @file HeatSolver.h
 * @brief Elmer FEM热传导求解器
 * 
 * 实现热传导方程的有限元求解，支持稳态和瞬态热分析
 */

#pragma once

#include "../../core/base/SolverBase.h"
#include "../core/base/Types.h"
#include <memory>
#include <vector>

namespace elmer {

/**
 * @brief 热传导求解器参数结构�?
 */
struct HeatSolverParameters {
    double thermalConductivity = 1.0;        ///< 热导�?[W/(m·K)]
    double density = 1.0;                    ///< 密度 [kg/m³]
    double specificHeat = 1.0;               ///< 比热�?[J/(kg·K)]
    double heatSource = 0.0;                 ///< 热源�?[W/m³]
    double initialTemperature = 293.15;      ///< 初始温度 [K]
    double ambientTemperature = 293.15;      ///< 环境温度 [K]
    double heatTransferCoefficient = 0.0;    ///< 热传导系�?[W/(m²·K)]
    
    // 边界条件类型
    enum BoundaryType {
        DIRICHLET,      ///< 狄利克雷边界条件（固定温度）
        NEUMANN,        ///< 诺伊曼边界条件（热通量�?
        ROBIN           ///< 罗宾边界条件（对流换热）
    };
    
    HeatSolverParameters() = default;
};

/**
 * @brief 热传导求解器�?
 * 
 * 实现热传导方程的有限元求解，支持稳态和瞬态分�?
 */
class HeatSolver : public LinearSolverBase {
private:
    HeatSolverParameters heatParams_;        ///< 热传导求解器参数
    std::vector<double> temperatureField_;   ///< 温度�?
    std::vector<double> heatFluxField_;      ///< 热通量�?
    
    // 边界条件数据
    std::vector<int> dirichletNodes_;        ///< 狄利克雷边界节点
    std::vector<double> dirichletValues_;    ///< 狄利克雷边界�?
    std::vector<int> neumannEdges_;          ///< 诺伊曼边界边
    std::vector<double> neumannValues_;      ///< 诺伊曼边界�?
    std::vector<int> robinEdges_;            ///< 罗宾边界�?
    std::vector<double> robinCoefficients_;  ///< 罗宾边界系数
    std::vector<double> robinAmbientTemps_;  ///< 罗宾边界环境温度
    
    // 瞬态分析相�?
    std::vector<double> prevTemperature_;    ///< 上一时间步温度场
    double timeIntegrationFactor_ = 1.0;     ///< 时间积分因子
    
public:
    /**
     * @brief 构造函�?
     */
    HeatSolver();
    
    /**
     * @brief 析构函数
     */
    virtual ~HeatSolver() = default;
    
    /**
     * @brief 设置热传导求解器参数
     */
    void setHeatParameters(const HeatSolverParameters& params);
    
    /**
     * @brief 获取热传导求解器参数
     */
    HeatSolverParameters getHeatParameters() const;
    
    /**
     * @brief 设置狄利克雷边界条件
     */
    void setDirichletBoundary(const std::vector<int>& nodes, const std::vector<double>& values);
    
    /**
     * @brief 设置诺伊曼边界条�?
     */
    void setNeumannBoundary(const std::vector<int>& edges, const std::vector<double>& values);
    
    /**
     * @brief 设置罗宾边界条件
     */
    void setRobinBoundary(const std::vector<int>& edges, const std::vector<double>& coefficients, 
                         const std::vector<double>& ambientTemps);
    
    /**
     * @brief 初始化求解器
     */
    bool initialize() override;
    
    /**
     * @brief 组装系统矩阵
     */
    bool assemble() override;
    
    /**
     * @brief 求解系统
     */
    bool solve() override;
    
    /**
     * @brief 获取求解结果（温度场�?
     */
    std::vector<double> getSolution() const override;
    
    /**
     * @brief 获取热通量�?
     */
    std::vector<double> getHeatFlux() const;
    
    /**
     * @brief 获取最大温�?
     */
    double getMaxTemperature() const;
    
    /**
     * @brief 获取最小温�?
     */
    double getMinTemperature() const;
    
    /**
     * @brief 获取平均温度
     */
    double getAverageTemperature() const;
    
    /**
     * @brief 检查收敛�?
     */
    bool checkConvergence() const;
    
    /**
     * @brief 获取残差
     */
    double getResidual() const;
    
    /**
     * @brief 支持瞬态计�?
     */
    bool supportsTransient() const override { return true; }
    
    /**
     * @brief 执行时间步进
     */
    bool executeTimeStep(int timeStepIndex, double currentTime);
    
    /**
     * @brief 保存求解器状�?
     */
    bool saveState(const std::string& filename) const;
    
    /**
     * @brief 加载求解器状�?
     */
    bool loadState(const std::string& filename);
    
private:
    /**
     * @brief 组装刚度矩阵
     */
    bool assembleStiffnessMatrix();
    
    /**
     * @brief 组装质量矩阵（用于瞬态分析）
     */
    bool assembleMassMatrix();
    
    /**
     * @brief 组装右端向量
     */
    bool assembleRhsVector();
    
    /**
     * @brief 应用边界条件
     */
    bool applyBoundaryConditions();
    
    /**
     * @brief 计算热通量
     */
    void computeHeatFlux();
    
    /**
     * @brief 计算单元热传导矩�?
     */
    void computeElementMatrix(int elementId, std::vector<std::vector<double>>& elementMatrix) const;
    
    /**
     * @brief 计算单元质量矩阵
     */
    void computeElementMassMatrix(int elementId, std::vector<std::vector<double>>& elementMatrix) const;
    
    /**
     * @brief 计算单元右端向量
     */
    void computeElementRhsVector(int elementId, std::vector<double>& elementVector) const;
    
    /**
     * @brief 计算形函数和导数
     */
    void computeShapeFunctions(double xi, double eta, std::vector<double>& N, 
                              std::vector<std::vector<double>>& dN) const;
};

} // namespace elmer

