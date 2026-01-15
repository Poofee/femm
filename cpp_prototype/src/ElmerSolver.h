/**
 * @file ElmerSolver.h
 * @brief Elmer FEM主求解程序
 * 
 * 移植自Fortran版本的ElmerSolver.F90，提供完整的求解流程控制
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include "SolverBase.h"
#include "Mesh.h"
#include "MaterialDatabase.h"
#include "BoundaryConditions.h"
#include "MPICommunicator.h"
#include "MPIUtils.h"

// 前向声明
namespace elmer {
    class SolverManager;
}

namespace elmer {

/**
 * @brief 仿真类型枚举
 */
enum class SimulationType {
    STEADY_STATE,   ///< 稳态仿真
    TRANSIENT,      ///< 瞬态仿真
    SCANNING,       ///< 参数扫描
    OPTIMIZATION    ///< 优化
};

/**
 * @brief 仿真参数结构体
 */
struct SimulationParameters {
    SimulationType type = SimulationType::STEADY_STATE; ///< 仿真类型
    std::string modelName;                              ///< 模型文件名
    std::string meshDir;                                ///< 网格目录
    std::string meshName;                               ///< 网格文件名
    
    // 时间步进参数
    double startTime = 0.0;                             ///< 开始时间
    double endTime = 1.0;                               ///< 结束时间
    double timeStep = 0.1;                              ///< 时间步长
    int numTimeSteps = 10;                              ///< 时间步数
    
    // 输出控制
    bool verbose = true;                                ///< 详细输出
    int outputInterval = 1;                             ///< 输出间隔
    std::string outputFormat = "vtk";                  ///< 输出格式
    
    // 并行参数
    bool useMPI = false;                                ///< 是否使用MPI
    bool useOpenMP = false;                             ///< 是否使用OpenMP
    int numThreads = 1;                                 ///< 线程数
    
    // 优化参数
    int maxOptimizationIterations = 100;                ///< 最大优化迭代次数
    double optimizationTolerance = 1.0e-6;              ///< 优化容差
    
    SimulationParameters() = default;
};

/**
 * @brief 仿真结果结构体
 */
struct SimulationResult {
    bool success = false;                               ///< 是否成功
    double cpuTime = 0.0;                               ///< CPU时间
    double realTime = 0.0;                              ///< 实际时间
    int numIterations = 0;                              ///< 迭代次数
    double finalResidual = 0.0;                         ///< 最终残差
    std::string errorMessage;                           ///< 错误信息
    
    // 求解器结果
    std::map<std::string, std::vector<double>> solutions; ///< 各求解器解
    
    SimulationResult() = default;
};

/**
 * @brief Elmer FEM主求解器类
 * 
 * 移植自Fortran版本的ElmerSolver.F90，提供完整的求解流程控制
 */
class ElmerSolver {
private:
    SimulationParameters parameters_;                   ///< 仿真参数
    SimulationResult result_;                           ///< 仿真结果
    
    std::shared_ptr<Mesh> mesh_;                        ///< 网格数据
    std::shared_ptr<MaterialDatabase> materialDB_;      ///< 材料数据库
    std::shared_ptr<BoundaryConditionManager> bc_;      ///< 边界条件管理器
    std::shared_ptr<MPICommunicator> comm_;             ///< MPI通信器
    
    std::unique_ptr<SolverManager> solverManager_;      ///< 求解器管理器
    
    // 时间控制
    std::chrono::steady_clock::time_point startTime_;   ///< 开始时间
    std::chrono::steady_clock::time_point endTime_;     ///< 结束时间
    
    // 状态标志
    bool initialized_ = false;                          ///< 是否已初始化
    bool meshLoaded_ = false;                           ///< 网格是否已加载
    
public:
    /**
     * @brief 构造函数
     */
    ElmerSolver();
    
    /**
     * @brief 析构函数
     */
    ~ElmerSolver();
    
    /**
     * @brief 设置仿真参数
     */
    void setParameters(const SimulationParameters& params);
    
    /**
     * @brief 获取仿真参数
     */
    SimulationParameters getParameters() const;
    
    /**
     * @brief 初始化求解器
     */
    bool initialize();
    
    /**
     * @brief 加载模型和网格
     */
    bool loadModel();
    
    /**
     * @brief 添加求解器
     */
    void addSolver(const std::string& solverName);
    
    /**
     * @brief 执行稳态仿真
     */
    bool executeSteadyState();
    
    /**
     * @brief 执行瞬态仿真
     */
    bool executeTransient();
    
    /**
     * @brief 执行参数扫描
     */
    bool executeScanning();
    
    /**
     * @brief 执行优化
     */
    bool executeOptimization();
    
    /**
     * @brief 执行仿真
     */
    SimulationResult execute();
    
    /**
     * @brief 获取仿真结果
     */
    SimulationResult getResult() const;
    
    /**
     * @brief 清理资源
     */
    void cleanup();
    
public:
    /**
     * @brief 处理命令行参数
     */
    void processCommandLineArguments(int argc, char** argv);
    
private:
    /**
     * @brief 打印横幅信息
     */
    void printBanner();
    
    /**
     * @brief 初始化并行环境
     */
    bool initializeParallelEnvironment();
    
    /**
     * @brief 初始化OpenMP线程
     */
    bool initializeOpenMP();
    
    /**
     * @brief 读取输入文件
     */
    bool readInputFile();
    
    /**
     * @brief 设置初始条件
     */
    bool setInitialConditions();
    
    /**
     * @brief 执行时间步进
     */
    bool executeTimeStep(int timeStepIndex, double currentTime);
    
    /**
     * @brief 保存结果
     */
    bool saveResults(int timeStepIndex, double currentTime);
    
    /**
     * @brief 检查收敛性
     */
    bool checkConvergence();
    
    /**
     * @brief 获取CPU时间
     */
    double getCPUTime() const;
    
    /**
     * @brief 获取实际时间
     */
    double getRealTime() const;
};

/**
 * @brief ElmerSolver主函数
 * 
 * 对应Fortran版本的ElmerSolver主程序入口
 */
int ElmerSolverMain(int argc, char** argv);

} // namespace elmer