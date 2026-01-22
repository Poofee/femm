#pragma once

#include "LinearAlgebra.h"
#include "Mesh.h"
#include "Types.h"
#include "SolverRegistry.h"
#include "BoundaryConditions.h"
#include "InputFileParser.h"
#include "LoggerInterface.h"
#include <algorithm>
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace elmer {

/**
 * @brief 仿真参数结构体
 */
struct SimulationParameters {
    double startTime = 0.0;          // 开始时间
    double endTime = 1.0;            // 结束时间
    double timeStep = 0.01;          // 时间步长
    std::string meshName;            // 网格文件名
    std::string meshDir;             // 网格目录
    std::string outputDir;           // 输出目录
    int maxIterations = 100;         // 最大迭代次数
    double tolerance = 1e-6;         // 收敛容差
    bool transient = false;          // 瞬态仿真
    bool nonlinear = false;          // 非线性求解
    int numThreads = 1;              // 线程数
    int numTimeSteps = 100;          // 时间步数
    int outputInterval = 10;         // 输出间隔
    bool verbose = true;             // 详细输出控制
    bool checkOnly = false;          // 只检查模式
};

/**
 * @brief 仿真结果结构体
 */
struct SimulationResult {
    bool success = false;            // 仿真是否成功
    double executionTime = 0.0;      // 执行时间
    int iterations = 0;              // 迭代次数
    double finalResidual = 0.0;      // 最终残差
    std::string message;             // 结果消息
    std::vector<double> solution;    // 解向量
};

/**
 * @brief Elmer求解器基类
 * 
 * 提供求解器的基本功能，包括网格管理、边界条件处理和求解过程
 */
class ElmerSolver {
public:
    ElmerSolver();
    virtual ~ElmerSolver();
    
    /**
     * @brief 初始化求解器
     */
    virtual bool initialize();
    
    /**
     * @brief 加载网格数据
     */
    virtual bool loadMesh(std::shared_ptr<Mesh> mesh);
    
    /**
     * @brief 添加边界条件
     */
    virtual void addBoundaryCondition(std::shared_ptr<BoundaryCondition> bc);
    
    /**
     * @brief 添加求解器
     */
    virtual void addSolver(const std::string& solverName);
    
    /**
     * @brief 组装系统矩阵
     */
    virtual bool assembleSystem();
    
    /**
     * @brief 求解系统
     */
    virtual bool solve();
    
    /**
     * @brief 验证仿真参数
     */
    virtual bool validateParameters(const SimulationParameters& params);
    
    /**
     * @brief 设置仿真参数
     */
    virtual void setParameters(const SimulationParameters& params);
    
    /**
     * @brief 执行稳态仿真
     */
    virtual bool executeSteadyState();
    
    /**
     * @brief 执行瞬态仿真
     */
    virtual bool executeTransient();
    
    /**
     * @brief 执行时间步进
     */
    virtual bool executeTimeStep(int timeStepIndex, double currentTime);
    
    /**
     * @brief 更新时变边界条件
     */
    virtual bool updateTimeDependentBoundaryConditions(double currentTime);
    
    /**
     * @brief 检查收敛性
     */
    virtual bool checkConvergence();
    
    /**
     * @brief 执行参数扫描仿真
     */
    virtual bool executeScanning();
    
    /**
     * @brief 执行仿真
     */
    virtual SimulationResult execute();
    
    /**
     * @brief 执行优化
     */
    virtual bool executeOptimization();
    
    /**
     * @brief 获取仿真参数
     */
    virtual SimulationParameters getParameters() const;
    
    /**
     * @brief 获取求解结果
     */
    virtual std::vector<double> getSolution() const;
    
    /**
     * @brief 保存结果
     */
    virtual bool saveResults(int step, double currentTime);
    
    /**
     * @brief 保存VTK格式结果
     */
    virtual bool saveResultsVTK(const std::string& filename, int timeStepIndex, double currentTime);
    
    /**
     * @brief 保存Gmsh格式结果
     */
    virtual bool saveResultsGmsh(const std::string& filename, int timeStepIndex, double currentTime);
    
    /**
     * @brief 保存CSV格式结果
     */
    virtual bool saveResultsCSV(const std::string& filename, int timeStepIndex, double currentTime);
    
    /**
     * @brief 获取求解器名称
     */
    virtual std::string getName() const;
    
    /**
     * @brief 获取求解器类型
     */
    virtual std::string getType() const;
    
    /**
     * @brief 检查求解器是否已初始化
     */
    virtual bool isInitialized() const;
    
    /**
     * @brief 检查网格是否已加载
     */
    virtual bool isMeshLoaded() const;

    /**
     * @brief 获取网格对象
     */
    virtual std::shared_ptr<Mesh> getMesh() const;
    
    /**
     * @brief 获取系统矩阵
     */
    virtual std::shared_ptr<Matrix> getSystemMatrix() const;
    
    /**
     * @brief 打印求解器横幅信息
     */
    virtual void printBanner();
    
    /**
     * @brief 处理命令行参数
     */
    virtual void processCommandLineArguments(int argc, char* argv[]);
    
    /**
     * @brief 获取CPU时间
     */
    virtual double getCPUTime() const;
    
    /**
     * @brief 获取实时时间
     */
    virtual double getRealTime() const;
    
    /**
     * @brief 启动计时器
     */
    virtual void startTimer();
    
    /**
     * @brief 停止计时器
     */
    virtual void stopTimer();
    
    /**
     * @brief 打印性能统计信息
     */
    virtual void printPerformanceStats() const;
    
    /**
     * @brief 初始化并行环境
     */
    virtual bool initializeParallelEnvironment();
    
    /**
     * @brief 注册电磁相关求解器
     */
    virtual bool registerElectromagneticSolvers();
    
    /**
     * @brief 初始化OpenMP并行环境
     */
    virtual bool initializeOpenMP();
    
    /**
     * @brief 读取输入文件
     */
    virtual bool readInputFile();
    
    /**
     * @brief 设置初始条件
     */
    virtual bool setInitialConditions();
    
    /**
     * @brief 加载模型
     */
    virtual bool loadModel();
    
    /**
     * @brief 获取右端向量
     */
    virtual std::vector<double> getRHSVector() const;
    
    /**
     * @brief 获取求解器状态信息
     */
    virtual std::string getStatus() const;
    
    /**
     * @brief 获取求解器统计信息
     */
    virtual std::unordered_map<std::string, double> getStatistics() const;
    
    /**
     * @brief 获取仿真结果
     */
    virtual SimulationResult getResult() const;
    
    /**
     * @brief 重置求解器状态
     */
    virtual void reset();
    
    /**
     * @brief 清理求解器资源
     */
    virtual void cleanup();

    /**
     * @brief 执行主程序循环
     */
    virtual bool executeMainLoop();
    
    /**
     * @brief Elmer求解器主程序
     */
    virtual bool elmerSolverMain(int initialize = 0);

protected:
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<BoundaryConditionManager> bcManager_;
    std::unique_ptr<SolverManager> solverManager_;
    std::unique_ptr<InputFileParser> inputParser_;
    std::shared_ptr<MaterialDatabase> materialDB_;
    std::shared_ptr<BoundaryConditionManager> bc_;
    SimulationParameters parameters_;
    SimulationResult result_;                          // 仿真结果
    bool initialized_ = false;
    bool meshLoaded_ = false;
    bool inputFileParsed_ = false;
    int timeStepIndex_ = 0;
    double currentTime_ = 0.0;
    double startCPUTime_ = 0.0;                        // 开始CPU时间
    double startRealTime_ = 0.0;                       // 开始实时时间
    double endCPUTime_ = 0.0;                          // 结束CPU时间
    double endRealTime_ = 0.0;                         // 结束实时时间
};


} // namespace elmer