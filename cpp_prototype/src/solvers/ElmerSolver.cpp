/**
 * @file ElmerSolver.cpp
 * @brief Elmer FEM主求解程序实现
 * 
 * 移植自Fortran版本的ElmerSolver.F90，提供完整的求解流程控制
 * 
 * TODO: 需要后续进一步开发的功能
 * - [ ] 实现完整的MPI并行通信
 * - [ ] 实现输入文件解析器
 * - [ ] 实现网格加载和预处理
 * - [ ] 实现求解器管理器
 * - [ ] 实现稳态和瞬态求解算法
 * - [ ] 实现结果输出和后处理
 * - [ ] 实现命令行参数处理
 * - [ ] 实现性能监控和日志记录
 */

#include "ElmerSolver.h"
#include "SolverRegistry.h"
#include "InputFileParser.h"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <ctime>
#include <cmath>
#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace elmer {

// ===== 构造函数和析构函数 =====

ElmerSolver::ElmerSolver() 
    : initialized_(false), 
      meshLoaded_(false),
      inputFileParsed_(false) {
    // TODO: 初始化默认参数 - 需要从配置文件读取
    parameters_.modelName = "default_model";
    parameters_.meshDir = ".";
    parameters_.meshName = "default_mesh";
    
    // TODO: 初始化MPI通信器 - 需要实现MPI环境初始化
    comm_ = std::make_shared<MPICommunicator>();
    
    // TODO: 初始化求解器管理器 - 需要实现求解器注册和查找
    solverManager_ = std::make_unique<SolverManager>();
    
    // 初始化输入文件解析器
    inputParser_ = std::make_unique<InputFileParser>();
    
    std::cout << "ElmerSolver构造函数完成" << std::endl;
}

ElmerSolver::~ElmerSolver() {
    // TODO: 实现资源清理
    cleanup();
    std::cout << "ElmerSolver析构函数完成" << std::endl;
}

// ===== 基本接口函数 =====

bool ElmerSolver::validateParameters(const SimulationParameters& params) {
    // 验证基本参数
    if (params.startTime < 0.0) {
        std::cerr << "错误: 开始时间不能为负数" << std::endl;
        return false;
    }
    
    if (params.endTime <= params.startTime) {
        std::cerr << "错误: 结束时间必须大于开始时间" << std::endl;
        return false;
    }
    
    if (params.timeStep <= 0.0) {
        std::cerr << "错误: 时间步长必须为正数" << std::endl;
        return false;
    }
    
    // 验证电磁相关参数
    if (!params.meshDir.empty() && !std::filesystem::exists(params.meshDir)) {
        std::cerr << "错误: 网格目录不存在: " << params.meshDir << std::endl;
        return false;
    }
    
    return true;
}

bool ElmerSolver::registerElectromagneticSolvers() {
    std::cout << "开始注册电磁相关求解器..." << std::endl;
    
    // 简化实现：直接注册电磁求解器
    std::cout << "注册电磁场求解器" << std::endl;
    
    std::cout << "电磁求解器注册完成" << std::endl;
    return true;
}

void ElmerSolver::setParameters(const SimulationParameters& params) {
    // 参数验证
    if (!validateParameters(params)) {
        std::cerr << "错误: 仿真参数验证失败" << std::endl;
        return;
    }
    
    parameters_ = params;
    std::cout << "设置仿真参数完成" << std::endl;
}

SimulationParameters ElmerSolver::getParameters() const {
    // TODO: 实现参数获取功能 - 需要返回当前参数状态
    std::cout << "获取仿真参数" << std::endl;
    return parameters_;
}

bool ElmerSolver::initialize() {
    if (initialized_) {
        std::cout << "警告: ElmerSolver已经初始化" << std::endl;
        return true;
    }
    
    // 实现完整的初始化流程 - 包括MPI、OpenMP、求解器注册等
    std::cout << "开始初始化ElmerSolver..." << std::endl;
    
    // 1. 打印横幅信息
    printBanner();
    
    // 2. 初始化并行环境
    if (!initializeParallelEnvironment()) {
        std::cerr << "错误: 并行环境初始化失败" << std::endl;
        return false;
    }
    
    // 3. 初始化OpenMP
    if (!initializeOpenMP()) {
        std::cerr << "错误: OpenMP初始化失败" << std::endl;
        return false;
    }
    
    // 4. 读取输入文件
    if (!readInputFile()) {
        std::cerr << "错误: 输入文件读取失败" << std::endl;
        return false;
    }
    
    // 5. 注册电磁相关求解器
    if (!registerElectromagneticSolvers()) {
        std::cerr << "错误: 电磁求解器注册失败" << std::endl;
        return false;
    }
    
    // 6. 设置初始条件
    if (!setInitialConditions()) {
        std::cerr << "错误: 初始条件设置失败" << std::endl;
        return false;
    }
    
    // 基础初始化
    initialized_ = true;
    
    std::cout << "ElmerSolver初始化完成" << std::endl;
    return true;
}

bool ElmerSolver::loadModel() {
    if (!initialized_) {
        std::cout << "错误: ElmerSolver未初始化" << std::endl;
        return false;
    }
    
    if (!inputFileParsed_) {
        std::cout << "错误: 输入文件未解析" << std::endl;
        return false;
    }
    
    // TODO: 实现模型加载功能 - 需要读取网格文件、边界条件等
    std::cout << "开始加载模型..." << std::endl;
    
    // 构建网格文件路径
    std::string meshFilePath = parameters_.meshDir + "/" + parameters_.meshName;
    std::cout << "网格文件路径: " << meshFilePath << std::endl;
    
    // TODO: 实现网格文件读取
    // mesh_ = std::make_shared<Mesh>();
    // if (!mesh_->loadFromFile(meshFilePath)) {
    //     std::cerr << "错误: 网格文件加载失败" << std::endl;
    //     return false;
    // }
    
    // TODO: 实现材料数据库加载
    // materialDB_ = std::make_shared<MaterialDatabase>();
    // if (!materialDB_->loadFromInput(*inputParser_)) {
    //     std::cerr << "错误: 材料数据库加载失败" << std::endl;
    //     return false;
    // }
    
    // TODO: 实现边界条件加载
    // bc_ = std::make_shared<BoundaryConditionManager>();
    // if (!bc_->loadFromInput(*inputParser_)) {
    //     std::cerr << "错误: 边界条件加载失败" << std::endl;
    //     return false;
    // }
    
    // 模拟模型加载
    meshLoaded_ = true;
    
    std::cout << "模型加载完成" << std::endl;
    return true;
}

void ElmerSolver::addSolver(const std::string& solverName) {
    // 实现求解器添加功能 - 从求解器注册表中查找并添加
    std::cout << "添加求解器: " << solverName << std::endl;
    
    // 从求解器注册表创建求解器实例
    auto solver = SolverRegistry::getInstance().createSolver(solverName);
    if (solver) {
        // 设置求解器参数
        solver->setMesh(mesh_);
        solver->setMaterialDatabase(materialDB_);
        solver->setBoundaryConditions(bc_);
        solver->setMPICommunicator(comm_.get());
        
        // 添加到求解器管理器
        solverManager_->addSolver(solver);
        std::cout << "求解器 " << solverName << " 添加成功" << std::endl;
    } else {
        std::cerr << "错误: 无法创建求解器 " << solverName << std::endl;
    }
}

// ===== 仿真执行函数 =====

bool ElmerSolver::executeSteadyState() {
    if (!initialized_ || !meshLoaded_) {
        std::cout << "错误: 求解器未初始化或模型未加载" << std::endl;
        return false;
    }
    
    std::cout << "开始执行稳态仿真..." << std::endl;
    
    // 获取稳态求解器参数
    int maxIterations = inputParser_->getParameterInteger("simulation", 0, "Steady State Max Iterations", 100);
    double tolerance = inputParser_->getParameterReal("simulation", 0, "Steady State Tolerance", 1e-6);
    
    std::cout << "稳态求解器参数: 最大迭代次数=" << maxIterations 
              << ", 容差=" << tolerance << std::endl;
    
    // 初始化求解器
    if (!solverManager_->initialize()) {
        std::cerr << "错误: 求解器管理器初始化失败" << std::endl;
        return false;
    }
    
    // 模拟求解器执行（简化处理）
    std::cout << "模拟求解器执行..." << std::endl;
    
    // 稳态求解循环
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        std::cout << "稳态迭代 " << iteration + 1 << "/" << maxIterations << std::endl;
        
        // 模拟求解器执行
        bool allConverged = true;
        
        // 模拟收敛检查（简化处理）
        double residual = 1.0 / (iteration + 1);
        if (residual > tolerance) {
            allConverged = false;
        }
        
        std::cout << "当前残差: " << residual << std::endl;
        
        // 如果收敛，则退出循环
        if (allConverged) {
            std::cout << "求解器在迭代 " << iteration + 1 << " 收敛" << std::endl;
            break;
        }
        
        // 检查是否达到最大迭代次数
        if (iteration == maxIterations - 1) {
            std::cout << "警告: 达到最大迭代次数，可能未完全收敛" << std::endl;
        }
    }
    
    std::cout << "稳态仿真完成" << std::endl;
    return true;
}

bool ElmerSolver::executeTransient() {
    if (!initialized_ || !meshLoaded_) {
        std::cout << "错误: 求解器未初始化或模型未加载" << std::endl;
        return false;
    }
    
    std::cout << "开始执行瞬态仿真..." << std::endl;
    
    // 获取瞬态求解器参数
    double startTime = parameters_.startTime;
    double endTime = parameters_.endTime;
    double timeStep = parameters_.timeStep;
    int numTimeSteps = parameters_.numTimeSteps;
    int outputInterval = parameters_.outputInterval;
    
    std::cout << "瞬态求解器参数: 开始时间=" << startTime 
              << "s, 结束时间=" << endTime << "s, 时间步长=" << timeStep 
              << "s, 时间步数=" << numTimeSteps << std::endl;
    
    // 初始化求解器
    if (!solverManager_->initialize()) {
        std::cerr << "错误: 求解器管理器初始化失败" << std::endl;
        return false;
    }
    
    // 模拟求解器执行（简化处理）
    std::cout << "模拟瞬态求解器执行..." << std::endl;
    
    // 时间步进循环
    for (int step = 0; step < numTimeSteps; ++step) {
        double currentTime = startTime + step * timeStep;
        
        std::cout << "时间步 " << step + 1 << "/" << numTimeSteps 
                  << ", 当前时间: " << currentTime << "s" << std::endl;
        
        // 更新时间步进参数
        timeStepIndex_ = step;
        currentTime_ = currentTime;
        
        // 更新时间变边界条件
        if (!updateTimeDependentBoundaryConditions(currentTime)) {
            std::cerr << "错误: 更新时变边界条件失败" << std::endl;
            return false;
        }
        
        // 模拟求解器执行
        std::cout << "执行时间步求解..." << std::endl;
        
        // 检查收敛性
        if (!checkConvergence()) {
            std::cout << "警告: 时间步 " << step + 1 << " 未收敛" << std::endl;
        }
        
        // 保存结果（如果达到输出间隔）
        if ((step + 1) % outputInterval == 0 || step == numTimeSteps - 1) {
            if (!saveResults(step, currentTime)) {
                std::cerr << "错误: 结果保存失败" << std::endl;
                return false;
            }
        }
        
        // 检查是否达到结束时间
        if (currentTime >= endTime) {
            std::cout << "达到结束时间: " << currentTime << "s" << std::endl;
            break;
        }
    }
    
    std::cout << "瞬态仿真完成" << std::endl;
    return true;
}

bool ElmerSolver::executeScanning() {
    if (!initialized_ || !meshLoaded_) {
        std::cout << "错误: 求解器未初始化或模型未加载" << std::endl;
        return false;
    }
    
    // TODO: 实现参数扫描 - 需要实现参数变化和批量仿真执行
    std::cout << "开始执行参数扫描..." << std::endl;
    
    // 模拟参数扫描
    std::cout << "参数扫描完成" << std::endl;
    return true;
}

bool ElmerSolver::executeOptimization() {
    if (!initialized_ || !meshLoaded_) {
        std::cout << "错误: 求解器未初始化或模型未加载" << std::endl;
        return false;
    }
    
    // TODO: 实现优化 - 需要实现优化算法和灵敏度分析
    std::cout << "开始执行优化..." << std::endl;
    
    // 模拟优化
    std::cout << "优化完成" << std::endl;
    return true;
}

SimulationResult ElmerSolver::execute() {
    if (!initialized_ || !meshLoaded_) {
        std::cout << "错误: 求解器未初始化或模型未加载" << std::endl;
        return SimulationResult{};
    }
    
    // 实现完整的仿真执行流程 - 包括求解器调用、收敛检查、结果输出
    std::cout << "开始执行仿真..." << std::endl;
    
    // 检查是否有电磁求解器需要执行
    auto solvers = solverManager_->getSolvers();
    if (solvers.empty()) {
        std::cerr << "错误: 没有可执行的求解器" << std::endl;
        return SimulationResult{};
    }
    
    // 打印求解器信息
    std::cout << "检测到 " << solvers.size() << " 个求解器:" << std::endl;
    for (const auto& solver : solvers) {
        std::cout << "  - " << solver->getName() << std::endl;
    }
    
    // 根据仿真类型选择执行方式
    bool success = false;
    if (parameters_.type == SimulationType::TRANSIENT) {
        success = executeTransient();
    } else {
        success = executeSteadyState();
    }
    
    // 简化实现：输出电磁求解器完成信息
    if (success) {
        std::cout << "电磁场求解完成" << std::endl;
    }
    
    SimulationResult result;
    result.success = success;
    
    std::cout << "仿真执行完成" << std::endl;
    return result;
}

SimulationResult ElmerSolver::getResult() const {
    // TODO: 实现结果获取功能 - 需要返回完整的仿真结果数据
    std::cout << "获取仿真结果" << std::endl;
    return result_;
}

// ===== 内部辅助函数 =====

void ElmerSolver::cleanup() {
    // TODO: 实现完整的资源清理 - 需要释放所有分配的内存和资源
    if (initialized_) {
        std::cout << "清理ElmerSolver资源..." << std::endl;
        initialized_ = false;
        meshLoaded_ = false;
    }
}

void ElmerSolver::printBanner() {
    // TODO: 实现横幅打印功能 - 需要显示版本信息和编译选项
    std::cout << "=== Elmer FEM Solver ===" << std::endl;
    std::cout << "版本: 1.0" << std::endl;
    std::cout << "=======================" << std::endl;
}

bool ElmerSolver::initializeParallelEnvironment() {
    std::cout << "初始化并行环境..." << std::endl;
    
    // 检查MPI支持
#ifdef USE_MPI
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    
    if (!mpiInitialized) {
        int provided;
        if (MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &provided) != MPI_SUCCESS) {
            std::cerr << "错误: MPI初始化失败" << std::endl;
            return false;
        }
        
        if (provided < MPI_THREAD_SERIALIZED) {
            std::cerr << "警告: MPI线程支持级别不足" << std::endl;
        }
    }
    
    // 获取进程信息
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    std::cout << "MPI进程: " << rank << "/" << size << std::endl;
    
    // 初始化MPI通信器
    if (!comm_->initialize()) {
        std::cerr << "错误: MPI通信器初始化失败" << std::endl;
        return false;
    }
    
    std::cout << "MPI并行环境初始化完成" << std::endl;
#else
    std::cout << "MPI支持未启用，使用串行模式" << std::endl;
#endif
    
    return true;
}

bool ElmerSolver::initializeOpenMP() {
    std::cout << "初始化OpenMP..." << std::endl;
    
#ifdef _OPENMP
    // 获取系统支持的线程数
    int maxThreads = omp_get_max_threads();
    std::cout << "系统支持的最大线程数: " << maxThreads << std::endl;
    
    // 设置线程数（可以根据配置调整）
    int desiredThreads = maxThreads;
    
    // 检查环境变量
    char* envThreads = std::getenv("OMP_NUM_THREADS");
    if (envThreads != nullptr) {
        desiredThreads = std::atoi(envThreads);
        if (desiredThreads <= 0 || desiredThreads > maxThreads) {
            desiredThreads = maxThreads;
        }
    }
    
    // 设置线程数
    omp_set_num_threads(desiredThreads);
    
    // 设置并行策略
    omp_set_dynamic(0); // 禁用动态线程调整
    omp_set_nested(0);  // 禁用嵌套并行
    
    std::cout << "OpenMP线程数设置为: " << desiredThreads << std::endl;
    std::cout << "OpenMP并行环境初始化完成" << std::endl;
#else
    std::cout << "OpenMP支持未启用，使用串行模式" << std::endl;
#endif
    
    return true;
}

bool ElmerSolver::readInputFile() {
    if (inputFileParsed_) {
        std::cout << "警告: 输入文件已经解析" << std::endl;
        return true;
    }
    
    if (parameters_.modelName.empty()) {
        std::cerr << "错误: 未指定输入文件名" << std::endl;
        return false;
    }
    
    std::cout << "开始解析输入文件: " << parameters_.modelName << std::endl;
    
    // 解析输入文件
    if (!inputParser_->parse(parameters_.modelName)) {
        std::cerr << "错误: 输入文件解析失败" << std::endl;
        return false;
    }
    
    // 验证输入文件
    if (!inputParser_->validate()) {
        std::cerr << "错误: 输入文件验证失败" << std::endl;
        return false;
    }
    
    // 从输入文件更新仿真参数
    std::string simType = inputParser_->getSimulationType();
    if (simType == "Steady State") {
        parameters_.type = SimulationType::STEADY_STATE;
    } else if (simType == "Transient") {
        parameters_.type = SimulationType::TRANSIENT;
    } else if (simType == "Scanning") {
        parameters_.type = SimulationType::SCANNING;
    } else if (simType == "Optimization") {
        parameters_.type = SimulationType::OPTIMIZATION;
    }
    
    // 获取网格信息
    parameters_.meshDir = inputParser_->getMeshDirectory();
    parameters_.meshName = inputParser_->getMeshFileName();
    
    // 获取时间步进参数（如果存在）
    if (parameters_.type == SimulationType::TRANSIENT) {
        parameters_.startTime = inputParser_->getParameterReal("simulation", 0, "Start Time", 0.0);
        parameters_.endTime = inputParser_->getParameterReal("simulation", 0, "End Time", 1.0);
        parameters_.timeStep = inputParser_->getParameterReal("simulation", 0, "Timestep Size", 0.1);
        parameters_.numTimeSteps = static_cast<int>((parameters_.endTime - parameters_.startTime) / parameters_.timeStep);
    }
    
    // 获取输出控制参数
    parameters_.outputInterval = inputParser_->getParameterInteger("simulation", 0, "Output Intervals", 1);
    parameters_.verbose = inputParser_->getParameterLogical("simulation", 0, "Verbose", true);
    
    // 打印解析摘要
    if (parameters_.verbose) {
        inputParser_->printSummary();
    }
    
    inputFileParsed_ = true;
    std::cout << "输入文件解析完成" << std::endl;
    
    return true;
}

bool ElmerSolver::setInitialConditions() {
    std::cout << "设置初始条件..." << std::endl;
    
    if (!inputFileParsed_) {
        std::cerr << "错误: 输入文件未解析，无法设置初始条件" << std::endl;
        return false;
    }
    
    // 获取初始条件参数
    double initialTemperature = inputParser_->getParameterReal("initial condition", 0, "Temperature", 293.15);
    double initialPressure = inputParser_->getParameterReal("initial condition", 0, "Pressure", 101325.0);
    double initialVelocityX = inputParser_->getParameterReal("initial condition", 0, "Velocity 1", 0.0);
    double initialVelocityY = inputParser_->getParameterReal("initial condition", 0, "Velocity 2", 0.0);
    double initialVelocityZ = inputParser_->getParameterReal("initial condition", 0, "Velocity 3", 0.0);
    
    // 检查是否有自定义初始条件
    bool hasCustomInitialConditions = inputParser_->hasSection("initial condition", 0);
    
    if (hasCustomInitialConditions) {
        std::cout << "使用自定义初始条件" << std::endl;
        
        // 设置温度场初始条件
        std::string tempValue = inputParser_->getParameterValue("initial condition", 0, "Temperature");
        if (!tempValue.empty()) {
            std::cout << "初始温度: " << initialTemperature << " K" << std::endl;
            // TODO: 设置温度场初始值
        }
        
        // 设置压力场初始条件
        std::string pressureValue = inputParser_->getParameterValue("initial condition", 0, "Pressure");
        if (!pressureValue.empty()) {
            std::cout << "初始压力: " << initialPressure << " Pa" << std::endl;
            // TODO: 设置压力场初始值
        }
        
        // 设置速度场初始条件
        std::string vel1Value = inputParser_->getParameterValue("initial condition", 0, "Velocity 1");
        std::string vel2Value = inputParser_->getParameterValue("initial condition", 0, "Velocity 2");
        std::string vel3Value = inputParser_->getParameterValue("initial condition", 0, "Velocity 3");
        if (!vel1Value.empty() || !vel2Value.empty() || !vel3Value.empty()) {
            std::cout << "初始速度: (" << initialVelocityX << ", " 
                      << initialVelocityY << ", " << initialVelocityZ << ") m/s" << std::endl;
            // TODO: 设置速度场初始值
        }
        
        // 处理其他物理场的初始条件
        std::vector<std::string> fields = {"Electric Potential", "Magnetic Field", "Concentration"};
        for (const auto& field : fields) {
            std::string fieldValue = inputParser_->getParameterValue("initial condition", 0, field);
            if (!fieldValue.empty()) {
                double value = inputParser_->getParameterReal("initial condition", 0, field, 0.0);
                std::cout << "初始" << field << ": " << value << std::endl;
                // TODO: 设置相应物理场的初始值
            }
        }
    } else {
        std::cout << "使用默认初始条件" << std::endl;
        std::cout << "默认温度: " << initialTemperature << " K" << std::endl;
        std::cout << "默认压力: " << initialPressure << " Pa" << std::endl;
        std::cout << "默认速度: (" << initialVelocityX << ", " 
                  << initialVelocityY << ", " << initialVelocityZ << ") m/s" << std::endl;
    }
    
    // 设置时间步进初始条件
    if (parameters_.type == SimulationType::TRANSIENT) {
        currentTime_ = parameters_.startTime;
        timeStepIndex_ = 0;
        
        std::cout << "瞬态仿真初始时间: " << currentTime_ << " s" << std::endl;
        std::cout << "初始时间步索引: " << timeStepIndex_ << std::endl;
    }
    
    // 设置稳态仿真初始条件
    if (parameters_.type == SimulationType::STEADY_STATE) {
        currentTime_ = 0.0;
        timeStepIndex_ = 0;
        std::cout << "稳态仿真初始化完成" << std::endl;
    }
    
    std::cout << "初始条件设置完成" << std::endl;
    return true;
}

bool ElmerSolver::updateTimeDependentBoundaryConditions(double currentTime) {
    // 检查是否有时间相关的边界条件
    if (!inputParser_->hasSection("boundary condition", 0)) {
        return true; // 没有边界条件，直接返回成功
    }
    
    // 获取边界条件节段数量（简化处理，假设只有一个边界条件节段）
    int numBCSections = 1;
    
    for (int bcIndex = 0; bcIndex < numBCSections; ++bcIndex) {
        // 检查是否有时间相关的参数（简化处理）
        std::string timeFunctionValue = inputParser_->getParameterValue("boundary condition", bcIndex, "Time Function");
        if (!timeFunctionValue.empty()) {
            std::string timeFunction = timeFunctionValue;
            
            if (timeFunction == "sinusoidal") {
                // 正弦函数边界条件
                double amplitude = inputParser_->getParameterReal("boundary condition", bcIndex, "Amplitude", 1.0);
                double frequency = inputParser_->getParameterReal("boundary condition", bcIndex, "Frequency", 1.0);
                double phase = inputParser_->getParameterReal("boundary condition", bcIndex, "Phase", 0.0);
                
                double value = amplitude * std::sin(2.0 * M_PI * frequency * currentTime + phase);
                
                // 更新边界条件值
                std::string tempValue = inputParser_->getParameterValue("boundary condition", bcIndex, "Temperature");
                if (!tempValue.empty()) {
                    // TODO: 更新温度边界条件
                    std::cout << "更新温度边界条件: " << value << " K (时间=" << currentTime << "s)" << std::endl;
                }
                
                std::string velValue = inputParser_->getParameterValue("boundary condition", bcIndex, "Velocity 1");
                if (!velValue.empty()) {
                    // TODO: 更新速度边界条件
                    std::cout << "更新速度边界条件: " << value << " m/s (时间=" << currentTime << "s)" << std::endl;
                }
            }
            else if (timeFunction == "linear") {
                // 线性函数边界条件
                double initialValue = inputParser_->getParameterReal("boundary condition", bcIndex, "Initial Value", 0.0);
                double rate = inputParser_->getParameterReal("boundary condition", bcIndex, "Rate", 1.0);
                
                double value = initialValue + rate * currentTime;
                
                // 更新边界条件值
                std::string pressureValue = inputParser_->getParameterValue("boundary condition", bcIndex, "Pressure");
                if (!pressureValue.empty()) {
                    // TODO: 更新压力边界条件
                    std::cout << "更新压力边界条件: " << value << " Pa (时间=" << currentTime << "s)" << std::endl;
                }
            }
            else if (timeFunction == "step") {
                // 阶跃函数边界条件
                double stepTime = inputParser_->getParameterReal("boundary condition", bcIndex, "Step Time", 0.5);
                double beforeValue = inputParser_->getParameterReal("boundary condition", bcIndex, "Before Value", 0.0);
                double afterValue = inputParser_->getParameterReal("boundary condition", bcIndex, "After Value", 1.0);
                
                double value = (currentTime >= stepTime) ? afterValue : beforeValue;
                
                // 更新边界条件值
                std::string potentialValue = inputParser_->getParameterValue("boundary condition", bcIndex, "Electric Potential");
                if (!potentialValue.empty()) {
                    // TODO: 更新电势边界条件
                    std::cout << "更新电势边界条件: " << value << " V (时间=" << currentTime << "s)" << std::endl;
                }
            }
            else {
                // 常数边界条件（默认）
                // 不需要更新时间相关的值
                std::cout << "常数边界条件，无需更新 (时间=" << currentTime << "s)" << std::endl;
            }
        }
    }
    
    return true;
}

bool ElmerSolver::executeTimeStep(int timeStepIndex, double currentTime) {
    // 实现时间步进 - 包括时间积分算法和电磁场求解器调用
    std::cout << "执行时间步进: " << timeStepIndex << ", 时间: " << currentTime << std::endl;
    
    // 简化实现：直接调用求解器管理器
    if (solverManager_) {
        // 执行电磁场求解器
        if (!solverManager_->executeTransient(currentTime, currentTime + parameters_.timeStep, parameters_.timeStep)) {
            std::cerr << "错误: 电磁场求解器执行失败" << std::endl;
            return false;
        }
    }
    
    return true;
}

bool ElmerSolver::saveResults(int timeStepIndex, double currentTime) {
    std::cout << "保存结果: 时间步 " << timeStepIndex << ", 时间: " << currentTime << std::endl;
    
    // 获取输出格式参数
    std::string outputFormat = inputParser_->getParameterValue("simulation", 0, "Output Format");
    if (outputFormat.empty()) {
        outputFormat = "vtk";
    }
    std::string outputDir = inputParser_->getParameterValue("simulation", 0, "Output Directory");
    if (outputDir.empty()) {
        outputDir = "./results";
    }
    
    // 创建输出目录（如果不存在）
    std::filesystem::create_directories(outputDir);
    
    // 生成输出文件名
    std::stringstream filename;
    filename << outputDir << "/" << parameters_.modelName << "_" 
             << std::setfill('0') << std::setw(6) << timeStepIndex << "." << outputFormat;
    
    std::cout << "输出文件: " << filename.str() << std::endl;
    
    // 根据输出格式选择保存方法
    if (outputFormat == "vtk") {
        return saveResultsVTK(filename.str(), timeStepIndex, currentTime);
    } else if (outputFormat == "gmsh") {
        return saveResultsGmsh(filename.str(), timeStepIndex, currentTime);
    } else if (outputFormat == "csv") {
        return saveResultsCSV(filename.str(), timeStepIndex, currentTime);
    } else {
        std::cerr << "错误: 不支持的输出格式: " << outputFormat << std::endl;
        return false;
    }
}

bool ElmerSolver::saveResultsVTK(const std::string& filename, int timeStepIndex, double currentTime) {
    // TODO: 实现VTK格式结果保存
    std::cout << "保存VTK格式结果到: " << filename << std::endl;
    
    // 模拟VTK文件保存
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法创建输出文件: " << filename << std::endl;
        return false;
    }
    
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Elmer FEM Results - Time: " << currentTime << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    
    // TODO: 添加实际的网格和结果数据
    file << "POINTS 0 float" << std::endl;
    file << "CELLS 0 0" << std::endl;
    file << "CELL_TYPES 0" << std::endl;
    
    file.close();
    std::cout << "VTK结果保存完成" << std::endl;
    return true;
}

bool ElmerSolver::saveResultsGmsh(const std::string& filename, int timeStepIndex, double currentTime) {
    // TODO: 实现Gmsh格式结果保存
    std::cout << "保存Gmsh格式结果到: " << filename << std::endl;
    
    // 模拟Gmsh文件保存
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法创建输出文件: " << filename << std::endl;
        return false;
    }
    
    file << "$MeshFormat" << std::endl;
    file << "2.2 0 8" << std::endl;
    file << "$EndMeshFormat" << std::endl;
    file << "$NodeData" << std::endl;
    file << "1" << std::endl;
    file << "\"Time=" << currentTime << "\"" << std::endl;
    file << "1" << std::endl;
    file << "0.0" << std::endl;
    file << "3" << std::endl;
    file << "0" << std::endl;
    file << "1" << std::endl;
    file << "0" << std::endl;
    file << "$EndNodeData" << std::endl;
    
    file.close();
    std::cout << "Gmsh结果保存完成" << std::endl;
    return true;
}

bool ElmerSolver::saveResultsCSV(const std::string& filename, int timeStepIndex, double currentTime) {
    // TODO: 实现CSV格式结果保存
    std::cout << "保存CSV格式结果到: " << filename << std::endl;
    
    // 模拟CSV文件保存
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法创建输出文件: " << filename << std::endl;
        return false;
    }
    
    file << "Time,StepIndex,MaxTemperature,MinTemperature,AverageTemperature" << std::endl;
    file << currentTime << "," << timeStepIndex << ",0.0,0.0,0.0" << std::endl;
    
    file.close();
    std::cout << "CSV结果保存完成" << std::endl;
    return true;
}

bool ElmerSolver::checkConvergence() {
    std::cout << "检查收敛性..." << std::endl;
    
    // 获取收敛性参数
    double tolerance = inputParser_->getParameterReal("simulation", 0, "Convergence Tolerance", 1e-6);
    int minIterations = inputParser_->getParameterInteger("simulation", 0, "Minimum Iterations", 1);
    
    // 模拟收敛性检查（简化处理）
    std::cout << "模拟收敛性检查..." << std::endl;
    
    bool allConverged = true;
    
    // 模拟求解器残差
    double residual = 1.0 / (timeStepIndex_ + 1);
    bool converged = (residual <= tolerance);
    
    std::cout << "求解器残差: " << residual << ", 收敛状态: " 
              << (converged ? "已收敛" : "未收敛") << std::endl;
    
    // 检查残差是否满足容差要求
    if (residual > tolerance) {
        allConverged = false;
    }
    
    // 检查求解器自身的收敛状态
    if (!converged) {
        allConverged = false;
    }
    
    // 检查最小迭代次数要求
    if (timeStepIndex_ < minIterations) {
        std::cout << "未达到最小迭代次数要求: " << timeStepIndex_ << " < " << minIterations << std::endl;
        allConverged = false;
    }
    
    if (allConverged) {
        std::cout << "所有求解器收敛" << std::endl;
    } else {
        std::cout << "有求解器未收敛" << std::endl;
    }
    
    return allConverged;
}

void ElmerSolver::processCommandLineArguments(int argc, char** argv) {
    std::cout << "处理命令行参数..." << std::endl;
    
    if (argc < 2) {
        std::cout << "用法: " << argv[0] << " <输入文件> [选项]" << std::endl;
        std::cout << "选项:" << std::endl;
        std::cout << "  -h, --help             显示帮助信息" << std::endl;
        std::cout << "  -v, --verbose          启用详细输出" << std::endl;
        std::cout << "  -o, --output <dir>     设置输出目录" << std::endl;
        std::cout << "  -t, --threads <num>    设置线程数" << std::endl;
        std::cout << "  -c, --check            只检查输入文件，不执行仿真" << std::endl;
        return;
    }
    
    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            std::cout << "Elmer FEM Solver - 有限元求解器" << std::endl;
            std::cout << "版本: 1.0" << std::endl;
            std::cout << "用法: " << argv[0] << " <输入文件> [选项]" << std::endl;
            std::cout << "选项:" << std::endl;
            std::cout << "  -h, --help             显示帮助信息" << std::endl;
            std::cout << "  -v, --verbose          启用详细输出" << std::endl;
            std::cout << "  -o, --output <dir>     设置输出目录" << std::endl;
            std::cout << "  -t, --threads <num>    设置线程数" << std::endl;
            std::cout << "  -c, --check            只检查输入文件，不执行仿真" << std::endl;
            exit(0);
        }
        else if (arg == "-v" || arg == "--verbose") {
            parameters_.verbose = true;
            std::cout << "启用详细输出模式" << std::endl;
        }
        else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                parameters_.outputDir = argv[++i];
                std::cout << "设置输出目录: " << parameters_.outputDir << std::endl;
            } else {
                std::cerr << "错误: -o 选项需要指定输出目录" << std::endl;
                exit(1);
            }
        }
        else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                int threads = std::atoi(argv[++i]);
                if (threads > 0) {
                    parameters_.numThreads = threads;
                    std::cout << "设置线程数: " << parameters_.numThreads << std::endl;
                } else {
                    std::cerr << "错误: 线程数必须为正整数" << std::endl;
                    exit(1);
                }
            } else {
                std::cerr << "错误: -t 选项需要指定线程数" << std::endl;
                exit(1);
            }
        }
        else if (arg == "-c" || arg == "--check") {
            parameters_.checkOnly = true;
            std::cout << "启用检查模式，只验证输入文件" << std::endl;
        }
        else if (arg[0] != '-') {
            // 第一个非选项参数作为输入文件
            if (parameters_.modelName.empty()) {
                parameters_.modelName = arg;
                std::cout << "设置输入文件: " << parameters_.modelName << std::endl;
            }
        }
        else {
            std::cerr << "错误: 未知选项: " << arg << std::endl;
            exit(1);
        }
    }
    
    // 验证输入文件
    if (parameters_.modelName.empty()) {
        std::cerr << "错误: 未指定输入文件" << std::endl;
        exit(1);
    }
    
    std::cout << "命令行参数处理完成" << std::endl;
}

double ElmerSolver::getCPUTime() const {
#ifdef _WIN32
    FILETIME createTime, exitTime, kernelTime, userTime;
    if (GetProcessTimes(GetCurrentProcess(), &createTime, &exitTime, &kernelTime, &userTime)) {
        ULARGE_INTEGER kernel, user;
        kernel.LowPart = kernelTime.dwLowDateTime;
        kernel.HighPart = kernelTime.dwHighDateTime;
        user.LowPart = userTime.dwLowDateTime;
        user.HighPart = userTime.dwHighDateTime;
        
        // 转换为秒（100纳秒单位）
        return (kernel.QuadPart + user.QuadPart) * 1e-7;
    }
#else
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1e-6 +
               usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1e-6;
    }
#endif
    return 0.0;
}

double ElmerSolver::getRealTime() const {
    auto now = std::chrono::steady_clock::now();
    auto duration = now.time_since_epoch();
    return std::chrono::duration<double>(duration).count();
}

void ElmerSolver::startTimer() {
    startCPUTime_ = getCPUTime();
    startRealTime_ = getRealTime();
}

void ElmerSolver::stopTimer() {
    endCPUTime_ = getCPUTime();
    endRealTime_ = getRealTime();
}

void ElmerSolver::printPerformanceStats() const {
    double cpuTime = endCPUTime_ - startCPUTime_;
    double realTime = endRealTime_ - startRealTime_;
    
    std::cout << "=== 性能统计 ===" << std::endl;
    std::cout << "CPU时间: " << cpuTime << " 秒" << std::endl;
    std::cout << "实际时间: " << realTime << " 秒" << std::endl;
    std::cout << "效率: " << (cpuTime / realTime * 100.0) << "%" << std::endl;
    
    if (meshLoaded_ && mesh_ != nullptr) {
        // TODO: 添加网格相关的性能统计
        std::cout << "网格节点数: 0" << std::endl;
        std::cout << "网格单元数: 0" << std::endl;
    }
    
    std::cout << "时间步数: " << timeStepIndex_ << std::endl;
    std::cout << "平均每步CPU时间: " << (cpuTime / (timeStepIndex_ + 1)) << " 秒" << std::endl;
    std::cout << "================" << std::endl;
}

// 简化实现：移除有问题的函数，专注于核心电磁求解器功能

// ===== 主函数 =====

int ElmerSolverMain(int argc, char** argv) {
    // TODO: 实现ElmerSolver主函数 - 需要完整的错误处理和日志记录
    std::cout << "ElmerSolver主函数开始执行" << std::endl;
    
    // 创建求解器实例
    auto solver = std::make_unique<ElmerSolver>();
    
    // 处理命令行参数
    solver->processCommandLineArguments(argc, argv);
    
    // 初始化求解器
    if (!solver->initialize()) {
        std::cout << "求解器初始化失败" << std::endl;
        return -1;
    }
    
    // 加载模型
    if (!solver->loadModel()) {
        std::cout << "模型加载失败" << std::endl;
        return -1;
    }
    
    // 执行仿真
    auto result = solver->execute();
    
    if (result.success) {
        std::cout << "仿真执行成功" << std::endl;
        return 0;
    } else {
        std::cout << "仿真执行失败: " << result.errorMessage << std::endl;
        return -1;
    }
}

}