#pragma once

#include "Types.h"
#include "LoggerInterface.h"
#include <memory>
#include <string>
#include <vector>

namespace elmer {

// 前向声明
class Solver;
class Mesh;
class Model;

/**
 * @brief 主程序工具类
 * 
 * 移植自Fortran的MainUtils模块，提供Elmer主程序所需的核心工具函数
 */
class MainUtils {
public:
    /**
     * @brief 检查线性求解器选项的可行性
     */
    static bool checkLinearSolverOptions(std::shared_ptr<Solver> solver);
    
    /**
     * @brief 设置归一化关键字
     */
    static void setNormalizedKeywords(std::shared_ptr<Model> model, std::shared_ptr<Mesh> mesh);
    
    /**
     * @brief 设置旋转属性
     */
    static void setRotatedProperties(std::shared_ptr<Model> model, std::shared_ptr<Mesh> mesh);
    
    /**
     * @brief 角度转旋转矩阵
     */
    static Eigen::Matrix3d anglesToRotationMatrix(const std::vector<double>& angles, 
                                                 const std::string& rotateOrder = "XYZ");
    
    /**
     * @brief 添加求解器过程
     */
    static void addSolverProcedure(std::shared_ptr<Solver> solver, const std::string& procedure);
    
    /**
     * @brief 交换网格
     */
    static bool swapMesh(std::shared_ptr<Model> model, std::shared_ptr<Mesh> mesh, 
                        const std::string& name, std::shared_ptr<Mesh> extMesh);
    
    /**
     * @brief 检查并创建DG索引
     */
    static bool checkAndCreateDGIndexes(std::shared_ptr<Mesh> mesh, bool activeElem = true);
    
    /**
     * @brief 创建DG排列
     */
    static std::vector<int> createDGPerm(std::shared_ptr<Solver> solver, 
                                        const std::string& maskName = "", 
                                        const std::string& secName = "");
    
    /**
     * @brief 创建节点排列
     */
    static std::vector<int> createNodalPerm(std::shared_ptr<Solver> solver, int nSize);
    
    /**
     * @brief 创建单元排列
     */
    static std::vector<int> createElementsPerm(std::shared_ptr<Solver> solver, int nsize,
                                              const std::string& maskName = "", 
                                              const std::string& secName = "");
    
    /**
     * @brief 添加执行时机标志
     */
    static void addExecWhenFlag(std::shared_ptr<Solver> solver);
    
    /**
     * @brief 添加方程基础
     */
    static bool addEquationBasics(std::shared_ptr<Solver> solver, const std::string& name, 
                                 bool transient = false);
    
    /**
     * @brief 创建时间导数变量
     */
    static bool createTimeDerivativeVariables(std::shared_ptr<Solver> solver, 
                                             const std::string& varName);
    
    /**
     * @brief 添加方程解
     */
    static bool addEquationSolution(std::shared_ptr<Solver> solver, bool transient = false);
    
    /**
     * @brief 创建子求解器
     */
    static std::shared_ptr<Solver> createChildSolver(std::shared_ptr<Solver> parentSolver,
                                                    const std::string& childVarName,
                                                    int childDofs,
                                                    const std::string& childPrefix,
                                                    bool noReuse = false);
    
    /**
     * @brief 求解方程系统
     */
    static bool solveEquations(std::shared_ptr<Model> model, double dt, 
                              bool transientSimulation = false);
    
    /**
     * @brief 耦合求解器
     */
    static bool coupledSolver(std::shared_ptr<Model> model, std::shared_ptr<Solver> solver,
                             double dt, bool transient = false);
    
    /**
     * @brief 块求解器
     */
    static bool blockSolver(std::shared_ptr<Model> model, std::shared_ptr<Solver> solver,
                           double dt, bool transient = false);
    
    /**
     * @brief 单求解器
     */
    static bool singleSolver(std::shared_ptr<Model> model, std::shared_ptr<Solver> solver,
                            double dt, bool transientSimulation = false);
    
    /**
     * @brief 求解器激活
     */
    static bool solverActivate(std::shared_ptr<Model> model, std::shared_ptr<Solver> solver,
                              double dt, bool transientSimulation = false);
    
    /**
     * @brief 预测-校正控制
     */
    static bool predictorCorrectorControl(std::shared_ptr<Model> model, double dt, 
                                         double realTimestep);
    
    /**
     * @brief 时间步长控制器
     */
    static double timeStepController(double dt, double dtOld, double eta, double etaOld,
                                    double epsilon, double beta1, double beta2);
    
    /**
     * @brief 时间步长限制器
     */
    static double timeStepLimiter(double dtOld, double dt, double gfactor, double k);
    
    /**
     * @brief 预测-校正误差估计
     */
    static double predCorrErrorEstimate(double eta, double dt, int predCorrOrder, 
                                       double timeError, double zeta);
    
    /**
     * @brief 读取预测-校正参数
     */
    static bool readPredCorrParams(std::shared_ptr<Model> model, 
                                  const std::vector<std::string>& solverParams,
                                  int& outOrder, double& outEps, double& outB1, double& outB2);
    
    /**
     * @brief 获取版本信息
     */
    static std::string getVersion();
    
    /**
     * @brief 获取修订信息
     */
    static std::string getRevision();
    
    /**
     * @brief 获取编译日期
     */
    static std::string getCompilationDate();
    
    /**
     * @brief 获取当前时间
     */
    static double getCurrentTime();
    
    /**
     * @brief 获取CPU时间
     */
    static double getCPUTime();
    
    /**
     * @brief 初始化并行环境
     */
    static bool initializeParallelEnvironment();
    
    /**
     * @brief 初始化元素描述
     */
    static bool initializeElementDescriptions();
    
    /**
     * @brief 设置OpenMP线程数
     */
    static void setOpenMPThreads(int numThreads);
    
    /**
     * @brief 获取环境变量
     */
    static std::string getEnvironmentVariable(const std::string& name);
    
    /**
     * @brief 信息输出
     */
    static void info(const std::string& module, const std::string& message, int level = 1);
    
    /**
     * @brief 警告输出
     */
    static void warn(const std::string& module, const std::string& message);
    
    /**
     * @brief 错误输出
     */
    static void error(const std::string& module, const std::string& message);
    
private:
    /**
     * @brief 检查线性求解器选项的私有实现
     */
    static bool checkLinearSolverOptionsImpl(std::shared_ptr<Solver> solver);
    
    /**
     * @brief 求解耦合系统的私有实现
     */
    static bool solveCoupledImpl(std::shared_ptr<Model> model);
    
    /**
     * @brief 边界残差计算
     */
    static double boundaryResidual(std::shared_ptr<Model> model, int edge, 
                                  std::shared_ptr<Mesh> mesh, 
                                  const std::vector<double>& quant, 
                                  const std::vector<int>& perm, double gnorm);
    
    /**
     * @brief 边残差计算
     */
    static double edgeResidual(std::shared_ptr<Model> model, int edge, 
                              std::shared_ptr<Mesh> mesh, 
                              const std::vector<double>& quant, 
                              const std::vector<int>& perm);
    
    /**
     * @brief 内部残差计算
     */
    static double insideResidual(std::shared_ptr<Model> model, int element, 
                                std::shared_ptr<Mesh> mesh, 
                                const std::vector<double>& quant, 
                                const std::vector<int>& perm, double fnorm);
};

} // namespace elmer