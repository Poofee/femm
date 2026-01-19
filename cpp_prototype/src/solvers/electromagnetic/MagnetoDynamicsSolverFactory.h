#pragma once

#include "MagnetoDynamicsSolverBase.h"
#include "MagnetoDynamics2DSolver.h"
#include "MagnetoDynamics3DSolver.h"
#include <memory>
#include <string>
#include <unordered_map>

namespace elmer {

/**
 * @brief 低频电磁求解器工厂类
 * 
 * 提供统一的接口来创建和管理2D和3D低频电磁求解器
 */
class MagnetoDynamicsSolverFactory {
public:
    MagnetoDynamicsSolverFactory() = default;
    ~MagnetoDynamicsSolverFactory() = default;
    
    /**
     * @brief 创建低频电磁求解器
     * 
     * @param dimension 求解器维度（2D或3D）
     * @param coordSystem 坐标系类型
     * @return std::unique_ptr<MagnetoDynamicsSolverBase> 求解器指针
     */
    static std::unique_ptr<MagnetoDynamicsSolverBase> createSolver(
        MagnetoDynamicsDimension dimension,
        CoordinateSystemType coordSystem = CoordinateSystemType::CARTESIAN);
    
    /**
     * @brief 创建2D磁动力学求解器
     * 
     * @param coordSystem 坐标系类型
     * @return std::unique_ptr<MagnetoDynamics2DSolver> 2D求解器指针
     */
    static std::unique_ptr<MagnetoDynamics2DSolver> create2DSolver(
        CoordinateSystemType coordSystem = CoordinateSystemType::CARTESIAN);
    
    /**
     * @brief 创建3D磁动力学求解器
     * 
     * @return std::unique_ptr<MagnetoDynamics3DSolver> 3D求解器指针
     */
    static std::unique_ptr<MagnetoDynamics3DSolver> create3DSolver();
    
    /**
     * @brief 从配置文件创建求解器
     * 
     * @param configFile 配置文件路径
     * @return std::unique_ptr<MagnetoDynamicsSolverBase> 求解器指针
     */
    static std::unique_ptr<MagnetoDynamicsSolverBase> createSolverFromConfig(
        const std::string& configFile);
    
    /**
     * @brief 从求解器名称创建求解器
     * 
     * @param solverName 求解器名称（"MagnetoDynamics2D" 或 "MagnetoDynamics3D"）
     * @param coordSystem 坐标系类型
     * @return std::unique_ptr<MagnetoDynamicsSolverBase> 求解器指针
     */
    static std::unique_ptr<MagnetoDynamicsSolverBase> createSolverByName(
        const std::string& solverName,
        CoordinateSystemType coordSystem = CoordinateSystemType::CARTESIAN);
    
    /**
     * @brief 获取支持的求解器列表
     * 
     * @return std::vector<std::string> 求解器名称列表
     */
    static std::vector<std::string> getSupportedSolvers();
    
    /**
     * @brief 检查求解器是否支持
     * 
     * @param solverName 求解器名称
     * @return bool 是否支持
     */
    static bool isSolverSupported(const std::string& solverName);
    
    /**
     * @brief 获取求解器的默认参数
     * 
     * @param dimension 求解器维度
     * @return MagnetoDynamicsParameters 默认参数
     */
    static MagnetoDynamicsParameters getDefaultParameters(MagnetoDynamicsDimension dimension);
    
    /**
     * @brief 获取2D求解器的默认参数
     * 
     * @param coordSystem 坐标系类型
     * @return MagnetoDynamicsParameters 默认参数
     */
    static MagnetoDynamicsParameters getDefault2DParameters(
        CoordinateSystemType coordSystem = CoordinateSystemType::CARTESIAN);
    
    /**
     * @brief 获取3D求解器的默认参数
     * 
     * @return MagnetoDynamicsParameters 默认参数
     */
    static MagnetoDynamicsParameters getDefault3DParameters();
    
    /**
     * @brief 验证求解器参数
     * 
     * @param params 求解器参数
     * @param dimension 求解器维度
     * @return std::string 验证错误信息，空字符串表示验证通过
     */
    static std::string validateParameters(const MagnetoDynamicsParameters& params, 
                                         MagnetoDynamicsDimension dimension);
    
    /**
     * @brief 获取求解器描述信息
     * 
     * @param dimension 求解器维度
     * @return std::string 描述信息
     */
    static std::string getSolverDescription(MagnetoDynamicsDimension dimension);
    
    /**
     * @brief 获取求解器版本信息
     * 
     * @param dimension 求解器维度
     * @return std::string 版本信息
     */
    static std::string getSolverVersion(MagnetoDynamicsDimension dimension);
    
    /**
     * @brief 获取求解器能力信息
     * 
     * @param dimension 求解器维度
     * @return std::unordered_map<std::string, bool> 能力映射表
     */
    static std::unordered_map<std::string, bool> getSolverCapabilities(
        MagnetoDynamicsDimension dimension);
    
private:
    // 求解器注册表
    static std::unordered_map<std::string, MagnetoDynamicsDimension> solverRegistry_;
    
    // 初始化求解器注册表
    static void initializeSolverRegistry();
    
    // 配置文件解析
    static MagnetoDynamicsParameters parseConfigFile(const std::string& configFile);
    static MagnetoDynamicsDimension parseDimensionFromConfig(const std::string& configFile);
    static CoordinateSystemType parseCoordinateSystemFromConfig(const std::string& configFile);
};

} // namespace elmer