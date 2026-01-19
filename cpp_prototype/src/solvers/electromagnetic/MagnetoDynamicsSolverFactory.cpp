#include "MagnetoDynamicsSolverFactory.h"
#include "LoggerFactory.h"
#include <fstream>
#include <sstream>
#include <algorithm>

namespace elmer {

// 静态成员变量定义
std::unordered_map<std::string, MagnetoDynamicsDimension> 
    MagnetoDynamicsSolverFactory::solverRegistry_;

std::unique_ptr<MagnetoDynamicsSolverBase> 
MagnetoDynamicsSolverFactory::createSolver(MagnetoDynamicsDimension dimension,
                                          CoordinateSystemType coordSystem) {
    
    // 确保注册表已初始化
    if (solverRegistry_.empty()) {
        initializeSolverRegistry();
    }
    
    switch (dimension) {
        case MagnetoDynamicsDimension::DIM_2D:
            return std::make_unique<MagnetoDynamics2DSolver>(coordSystem);
            
        case MagnetoDynamicsDimension::DIM_3D:
            return std::make_unique<MagnetoDynamics3DSolver>();
            
        default:
            ELMER_ERROR("错误: 不支持的求解器维度");
            return nullptr;
    }
}

std::unique_ptr<MagnetoDynamics2DSolver> 
MagnetoDynamicsSolverFactory::create2DSolver(CoordinateSystemType coordSystem) {
    ELMER_INFO("创建2D磁动力学求解器");
    return std::make_unique<MagnetoDynamics2DSolver>(coordSystem);
}

std::unique_ptr<MagnetoDynamics3DSolver> 
MagnetoDynamicsSolverFactory::create3DSolver() {
    ELMER_INFO("创建3D磁动力学求解器");
    return std::make_unique<MagnetoDynamics3DSolver>();
}

std::unique_ptr<MagnetoDynamicsSolverBase> 
MagnetoDynamicsSolverFactory::createSolverFromConfig(const std::string& configFile) {
    ELMER_INFO("从配置文件创建求解器: {}", configFile);
    
    // 解析配置文件
    MagnetoDynamicsParameters params = parseConfigFile(configFile);
    MagnetoDynamicsDimension dimension = parseDimensionFromConfig(configFile);
    CoordinateSystemType coordSystem = parseCoordinateSystemFromConfig(configFile);
    
    // 创建求解器
    auto solver = createSolver(dimension, coordSystem);
    if (solver) {
        solver->setParameters(params);
    }
    
    return solver;
}

std::unique_ptr<MagnetoDynamicsSolverBase> 
MagnetoDynamicsSolverFactory::createSolverByName(const std::string& solverName,
                                                CoordinateSystemType coordSystem) {
    ELMER_INFO("按名称创建求解器: {}", solverName);
    
    // 确保注册表已初始化
    if (solverRegistry_.empty()) {
        initializeSolverRegistry();
    }
    
    // 查找求解器
    auto it = solverRegistry_.find(solverName);
    if (it == solverRegistry_.end()) {
        ELMER_ERROR("错误: 不支持的求解器名称: {}", solverName);
        return nullptr;
    }
    
    // 创建求解器
    return createSolver(it->second, coordSystem);
}

std::vector<std::string> MagnetoDynamicsSolverFactory::getSupportedSolvers() {
    // 确保注册表已初始化
    if (solverRegistry_.empty()) {
        initializeSolverRegistry();
    }
    
    std::vector<std::string> solvers;
    for (const auto& entry : solverRegistry_) {
        solvers.push_back(entry.first);
    }
    
    return solvers;
}

bool MagnetoDynamicsSolverFactory::isSolverSupported(const std::string& solverName) {
    // 确保注册表已初始化
    if (solverRegistry_.empty()) {
        initializeSolverRegistry();
    }
    
    return solverRegistry_.find(solverName) != solverRegistry_.end();
}

MagnetoDynamicsParameters 
MagnetoDynamicsSolverFactory::getDefaultParameters(MagnetoDynamicsDimension dimension) {
    MagnetoDynamicsParameters params;
    
    switch (dimension) {
        case MagnetoDynamicsDimension::DIM_2D:
            return getDefault2DParameters();
            
        case MagnetoDynamicsDimension::DIM_3D:
            return getDefault3DParameters();
            
        default:
            ELMER_WARN("警告: 未知维度，返回默认参数");
            return params;
    }
}

MagnetoDynamicsParameters 
MagnetoDynamicsSolverFactory::getDefault2DParameters(CoordinateSystemType coordSystem) {
    MagnetoDynamicsParameters params;
    
    // 2D求解器默认参数
    params.tolerance = 1.0e-6;
    params.maxIterations = 1000;
    params.maxNonlinearIterations = 50;
    params.includeEddyCurrents = true;
    params.includeConvection = false;
    params.isTransient = false;
    params.isHarmonic = false;
    params.frequency = 0.0;
    params.timeStep = 0.01;
    params.endTime = 1.0;
    
    // 根据坐标系调整参数
    if (coordSystem == CoordinateSystemType::AXISYMMETRIC) {
        params.includeEddyCurrents = true; // 轴对称通常包含涡流
    }
    
    return params;
}

MagnetoDynamicsParameters MagnetoDynamicsSolverFactory::getDefault3DParameters() {
    MagnetoDynamicsParameters params;
    
    // 3D求解器默认参数
    params.tolerance = 1.0e-6;
    params.maxIterations = 2000;
    params.maxNonlinearIterations = 100;
    params.includeEddyCurrents = true;
    params.includeConvection = false;
    params.isTransient = false;
    params.isHarmonic = false;
    params.frequency = 0.0;
    params.timeStep = 0.01;
    params.endTime = 1.0;
    
    return params;
}

std::string MagnetoDynamicsSolverFactory::validateParameters(
    const MagnetoDynamicsParameters& params, MagnetoDynamicsDimension dimension) {
    
    std::stringstream errors;
    
    // 验证容差
    if (params.tolerance <= 0.0) {
        errors << "容差必须为正数\n";
    }
    
    // 验证最大迭代次数
    if (params.maxIterations <= 0) {
        errors << "最大迭代次数必须为正数\n";
    }
    
    // 验证非线性迭代次数
    if (params.maxNonlinearIterations <= 0) {
        errors << "非线性迭代次数必须为正数\n";
    }
    
    // 验证频率
    if (params.isHarmonic && params.frequency <= 0.0) {
        errors << "谐波分析时频率必须为正数\n";
    }
    
    // 验证时间步长
    if (params.isTransient && params.timeStep <= 0.0) {
        errors << "瞬态分析时时间步长必须为正数\n";
    }
    
    // 验证结束时间
    if (params.isTransient && params.endTime <= 0.0) {
        errors << "瞬态分析时结束时间必须为正数\n";
    }
    
    // 维度特定验证
    if (dimension == MagnetoDynamicsDimension::DIM_3D) {
        if (params.maxIterations < 500) {
            errors << "3D求解器建议使用更大的最大迭代次数\n";
        }
    }
    
    return errors.str();
}

std::string MagnetoDynamicsSolverFactory::getSolverDescription(MagnetoDynamicsDimension dimension) {
    switch (dimension) {
        case MagnetoDynamicsDimension::DIM_2D:
            return "2D磁动力学求解器 - 用于二维电磁场分析，支持笛卡尔、轴对称和柱对称坐标系";
            
        case MagnetoDynamicsDimension::DIM_3D:
            return "3D磁动力学求解器 - 用于三维电磁场分析，支持Whitney边元和Lagrange单元";
            
        default:
            return "未知求解器";
    }
}

std::string MagnetoDynamicsSolverFactory::getSolverVersion(MagnetoDynamicsDimension dimension) {
    switch (dimension) {
        case MagnetoDynamicsDimension::DIM_2D:
            return "2.0.0";
            
        case MagnetoDynamicsDimension::DIM_3D:
            return "3.0.0";
            
        default:
            return "1.0.0";
    }
}

std::unordered_map<std::string, bool> 
MagnetoDynamicsSolverFactory::getSolverCapabilities(MagnetoDynamicsDimension dimension) {
    std::unordered_map<std::string, bool> capabilities;
    
    switch (dimension) {
        case MagnetoDynamicsDimension::DIM_2D:
            capabilities["稳态分析"] = true;
            capabilities["瞬态分析"] = true;
            capabilities["谐波分析"] = true;
            capabilities["涡流分析"] = true;
            capabilities["非线性材料"] = true;
            capabilities["轴对称"] = true;
            capabilities["柱对称"] = true;
            break;
            
        case MagnetoDynamicsDimension::DIM_3D:
            capabilities["稳态分析"] = true;
            capabilities["瞬态分析"] = true;
            capabilities["谐波分析"] = true;
            capabilities["涡流分析"] = true;
            capabilities["非线性材料"] = true;
            capabilities["Whitney边元"] = true;
            capabilities["Lagrange单元"] = true;
            capabilities["Piola变换"] = true;
            break;
            
        default:
            break;
    }
    
    return capabilities;
}

// 私有方法实现
void MagnetoDynamicsSolverFactory::initializeSolverRegistry() {
    solverRegistry_["MagnetoDynamics2D"] = MagnetoDynamicsDimension::DIM_2D;
    solverRegistry_["MagnetoDynamics3D"] = MagnetoDynamicsDimension::DIM_3D;
    solverRegistry_["MagnetoDynamics2DSolver"] = MagnetoDynamicsDimension::DIM_2D;
    solverRegistry_["MagnetoDynamics3DSolver"] = MagnetoDynamicsDimension::DIM_3D;
    
    ELMER_DEBUG("求解器注册表初始化完成");
}

MagnetoDynamicsParameters 
MagnetoDynamicsSolverFactory::parseConfigFile(const std::string& configFile) {
    MagnetoDynamicsParameters params;
    
    // 简单的配置文件解析
    std::ifstream file(configFile);
    if (!file.is_open()) {
        ELMER_ERROR("错误: 无法打开配置文件: {}", configFile);
        return params;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // 跳过注释和空行
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            // 去除空格
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);
            
            // 解析参数
            if (key == "tolerance") {
                params.tolerance = std::stod(value);
            } else if (key == "maxIterations") {
                params.maxIterations = std::stoi(value);
            } else if (key == "maxNonlinearIterations") {
                params.maxNonlinearIterations = std::stoi(value);
            } else if (key == "includeEddyCurrents") {
                params.includeEddyCurrents = (value == "true" || value == "1");
            } else if (key == "includeConvection") {
                params.includeConvection = (value == "true" || value == "1");
            } else if (key == "isTransient") {
                params.isTransient = (value == "true" || value == "1");
            } else if (key == "isHarmonic") {
                params.isHarmonic = (value == "true" || value == "1");
            } else if (key == "frequency") {
                params.frequency = std::stod(value);
            } else if (key == "timeStep") {
                params.timeStep = std::stod(value);
            } else if (key == "endTime") {
                params.endTime = std::stod(value);
            }
        }
    }
    
    ELMER_INFO("从配置文件解析参数完成");
    return params;
}

MagnetoDynamicsDimension 
MagnetoDynamicsSolverFactory::parseDimensionFromConfig(const std::string& configFile) {
    // 从配置文件名或内容推断维度
    std::string filename = configFile;
    
    // 从文件名推断
    if (filename.find("2D") != std::string::npos || filename.find("2d") != std::string::npos) {
        return MagnetoDynamicsDimension::DIM_2D;
    } else if (filename.find("3D") != std::string::npos || filename.find("3d") != std::string::npos) {
        return MagnetoDynamicsDimension::DIM_3D;
    }
    
    // 默认返回2D
    ELMER_WARN("警告: 无法从配置文件名推断维度，默认使用2D");
    return MagnetoDynamicsDimension::DIM_2D;
}

CoordinateSystemType 
MagnetoDynamicsSolverFactory::parseCoordinateSystemFromConfig(const std::string& configFile) {
    // 从配置文件名或内容推断坐标系
    std::string filename = configFile;
    
    // 从文件名推断
    if (filename.find("axisymmetric") != std::string::npos || 
        filename.find("axisym") != std::string::npos) {
        return CoordinateSystemType::AXISYMMETRIC;
    } else if (filename.find("cylindric") != std::string::npos || 
               filename.find("cyl") != std::string::npos) {
        return CoordinateSystemType::CYLINDRIC_SYMMETRIC;
    }
    
    // 默认返回笛卡尔坐标系
    return CoordinateSystemType::CARTESIAN;
}

} // namespace elmer