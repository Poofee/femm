/**
 * @file test_input_parser.cpp
 * @brief 输入文件解析器测试程序
 */

#include "InputFileParser.h"
#include <iostream>
#include <string>

using namespace elmer;

int main() {
    std::cout << "=== Elmer输入文件解析器测试 ===" << std::endl;
    
    // 创建解析器实例
    InputFileParser parser;
    
    // 测试文件路径
    std::string testFile = "..\\..\\..\\tests\\1dtests\\1d.sif";
    
    std::cout << "测试文件: " << testFile << std::endl;
    
    // 解析输入文件
    if (!parser.parse(testFile)) {
        std::cerr << "错误: 文件解析失败" << std::endl;
        return -1;
    }
    
    // 验证输入文件
    if (!parser.validate()) {
        std::cerr << "错误: 文件验证失败" << std::endl;
        return -1;
    }
    
    // 打印解析摘要
    parser.printSummary();
    
    // 测试参数获取功能
    std::cout << "\n=== 参数获取测试 ===" << std::endl;
    
    // 获取仿真类型
    std::string simType = parser.getSimulationType();
    std::cout << "仿真类型: " << simType << std::endl;
    
    // 获取网格信息
    std::string meshDir = parser.getMeshDirectory();
    std::string meshFile = parser.getMeshFileName();
    std::cout << "网格目录: " << meshDir << std::endl;
    std::cout << "网格文件: " << meshFile << std::endl;
    
    // 获取具体参数值
    double maxOutputLevel = parser.getParameterReal("simulation", 0, "Max Output Level", 0.0);
    std::cout << "最大输出级别: " << maxOutputLevel << std::endl;
    
    std::string coordSystem = parser.getParameterValue("simulation", 0, "Coordinate System");
    std::cout << "坐标系: " << coordSystem << std::endl;
    
    // 获取激活的求解器
    auto activeSolvers = parser.getActiveSolvers();
    std::cout << "激活求解器: ";
    for (int solverId : activeSolvers) {
        std::cout << solverId << " ";
    }
    std::cout << std::endl;
    
    // 测试特定节段参数
    std::cout << "\n=== 节段参数测试 ===" << std::endl;
    
    // 获取Solver节段参数
    std::string solverEquation = parser.getParameterValue("solver", 1, "Equation");
    std::cout << "求解器1方程: " << solverEquation << std::endl;
    
    std::string solverVariable = parser.getParameterValue("solver", 1, "Variable");
    std::cout << "求解器1变量: " << solverVariable << std::endl;
    
    double tolerance = parser.getParameterReal("solver", 1, "Steady State Convergence Tolerance", 1e-6);
    std::cout << "收敛容差: " << tolerance << std::endl;
    
    // 获取边界条件参数
    std::string targetBoundaries = parser.getParameterValue("boundary condition", 1, "Target Boundaries");
    std::cout << "边界条件1目标边界: " << targetBoundaries << std::endl;
    
    double potentialValue = parser.getParameterReal("boundary condition", 1, "Potential", 0.0);
    std::cout << "边界条件1电势值: " << potentialValue << std::endl;
    
    std::cout << "\n=== 测试完成 ===" << std::endl;
    
    return 0;
}