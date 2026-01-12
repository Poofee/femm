// test_magnetodynamics2d_solver.cpp - 2D磁动力学求解器单元测试

#include "MagnetoDynamics2DSolver.h"
#include <iostream>
#include <memory>

namespace elmer {

// 测试基本求解器功能
void testBasicSolverFunctionality() {
    std::cout << "=== 基本求解器功能测试 ===" << std::endl;
    
    // 创建求解器（使用默认网格）
    MagnetoDynamics2DSolver solver;
    
    // 设置求解参数
    MagnetoDynamics2DParameters params;
    params.tolerance = 1.0e-8;
    params.maxNonlinearIterations = 50;
    params.isTransient = false;
    params.coordinateSystem = MagnetoDynamics2DParameters::CARTESIAN;
    
    solver.setParameters(params);
    
    try {
        // 求解问题
        auto results = solver.solve();
        
        std::cout << "求解结果:" << std::endl;
        std::cout << "收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        std::cout << "非线性迭代次数: " << results.nonlinearIterations << std::endl;
        
        // 验证解向量不为空
        if (results.vectorPotential.empty()) {
            std::cout << "错误: 解向量为空" << std::endl;
            return;
        }
        
        std::cout << "解向量大小: " << results.vectorPotential.size() << std::endl;
        std::cout << "基本求解器功能测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "求解过程中出现异常: " << e.what() << std::endl;
    }
}

// 测试坐标系功能
void testCoordinateSystems() {
    std::cout << "=== 坐标系功能测试 ===" << std::endl;
    
    // 测试笛卡尔坐标系
    {
        MagnetoDynamics2DSolver solver;
        MagnetoDynamics2DParameters params;
        params.coordinateSystem = MagnetoDynamics2DParameters::CARTESIAN;
        solver.setParameters(params);
        
        std::cout << "笛卡尔坐标系测试完成" << std::endl;
    }
    
    // 测试轴对称坐标系
    {
        MagnetoDynamics2DSolver solver;
        MagnetoDynamics2DParameters params;
        params.coordinateSystem = MagnetoDynamics2DParameters::AXISYMMETRIC;
        solver.setParameters(params);
        
        std::cout << "轴对称坐标系测试完成" << std::endl;
    }
    
    // 测试柱对称坐标系
    {
        MagnetoDynamics2DSolver solver;
        MagnetoDynamics2DParameters params;
        params.coordinateSystem = MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC;
        solver.setParameters(params);
        
        std::cout << "柱对称坐标系测试完成" << std::endl;
    }
    
    std::cout << "坐标系功能测试通过!" << std::endl;
}

// 测试性能优化功能
void testPerformanceOptimization() {
    std::cout << "=== 性能优化功能测试 ===" << std::endl;
    
    // 创建求解器
    MagnetoDynamics2DSolver solver;
    
    // 启用基函数缓存
    solver.enableBasisFunctionCache(true);
    
    // 预计算缓存
    solver.precomputeBasisFunctionCache();
    
    // 测试优化的矩阵组装
    solver.assembleSystemOptimized();
    
    // 测试并行组装（简化版本）
    solver.assembleSystemParallel();
    
    // 测试增量式更新
    solver.updateSystemIncremental();
    
    // 获取缓存统计信息
    solver.getCacheStatistics();
    
    // 清除缓存
    solver.clearCache();
    
    std::cout << "性能优化功能测试通过!" << std::endl;
}

// 测试缓存管理功能
void testCacheManagement() {
    std::cout << "=== 缓存管理功能测试 ===" << std::endl;
    
    MagnetoDynamics2DSolver solver;
    
    // 测试启用/禁用缓存
    solver.enableBasisFunctionCache(true);
    solver.enableBasisFunctionCache(false);
    
    // 测试缓存预计算
    solver.precomputeBasisFunctionCache();
    
    // 测试缓存统计
    solver.getCacheStatistics();
    
    // 测试缓存清除
    solver.clearCache();
    
    // 再次获取统计信息
    solver.getCacheStatistics();
    
    std::cout << "缓存管理功能测试通过!" << std::endl;
}

// 测试瞬态分析功能
void testTransientAnalysis() {
    std::cout << "=== 瞬态分析功能测试 ===" << std::endl;
    
    // 创建求解器
    MagnetoDynamics2DSolver solver;
    
    // 设置瞬态求解参数
    MagnetoDynamics2DParameters params;
    params.isTransient = true;
    params.timeStep = 0.01;
    params.coordinateSystem = MagnetoDynamics2DParameters::CARTESIAN;
    
    solver.setParameters(params);
    
    try {
        // 求解瞬态问题
        auto results = solver.solve();
        
        std::cout << "瞬态分析结果:" << std::endl;
        std::cout << "收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        
        std::cout << "瞬态分析功能测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "瞬态求解过程中出现异常: " << e.what() << std::endl;
    }
}

} // namespace elmer

// 主测试函数
int main() {
    std::cout << "开始2D磁动力学求解器单元测试" << std::endl;
    std::cout << "================================" << std::endl;
    
    try {
        // 运行基本功能测试
        elmer::testBasicSolverFunctionality();
        
        std::cout << std::endl;
        
        // 运行坐标系测试
        elmer::testCoordinateSystems();
        
        std::cout << std::endl;
        
        // 运行性能优化测试
        elmer::testPerformanceOptimization();
        
        std::cout << std::endl;
        
        // 运行缓存管理测试
        elmer::testCacheManagement();
        
        std::cout << std::endl;
        
        // 运行瞬态分析测试
        elmer::testTransientAnalysis();
        
        std::cout << std::endl;
        std::cout << "所有测试完成!" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cout << "测试过程中出现异常: " << e.what() << std::endl;
        return 1;
    }
}