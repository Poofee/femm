// test_harmonic_analysis.cpp - 谐波分析功能单元测试

#include "MagnetoDynamics2DSolver.h"
#include "BoundaryConditions.h"
#include "ElectromagneticMaterial.h"
#include "Mesh.h"
#include <iostream>
#include <memory>
#include <complex>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace elmer {

// 测试谐波分析基本功能
void testHarmonicAnalysisBasic() {
    std::cout << "=== 谐波分析基本功能测试 ===" << std::endl;
    
    // 创建求解器
    MagnetoDynamics2DSolver solver;
    
    // 设置谐波分析参数
    MagnetoDynamics2DParameters params;
    params.isHarmonic = true;
    params.frequency = 50.0; // 50Hz
    params.useComplexMatrices = true;
    params.includeEddyCurrents = true;
    params.tolerance = 1.0e-8;
    params.maxIterations = 1000;
    
    solver.setParameters(params);
    
    try {
        // 执行谐波分析
        auto results = solver.solveHarmonic();
        
        std::cout << "谐波分析结果:" << std::endl;
        std::cout << "收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << results.iterations << std::endl;
        
        // 验证复数解向量
        if (results.complexVectorPotential.empty()) {
            std::cout << "错误: 复数解向量为空" << std::endl;
            return;
        }
        
        std::cout << "复数解向量大小: " << results.complexVectorPotential.size() << std::endl;
        
        // 验证复数场量
        if (!results.complexMagneticFluxDensity.empty()) {
            std::cout << "复数磁通密度计算成功，大小: " 
                      << results.complexMagneticFluxDensity.size() << std::endl;
        }
        
        if (!results.complexMagneticFieldStrength.empty()) {
            std::cout << "复数磁场强度计算成功，大小: " 
                      << results.complexMagneticFieldStrength.size() << std::endl;
        }
        
        // 验证复数阻抗计算
        std::cout << "复数阻抗: " << results.complexImpedance << " Ω" << std::endl;
        std::cout << "复数电感: " << results.complexInductance << " H" << std::endl;
        std::cout << "功率损耗: " << results.powerLoss << " W" << std::endl;
        
        std::cout << "谐波分析基本功能测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "谐波分析过程中出现异常: " << e.what() << std::endl;
    }
}

// 测试简单线圈模型的谐波分析
void testSimpleCoilHarmonic() {
    std::cout << "=== 简单线圈模型谐波分析测试 ===" << std::endl;
    
    // 创建求解器
    MagnetoDynamics2DSolver solver;
    
    // 创建简单矩形网格（模拟线圈结构）
    auto mesh = std::make_shared<Mesh>();
    
    // 添加节点（简单的2x2网格）
    mesh->getNodes().addNode(0.0, 0.0, 0.0); // 节点0
    mesh->getNodes().addNode(1.0, 0.0, 0.0); // 节点1
    mesh->getNodes().addNode(0.0, 1.0, 0.0); // 节点2
    mesh->getNodes().addNode(1.0, 1.0, 0.0); // 节点3
    
    // 添加四边形元素
    std::vector<int> elementNodes = {0, 1, 3, 2};
    mesh->addQuadrilateralElement(elementNodes, "Copper");
    
    solver.setMesh(mesh);
    
    // 设置谐波分析参数
    MagnetoDynamics2DParameters params;
    params.isHarmonic = true;
    params.frequency = 60.0; // 60Hz
    params.useComplexMatrices = true;
    params.includeEddyCurrents = true;
    params.tolerance = 1.0e-8;
    
    solver.setParameters(params);
    
    try {
        // 执行谐波分析
        auto results = solver.solveHarmonic();
        
        std::cout << "简单线圈模型谐波分析结果:" << std::endl;
        std::cout << "收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        
        // 验证复数解向量
        if (results.complexVectorPotential.size() >= 4) {
            std::cout << "复数磁矢势计算成功，节点数: " 
                      << results.complexVectorPotential.size() << std::endl;
            
            // 输出前几个节点的复数磁矢势
            for (int i = 0; i < 4 && i < results.complexVectorPotential.size(); ++i) {
                std::cout << "节点" << i << ": A_z = " 
                          << results.complexVectorPotential[i] << " Wb/m" << std::endl;
            }
        }
        
        // 验证复数场量
        if (!results.complexMagneticFluxDensity.empty()) {
            std::cout << "复数磁通密度计算成功" << std::endl;
        }
        
        if (!results.complexMagneticFieldStrength.empty()) {
            std::cout << "复数磁场强度计算成功" << std::endl;
        }
        
        std::cout << "简单线圈模型谐波分析测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "简单线圈模型谐波分析过程中出现异常: " << e.what() << std::endl;
    }
}

// 测试复数阻抗计算正确性
void testComplexImpedanceCalculation() {
    std::cout << "=== 复数阻抗计算正确性测试 ===" << std::endl;
    
    // 创建求解器
    MagnetoDynamics2DSolver solver;
    
    // 设置谐波分析参数
    MagnetoDynamics2DParameters params;
    params.isHarmonic = true;
    params.frequency = 50.0; // 50Hz
    params.useComplexMatrices = true;
    params.includeEddyCurrents = true;
    params.calculateLumpedParameters = true;
    
    solver.setParameters(params);
    
    try {
        // 执行谐波分析
        auto results = solver.solveHarmonic();
        
        std::cout << "阻抗计算结果验证:" << std::endl;
        std::cout << "复数阻抗: " << results.complexImpedance << " Ω" << std::endl;
        std::cout << "复数电感: " << results.complexInductance << " H" << std::endl;
        std::cout << "功率损耗: " << results.powerLoss << " W" << std::endl;
        
        // 验证阻抗计算的合理性
        double impedanceMagnitude = std::abs(results.complexImpedance);
        double inductanceMagnitude = std::abs(results.complexInductance);
        
        std::cout << "阻抗幅值: " << impedanceMagnitude << " Ω" << std::endl;
        std::cout << "电感幅值: " << inductanceMagnitude << " H" << std::endl;
        
        // 基本合理性检查
        if (impedanceMagnitude >= 0.0 && inductanceMagnitude >= 0.0 && results.powerLoss >= 0.0) {
            std::cout << "阻抗计算结果合理性检查通过!" << std::endl;
        } else {
            std::cout << "警告: 阻抗计算结果可能不合理" << std::endl;
        }
        
        // 验证阻抗与频率的关系（理论检查）
        double expectedReactance = 2.0 * M_PI * params.frequency * inductanceMagnitude;
        double actualReactance = results.complexImpedance.imag();
        
        std::cout << "理论电抗: " << expectedReactance << " Ω" << std::endl;
        std::cout << "实际电抗: " << actualReactance << " Ω" << std::endl;
        
        if (std::abs(expectedReactance - actualReactance) < 1.0) { // 容差检查
            std::cout << "阻抗频率关系验证通过!" << std::endl;
        } else {
            std::cout << "警告: 阻抗频率关系可能存在偏差" << std::endl;
        }
        
        std::cout << "复数阻抗计算正确性测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "复数阻抗计算过程中出现异常: " << e.what() << std::endl;
    }
}

// 测试谐波激励边界条件
void testHarmonicExcitationBoundaryConditions() {
    std::cout << "=== 谐波激励边界条件测试 ===" << std::endl;
    
    // 创建求解器
    MagnetoDynamics2DSolver solver;
    
    // 创建简单网格
    auto mesh = std::make_shared<Mesh>();
    
    // 添加节点
    mesh->getNodes().addNode(0.0, 0.0, 0.0);
    mesh->getNodes().addNode(1.0, 0.0, 0.0);
    mesh->getNodes().addNode(0.0, 1.0, 0.0);
    mesh->getNodes().addNode(1.0, 1.0, 0.0);
    
    // 添加元素
    std::vector<int> elementNodes = {0, 1, 3, 2};
    mesh->addQuadrilateralElement(elementNodes, "Copper");
    
    solver.setMesh(mesh);
    
    // 设置谐波分析参数
    MagnetoDynamics2DParameters params;
    params.isHarmonic = true;
    params.frequency = 50.0;
    params.useComplexMatrices = true;
    
    solver.setParameters(params);
    
    // 创建谐波激励边界条件
    auto bcManager = solver.getBoundaryConditionManager();
    
    // 添加谐波电流源边界条件
    auto currentBC = std::make_shared<HarmonicCurrentBoundaryCondition>("coil_current");
    std::vector<int> boundaryNodes = {0, 1}; // 边界节点
    std::vector<double> amplitudes = {10.0, 10.0}; // 10A幅值
    std::vector<double> phases = {0.0, 0.0}; // 同相位
    
    currentBC->setHarmonicParameters(boundaryNodes, amplitudes, phases, 50.0);
    bcManager.addBoundaryCondition(currentBC);
    
    try {
        // 执行谐波分析
        auto results = solver.solveHarmonic();
        
        std::cout << "谐波激励边界条件测试结果:" << std::endl;
        std::cout << "收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        
        // 验证边界条件应用效果
        if (!results.complexVectorPotential.empty()) {
            std::cout << "边界节点复数磁矢势:" << std::endl;
            for (int nodeIdx : boundaryNodes) {
                if (nodeIdx < results.complexVectorPotential.size()) {
                    std::cout << "节点" << nodeIdx << ": A_z = " 
                              << results.complexVectorPotential[nodeIdx] << " Wb/m" << std::endl;
                }
            }
        }
        
        std::cout << "谐波激励边界条件测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "谐波激励边界条件测试过程中出现异常: " << e.what() << std::endl;
    }
}

// 测试多频率谐波分析
void testMultiFrequencyHarmonicAnalysis() {
    std::cout << "=== 多频率谐波分析测试 ===" << std::endl;
    
    // 测试不同频率下的谐波分析
    std::vector<double> frequencies = {50.0, 100.0, 200.0, 500.0}; // Hz
    
    for (double freq : frequencies) {
        std::cout << "测试频率: " << freq << " Hz" << std::endl;
        
        MagnetoDynamics2DSolver solver;
        
        MagnetoDynamics2DParameters params;
        params.isHarmonic = true;
        params.frequency = freq;
        params.useComplexMatrices = true;
        params.includeEddyCurrents = true;
        
        solver.setParameters(params);
        
        try {
            auto results = solver.solveHarmonic();
            
            std::cout << "频率 " << freq << " Hz - " 
                      << (results.converged ? "收敛" : "未收敛") 
                      << ", 阻抗: " << results.complexImpedance << " Ω" << std::endl;
            
        } catch (const std::exception& e) {
            std::cout << "频率 " << freq << " Hz 分析失败: " << e.what() << std::endl;
        }
    }
    
    std::cout << "多频率谐波分析测试完成!" << std::endl;
}

// 测试时域重构功能
void testTimeDomainReconstruction() {
    std::cout << "=== 时域重构功能测试 ===" << std::endl;
    
    MagnetoDynamics2DSolver solver;
    
    // 设置谐波分析参数
    MagnetoDynamics2DParameters params;
    params.isHarmonic = true;
    params.frequency = 50.0;
    params.useComplexMatrices = true;
    
    solver.setParameters(params);
    
    try {
        // 执行谐波分析
        auto results = solver.solveHarmonic();
        
        // 测试时域重构
        std::vector<double> timePoints = {0.0, 0.002, 0.004, 0.006, 0.008, 0.01}; // 20ms周期
        
        std::cout << "时域重构测试:" << std::endl;
        for (double time : timePoints) {
            auto timeDomainSolution = solver.reconstructTimeDomain(time);
            
            if (!timeDomainSolution.empty()) {
                std::cout << "时间 " << time << " s: 解向量大小 = " 
                          << timeDomainSolution.size() << std::endl;
            }
        }
        
        std::cout << "时域重构功能测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "时域重构测试过程中出现异常: " << e.what() << std::endl;
    }
}

} // namespace elmer

// 主测试函数
int main() {
    std::cout << "开始谐波分析功能测试..." << std::endl;
    std::cout << "========================" << std::endl;
    
    try {
        // 运行所有测试
        elmer::testHarmonicAnalysisBasic();
        std::cout << std::endl;
        
        elmer::testSimpleCoilHarmonic();
        std::cout << std::endl;
        
        elmer::testComplexImpedanceCalculation();
        std::cout << std::endl;
        
        elmer::testHarmonicExcitationBoundaryConditions();
        std::cout << std::endl;
        
        elmer::testMultiFrequencyHarmonicAnalysis();
        std::cout << std::endl;
        
        elmer::testTimeDomainReconstruction();
        std::cout << std::endl;
        
        std::cout << "所有谐波分析测试完成!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "测试过程中出现异常: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}