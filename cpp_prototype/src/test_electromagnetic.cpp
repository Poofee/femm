// test_electromagnetic.cpp - 电磁场求解器功能验证测试
#include "MagnetodynamicsSolver.h"
#include "ElectromagneticMaterial.h"
#include "Mesh.h"
#include "MeshIO.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace elmer;

void testMaterialProperties() {
    std::cout << "=== 电磁材料属性测试 ===" << std::endl;
    
    // 测试基本材料属性
    ElectromagneticMaterial copper(1.0, 1.0, 5.96e7, "Copper");
    
    std::cout << "铜材料属性:" << std::endl;
    std::cout << "相对介电常数: " << copper.relativePermittivity << std::endl;
    std::cout << "相对磁导率: " << copper.relativePermeability << std::endl;
    std::cout << "电导率: " << copper.conductivity << " S/m" << std::endl;
    std::cout << "介电常数: " << copper.permittivity() << " F/m" << std::endl;
    std::cout << "磁导率: " << copper.permeability() << " H/m" << std::endl;
    
    // 测试频率相关属性
    copper.frequency = 50.0; // 50 Hz
    std::cout << "角频率: " << copper.angularFrequency() << " rad/s" << std::endl;
    std::cout << "趋肤深度: " << copper.skinDepth() << " m" << std::endl;
    std::cout << "波数: " << copper.waveNumber() << " rad/m" << std::endl;
    
    std::cout << std::endl;
}

void testMaterialDatabase() {
    std::cout << "=== 材料数据库测试 ===" << std::endl;
    
    MaterialDatabase db;
    db.createPredefinedMaterials();
    
    std::cout << "预定义材料:" << std::endl;
    auto materialNames = db.getMaterialNames();
    for (const auto& name : materialNames) {
        auto material = db.getMaterial(name);
        std::cout << "- " << name << ": ε_r=" << material.relativePermittivity 
                  << ", μ_r=" << material.relativePermeability 
                  << ", σ=" << material.conductivity << std::endl;
    }
    
    std::cout << std::endl;
}

void testSimpleElectromagneticProblem() {
    std::cout << "=== 简单电磁问题测试 ===" << std::endl;
    
    // 创建简单网格（1D问题）
    auto mesh = std::make_shared<Mesh>();
    
    // 添加节点
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    mesh->addNode(Node(2.0, 0.0, 0.0));
    
    // 添加单元
    std::vector<int> element1Nodes = {0, 1};
    std::vector<int> element2Nodes = {1, 2};
    mesh->addElement(Element(ElementType::LINEAR, element1Nodes, "Copper"));
    mesh->addElement(Element(ElementType::LINEAR, element2Nodes, "Copper"));
    
    // 创建求解器
    MagnetodynamicsSolver solver(mesh);
    
    // 设置求解参数
    MagnetodynamicsParameters params;
    params.tolerance = 1.0e-8;
    params.maxIterations = 100;
    params.isHarmonic = false;
    
    solver.setParameters(params);
    
    // 添加边界条件
    EMBoundaryCondition bc1(EM_BoundaryType::DIRICHLET, "LeftBoundary");
    bc1.nodeIndices = {0};
    bc1.values = {0.0}; // 零电位
    
    EMBoundaryCondition bc2(EM_BoundaryType::DIRICHLET, "RightBoundary");
    bc2.nodeIndices = {2};
    bc2.values = {10.0}; // 10V电位
    
    solver.addBoundaryCondition(bc1);
    solver.addBoundaryCondition(bc2);
    
    try {
        // 求解问题
        auto results = solver.solve();
        
        std::cout << "求解结果:" << std::endl;
        std::cout << "收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << results.iterations << std::endl;
        std::cout << "最终残差: " << results.residual << std::endl;
        
        std::cout << "节点电位:" << std::endl;
        for (size_t i = 0; i < results.potentialReal.size(); ++i) {
            std::cout << "节点 " << i << ": " << results.potentialReal[i] << " V" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "求解失败: " << e.what() << std::endl;
    }
    
    std::cout << std::endl;
}

void testHarmonicAnalysis() {
    std::cout << "=== 谐波分析测试 ===" << std::endl;
    
    // 创建简单网格
    auto mesh = std::make_shared<Mesh>();
    
    // 添加节点
    mesh->addNode(Node(0.0, 0.0, 0.0));
    mesh->addNode(Node(0.5, 0.0, 0.0));
    mesh->addNode(Node(1.0, 0.0, 0.0));
    
    // 添加单元
    std::vector<int> elementNodes = {0, 1, 2};
    mesh->addElement(Element(ElementType::LINEAR, elementNodes, "Copper"));
    
    // 创建求解器
    MagnetodynamicsSolver solver(mesh);
    
    // 设置谐波分析参数
    MagnetodynamicsParameters params;
    params.isHarmonic = true;
    params.frequency = 50.0; // 50 Hz
    params.tolerance = 1.0e-6;
    params.maxIterations = 200;
    
    solver.setParameters(params);
    
    // 添加边界条件
    EMBoundaryCondition bc1(EM_BoundaryType::DIRICHLET, "Source");
    bc1.nodeIndices = {0};
    bc1.values = {1.0};      // 实部 1V
    bc1.valuesImag = {0.0};  // 虚部 0V
    
    EMBoundaryCondition bc2(EM_BoundaryType::DIRICHLET, "Ground");
    bc2.nodeIndices = {2};
    bc2.values = {0.0};
    bc2.valuesImag = {0.0};
    
    solver.addBoundaryCondition(bc1);
    solver.addBoundaryCondition(bc2);
    
    try {
        auto results = solver.solve();
        
        std::cout << "谐波分析结果:" << std::endl;
        std::cout << "收敛: " << (results.converged ? "是" : "否") << std::endl;
        std::cout << "迭代: " << results.iterations << std::endl;
        
        if (!results.potentialImag.empty()) {
            std::cout << "复数电位:" << std::endl;
            for (size_t i = 0; i < results.potentialReal.size(); ++i) {
                std::cout << "节点 " << i << ": " 
                          << results.potentialReal[i] << " + j" 
                          << results.potentialImag[i] << " V" << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "谐波分析失败: " << e.what() << std::endl;
    }
    
    std::cout << std::endl;
}

void testNonlinearMaterial() {
    std::cout << "=== 非线性材料测试 ===" << std::endl;
    
    // 创建非线性材料（铁磁材料）
    ElectromagneticMaterial iron(1.0, 1000.0, 1.0e7, "Iron");
    
    // 设置B-H曲线
    std::vector<double> H_values = {0, 100, 500, 1000, 5000, 10000}; // A/m
    std::vector<double> B_values = {0, 0.1, 0.5, 1.0, 1.5, 1.6};     // T
    
    iron.setBHCurve(H_values, B_values);
    
    std::cout << "非线性铁材料测试:" << std::endl;
    std::cout << "H=100 A/m时的磁导率: " << iron.getNonlinearPermeability(100) << " H/m" << std::endl;
    std::cout << "H=1000 A/m时的磁导率: " << iron.getNonlinearPermeability(1000) << " H/m" << std::endl;
    std::cout << "H=5000 A/m时的磁导率: " << iron.getNonlinearPermeability(5000) << " H/m" << std::endl;
    
    std::cout << std::endl;
}

int main() {
    std::cout << "Elmer FEM C++ 电磁场求解器功能验证测试" << std::endl;
    std::cout << "======================================" << std::endl;
    
    try {
        testMaterialProperties();
        testMaterialDatabase();
        testSimpleElectromagneticProblem();
        testHarmonicAnalysis();
        testNonlinearMaterial();
        
        std::cout << "✅ 所有电磁场功能测试完成！" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}