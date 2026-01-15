#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include "MultiphysicsSolver.h"
#include "MagnetodynamicsSolver.h"
#include "HeatTransferSolver.h"
#include "StructuralSolver.h"
#include "Mesh.h"
#include "Material.h"
#include "LinearAlgebra.h"
#include "Types.h"

using namespace elmer;
using namespace ElmerCpp;

/**
 * @brief 创建简单的测试网格用于多物理场模拟
 */
std::shared_ptr<Mesh> createSimpleTestMesh() {
    auto mesh = std::make_shared<Mesh>();
    
    // 创建4个节点
    std::vector<Node> nodes = {
        Node(0.0, 0.0, 0.0),
        Node(1.0, 0.0, 0.0),
        Node(1.0, 1.0, 0.0),
        Node(0.0, 1.0, 0.0)
    };
    
    // 创建单个四边形单元
    std::vector<size_t> elementNodes = {0, 1, 2, 3};
    Element element(ElementType::QUADRATIC, 0);
    element.setNodeIndices(elementNodes);
    
    // 添加节点到网格
    for (const auto& node : nodes) {
        mesh->getNodes().addNode(node);
    }
    
    // 添加单元到网格
    mesh->addBulkElement(element);
    
    return mesh;
}

/**
 * @brief 创建材料数据库
 */
elmer::MaterialDatabase createSimpleMaterialDatabase() {
    elmer::MaterialDatabase matDB;
    
    // 铜材料属性
    Material copper;
    copper.name = "Copper";
    copper.density = 8960.0;           // kg/m³
    copper.specificHeat = 385.0;       // J/(kg·K)
    copper.thermalConductivity = 400.0; // W/(m·K)
    copper.conductivity = 5.96e7;      // S/m
    copper.relativePermeability = 1.0;
    copper.youngsModulus = 110e9;      // Pa
    copper.poissonsRatio = 0.34;
    copper.thermalExpansion = 16.5e-6; // 1/K
    
    matDB.addMaterial(copper.name, copper);
    
    return matDB;
}

/**
 * @brief 测试单个物理场求解器
 */
void testSinglePhysicsSolvers() {
    std::cout << "=== 测试单个物理场求解器 ===" << std::endl;
    
    auto mesh = createSimpleTestMesh();
    auto matDB = createSimpleMaterialDatabase();
    
    // 测试电磁求解器
    {
        std::cout << "\n1. 测试电磁求解器" << std::endl;
        
        MagnetodynamicsSolver emSolver(mesh);
        emSolver.setMaterialDatabase(matDB);
        
        // 设置求解器参数
        MagnetodynamicsParameters params;
        params.tolerance = 1.0e-8;
        params.maxIterations = 100;
        emSolver.setParameters(params);
        
        // 求解电磁问题
        auto results = emSolver.solve();
        
        std::cout << "  ✓ 电磁求解器在 " << results.iterations << " 次迭代后收敛" << std::endl;
        std::cout << "  ✓ 最终残差: " << results.residual << std::endl;
    }
    
    // 测试热求解器
    {
        std::cout << "\n2. 测试热求解器" << std::endl;
        
        HeatTransferSolver thermalSolver(mesh);
        thermalSolver.setMaterialDatabase(matDB);
        
        // 设置求解器参数
        HeatTransferParameters params;
        params.tolerance = 1.0e-8;
        params.maxIterations = 100;
        thermalSolver.setParameters(params);
        
        // 求解热问题
        auto results = thermalSolver.solve();
        
        std::cout << "  ✓ 热求解器在 " << results.iterations << " 次迭代后收敛" << std::endl;
        std::cout << "  ✓ 最终残差: " << results.residual << std::endl;
    }
    
    // 测试结构求解器
    {
        std::cout << "\n3. 测试结构求解器" << std::endl;
        
        StructuralSolver structuralSolver(mesh);
        structuralSolver.setMaterialDatabase(matDB);
        
        // 设置求解器参数
        StructuralParameters params;
        params.tolerance = 1.0e-8;
        params.maxIterations = 100;
        structuralSolver.setParameters(params);
        
        // 求解结构问题
        auto results = structuralSolver.solve();
        
        std::cout << "  ✓ 结构求解器在 " << results.iterations << " 次迭代后收敛" << std::endl;
        std::cout << "  ✓ 最终残差: " << results.residual << std::endl;
    }
}

/**
 * @brief 测试弱耦合多物理场求解器
 */
void testWeakCouplingMultiphysics() {
    std::cout << "\n=== 测试弱耦合多物理场求解器 ===" << std::endl;
    
    auto mesh = createSimpleTestMesh();
    auto matDB = createSimpleMaterialDatabase();
    
    // 创建多物理场求解器
    MultiphysicsSolver multiphysicsSolver;
    multiphysicsSolver.setMesh(mesh);
    multiphysicsSolver.setMaterialDatabase(matDB);
    multiphysicsSolver.setCouplingType(CouplingType::WEAK);
    
    // 设置耦合参数
    CouplingParameters params;
    params.couplingTolerance = 1.0e-6;
    params.maxCouplingIterations = 3;
    params.enableElectroThermalCoupling = true;
    params.enableElectroStructuralCoupling = false;
    params.enableThermoStructuralCoupling = false;
    multiphysicsSolver.setCouplingParameters(params);
    
    // 使用弱耦合求解
    multiphysicsSolver.solveWeakCoupling();
    
    // 获取解场
    auto potentialField = multiphysicsSolver.getPotentialField();
    auto temperatureField = multiphysicsSolver.getTemperatureField();
    
    std::cout << "✓ 电磁-热耦合测试完成" << std::endl;
    std::cout << "  电势场大小: " << potentialField.size() << std::endl;
    std::cout << "  温度场大小: " << temperatureField.size() << std::endl;
}

/**
 * @brief 测试迭代耦合多物理场求解器
 */
void testIterativeCouplingMultiphysics() {
    std::cout << "\n=== 测试迭代耦合多物理场求解器 ===" << std::endl;
    
    auto mesh = createSimpleTestMesh();
    auto matDB = createSimpleMaterialDatabase();
    
    // 创建多物理场求解器
    MultiphysicsSolver multiphysicsSolver;
    multiphysicsSolver.setMesh(mesh);
    multiphysicsSolver.setMaterialDatabase(matDB);
    multiphysicsSolver.setCouplingType(CouplingType::ITERATIVE);
    
    // 设置耦合参数
    CouplingParameters params;
    params.couplingTolerance = 1.0e-6;
    params.maxCouplingIterations = 5;
    multiphysicsSolver.setCouplingParameters(params);
    
    // 使用迭代耦合求解
    multiphysicsSolver.solveIterativeCoupling();
    
    // 获取解场
    auto potentialField = multiphysicsSolver.getPotentialField();
    auto temperatureField = multiphysicsSolver.getTemperatureField();
    auto displacementField = multiphysicsSolver.getDisplacementField();
    
    std::cout << "✓ 迭代耦合多物理场求解器测试完成" << std::endl;
    std::cout << "  电势场大小: " << potentialField.size() << std::endl;
    std::cout << "  温度场大小: " << temperatureField.size() << std::endl;
    std::cout << "  位移场大小: " << displacementField.size() << std::endl;
}

/**
 * @brief 主函数 - 运行所有测试
 */
int main() {
    std::cout << "开始多物理场求解器测试..." << std::endl;
    
    try {
        // 测试单个物理场求解器
        testSinglePhysicsSolvers();
        
        // 测试弱耦合多物理场求解器
        testWeakCouplingMultiphysics();
        
        // 测试迭代耦合多物理场求解器
        testIterativeCouplingMultiphysics();
        
        std::cout << "\n✓ 所有多物理场求解器测试完成" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}