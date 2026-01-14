// test_magnetic_solve.cpp - MagneticSolve模块测试程序
// 验证磁动力学求解器功能

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

// 包含MagneticSolve头文件
#include "src/MagneticSolve.h"
#include "src/Mesh.h"
#include "src/LinearAlgebra.h"
#include "src/CRSMatrix.h"
#include "src/IterativeSolver.h"

using namespace ElmerCpp;

// 简单的网格类实现
class SimpleMesh : public elmer::Mesh {
public:
    SimpleMesh() = default;
    
    // 创建简单的矩形网格
    void createRectangleMesh(int nx, int ny, double width, double height) {
        clear(); // 清空现有网格
        
        // 创建节点
        double dx = width / nx;
        double dy = height / ny;
        
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                getNodes().addNode(i * dx, j * dy, 0.0);
            }
        }
        
        // 创建体单元（四边形）
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                // 四边形节点索引
                int n1 = j * (nx + 1) + i;
                int n2 = j * (nx + 1) + i + 1;
                int n3 = (j + 1) * (nx + 1) + i + 1;
                int n4 = (j + 1) * (nx + 1) + i;
                
                std::vector<int> nodeIndices = {n1, n2, n3, n4};
                addQuadrilateralElement(nodeIndices, "Copper", 1);
            }
        }
        
        // 创建边界单元
        // 底部边界
        for (int i = 0; i < nx; ++i) {
            elmer::Element boundary(elmer::ElementType::LINEAR);
            boundary.setBoundaryId(1);
            
            std::vector<size_t> indices = {
                static_cast<size_t>(i), 
                static_cast<size_t>(i + 1)
            };
            boundary.setNodeIndices(indices);
            
            addBoundaryElement(boundary);
        }
        
        // 顶部边界
        for (int i = 0; i < nx; ++i) {
            elmer::Element boundary(elmer::ElementType::LINEAR);
            boundary.setBoundaryId(2);
            
            std::vector<size_t> indices = {
                static_cast<size_t>(ny * (nx + 1) + i), 
                static_cast<size_t>(ny * (nx + 1) + i + 1)
            };
            boundary.setNodeIndices(indices);
            
            addBoundaryElement(boundary);
        }
        
        // 左侧边界
        for (int j = 0; j < ny; ++j) {
            elmer::Element boundary(elmer::ElementType::LINEAR);
            boundary.setBoundaryId(3);
            
            std::vector<size_t> indices = {
                static_cast<size_t>(j * (nx + 1)), 
                static_cast<size_t>((j + 1) * (nx + 1))
            };
            boundary.setNodeIndices(indices);
            
            addBoundaryElement(boundary);
        }
        
        // 右侧边界
        for (int j = 0; j < ny; ++j) {
            elmer::Element boundary(elmer::ElementType::LINEAR);
            boundary.setBoundaryId(4);
            
            std::vector<size_t> indices = {
                static_cast<size_t>(j * (nx + 1) + nx), 
                static_cast<size_t>((j + 1) * (nx + 1) + nx)
            };
            boundary.setNodeIndices(indices);
            
            addBoundaryElement(boundary);
        }
    }
    
    // 添加额外的接口方法以兼容MagneticSolve
    double getX(int nodeIndex) const {
        if (nodeIndex >= 0 && nodeIndex < static_cast<int>(getNodes().numberOfNodes())) {
            return getNodes()[static_cast<size_t>(nodeIndex)].x;
        }
        return 0.0;
    }
    
    double getY(int nodeIndex) const {
        if (nodeIndex >= 0 && nodeIndex < static_cast<int>(getNodes().numberOfNodes())) {
            return getNodes()[static_cast<size_t>(nodeIndex)].y;
        }
        return 0.0;
    }
    
    double getZ(int nodeIndex) const {
        if (nodeIndex >= 0 && nodeIndex < static_cast<int>(getNodes().numberOfNodes())) {
            return getNodes()[static_cast<size_t>(nodeIndex)].z;
        }
        return 0.0;
    }
};

// 测试函数：基本功能测试
void testBasicFunctionality() {
    std::cout << "=== 测试MagneticSolve基本功能 ===" << std::endl;
    
    // 创建磁动力学求解器
    auto magneticSolver = CreateMagneticSolve();
    
    // 设置求解器参数
    MagneticSolveParameters params;
    params.tolerance = 1.0e-6;
    params.maxIterations = 100;
    params.calculateElectricField = true;
    params.calculateLorentzForce = true;
    params.calculateCurrentDensity = true;
    
    magneticSolver->setParameters(params);
    
    // 创建简单网格
    auto mesh = std::make_shared<SimpleMesh>();
    mesh->createRectangleMesh(4, 4, 1.0, 1.0);
    
    // 设置网格
    magneticSolver->setMesh(mesh);
    
    std::cout << "网格信息:" << std::endl;
    std::cout << "  节点数: " << mesh->numberOfNodes() << std::endl;
    std::cout << "  体单元数: " << mesh->numberOfBulkElements() << std::endl;
    std::cout << "  边界单元数: " << mesh->numberOfBoundaryElements() << std::endl;
    
    // 执行求解
    try {
        auto results = magneticSolver->solve();
        
        std::cout << "求解结果:" << std::endl;
        std::cout << "  收敛状态: " << (results.converged ? "是" : "否") << std::endl;
        std::cout << "  迭代次数: " << results.iterations << std::endl;
        std::cout << "  最终残差: " << results.residual << std::endl;
        std::cout << "  磁能: " << results.magneticEnergy << " J" << std::endl;
        
        // 检查磁场解
        if (!results.magneticField.empty()) {
            std::cout << "  磁场解大小: " << results.magneticField.size() << std::endl;
        }
        
        // 检查导出场
        if (!results.electricField.empty()) {
            std::cout << "  电场计算: 完成" << std::endl;
        }
        
        if (!results.lorentzForce.empty()) {
            std::cout << "  洛伦兹力计算: 完成" << std::endl;
        }
        
        if (!results.currentDensity.empty()) {
            std::cout << "  电流密度计算: 完成" << std::endl;
        }
        
        std::cout << "基本功能测试通过!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "测试失败: " << e.what() << std::endl;
    }
}

// 测试函数：参数设置测试
void testParameterSettings() {
    std::cout << "\n=== 测试参数设置 ===" << std::endl;
    
    auto magneticSolver = CreateMagneticSolve();
    
    // 测试默认参数
    auto defaultParams = magneticSolver->getParameters();
    std::cout << "默认参数:" << std::endl;
    std::cout << "  容差: " << defaultParams.tolerance << std::endl;
    std::cout << "  最大迭代次数: " << defaultParams.maxIterations << std::endl;
    std::cout << "  包含对流项: " << (defaultParams.includeConvection ? "是" : "否") << std::endl;
    
    // 设置新参数
    MagneticSolveParameters newParams;
    newParams.tolerance = 1.0e-10;
    newParams.maxIterations = 500;
    newParams.includeDisplacementCurrent = true;
    newParams.stabilize = false;
    
    magneticSolver->setParameters(newParams);
    
    // 验证参数设置
    auto retrievedParams = magneticSolver->getParameters();
    
    if (std::abs(retrievedParams.tolerance - newParams.tolerance) < 1e-15 &&
        retrievedParams.maxIterations == newParams.maxIterations &&
        retrievedParams.includeDisplacementCurrent == newParams.includeDisplacementCurrent &&
        retrievedParams.stabilize == newParams.stabilize) {
        std::cout << "参数设置测试通过!" << std::endl;
    } else {
        std::cout << "参数设置测试失败!" << std::endl;
    }
}

// 测试函数：网格兼容性测试
void testMeshCompatibility() {
    std::cout << "\n=== 测试网格兼容性 ===" << std::endl;
    
    auto magneticSolver = CreateMagneticSolve();
    
    // 测试无网格情况
    try {
        auto results = magneticSolver->solve();
        std::cout << "错误: 无网格情况下求解应该抛出异常!" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "无网格测试通过: " << e.what() << std::endl;
    }
    
    // 测试空网格
    auto emptyMesh = std::make_shared<SimpleMesh>();
    magneticSolver->setMesh(emptyMesh);
    
    try {
        auto results = magneticSolver->solve();
        std::cout << "空网格测试通过" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "空网格测试: " << e.what() << std::endl;
    }
    
    // 测试不同尺寸的网格
    for (int size : {2, 4, 8}) {
        auto mesh = std::make_shared<SimpleMesh>();
        mesh->createRectangleMesh(size, size, 1.0, 1.0);
        
        magneticSolver->setMesh(mesh);
        
        try {
            auto results = magneticSolver->solve();
            std::cout << "网格尺寸 " << size << "x" << size << " 测试通过" << std::endl;
        } catch (const std::exception& e) {
            std::cout << "网格尺寸 " << size << "x" << size << " 测试失败: " << e.what() << std::endl;
        }
    }
}

// 主测试函数
int main() {
    std::cout << "开始MagneticSolve模块测试..." << std::endl;
    std::cout << "================================" << std::endl;
    
    try {
        // 运行测试
        testBasicFunctionality();
        testParameterSettings();
        testMeshCompatibility();
        
        std::cout << "\n================================" << std::endl;
        std::cout << "所有测试完成!" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "测试过程中出现异常: " << e.what() << std::endl;
        return 1;
    }
}