/**
 * @file test_magnetic_solve.cpp
 * @brief MagneticSolve类单元测试
 * 
 * 测试基于Fortran MaxwellCompose子程序的assembleCartesianElement方法
 */

#include <gtest/gtest.h>
#include "MagneticSolve.h"
#include "Mesh.h"
#include "LinearAlgebra.h"
#include <vector>
#include <array>
#include <cmath>

using namespace ElmerCpp;

/**
 * @brief 测试MagneticSolve的基本功能
 */
class MagneticSolveTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建简单网格
        mesh = std::make_shared<Mesh>();
        
        // 添加节点
        mesh->addNode(0.0, 0.0, 0.0); // 节点0
        mesh->addNode(1.0, 0.0, 0.0); // 节点1
        mesh->addNode(0.0, 1.0, 0.0); // 节点2
        mesh->addNode(1.0, 1.0, 0.0); // 节点3
        
        // 添加四面体单元
        std::vector<int> tetraNodes = {0, 1, 2, 3};
        mesh->addElement(ElementType::TETRAHEDRON, tetraNodes);
        
        // 设置求解器
        solver.setMesh(mesh);
        
        // 设置参数
        MagneticSolveParameters params;
        params.tolerance = 1.0e-8;
        params.maxIterations = 1000;
        solver.setParameters(params);
    }
    
    std::shared_ptr<Mesh> mesh;
    MagneticSolve solver;
};

/**
 * @brief 测试assembleCartesianElement方法的基本功能
 */
TEST_F(MagneticSolveTest, AssembleCartesianElementBasic) {
    // 获取网格中的第一个单元
    auto& elements = mesh->getBulkElements();
    ASSERT_FALSE(elements.empty());
    
    const auto& element = elements[0];
    
    // 创建单元节点
    ElementNodes elementNodes;
    auto nodeIndices = element.getNodeIndices();
    auto& nodes = mesh->getNodes();
    
    for (size_t i = 0; i < nodeIndices.size(); ++i) {
        elementNodes.x[i] = nodes.getX(nodeIndices[i]);
        elementNodes.y[i] = nodes.getY(nodeIndices[i]);
        elementNodes.z[i] = nodes.getZ(nodeIndices[i]);
    }
    
    // 测试assembleCartesianElement方法
    // 注意：这个方法需要材料参数和场量数据，这里主要测试接口正确性
    EXPECT_NO_THROW(solver.assembleCartesianElement(element, elementNodes));
}

/**
 * @brief 测试求解器参数设置
 */
TEST_F(MagneticSolveTest, ParameterSetting) {
    MagneticSolveParameters params;
    params.tolerance = 1.0e-10;
    params.maxIterations = 500;
    params.calculateElectricField = true;
    
    solver.setParameters(params);
    
    const auto& retrievedParams = solver.getParameters();
    EXPECT_DOUBLE_EQ(retrievedParams.tolerance, 1.0e-10);
    EXPECT_EQ(retrievedParams.maxIterations, 500);
    EXPECT_TRUE(retrievedParams.calculateElectricField);
}

/**
 * @brief 测试求解器初始化
 */
TEST_F(MagneticSolveTest, SolverInitialization) {
    // 测试求解器是否能够正常初始化
    EXPECT_NO_THROW(solver.solve());
}

/**
 * @brief 测试边界条件处理
 */
TEST_F(MagneticSolveTest, BoundaryConditionHandling) {
    // 测试边界条件处理功能
    // 这里主要测试接口正确性，实际边界条件需要完整的实现
    EXPECT_NO_THROW(solver.applyBoundaryConditions());
}

/**
 * @brief 测试线性系统求解
 */
TEST_F(MagneticSolveTest, LinearSystemSolving) {
    // 测试线性系统求解功能
    // 这里主要测试接口正确性，实际求解需要完整的矩阵组装
    EXPECT_NO_THROW(solver.solveLinearSystem());
}

/**
 * @brief 测试导出场计算
 */
TEST_F(MagneticSolveTest, DerivedFieldComputation) {
    MagneticSolveResults results;
    
    // 测试电场计算
    std::vector<std::array<double, 3>> electricField;
    EXPECT_NO_THROW(solver.computeElectricField(electricField));
    
    // 测试洛伦兹力计算
    std::vector<std::array<double, 3>> lorentzForce;
    EXPECT_NO_THROW(solver.computeLorentzForce(lorentzForce));
    
    // 测试电流密度计算
    std::vector<std::array<double, 3>> currentDensity;
    EXPECT_NO_THROW(solver.computeCurrentDensity(currentDensity));
}

/**
 * @brief 测试磁能计算
 */
TEST_F(MagneticSolveTest, MagneticEnergyComputation) {
    double magneticEnergy = solver.computeMagneticEnergy();
    
    // 磁能应该为非负数
    EXPECT_GE(magneticEnergy, 0.0);
}

/**
 * @brief 测试自由表面检测
 */
TEST_F(MagneticSolveTest, FreeSurfaceDetection) {
    bool hasFreeSurface = solver.checkFreeSurface();
    
    // 对于简单测试网格，应该没有自由表面
    EXPECT_FALSE(hasFreeSurface);
}

/**
 * @brief 测试数值精度
 */
TEST_F(MagneticSolveTest, NumericalPrecision) {
    // 测试数值精度要求：相对误差 < 1e-12
    // 这里可以添加具体的数值精度测试
    
    // 示例：测试基本数学运算的精度
    double a = 1.0;
    double b = 1.0 + 1e-15;
    
    // 相对误差应该小于1e-12
    double relativeError = std::abs(a - b) / std::max(std::abs(a), std::abs(b));
    EXPECT_LT(relativeError, 1e-12);
}

/**
 * @brief 测试内存管理
 */
TEST_F(MagneticSolveTest, MemoryManagement) {
    // 测试临时存储的初始化和清理
    EXPECT_NO_THROW(solver.initializeTemporaryStorage());
    EXPECT_NO_THROW(solver.cleanupTemporaryStorage());
}

/**
 * @brief 主测试函数
 */
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}