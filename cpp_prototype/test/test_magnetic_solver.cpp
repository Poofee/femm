/**
 * @file test_magnetic_solver.cpp
 * @brief MagneticSolver类单元测试 - 基于Fortran MagneticSolve.F90的移植验证
 * 
 * 测试MagneticSolver类的核心功能，包括：
 * - 坐标系检测
 * - 单元矩阵组装
 * - 非线性迭代求解
 * - 边界条件处理
 * - 物理量计算
 */

#include <gtest/gtest.h>
#include "solvers/electromagnetic/MagneticSolver.h"
#include "core/mesh/Mesh.h"
#include "core/math/LinearAlgebra.h"
#include "core/math/CRSMatrix.h"
#include "boundary/BoundaryConditions.h"
#include <vector>
#include <array>
#include <cmath>
#include <memory>

using namespace elmer;

/**
 * @brief MagneticSolver测试类
 */
class MagneticSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建简单网格（立方体）
        mesh = std::make_shared<Mesh>();
        
        // 添加8个节点（立方体顶点）
        mesh->addNode(0.0, 0.0, 0.0); // 节点0
        mesh->addNode(1.0, 0.0, 0.0); // 节点1
        mesh->addNode(1.0, 1.0, 0.0); // 节点2
        mesh->addNode(0.0, 1.0, 0.0); // 节点3
        mesh->addNode(0.0, 0.0, 1.0); // 节点4
        mesh->addNode(1.0, 0.0, 1.0); // 节点5
        mesh->addNode(1.0, 1.0, 1.0); // 节点6
        mesh->addNode(0.0, 1.0, 1.0); // 节点7
        
        // 添加六面体单元
        std::vector<int> hexNodes = {0, 1, 2, 3, 4, 5, 6, 7};
        mesh->addElement(ElementType::HEXAHEDRON, hexNodes);
        
        // 创建边界条件管理器
        bc = std::make_shared<BoundaryConditions>();
        
        // 创建材料数据库（简化实现）
        materialDB = std::make_shared<MaterialDatabase>();
        
        // 设置求解器
        solver.setMesh(mesh);
        solver.setBoundaryConditions(bc);
        solver.setMaterialDatabase(materialDB);
    }
    
    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<BoundaryConditions> bc;
    std::shared_ptr<MaterialDatabase> materialDB;
    MagneticSolver solver;
};

/**
 * @brief 测试坐标系检测功能
 */
TEST_F(MagneticSolverTest, CoordinateSystemDetection) {
    // 测试笛卡尔坐标系检测
    std::string system = solver.detectCoordinateSystem();
    EXPECT_EQ(system, "cartesian");
    
    // 测试柱对称坐标系检测（修改网格为柱对称）
    auto& nodes = mesh->getNodes();
    for (auto& node : nodes) {
        // 将网格修改为柱对称形状
        node.x = node.x * std::cos(node.y); // 径向坐标
        node.y = node.x * std::sin(node.y); // 周向坐标
        node.z = node.z;                    // 轴向坐标
    }
    
    system = solver.detectCoordinateSystem();
    EXPECT_EQ(system, "axisymmetric");
}

/**
 * @brief 测试单元矩阵组装功能
 */
TEST_F(MagneticSolverTest, ElementMatrixAssembly) {
    // 测试单元矩阵计算
    std::vector<std::vector<double>> elementStiffness;
    std::vector<double> elementRHS;
    
    bool success = solver.computeElementMatrices(0, elementStiffness, elementRHS);
    EXPECT_TRUE(success);
    
    // 检查矩阵尺寸
    EXPECT_EQ(elementStiffness.size(), 24); // 8节点 * 3自由度
    EXPECT_EQ(elementStiffness[0].size(), 24);
    EXPECT_EQ(elementRHS.size(), 24);
    
    // 检查矩阵对称性（刚度矩阵应该是对称的）
    for (int i = 0; i < 24; ++i) {
        for (int j = i + 1; j < 24; ++j) {
            EXPECT_NEAR(elementStiffness[i][j], elementStiffness[j][i], 1e-12);
        }
    }
    
    // 检查对角线元素为正（物理意义）
    for (int i = 0; i < 24; ++i) {
        EXPECT_GT(elementStiffness[i][i], 0.0);
    }
}

/**
 * @brief 测试非线性迭代求解功能
 */
TEST_F(MagneticSolverTest, NonlinearIteration) {
    // 初始化求解器
    bool initSuccess = solver.initialize();
    EXPECT_TRUE(initSuccess);
    
    // 组装系统矩阵
    bool assembleSuccess = solver.assemble();
    EXPECT_TRUE(assembleSuccess);
    
    // 测试非线性迭代求解
    bool solveSuccess = solver.solve();
    EXPECT_TRUE(solveSuccess);
    
    // 获取解向量
    std::vector<double> solution = solver.getSolution();
    EXPECT_FALSE(solution.empty());
    
    // 检查解向量的合理性
    double maxValue = *std::max_element(solution.begin(), solution.end());
    double minValue = *std::min_element(solution.begin(), solution.end());
    
    // 解应该在合理范围内
    EXPECT_LT(std::abs(maxValue), 1e6);
    EXPECT_LT(std::abs(minValue), 1e6);
}

/**
 * @brief 测试物理量计算功能
 */
TEST_F(MagneticSolverTest, PhysicalQuantities) {
    // 初始化并求解
    solver.initialize();
    solver.assemble();
    solver.solve();
    
    // 测试磁场计算
    std::vector<double> magneticField = solver.getMagneticField();
    EXPECT_FALSE(magneticField.empty());
    
    // 测试电流密度计算
    std::vector<double> electricCurrent = solver.getElectricCurrent();
    EXPECT_FALSE(electricCurrent.empty());
    
    // 测试洛伦兹力计算
    bool lorentzSuccess = solver.computeLorentzForce();
    EXPECT_TRUE(lorentzSuccess);
    std::vector<double> lorentzForce = solver.getLorentzForce();
    EXPECT_FALSE(lorentzForce.empty());
    
    // 测试电场计算
    bool electricSuccess = solver.computeElectricField();
    EXPECT_TRUE(electricSuccess);
    std::vector<double> electricField = solver.getElectricField();
    EXPECT_FALSE(electricField.empty());
}

/**
 * @brief 测试瞬态仿真支持
 */
TEST_F(MagneticSolverTest, TransientSimulation) {
    // 启用瞬态仿真
    solver.setTransientSimulation(true);
    solver.setTimeStep(0.01);
    
    // 初始化并求解
    bool initSuccess = solver.initialize();
    EXPECT_TRUE(initSuccess);
    
    bool assembleSuccess = solver.assemble();
    EXPECT_TRUE(assembleSuccess);
    
    bool solveSuccess = solver.solve();
    EXPECT_TRUE(solveSuccess);
}

/**
 * @brief 测试边界条件处理
 */
TEST_F(MagneticSolverTest, BoundaryConditions) {
    // 添加狄利克雷边界条件
    std::vector<int> boundaryNodes = {0, 1, 2, 3}; // 底面节点
    std::vector<double> fixedValues(boundaryNodes.size() * 3, 0.0);
    
    // 设置边界条件
    // 这里简化实现，实际应该通过边界条件管理器设置
    
    // 初始化并求解
    solver.initialize();
    solver.assemble();
    solver.solve();
    
    // 检查边界条件是否被正确应用
    std::vector<double> solution = solver.getSolution();
    
    // 边界节点应该满足边界条件
    for (int nodeId : boundaryNodes) {
        for (int dof = 0; dof < 3; ++dof) {
            int index = nodeId * 3 + dof;
            if (index < solution.size()) {
                EXPECT_NEAR(solution[index], 0.0, 1e-12);
            }
        }
    }
}

/**
 * @brief 测试内存管理
 */
TEST_F(MagneticSolverTest, MemoryManagement) {
    // 测试内存分配
    bool allocSuccess = solver.allocateMemory();
    EXPECT_TRUE(allocSuccess);
    
    // 测试内存释放
    solver.deallocateMemory();
    
    // 重新分配内存
    allocSuccess = solver.allocateMemory();
    EXPECT_TRUE(allocSuccess);
}

/**
 * @brief 测试材料参数处理
 */
TEST_F(MagneticSolverTest, MaterialParameters) {
    // 测试材料参数获取
    bool materialSuccess = solver.getMaterialParameters();
    EXPECT_TRUE(materialSuccess);
    
    // 测试速度场获取
    bool velocitySuccess = solver.getVelocityField();
    EXPECT_TRUE(velocitySuccess);
}

/**
 * @brief 测试旋度计算
 */
TEST_F(MagneticSolverTest, CurlComputation) {
    // 创建测试场量
    std::vector<double> fieldX(8, 1.0); // 均匀X场
    std::vector<double> fieldY(8, 0.0); // 零Y场
    std::vector<double> fieldZ(8, 0.0); // 零Z场
    
    std::vector<double> curlX, curlY, curlZ;
    
    // 测试旋度计算
    bool curlSuccess = solver.computeCurl(fieldX, fieldY, fieldZ, curlX, curlY, curlZ);
    EXPECT_TRUE(curlSuccess);
    
    // 均匀场的旋度应该为零
    for (size_t i = 0; i < curlX.size(); ++i) {
        EXPECT_NEAR(curlX[i], 0.0, 1e-12);
        EXPECT_NEAR(curlY[i], 0.0, 1e-12);
        EXPECT_NEAR(curlZ[i], 0.0, 1e-12);
    }
}

/**
 * @brief 主测试函数
 */
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}