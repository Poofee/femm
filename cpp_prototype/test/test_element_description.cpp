/**
 * @file test_element_description.cpp
 * @brief ElementInfo和EdgeElementInfo接口单元测试
 * @author Elmer FEM C++移植项目组
 * @date 2026-01-21
 */

#include <gtest/gtest.h>
#include "ElementDescription.h"
#include <vector>
#include <cmath>

using namespace elmer;

/**
 * @brief ElementInfo接口测试类
 */
class ElementInfoTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 设置测试用的元素和节点数据
        setupTriangleElement();
        setupQuadElement();
        setupTetraElement();
    }

    void setupTriangleElement() {
        // 三角形元素设置
        triangleElement.type.elementCode = 303; // 三角形元素代码
        triangleElement.type.numberOfNodes = 3;
        triangleElement.type.dimension = 2;
        
        triangleNodes.x = {0.0, 1.0, 0.0};
        triangleNodes.y = {0.0, 0.0, 1.0};
        triangleNodes.z = {0.0, 0.0, 0.0};
    }

    void setupQuadElement() {
        // 四边形元素设置
        quadElement.type.elementCode = 404; // 四边形元素代码
        quadElement.type.numberOfNodes = 4;
        quadElement.type.dimension = 2;
        
        quadNodes.x = {0.0, 1.0, 1.0, 0.0};
        quadNodes.y = {0.0, 0.0, 1.0, 1.0};
        quadNodes.z = {0.0, 0.0, 0.0, 0.0};
    }

    void setupTetraElement() {
        // 四面体元素设置
        tetraElement.type.elementCode = 504; // 四面体元素代码
        tetraElement.type.numberOfNodes = 4;
        tetraElement.type.dimension = 3;
        
        tetraNodes.x = {0.0, 1.0, 0.0, 0.0};
        tetraNodes.y = {0.0, 0.0, 1.0, 0.0};
        tetraNodes.z = {0.0, 0.0, 0.0, 1.0};
    }

    Element triangleElement;
    Element quadElement;
    Element tetraElement;
    Nodes triangleNodes;
    Nodes quadNodes;
    Nodes tetraNodes;
};

/**
 * @brief 测试三角形元素的ElementInfo计算
 */
TEST_F(ElementInfoTest, TriangleElementInfo) {
    // 测试参数
    double u = 0.25, v = 0.25, w = 0.0;
    double detJ = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> dBasisdx;
    std::vector<std::vector<std::vector<double>>> ddBasisddx;
    bool compute2ndDerivatives = false;

    // 调用ElementInfo函数
    bool result = ElementInfo(triangleElement, triangleNodes, u, v, w, detJ,
                             basis, &dBasisdx, &ddBasisddx, &compute2ndDerivatives);

    // 验证结果
    EXPECT_TRUE(result);
    EXPECT_GT(detJ, 0.0); // 雅可比行列式应为正
    
    // 验证基函数数量
    EXPECT_EQ(basis.size(), static_cast<size_t>(3));
    
    // 验证基函数值（应在[0,1]范围内）
    for (const auto& val : basis) {
        EXPECT_GE(val, 0.0);
        EXPECT_LE(val, 1.0);
    }
    
    // 验证基函数和为1
    double sum = 0.0;
    for (const auto& val : basis) {
        sum += val;
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);
    
    // 验证导数矩阵维度
    EXPECT_EQ(dBasisdx.size(), static_cast<size_t>(3));
    if (dBasisdx.size() > 0) {
        EXPECT_EQ(dBasisdx[0].size(), static_cast<size_t>(2)); // 2D导数
    }
}

/**
 * @brief 测试四边形元素的ElementInfo计算
 */
TEST_F(ElementInfoTest, QuadElementInfo) {
    // 测试参数
    double u = 0.0, v = 0.0, w = 0.0; // 中心点
    double detJ = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> dBasisdx;
    std::vector<std::vector<std::vector<double>>> ddBasisddx;
    bool compute2ndDerivatives = false;

    // 调用ElementInfo函数
    bool result = ElementInfo(quadElement, quadNodes, u, v, w, detJ,
                             basis, &dBasisdx, &ddBasisddx, &compute2ndDerivatives);

    // 验证结果
    EXPECT_TRUE(result);
    EXPECT_GT(detJ, 0.0);
    
    // 验证基函数数量
    EXPECT_EQ(basis.size(), 4);
    
    // 验证基函数值
    for (const auto& val : basis) {
        EXPECT_GE(val, 0.0);
        EXPECT_LE(val, 1.0);
    }
    
    // 验证基函数和为1
    double sum = 0.0;
    for (const auto& val : basis) {
        sum += val;
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);
    
    // 验证导数矩阵维度
    EXPECT_EQ(dBasisdx.size(), 4);
    if (dBasisdx.size() > 0) {
        EXPECT_EQ(dBasisdx[0].size(), 2);
    }
}

/**
 * @brief 测试四面体元素的ElementInfo计算
 */
TEST_F(ElementInfoTest, TetraElementInfo) {
    // 测试参数
    double u = 0.25, v = 0.25, w = 0.25; // 重心坐标
    double detJ = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> dBasisdx;
    std::vector<std::vector<std::vector<double>>> ddBasisddx;
    bool compute2ndDerivatives = false;

    // 调用ElementInfo函数
    bool result = ElementInfo(tetraElement, tetraNodes, u, v, w, detJ,
                             basis, &dBasisdx, &ddBasisddx, &compute2ndDerivatives);

    // 验证结果
    EXPECT_TRUE(result);
    EXPECT_GT(detJ, 0.0);
    
    // 验证基函数数量
    EXPECT_EQ(basis.size(), 4);
    
    // 验证基函数值
    for (const auto& val : basis) {
        EXPECT_GE(val, 0.0);
        EXPECT_LE(val, 1.0);
    }
    
    // 验证基函数和为1
    double sum = 0.0;
    for (const auto& val : basis) {
        sum += val;
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);
    
    // 验证导数矩阵维度
    EXPECT_EQ(dBasisdx.size(), 4);
    if (dBasisdx.size() > 0) {
        EXPECT_EQ(dBasisdx[0].size(), 3); // 3D导数
    }
}

/**
 * @brief 测试ElementInfo在边界点的情况
 */
TEST_F(ElementInfoTest, BoundaryPoints) {
    // 测试三角形边界点
    std::vector<std::vector<double>> boundaryPoints = {
        {0.0, 0.0, 0.0}, // 顶点1
        {1.0, 0.0, 0.0}, // 顶点2
        {0.0, 1.0, 0.0}, // 顶点3
        {0.5, 0.0, 0.0}, // 边中点
        {0.0, 0.5, 0.0}, // 边中点
        {0.5, 0.5, 0.0}  // 内部点
    };

    for (const auto& point : boundaryPoints) {
        double u = point[0], v = point[1], w = point[2];
        double detJ = 0.0;
        std::vector<double> basis;
        std::vector<std::vector<double>> dBasisdx;
        std::vector<std::vector<std::vector<double>>> ddBasisddx;
        bool compute2ndDerivatives = false;

        bool result = ElementInfo(triangleElement, triangleNodes, u, v, w, detJ,
                                 basis, &dBasisdx, &ddBasisddx, &compute2ndDerivatives);

        EXPECT_TRUE(result);
        EXPECT_GT(detJ, 0.0);
        EXPECT_EQ(basis.size(), 3);
        
        // 验证基函数和为1
        double sum = 0.0;
        for (const auto& val : basis) {
            sum += val;
        }
        EXPECT_NEAR(sum, 1.0, 1e-12);
    }
}

/**
 * @brief 测试ElementInfo二阶导数计算
 */
TEST_F(ElementInfoTest, SecondDerivatives) {
    double u = 0.25, v = 0.25, w = 0.0;
    double detJ = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> dBasisdx;
    std::vector<std::vector<std::vector<double>>> ddBasisddx;
    bool compute2ndDerivatives = true;

    // 调用ElementInfo函数，启用二阶导数计算
    bool result = ElementInfo(triangleElement, triangleNodes, u, v, w, detJ,
                             basis, &dBasisdx, &ddBasisddx, &compute2ndDerivatives,
                             nullptr, nullptr, nullptr, nullptr, nullptr);

    EXPECT_TRUE(result);
    
    // 验证二阶导数矩阵维度
    EXPECT_EQ(ddBasisddx.size(), 3u);
    if (ddBasisddx.size() > 0) {
        EXPECT_EQ(ddBasisddx[0].size(), 2u); // 2x2二阶导数矩阵
        if (ddBasisddx[0].size() > 0) {
            EXPECT_EQ(ddBasisddx[0][0].size(), 2u);
        }
    }
}

/**
 * @brief EdgeElementInfo接口测试类
 */
class EdgeElementInfoTest : public ::testing::Test {
protected:
    void SetUp() override {
        setupTriangleElement();
        setupQuadElement();
        setupTetraElement();
    }

    void setupTriangleElement() {
        triangleElement.type.elementCode = 303;
        triangleElement.type.numberOfNodes = 3;
        triangleElement.type.dimension = 2;
        
        triangleNodes.x = {0.0, 1.0, 0.0};
        triangleNodes.y = {0.0, 0.0, 1.0};
        triangleNodes.z = {0.0, 0.0, 0.0};
    }

    void setupQuadElement() {
        quadElement.type.elementCode = 404;
        quadElement.type.numberOfNodes = 4;
        quadElement.type.dimension = 2;
        
        quadNodes.x = {0.0, 1.0, 1.0, 0.0};
        quadNodes.y = {0.0, 0.0, 1.0, 1.0};
        quadNodes.z = {0.0, 0.0, 0.0, 0.0};
    }

    void setupTetraElement() {
        tetraElement.type.elementCode = 504;
        tetraElement.type.numberOfNodes = 4;
        tetraElement.type.dimension = 3;
        
        tetraNodes.x = {0.0, 1.0, 0.0, 0.0};
        tetraNodes.y = {0.0, 0.0, 1.0, 0.0};
        tetraNodes.z = {0.0, 0.0, 0.0, 1.0};
    }

    Element triangleElement;
    Element quadElement;
    Element tetraElement;
    Nodes triangleNodes;
    Nodes quadNodes;
    Nodes tetraNodes;
};

/**
 * @brief 测试三角形元素的EdgeElementInfo计算
 */
TEST_F(EdgeElementInfoTest, TriangleEdgeElementInfo) {
    double u = 0.25, v = 0.25, w = 0.0;
    double detF = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> edgeBasis;
    
    // 调用EdgeElementInfo函数
    bool result = EdgeElementInfo(triangleElement, triangleNodes, u, v, w,
                                 nullptr, nullptr, detF, basis, edgeBasis,
                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                 nullptr, nullptr, nullptr);

    EXPECT_TRUE(result);
    EXPECT_GT(detF, 0.0);
    
    // 验证基函数数量
    EXPECT_EQ(basis.size(), 3);
    
    // 验证边基函数维度
    EXPECT_GT(edgeBasis.size(), 0);
    if (edgeBasis.size() > 0) {
        EXPECT_EQ(edgeBasis[0].size(), 3); // 3D向量
    }
    
    // 验证基函数和为1
    double sum = 0.0;
    for (const auto& val : basis) {
        sum += val;
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);
}

/**
 * @brief 测试四边形元素的EdgeElementInfo计算
 */
TEST_F(EdgeElementInfoTest, QuadEdgeElementInfo) {
    double u = 0.0, v = 0.0, w = 0.0;
    double detF = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> edgeBasis;
    
    bool result = EdgeElementInfo(quadElement, quadNodes, u, v, w,
                                 nullptr, nullptr, detF, basis, edgeBasis,
                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                 nullptr, nullptr, nullptr);

    EXPECT_TRUE(result);
    EXPECT_GT(detF, 0.0);
    EXPECT_EQ(basis.size(), 4);
    EXPECT_GT(edgeBasis.size(), 0);
}

/**
 * @brief 测试四面体元素的EdgeElementInfo计算
 */
TEST_F(EdgeElementInfoTest, TetraEdgeElementInfo) {
    double u = 0.25, v = 0.25, w = 0.25;
    double detF = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> edgeBasis;
    
    bool result = EdgeElementInfo(tetraElement, tetraNodes, u, v, w,
                                 nullptr, nullptr, detF, basis, edgeBasis,
                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                 nullptr, nullptr, nullptr);

    EXPECT_TRUE(result);
    EXPECT_GT(detF, 0.0);
    EXPECT_EQ(basis.size(), 4);
    EXPECT_GT(edgeBasis.size(), 0);
}

/**
 * @brief 测试EdgeElementInfo的旋度计算
 */
TEST_F(EdgeElementInfoTest, RotBasisCalculation) {
    double u = 0.25, v = 0.25, w = 0.0;
    double detF = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> edgeBasis;
    std::vector<std::vector<double>> rotBasis;
    
    // 调用EdgeElementInfo函数，启用旋度计算
    bool result = EdgeElementInfo(triangleElement, triangleNodes, u, v, w,
                                 nullptr, nullptr, detF, basis, edgeBasis,
                                 &rotBasis, nullptr, nullptr, nullptr, nullptr,
                                 nullptr, nullptr, nullptr);

    EXPECT_TRUE(result);
    
    // 验证旋度矩阵维度
    EXPECT_EQ(rotBasis.size(), edgeBasis.size());
    if (rotBasis.size() > 0) {
        EXPECT_EQ(rotBasis[0].size(), 3); // 3D旋度向量
    }
}

/**
 * @brief 测试EdgeElementInfo的Piola变换
 */
TEST_F(EdgeElementInfoTest, PiolaTransform) {
    double u = 0.25, v = 0.25, w = 0.0;
    double detF = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> edgeBasis;
    bool applyPiolaTransform = true;
    
    // 调用EdgeElementInfo函数，启用Piola变换
    bool result = EdgeElementInfo(triangleElement, triangleNodes, u, v, w,
                                 nullptr, nullptr, detF, basis, edgeBasis,
                                 nullptr, nullptr, nullptr, nullptr, &applyPiolaTransform,
                                 nullptr, nullptr, nullptr);

    EXPECT_TRUE(result);
    EXPECT_GT(std::abs(detF), 0.0); // Piola变换后detF应为绝对值
}

/**
 * @brief 测试EdgeElementInfo的预计算基函数
 */
TEST_F(EdgeElementInfoTest, ReadyBasisFunctions) {
    double u = 0.25, v = 0.25, w = 0.0;
    double detF = 0.0;
    std::vector<double> basis;
    std::vector<std::vector<double>> edgeBasis;
    
    // 创建预计算的基函数
    std::vector<std::vector<double>> readyEdgeBasis = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    std::vector<std::vector<double>> readyRotBasis = {{0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}};
    
    bool result = EdgeElementInfo(triangleElement, triangleNodes, u, v, w,
                                 nullptr, nullptr, detF, basis, edgeBasis,
                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                 &readyEdgeBasis, &readyRotBasis, nullptr);

    EXPECT_TRUE(result);
    
    // 验证使用了预计算的基函数
    EXPECT_EQ(edgeBasis.size(), readyEdgeBasis.size());
}

/**
 * @brief 测试EdgeElementInfo在边界点的情况
 */
TEST_F(EdgeElementInfoTest, BoundaryPoints) {
    std::vector<std::vector<double>> boundaryPoints = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.5, 0.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.5, 0.5, 0.0}
    };

    for (const auto& point : boundaryPoints) {
        double u = point[0], v = point[1], w = point[2];
        double detF = 0.0;
        std::vector<double> basis;
        std::vector<std::vector<double>> edgeBasis;
        
        bool result = EdgeElementInfo(triangleElement, triangleNodes, u, v, w,
                                     nullptr, nullptr, detF, basis, edgeBasis,
                                     nullptr, nullptr, nullptr, nullptr, nullptr,
                                     nullptr, nullptr, nullptr);

        EXPECT_TRUE(result);
        EXPECT_GT(detF, 0.0);
        EXPECT_EQ(basis.size(), 3);
        EXPECT_GT(edgeBasis.size(), 0);
    }
}

/**
 * @brief 主测试入口点
 */
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}