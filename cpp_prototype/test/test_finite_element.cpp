// test_finite_element.cpp - 有限元功能验证测试
#include "ElementMatrix.h"
#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include <iostream>
#include <iomanip>

using namespace elmer;

void test1DLinearElement() {
    std::cout << "=== 1D线性单元测试 ===" << std::endl;
    
    // 创建1D线性单元（2节点）
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0},  // 节点1
        {1.0, 0.0, 0.0}   // 节点2
    };
    
    // 材料属性
    ElementMatrix::MaterialProperties material;
    material.youngsModulus = 1.0e6;  // 杨氏模量
    material.density = 1.0;          // 密度
    
    // 计算刚度矩阵
    auto stiffnessMatrix = ElementMatrix::computeStiffnessMatrix(
        ElementType::LINEAR, nodes, material, 2);
    
    std::cout << "1D刚度矩阵:" << std::endl;
    for (const auto& row : stiffnessMatrix) {
        for (double val : row) {
            std::cout << std::setw(12) << val << " ";
        }
        std::cout << std::endl;
    }
    
    // 计算质量矩阵
    auto massMatrix = ElementMatrix::computeMassMatrix(
        ElementType::LINEAR, nodes, material, 2);
    
    std::cout << "1D质量矩阵:" << std::endl;
    for (const auto& row : massMatrix) {
        for (double val : row) {
            std::cout << std::setw(12) << val << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
}

void test2DQuadrilateralElement() {
    std::cout << "=== 2D四边形单元测试 ===" << std::endl;
    
    // 创建2D四边形单元（4节点）
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0},  // 节点1
        {1.0, 0.0, 0.0},  // 节点2
        {1.0, 1.0, 0.0},  // 节点3
        {0.0, 1.0, 0.0}   // 节点4
    };
    
    // 材料属性
    ElementMatrix::MaterialProperties material;
    material.youngsModulus = 1.0e6;
    material.poissonsRatio = 0.3;
    material.density = 1.0;
    
    // 计算刚度矩阵
    auto stiffnessMatrix = ElementMatrix::computeStiffnessMatrix(
        ElementType::QUADRILATERAL, nodes, material, 2);
    
    std::cout << "2D刚度矩阵大小: " << stiffnessMatrix.size() << "x" << stiffnessMatrix[0].size() << std::endl;
    
    // 计算质量矩阵
    auto massMatrix = ElementMatrix::computeMassMatrix(
        ElementType::QUADRILATERAL, nodes, material, 2);
    
    std::cout << "2D质量矩阵大小: " << massMatrix.size() << "x" << massMatrix[0].size() << std::endl;
    
    // 计算载荷向量
    std::array<double, 3> bodyForce = {0.0, -9.81, 0.0};  // 重力载荷
    auto loadVector = ElementMatrix::computeLoadVector(
        ElementType::QUADRILATERAL, nodes, bodyForce, 2);
    
    std::cout << "2D载荷向量大小: " << loadVector.size() << std::endl;
    
    std::cout << std::endl;
}

void testShapeFunctions() {
    std::cout << "=== 形状函数测试 ===" << std::endl;
    
    // 测试1D形状函数
    std::vector<Node> nodes1D = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    auto shapeResult1D = ShapeFunctions::computeShapeFunctions(
        ElementType::LINEAR, nodes1D, 0.0, 0.0, 0.0);
    
    std::cout << "1D形状函数值: ";
    for (double val : shapeResult1D.values) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Jacobian行列式: " << shapeResult1D.detJ << std::endl;
    
    std::cout << std::endl;
}

void testGaussIntegration() {
    std::cout << "=== 高斯积分测试 ===" << std::endl;
    
    // 测试1D高斯积分
    auto points1D = GaussIntegration::get1DPoints(2);
    std::cout << "1D高斯积分点数量: " << points1D.size() << std::endl;
    
    // 测试函数 f(x) = x^2
    auto result = GaussIntegration::integrate(points1D, 
        [](double xi, double eta, double zeta) { return xi * xi; });
    
    std::cout << "∫x²dx在[-1,1]上的积分结果: " << result << std::endl;
    std::cout << "理论值: " << 2.0/3.0 << std::endl;
    
    std::cout << std::endl;
}

int main() {
    std::cout << "Elmer FEM C++ 有限元功能验证测试" << std::endl;
    std::cout << "================================" << std::endl;
    
    try {
        test1DLinearElement();
        test2DQuadrilateralElement();
        testShapeFunctions();
        testGaussIntegration();
        
        std::cout << "✅ 所有测试完成！" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}