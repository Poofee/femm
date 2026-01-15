/**
 * @file test_shapefunctions.cpp
 * @brief ShapeFunctions模块测试文件
 * 
 * 测试ShapeFunctions模块的形函数计算功能，包括线性、二次形函数及其导数。
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "../src/ShapeFunctions.h"

using namespace elmer;

/**
 * @brief 自定义断言宏，用于测试
 */
#define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "断言失败: " << message << std::endl; \
            std::cerr << "位置: " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

/**
 * @brief 自定义近似相等断言宏
 */
#define ASSERT_NEAR(a, b, tolerance) \
    do { \
        if (std::abs((a) - (b)) > (tolerance)) { \
            std::cerr << "近似相等断言失败: " << (a) << " != " << (b) << std::endl; \
            std::cerr << "容差: " << (tolerance) << std::endl; \
            std::cerr << "位置: " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

/**
 * @brief 测试1D线性形函数
 */
bool Test1DLinearShapeFunctions() {
    std::cout << "测试1D线性形函数..." << std::endl;
    
    // 测试端点
    auto N_minus1 = ShapeFunctions::linear1D(-1.0, 2);
    ASSERT(N_minus1.size() == 2, "2节点形函数应有2个值");
    ASSERT_NEAR(N_minus1[0], 1.0, 1e-12);
    ASSERT_NEAR(N_minus1[1], 0.0, 1e-12);
    
    auto N_plus1 = ShapeFunctions::linear1D(1.0, 2);
    ASSERT_NEAR(N_plus1[0], 0.0, 1e-12);
    ASSERT_NEAR(N_plus1[1], 1.0, 1e-12);
    
    // 测试中点
    auto N_zero = ShapeFunctions::linear1D(0.0, 2);
    ASSERT_NEAR(N_zero[0], 0.5, 1e-12);
    ASSERT_NEAR(N_zero[1], 0.5, 1e-12);
    
    // 测试形函数和为1
    double sum = N_zero[0] + N_zero[1];
    ASSERT_NEAR(sum, 1.0, 1e-12);
    
    std::cout << "1D线性形函数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试1D二次形函数
 */
bool Test1DQuadraticShapeFunctions() {
    std::cout << "测试1D二次形函数..." << std::endl;
    
    // 测试端点
    auto N_minus1 = ShapeFunctions::linear1D(-1.0, 3);
    ASSERT(N_minus1.size() == 3, "3节点形函数应有3个值");
    ASSERT_NEAR(N_minus1[0], 1.0, 1e-12);
    ASSERT_NEAR(N_minus1[1], 0.0, 1e-12);
    ASSERT_NEAR(N_minus1[2], 0.0, 1e-12);
    
    auto N_plus1 = ShapeFunctions::linear1D(1.0, 3);
    ASSERT_NEAR(N_plus1[0], 0.0, 1e-12);
    ASSERT_NEAR(N_plus1[1], 0.0, 1e-12);
    ASSERT_NEAR(N_plus1[2], 1.0, 1e-12);
    
    // 测试中点
    auto N_zero = ShapeFunctions::linear1D(0.0, 3);
    ASSERT_NEAR(N_zero[0], 0.0, 1e-12);
    ASSERT_NEAR(N_zero[1], 1.0, 1e-12);
    ASSERT_NEAR(N_zero[2], 0.0, 1e-12);
    
    // 测试形函数和为1
    double sum = N_zero[0] + N_zero[1] + N_zero[2];
    ASSERT_NEAR(sum, 1.0, 1e-12);
    
    std::cout << "1D二次形函数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试1D形函数导数
 */
bool Test1DShapeFunctionDerivatives() {
    std::cout << "测试1D形函数导数..." << std::endl;
    
    // 测试线性形函数导数
    auto dN_linear = ShapeFunctions::linear1DDerivatives(0.0, 2);
    ASSERT(dN_linear.size() == 2, "线性形函数导数应有2个值");
    ASSERT_NEAR(dN_linear[0], -0.5, 1e-12);
    ASSERT_NEAR(dN_linear[1], 0.5, 1e-12);
    
    // 测试二次形函数导数
    auto dN_quadratic = ShapeFunctions::linear1DDerivatives(0.0, 3);
    ASSERT(dN_quadratic.size() == 3, "二次形函数导数应有3个值");
    ASSERT_NEAR(dN_quadratic[0], -0.5, 1e-12);
    ASSERT_NEAR(dN_quadratic[1], 0.0, 1e-12);
    ASSERT_NEAR(dN_quadratic[2], 0.5, 1e-12);
    
    // 测试导数在端点
    auto dN_end = ShapeFunctions::linear1DDerivatives(1.0, 3);
    ASSERT_NEAR(dN_end[0], 0.5, 1e-12);
    ASSERT_NEAR(dN_end[1], -2.0, 1e-12);
    ASSERT_NEAR(dN_end[2], 1.5, 1e-12);
    
    std::cout << "1D形函数导数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试四边形线性形函数
 */
bool TestQuadrilateralLinearShapeFunctions() {
    std::cout << "测试四边形线性形函数..." << std::endl;
    
    // 测试角点
    auto N_corner = ShapeFunctions::linearQuadrilateral(-1.0, -1.0, 4);
    ASSERT(N_corner.size() == 4, "四边形形函数应有4个值");
    ASSERT_NEAR(N_corner[0], 1.0, 1e-12);
    ASSERT_NEAR(N_corner[1], 0.0, 1e-12);
    ASSERT_NEAR(N_corner[2], 0.0, 1e-12);
    ASSERT_NEAR(N_corner[3], 0.0, 1e-12);
    
    // 测试中心
    auto N_center = ShapeFunctions::linearQuadrilateral(0.0, 0.0, 4);
    ASSERT_NEAR(N_center[0], 0.25, 1e-12);
    ASSERT_NEAR(N_center[1], 0.25, 1e-12);
    ASSERT_NEAR(N_center[2], 0.25, 1e-12);
    ASSERT_NEAR(N_center[3], 0.25, 1e-12);
    
    // 测试形函数和为1
    double sum = 0.0;
    for (double n : N_center) {
        sum += n;
    }
    ASSERT_NEAR(sum, 1.0, 1e-12);
    
    std::cout << "四边形线性形函数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试四边形形函数导数
 */
bool TestQuadrilateralShapeFunctionDerivatives() {
    std::cout << "测试四边形形函数导数..." << std::endl;
    
    // 测试中心点导数
    auto derivatives = ShapeFunctions::linearQuadrilateralDerivatives(0.0, 0.0, 4);
    auto& dNdxi = derivatives.first;
    auto& dNdeta = derivatives.second;
    
    ASSERT(dNdxi.size() == 4, "四边形形函数导数应有4个值");
    ASSERT(dNdeta.size() == 4, "四边形形函数导数应有4个值");
    
    // 检查dNdxi
    ASSERT_NEAR(dNdxi[0], -0.25, 1e-12);
    ASSERT_NEAR(dNdxi[1], 0.25, 1e-12);
    ASSERT_NEAR(dNdxi[2], 0.25, 1e-12);
    ASSERT_NEAR(dNdxi[3], -0.25, 1e-12);
    
    // 检查dNdeta
    ASSERT_NEAR(dNdeta[0], -0.25, 1e-12);
    ASSERT_NEAR(dNdeta[1], -0.25, 1e-12);
    ASSERT_NEAR(dNdeta[2], 0.25, 1e-12);
    ASSERT_NEAR(dNdeta[3], 0.25, 1e-12);
    
    std::cout << "四边形形函数导数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试六面体线性形函数
 */
bool TestHexahedronLinearShapeFunctions() {
    std::cout << "测试六面体线性形函数..." << std::endl;
    
    // 测试角点
    auto N_corner = ShapeFunctions::linearHexahedron(-1.0, -1.0, -1.0, 8);
    ASSERT(N_corner.size() == 8, "六面体形函数应有8个值");
    ASSERT_NEAR(N_corner[0], 1.0, 1e-12);
    for (int i = 1; i < 8; ++i) {
        ASSERT_NEAR(N_corner[i], 0.0, 1e-12);
    }
    
    // 测试中心
    auto N_center = ShapeFunctions::linearHexahedron(0.0, 0.0, 0.0, 8);
    for (int i = 0; i < 8; ++i) {
        ASSERT_NEAR(N_center[i], 0.125, 1e-12);
    }
    
    // 测试形函数和为1
    double sum = 0.0;
    for (double n : N_center) {
        sum += n;
    }
    ASSERT_NEAR(sum, 1.0, 1e-12);
    
    std::cout << "六面体线性形函数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试六面体形函数导数
 */
bool TestHexahedronShapeFunctionDerivatives() {
    std::cout << "测试六面体形函数导数..." << std::endl;
    
    // 测试中心点导数
    auto derivatives = ShapeFunctions::linearHexahedronDerivatives(0.0, 0.0, 0.0, 8);
    auto& dNdxi = std::get<0>(derivatives);
    auto& dNdeta = std::get<1>(derivatives);
    auto& dNdzeta = std::get<2>(derivatives);
    
    ASSERT(dNdxi.size() == 8, "六面体形函数导数应有8个值");
    ASSERT(dNdeta.size() == 8, "六面体形函数导数应有8个值");
    ASSERT(dNdzeta.size() == 8, "六面体形函数导数应有8个值");
    
    // 检查dNdxi
    ASSERT_NEAR(dNdxi[0], -0.125, 1e-12);
    ASSERT_NEAR(dNdxi[1], 0.125, 1e-12);
    ASSERT_NEAR(dNdxi[2], 0.125, 1e-12);
    ASSERT_NEAR(dNdxi[3], -0.125, 1e-12);
    ASSERT_NEAR(dNdxi[4], -0.125, 1e-12);
    ASSERT_NEAR(dNdxi[5], 0.125, 1e-12);
    ASSERT_NEAR(dNdxi[6], 0.125, 1e-12);
    ASSERT_NEAR(dNdxi[7], -0.125, 1e-12);
    
    std::cout << "六面体形函数导数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试形函数的插值性质
 */
bool TestShapeFunctionInterpolation() {
    std::cout << "测试形函数的插值性质..." << std::endl;
    
    // 测试1D线性插值
    std::vector<double> nodalValues = {2.0, 4.0}; // 节点值
    
    // 在端点处插值
    auto N_left = ShapeFunctions::linear1D(-1.0, 2);
    double interpolated_left = 0.0;
    for (size_t i = 0; i < N_left.size(); ++i) {
        interpolated_left += N_left[i] * nodalValues[i];
    }
    ASSERT_NEAR(interpolated_left, 2.0, 1e-12);
    
    // 在右端点处插值
    auto N_right = ShapeFunctions::linear1D(1.0, 2);
    double interpolated_right = 0.0;
    for (size_t i = 0; i < N_right.size(); ++i) {
        interpolated_right += N_right[i] * nodalValues[i];
    }
    ASSERT_NEAR(interpolated_right, 4.0, 1e-12);
    
    // 在中点处插值
    auto N_mid = ShapeFunctions::linear1D(0.0, 2);
    double interpolated_mid = 0.0;
    for (size_t i = 0; i < N_mid.size(); ++i) {
        interpolated_mid += N_mid[i] * nodalValues[i];
    }
    ASSERT_NEAR(interpolated_mid, 3.0, 1e-12);
    
    std::cout << "形函数的插值性质测试通过" << std::endl;
    return true;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始ShapeFunctions模块测试..." << std::endl;
    
    bool allPassed = true;
    
    try {
        allPassed &= Test1DLinearShapeFunctions();
        allPassed &= Test1DQuadraticShapeFunctions();
        allPassed &= Test1DShapeFunctionDerivatives();
        allPassed &= TestQuadrilateralLinearShapeFunctions();
        allPassed &= TestQuadrilateralShapeFunctionDerivatives();
        allPassed &= TestHexahedronLinearShapeFunctions();
        allPassed &= TestHexahedronShapeFunctionDerivatives();
        allPassed &= TestShapeFunctionInterpolation();
        
        if (allPassed) {
            std::cout << "\n✅ 所有ShapeFunctions模块测试通过！" << std::endl;
            return 0;
        } else {
            std::cout << "\n❌ 部分ShapeFunctions模块测试失败！" << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "测试异常: " << e.what() << std::endl;
        return 1;
    }
}