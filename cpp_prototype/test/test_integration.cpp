/**
 * @file test_integration.cpp
 * @brief Integration模块测试文件
 * 
 * 测试Integration模块的数值积分功能，包括Gauss积分点和权重计算。
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "../src/Integration.h"

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
 * @brief 测试1D Gauss积分点
 */
bool TestGaussPoints1D() {
    std::cout << "测试1D Gauss积分点..." << std::endl;
    
    // 测试1点积分
    auto points1 = Integration::getGaussPoints1D(1);
    ASSERT(points1.size() == 1, "1点积分应有1个积分点");
    ASSERT_NEAR(points1[0], 0.0, 1e-12);
    
    // 测试2点积分
    auto points2 = Integration::getGaussPoints1D(2);
    ASSERT(points2.size() == 2, "2点积分应有2个积分点");
    double expected2 = 1.0 / std::sqrt(3.0);
    ASSERT_NEAR(points2[0], -expected2, 1e-12);
    ASSERT_NEAR(points2[1], expected2, 1e-12);
    
    // 测试3点积分
    auto points3 = Integration::getGaussPoints1D(3);
    ASSERT(points3.size() == 3, "3点积分应有3个积分点");
    double expected3 = std::sqrt(3.0 / 5.0);
    ASSERT_NEAR(points3[0], -expected3, 1e-12);
    ASSERT_NEAR(points3[1], 0.0, 1e-12);
    ASSERT_NEAR(points3[2], expected3, 1e-12);
    
    std::cout << "1D Gauss积分点测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试1D Gauss积分权重
 */
bool TestGaussWeights1D() {
    std::cout << "测试1D Gauss积分权重..." << std::endl;
    
    // 测试1点权重
    auto weights1 = Integration::getGaussWeights1D(1);
    ASSERT(weights1.size() == 1, "1点积分应有1个权重");
    ASSERT_NEAR(weights1[0], 2.0, 1e-12);
    
    // 测试2点权重
    auto weights2 = Integration::getGaussWeights1D(2);
    ASSERT(weights2.size() == 2, "2点积分应有2个权重");
    ASSERT_NEAR(weights2[0], 1.0, 1e-12);
    ASSERT_NEAR(weights2[1], 1.0, 1e-12);
    
    // 测试3点权重
    auto weights3 = Integration::getGaussWeights1D(3);
    ASSERT(weights3.size() == 3, "3点积分应有3个权重");
    ASSERT_NEAR(weights3[0], 5.0/9.0, 1e-12);
    ASSERT_NEAR(weights3[1], 8.0/9.0, 1e-12);
    ASSERT_NEAR(weights3[2], 5.0/9.0, 1e-12);
    
    std::cout << "1D Gauss积分权重测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试三角形积分规则
 */
bool TestTriangleIntegration() {
    std::cout << "测试三角形积分规则..." << std::endl;
    
    // 测试1点三角形积分
    auto rule1 = Integration::getTriangleRule(1);
    ASSERT(rule1.u.size() == 1, "1点三角形积分应有1个积分点");
    ASSERT(rule1.v.size() == 1, "1点三角形积分应有1个积分点");
    ASSERT(rule1.weights.size() == 1, "1点三角形积分应有1个权重");
    ASSERT_NEAR(rule1.u[0], 1.0/3.0, 1e-12);
    ASSERT_NEAR(rule1.v[0], 1.0/3.0, 1e-12);
    ASSERT_NEAR(rule1.weights[0], 1.0, 1e-12);
    
    // 测试3点三角形积分
    auto rule3 = Integration::getTriangleRule(3);
    ASSERT(rule3.u.size() == 3, "3点三角形积分应有3个积分点");
    ASSERT(rule3.v.size() == 3, "3点三角形积分应有3个积分点");
    ASSERT(rule3.weights.size() == 3, "3点三角形积分应有3个权重");
    
    // 检查积分点坐标
    double one_sixth = 1.0/6.0;
    double two_thirds = 2.0/3.0;
    ASSERT_NEAR(rule3.u[0], one_sixth, 1e-12);
    ASSERT_NEAR(rule3.u[1], two_thirds, 1e-12);
    ASSERT_NEAR(rule3.u[2], one_sixth, 1e-12);
    
    ASSERT_NEAR(rule3.v[0], one_sixth, 1e-12);
    ASSERT_NEAR(rule3.v[1], one_sixth, 1e-12);
    ASSERT_NEAR(rule3.v[2], two_thirds, 1e-12);
    
    // 检查权重
    double weight = 1.0/3.0;
    ASSERT_NEAR(rule3.weights[0], weight, 1e-12);
    ASSERT_NEAR(rule3.weights[1], weight, 1e-12);
    ASSERT_NEAR(rule3.weights[2], weight, 1e-12);
    
    // 测试权重和为1（面积归一化）
    double sum_weights = 0.0;
    for (double w : rule3.weights) {
        sum_weights += w;
    }
    ASSERT_NEAR(sum_weights, 1.0, 1e-12);
    
    std::cout << "三角形积分规则测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试四边形积分规则
 */
bool TestQuadrilateralIntegration() {
    std::cout << "测试四边形积分规则..." << std::endl;
    
    // 测试1点四边形积分
    auto rule1 = Integration::getQuadrilateralRule(1);
    auto& points1 = rule1.first;
    auto& weights1 = rule1.second;
    
    ASSERT(points1.size() == 1, "1点四边形积分应有1个积分点");
    ASSERT(weights1.size() == 1, "1点四边形积分应有1个权重");
    ASSERT_NEAR(points1[0], 0.0, 1e-12);
    ASSERT_NEAR(weights1[0], 4.0, 1e-12); // 面积归一化
    
    // 测试2点四边形积分（张量积：2x2=4点）
    auto rule2 = Integration::getQuadrilateralRule(2);
    auto& points2 = rule2.first;
    auto& weights2 = rule2.second;
    
    ASSERT(points2.size() == 4, "2点四边形积分应有4个积分点（2x2张量积）");
    ASSERT(weights2.size() == 4, "2点四边形积分应有4个权重（2x2张量积）");
    
    // 检查权重和为4（面积归一化）
    double sum_weights = 0.0;
    for (double w : weights2) {
        sum_weights += w;
    }
    ASSERT_NEAR(sum_weights, 4.0, 1e-12);
    
    std::cout << "四边形积分规则测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试砖形单元积分规则
 */
bool TestBrickIntegration() {
    std::cout << "测试砖形单元积分规则..." << std::endl;
    
    // 测试1点砖形单元积分
    auto rule1 = Integration::getBrickRule(1);
    auto& points1 = std::get<0>(rule1);
    auto& weights1 = std::get<3>(rule1);
    
    ASSERT(points1.size() == 1, "1点砖形单元积分应有1个积分点");
    ASSERT(weights1.size() == 1, "1点砖形单元积分应有1个权重");
    ASSERT_NEAR(points1[0], 0.0, 1e-12);
    ASSERT_NEAR(weights1[0], 8.0, 1e-12); // 体积归一化
    
    // 测试2点砖形单元积分（张量积：2x2x2=8点）
    auto rule2 = Integration::getBrickRule(2);
    auto& points2 = std::get<0>(rule2);
    auto& weights2 = std::get<3>(rule2);
    
    ASSERT(points2.size() == 8, "2点砖形单元积分应有8个积分点（2x2x2张量积）");
    ASSERT(weights2.size() == 8, "2点砖形单元积分应有8个权重（2x2x2张量积）");
    
    // 检查权重和为8（体积归一化）
    double sum_weights = 0.0;
    for (double w : weights2) {
        sum_weights += w;
    }
    ASSERT_NEAR(sum_weights, 8.0, 1e-12);
    
    std::cout << "砖形单元积分规则测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试数值积分精度
 */
bool TestIntegrationAccuracy() {
    std::cout << "测试数值积分精度..." << std::endl;
    
    // 测试1D多项式积分精度
    auto points3 = Integration::getGaussPoints1D(3);
    auto weights3 = Integration::getGaussWeights1D(3);
    
    // 积分 x^2 在 [-1, 1] 区间，精确值为 2/3
    double integral_x2 = 0.0;
    for (size_t i = 0; i < points3.size(); ++i) {
        integral_x2 += weights3[i] * (points3[i] * points3[i]);
    }
    ASSERT_NEAR(integral_x2, 2.0/3.0, 1e-12);
    
    // 积分 x^4 在 [-1, 1] 区间，精确值为 2/5
    double integral_x4 = 0.0;
    for (size_t i = 0; i < points3.size(); ++i) {
        integral_x4 += weights3[i] * (points3[i] * points3[i] * points3[i] * points3[i]);
    }
    ASSERT_NEAR(integral_x4, 2.0/5.0, 1e-12);
    
    std::cout << "数值积分精度测试通过" << std::endl;
    return true;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始Integration模块测试..." << std::endl;
    
    bool allPassed = true;
    
    try {
        allPassed &= TestGaussPoints1D();
        allPassed &= TestGaussWeights1D();
        allPassed &= TestTriangleIntegration();
        allPassed &= TestQuadrilateralIntegration();
        allPassed &= TestBrickIntegration();
        allPassed &= TestIntegrationAccuracy();
        
        if (allPassed) {
            std::cout << "\n✅ 所有Integration模块测试通过！" << std::endl;
            return 0;
        } else {
            std::cout << "\n❌ 部分Integration模块测试失败！" << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "测试异常: " << e.what() << std::endl;
        return 1;
    }
}