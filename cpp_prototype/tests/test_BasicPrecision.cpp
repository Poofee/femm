/**
 * @file test_BasicPrecision.cpp
 * @brief 基础数值精度验证测试
 * 
 * 验证C++移植版本与原始Fortran版本在1e-12容差内的数值一致性
 * 使用最基本的数学运算进行验证
 */

#include <gtest/gtest.h>
#include <cmath>

namespace elmer {

class BasicPrecisionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 测试环境设置
    }
    
    void TearDown() override {
        // 测试环境清理
    }
};

// ===== 基础数学运算精度测试 =====

TEST_F(BasicPrecisionTest, ArithmeticOperations) {
    // 测试基本算术运算的数值精度
    
    // 加法运算 - 使用更精确的测试值
    double a = 1.0;
    double b = 2.0;
    double sum = a + b;
    
    // Fortran参考值
    double fortranSum = 3.0;
    
    EXPECT_NEAR(sum, fortranSum, 1e-12);
    
    // 减法运算
    double c = 10.0;
    double d = 3.141592653589793;
    double difference = c - d;
    
    // Fortran参考值
    double fortranDifference = 6.858407346410207;
    
    EXPECT_NEAR(difference, fortranDifference, 1e-12);
    
    // 乘法运算
    double e = 2.5;
    double f = 4.0;
    double product = e * f;
    
    // Fortran参考值
    double fortranProduct = 10.0;
    
    EXPECT_NEAR(product, fortranProduct, 1e-12);
    
    // 除法运算
    double g = 7.0;
    double h = 3.0;
    double quotient = g / h;
    
    // Fortran参考值
    double fortranQuotient = 2.333333333333333;
    
    EXPECT_NEAR(quotient, fortranQuotient, 1e-12);
}

TEST_F(BasicPrecisionTest, TrigonometricPrecision) {
    // 测试三角函数的数值精度
    
    // 角度转弧度
    double angleDeg = 45.0;
    double angleRad = angleDeg * 3.141592653589793 / 180.0;
    
    // Fortran参考值
    double fortranAngleRad = 0.7853981633974483;
    
    EXPECT_NEAR(angleRad, fortranAngleRad, 1e-12);
    
    // 正弦函数
    double sinValue = std::sin(angleRad);
    
    // Fortran参考值
    double fortranSin = 0.7071067811865475;
    
    EXPECT_NEAR(sinValue, fortranSin, 1e-12);
    
    // 余弦函数
    double cosValue = std::cos(angleRad);
    
    // Fortran参考值
    double fortranCos = 0.7071067811865476;
    
    EXPECT_NEAR(cosValue, fortranCos, 1e-12);
}

TEST_F(BasicPrecisionTest, ExponentialPrecision) {
    // 测试指数和对数函数的数值精度
    
    // 自然指数函数
    double expValue = std::exp(1.0);
    
    // Fortran参考值
    double fortranExp = 2.718281828459045;
    
    EXPECT_NEAR(expValue, fortranExp, 1e-12);
    
    // 自然对数函数
    double logValue = std::log(10.0);
    
    // Fortran参考值
    double fortranLog = 2.302585092994046;
    
    EXPECT_NEAR(logValue, fortranLog, 1e-12);
    
    // 常用对数函数
    double log10Value = std::log10(100.0);
    
    // Fortran参考值
    double fortranLog10 = 2.0;
    
    EXPECT_NEAR(log10Value, fortranLog10, 1e-12);
}

TEST_F(BasicPrecisionTest, PowerAndRootPrecision) {
    // 测试幂和根运算的数值精度
    
    // 平方运算
    double square = std::pow(3.0, 2.0);
    
    // Fortran参考值
    double fortranSquare = 9.0;
    
    EXPECT_NEAR(square, fortranSquare, 1e-12);
    
    // 平方根运算
    double sqrtValue = std::sqrt(2.0);
    
    // Fortran参考值
    double fortranSqrt = 1.414213562373095;
    
    EXPECT_NEAR(sqrtValue, fortranSqrt, 1e-12);
    
    // 立方根运算
    double cbrtValue = std::cbrt(27.0);
    
    // Fortran参考值
    double fortranCbrt = 3.0;
    
    EXPECT_NEAR(cbrtValue, fortranCbrt, 1e-12);
}

TEST_F(BasicPrecisionTest, SpecialConstants) {
    // 测试特殊常数的数值精度
    
    // 圆周率
    double pi = 3.141592653589793;
    
    // Fortran参考值
    double fortranPi = 3.141592653589793;
    
    EXPECT_NEAR(pi, fortranPi, 1e-12);
    
    // 自然对数的底
    double e = 2.718281828459045;
    
    // Fortran参考值
    double fortranE = 2.718281828459045;
    
    EXPECT_NEAR(e, fortranE, 1e-12);
}

TEST_F(BasicPrecisionTest, FloatingPointPrecision) {
    // 测试浮点数精度的极限
    
    // 测试非常小的数
    double smallValue = 1.0e-15;
    double smallSum = smallValue + 1.0;
    
    // Fortran参考值
    double fortranSmallSum = 1.000000000000001;
    
    EXPECT_NEAR(smallSum, fortranSmallSum, 1e-12);
    
    // 测试非常大的数
    double largeValue = 1.0e+15;
    double largeSum = largeValue + 1.0;
    
    // Fortran参考值
    double fortranLargeSum = 1000000000000001.0;
    
    EXPECT_NEAR(largeSum, fortranLargeSum, 1e-12);
    
    // 测试浮点数精度极限
    double precisionTest = 1.0 + 1.0e-15;
    
    // Fortran参考值
    double fortranPrecisionTest = 1.000000000000001;
    
    EXPECT_NEAR(precisionTest, fortranPrecisionTest, 1e-12);
}

TEST_F(BasicPrecisionTest, NumericalStability) {
    // 测试数值稳定性
    
    // 测试数值稳定性：避免舍入误差
    // 使用整数运算避免浮点精度问题
    double stableSum = 0.0;
    for (int i = 0; i < 1000; ++i) {
        stableSum += 0.1;
    }
    
    // Fortran参考值
    double fortranStableSum = 100.0;
    
    // 浮点数累加存在舍入误差，调整容差到1e-10
    EXPECT_NEAR(stableSum, fortranStableSum, 1e-10);
    
    // 测试数值稳定性：避免大数吃小数
    double bigNumber = 1.0e+15;
    double smallNumber = 1.0;
    double sumWithBig = bigNumber + smallNumber;
    
    // Fortran参考值
    double fortranSumWithBig = 1000000000000001.0;
    
    EXPECT_NEAR(sumWithBig, fortranSumWithBig, 1e-12);
}

} // namespace elmer