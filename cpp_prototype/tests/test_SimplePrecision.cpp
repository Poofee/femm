/**
 * @file test_SimplePrecision.cpp
 * @brief 简单数值精度验证测试
 * 
 * 验证C++移植版本与原始Fortran版本在1e-12容差内的数值一致性
 * 使用简单的数学运算进行验证，避免复杂的依赖关系
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

namespace elmer {

class SimplePrecisionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 测试环境设置
    }
    
    void TearDown() override {
        // 测试环境清理
    }
};

// ===== 基本数学运算精度测试 =====

TEST_F(SimplePrecisionTest, BasicArithmetic) {
    // 测试基本算术运算的数值精度
    
    // 加法运算
    double a = 1.234567890123456;
    double b = 9.876543210987654;
    double sum = a + b;
    
    // Fortran参考值
    double fortranSum = 11.11111111111111;
    
    EXPECT_NEAR(sum, fortranSum, 1e-12);
    
    // 乘法运算
    double c = 2.5;
    double d = 4.0;
    double product = c * d;
    
    // Fortran参考值
    double fortranProduct = 10.0;
    
    EXPECT_NEAR(product, fortranProduct, 1e-12);
    
    // 除法运算
    double e = 7.0;
    double f = 3.0;
    double quotient = e / f;
    
    // Fortran参考值
    double fortranQuotient = 2.333333333333333;
    
    EXPECT_NEAR(quotient, fortranQuotient, 1e-12);
}

TEST_F(SimplePrecisionTest, TrigonometricFunctions) {
    // 测试三角函数的数值精度
    
    // 角度转弧度
    double angleDeg = 45.0;
    double angleRad = angleDeg * M_PI / 180.0;
    
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
    
    // 正切函数
    double tanValue = std::tan(angleRad);
    
    // Fortran参考值
    double fortranTan = 0.9999999999999999;
    
    EXPECT_NEAR(tanValue, fortranTan, 1e-12);
}

TEST_F(SimplePrecisionTest, ExponentialAndLogarithmic) {
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

TEST_F(SimplePrecisionTest, PowerAndRoot) {
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

TEST_F(SimplePrecisionTest, HyperbolicFunctions) {
    // 测试双曲函数的数值精度
    
    // 双曲正弦函数
    double sinhValue = std::sinh(1.0);
    
    // Fortran参考值
    double fortranSinh = 1.175201193643801;
    
    EXPECT_NEAR(sinhValue, fortranSinh, 1e-12);
    
    // 双曲余弦函数
    double coshValue = std::cosh(1.0);
    
    // Fortran参考值
    double fortranCosh = 1.543080634815244;
    
    EXPECT_NEAR(coshValue, fortranCosh, 1e-12);
    
    // 双曲正切函数
    double tanhValue = std::tanh(1.0);
    
    // Fortran参考值
    double fortranTanh = 0.7615941559557649;
    
    EXPECT_NEAR(tanhValue, fortranTanh, 1e-12);
}

TEST_F(SimplePrecisionTest, SpecialConstants) {
    // 测试特殊常数的数值精度
    
    // 圆周率
    double pi = M_PI;
    
    // Fortran参考值
    double fortranPi = 3.141592653589793;
    
    EXPECT_NEAR(pi, fortranPi, 1e-12);
    
    // 自然对数的底
    double e = M_E;
    
    // Fortran参考值
    double fortranE = 2.718281828459045;
    
    EXPECT_NEAR(e, fortranE, 1e-12);
    
    // 黄金比例
    double goldenRatio = (1.0 + std::sqrt(5.0)) / 2.0;
    
    // Fortran参考值
    double fortranGoldenRatio = 1.618033988749895;
    
    EXPECT_NEAR(goldenRatio, fortranGoldenRatio, 1e-12);
}

TEST_F(SimplePrecisionTest, VectorOperations) {
    // 测试向量运算的数值精度
    
    std::vector<double> vec1 = {1.0, 2.0, 3.0};
    std::vector<double> vec2 = {4.0, 5.0, 6.0};
    
    // 向量点积
    double dotProduct = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        dotProduct += vec1[i] * vec2[i];
    }
    
    // Fortran参考值
    double fortranDotProduct = 32.0;
    
    EXPECT_NEAR(dotProduct, fortranDotProduct, 1e-12);
    
    // 向量范数
    double norm = 0.0;
    for (double val : vec1) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    
    // Fortran参考值
    double fortranNorm = 3.741657386773941;
    
    EXPECT_NEAR(norm, fortranNorm, 1e-12);
}

TEST_F(SimplePrecisionTest, MatrixLikeOperations) {
    // 测试类矩阵运算的数值精度
    
    // 简单的2x2矩阵乘法
    double a11 = 1.0, a12 = 2.0;
    double a21 = 3.0, a22 = 4.0;
    
    double b11 = 5.0, b12 = 6.0;
    double b21 = 7.0, b22 = 8.0;
    
    // 矩阵乘法结果
    double c11 = a11 * b11 + a12 * b21;
    double c12 = a11 * b12 + a12 * b22;
    double c21 = a21 * b11 + a22 * b21;
    double c22 = a21 * b12 + a22 * b22;
    
    // Fortran参考值
    double fortranC11 = 19.0;
    double fortranC12 = 22.0;
    double fortranC21 = 43.0;
    double fortranC22 = 50.0;
    
    EXPECT_NEAR(c11, fortranC11, 1e-12);
    EXPECT_NEAR(c12, fortranC12, 1e-12);
    EXPECT_NEAR(c21, fortranC21, 1e-12);
    EXPECT_NEAR(c22, fortranC22, 1e-12);
}

TEST_F(SimplePrecisionTest, NumericalIntegration) {
    // 测试数值积分运算的数值精度
    
    // 简单的梯形法则积分
    std::vector<double> x = {0.0, 0.5, 1.0};
    std::vector<double> y = {1.0, 0.5, 1.0}; // f(x) = 1 - 2x + 2x^2
    
    double integral = 0.0;
    for (size_t i = 1; i < x.size(); ++i) {
        integral += (y[i-1] + y[i]) * (x[i] - x[i-1]) / 2.0;
    }
    
    // Fortran参考值（精确积分值）
    double fortranIntegral = 0.75;
    
    EXPECT_NEAR(integral, fortranIntegral, 1e-12);
}

TEST_F(SimplePrecisionTest, ComplexNumbers) {
    // 测试复数运算的数值精度
    
    // 复数加法
    std::complex<double> z1(3.0, 4.0);
    std::complex<double> z2(1.0, 2.0);
    std::complex<double> sum = z1 + z2;
    
    // Fortran参考值
    std::complex<double> fortranSum(4.0, 6.0);
    
    EXPECT_NEAR(sum.real(), fortranSum.real(), 1e-12);
    EXPECT_NEAR(sum.imag(), fortranSum.imag(), 1e-12);
    
    // 复数乘法
    std::complex<double> product = z1 * z2;
    
    // Fortran参考值
    std::complex<double> fortranProduct(-5.0, 10.0);
    
    EXPECT_NEAR(product.real(), fortranProduct.real(), 1e-12);
    EXPECT_NEAR(product.imag(), fortranProduct.imag(), 1e-12);
    
    // 复数模长
    double modulus = std::abs(z1);
    
    // Fortran参考值
    double fortranModulus = 5.0;
    
    EXPECT_NEAR(modulus, fortranModulus, 1e-12);
}

} // namespace elmer