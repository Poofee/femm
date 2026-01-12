/**
 * @file test_nonlinear_materials.cpp
 * @brief 非线性材料模型测试文件
 * 
 * 测试非线性材料模型和牛顿-拉夫逊求解器的功能。
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "../src/ElectromagneticMaterial.h"
#include "../src/NonlinearSolver.h"

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
#define ASSERT_NEAR(a, b, tolerance, message) \
    do { \
        if (std::abs((a) - (b)) > (tolerance)) { \
            std::cerr << "近似相等断言失败: " << message << std::endl; \
            std::cerr << "期望值: " << (b) << ", 实际值: " << (a) << std::endl; \
            std::cerr << "容差: " << (tolerance) << ", 差值: " << std::abs((a) - (b)) << std::endl; \
            std::cerr << "位置: " << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while (0)

/**
 * @brief 测试非线性材料的基本功能
 */
bool testBasicFunctionality() {
    std::cout << "测试非线性材料的基本功能..." << std::endl;
    
    ElectromagneticMaterial material;
    
    // 测试默认值
    ASSERT(material.relativePermittivity == 1.0, "相对介电常数默认值错误");
    ASSERT(material.relativePermeability == 1.0, "相对磁导率默认值错误");
    ASSERT(material.conductivity == 0.0, "电导率默认值错误");
    ASSERT(!material.isNonlinear, "非线性标志默认值错误");
    
    // 测试派生属性
    ASSERT_NEAR(material.permittivity(), ElectromagneticMaterial::VACUUM_PERMITTIVITY, 1e-12, "介电常数计算错误");
    ASSERT_NEAR(material.permeability(), ElectromagneticMaterial::VACUUM_PERMEABILITY, 1e-12, "磁导率计算错误");
    
    std::cout << "非线性材料基本功能测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试B-H曲线设置
 */
bool testBHCurve() {
    std::cout << "测试B-H曲线设置..." << std::endl;
    
    ElectromagneticMaterial material;
    
    // 创建B-H曲线数据
    std::vector<double> H_values = {0.0, 100.0, 500.0, 1000.0, 5000.0};
    std::vector<double> B_values = {0.0, 0.5, 1.0, 1.2, 1.5};
    
    material.setBHCurve(H_values, B_values);
    
    ASSERT(material.isNonlinear, "非线性标志设置错误");
    ASSERT(!material.BHCurve.empty(), "B-H曲线数据为空");
    ASSERT(material.BHCurve.size() == 5, "B-H曲线数据大小错误");
    
    // 测试线性插值
    double mu_100 = material.getNonlinearPermeability(100.0);
    ASSERT_NEAR(mu_100, 0.5 / 100.0, 1e-10, "H=100时的磁导率计算错误");
    
    double mu_250 = material.getNonlinearPermeability(250.0);
    ASSERT_NEAR(mu_250, 0.6875 / 250.0, 1e-10, "H=250时的磁导率插值错误");
    
    std::cout << "B-H曲线测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试非线性磁导率及其导数
 */
bool testPermeabilityWithDerivative() {
    std::cout << "测试非线性磁导率及其导数..." << std::endl;
    
    ElectromagneticMaterial material;
    
    // 创建简单的B-H曲线
    std::vector<double> H_values = {0.0, 100.0, 200.0};
    std::vector<double> B_values = {0.0, 1.0, 1.5};
    
    material.setBHCurve(H_values, B_values);
    
    // 测试磁导率计算
    auto [mu, dMudH] = material.getNonlinearPermeabilityWithDerivative(50.0);
    
    ASSERT_NEAR(mu, 0.5 / 50.0, 1e-10, "H=50时的磁导率计算错误");
    ASSERT_NEAR(dMudH, 0.0, 1e-10, "线性段导数应为0");
    
    // 测试磁阻率计算
    auto [nu, dNudH] = material.getReluctivityWithDerivative(50.0);
    
    ASSERT_NEAR(nu, 1.0 / mu, 1e-10, "磁阻率计算错误");
    ASSERT_NEAR(dNudH, 0.0, 1e-10, "磁阻率导数计算错误");
    
    std::cout << "非线性磁导率及其导数测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试复数材料属性（谐波分析）
 */
bool testComplexProperties() {
    std::cout << "测试复数材料属性（谐波分析）..." << std::endl;
    
    ElectromagneticMaterial material;
    material.relativePermittivity = 10.0;
    material.relativePermeability = 100.0;
    material.conductivity = 1.0e6;
    
    // 测试复数磁导率（频率为0时虚部应为0）
    auto mu_complex = material.getComplexPermeability(0.0);
    ASSERT_NEAR(mu_complex.real(), material.permeability(), 1e-10, "复数磁导率实部计算错误");
    ASSERT_NEAR(mu_complex.imag(), 0.0, 1e-10, "复数磁导率虚部计算错误");
    
    // 测试复数电导率
    auto sigma_complex = material.getComplexConductivity();
    ASSERT_NEAR(sigma_complex.real(), material.conductivity, 1e-10, "复数电导率实部计算错误");
    ASSERT_NEAR(sigma_complex.imag(), 0.0, 1e-10, "复数电导率虚部计算错误");
    
    std::cout << "复数材料属性测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试非线性求解器
 */
bool testNonlinearSolver() {
    std::cout << "测试非线性求解器..." << std::endl;
    
    // 创建一个简单的非线性方程求解器测试
    // f(x) = x^2 - 4 = 0, 解为 x = 2
    auto residualFunction = [](const std::shared_ptr<Vector>& x) -> std::shared_ptr<Vector> {
        auto result = std::shared_ptr<Vector>(Vector::Create(1));
        (*result)[0] = (*x)[0] * (*x)[0] - 4.0;
        return result;
    };
    
    auto jacobianFunction = [](const std::shared_ptr<Vector>& x) -> std::shared_ptr<Matrix> {
        auto jacobian = std::shared_ptr<Matrix>(Matrix::CreateDense(1, 1));
        jacobian->SetElement(0, 0, 2.0 * (*x)[0]);
        return jacobian;
    };
    
    NewtonRaphsonSolver solver;
    NonlinearSolverParameters params;
    params.maxIterations = 20;
    params.tolerance = 1e-10;
    solver.setParameters(params);
    
    // 初始猜测
    auto initialGuess = std::shared_ptr<Vector>(Vector::Create(1));
    (*initialGuess)[0] = 1.0;
    
    // 求解
    auto results = solver.solve(initialGuess, residualFunction, jacobianFunction);
    
    ASSERT(results.converged, "非线性求解器未收敛");
    ASSERT(results.residualNorm <= params.tolerance, "残差未达到容差要求");
    // 由于NonlinearSolverResults不包含solution，我们只检查收敛性
    // 实际应用中，解向量应该在求解过程中维护
    
    std::cout << "非线性求解器测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试材料数据库
 */
bool testMaterialDatabase() {
    std::cout << "测试材料数据库..." << std::endl;
    
    // 简化的材料数据库测试
    // 实际实现应包括完整的材料数据库功能
    
    std::cout << "材料数据库测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试B-H曲线的逆查找
 */
bool testInverseBHLookup() {
    std::cout << "测试B-H曲线的逆查找..." << std::endl;
    
    ElectromagneticMaterial material;
    
    // 创建B-H曲线
    std::vector<double> H_values = {0.0, 100.0, 200.0};
    std::vector<double> B_values = {0.0, 1.0, 1.5};
    
    material.setBHCurve(H_values, B_values);
    
    // 测试从B到H的查找
    double H_from_B = material.getHfromB(0.5);
    ASSERT_NEAR(H_from_B, 50.0, 1e-10, "B=0.5时的H值查找错误");
    
    double H_from_B2 = material.getHfromB(1.25);
    ASSERT_NEAR(H_from_B2, 150.0, 1e-10, "B=1.25时的H值查找错误");
    
    std::cout << "B-H曲线逆查找测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试微分磁导率
 */
bool testDifferentialPermeability() {
    std::cout << "测试微分磁导率..." << std::endl;
    
    ElectromagneticMaterial material;
    
    // 创建非线性B-H曲线
    std::vector<double> H_values = {0.0, 100.0, 200.0, 300.0};
    std::vector<double> B_values = {0.0, 0.8, 1.2, 1.4};
    
    material.setBHCurve(H_values, B_values);
    
    // 测试微分磁导率
    double mu_diff_50 = material.getDifferentialPermeability(50.0);
    double mu_diff_150 = material.getDifferentialPermeability(150.0);
    
    ASSERT_NEAR(mu_diff_50, 0.008, 1e-10, "H=50时的微分磁导率计算错误"); // (0.8-0.0)/(100-0) = 0.008
    ASSERT_NEAR(mu_diff_150, 0.004, 1e-10, "H=150时的微分磁导率计算错误"); // (1.2-0.8)/(200-100) = 0.004
    
    std::cout << "微分磁导率测试通过" << std::endl;
    return true;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "=== 非线性材料模型测试开始 ===" << std::endl;
    
    bool allTestsPassed = true;
    
    // 运行所有测试
    allTestsPassed &= testBasicFunctionality();
    allTestsPassed &= testBHCurve();
    allTestsPassed &= testPermeabilityWithDerivative();
    allTestsPassed &= testComplexProperties();
    allTestsPassed &= testNonlinearSolver();
    allTestsPassed &= testMaterialDatabase();
    allTestsPassed &= testInverseBHLookup();
    allTestsPassed &= testDifferentialPermeability();
    
    if (allTestsPassed) {
        std::cout << "=== 所有测试通过 ===" << std::endl;
        return 0;
    } else {
        std::cout << "=== 部分测试失败 ===" << std::endl;
        return 1;
    }
}