#include <iostream>
#include <cmath>
#include <vector>
#include "Interpolation.h"
#include "ElementDescription.h"
#include "ElementUtils.h"

using namespace elmer;
using namespace std;

/**
 * @brief 断言宏（简化版）
 */
#define ASSERT_NEAR(a, b, tolerance) \
    if (std::abs((a) - (b)) > (tolerance)) { \
        std::cout << "断言失败: " << (a) << " != " << (b) << " (容差: " << (tolerance) << ")" << std::endl; \
        return false; \
    }

#define ASSERT_TRUE(condition) \
    if (!(condition)) { \
        std::cout << "断言失败: 期望为真" << std::endl; \
        return false; \
    }

#define ASSERT_FALSE(condition) \
    if ((condition)) { \
        std::cout << "断言失败: 期望为假" << std::endl; \
        return false; \
    }

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::cout << "断言失败: " << (a) << " != " << (b) << std::endl; \
        return false; \
    }

#define ASSERT_NE(a, b) \
    if ((a) == (b)) { \
        std::cout << "断言失败: " << (a) << " == " << (b) << std::endl; \
        return false; \
    }

/**
 * @brief 测试1D插值功能
 */
bool Test1DInterpolation() {
    std::cout << "测试1D插值功能..." << std::endl;
    
    // 创建1D元素
    ElementTypeStruct elementType1D;
    elementType1D.numberOfNodes = 2;
    elementType1D.dimension = 1;
    elementType1D.elementCode = 101; // 线单元
    
    // 设置线性基函数
    BasisFunction basis1, basis2;
    basis1.n = 2;
    basis1.p = {0, 1};
    basis1.q = {0, 0};
    basis1.r = {0, 0};
    basis1.coeff = {1.0, -1.0}; // 1 - u
    
    basis2.n = 2;
    basis2.p = {0, 1};
    basis2.q = {0, 0};
    basis2.r = {0, 0};
    basis2.coeff = {0.0, 1.0};  // u
    
    elementType1D.basisFunctions = {basis1, basis2};
    
    Element element1D;
    element1D.type = elementType1D;
    
    // 测试节点值
    std::vector<double> nodalValues = {1.0, 2.0}; // 节点1=1.0, 节点2=2.0
    
    // 在中间点插值
    double result = Interpolation::InterpolateInElement1D(element1D, nodalValues, 0.5);
    
    // 线性插值：1.0 * (1-0.5) + 2.0 * 0.5 = 1.5
    ASSERT_NEAR(result, 1.5, 1e-12);
    
    std::cout << "1D插值测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试2D插值功能
 */
bool Test2DInterpolation() {
    std::cout << "测试2D插值功能..." << std::endl;
    
    // 创建测试元素类型
    ElementTypeStruct testElementType;
    testElementType.numberOfNodes = 4;
    testElementType.dimension = 2;
    testElementType.elementCode = 404; // 四边形单元
    
    // 设置双线性基函数
    BasisFunction basis1, basis2, basis3, basis4;
    
    // 基函数1: (1-u)(1-v)
    basis1.n = 4;
    basis1.p = {0, 1, 0, 1};
    basis1.q = {0, 0, 1, 1};
    basis1.r = {0, 0, 0, 0};
    basis1.coeff = {1.0, -1.0, -1.0, 1.0};
    
    // 基函数2: u(1-v)
    basis2.n = 4;
    basis2.p = {0, 1, 0, 1};
    basis2.q = {0, 0, 1, 1};
    basis2.r = {0, 0, 0, 0};
    basis2.coeff = {0.0, 1.0, 0.0, -1.0};
    
    // 基函数3: (1-u)v
    basis3.n = 4;
    basis3.p = {0, 1, 0, 1};
    basis3.q = {0, 0, 1, 1};
    basis3.r = {0, 0, 0, 0};
    basis3.coeff = {0.0, 0.0, 1.0, -1.0};
    
    // 基函数4: uv
    basis4.n = 4;
    basis4.p = {0, 1, 0, 1};
    basis4.q = {0, 0, 1, 1};
    basis4.r = {0, 0, 0, 0};
    basis4.coeff = {0.0, 0.0, 0.0, 1.0};
    
    testElementType.basisFunctions = {basis1, basis2, basis3, basis4};
    
    Element testElement;
    testElement.type = testElementType;
    
    // 测试节点值
    std::vector<double> nodalValues = {1.0, 2.0, 3.0, 4.0}; // 四个角点的值
    
    // 在中心点插值
    double result = Interpolation::InterpolateInElement2D(testElement, nodalValues, 0.5, 0.5);
    
    // 双线性插值：平均值应为2.5
    ASSERT_NEAR(result, 2.5, 1e-12);
    
    std::cout << "2D插值测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试通用插值接口
 */
bool TestGeneralInterpolation() {
    std::cout << "测试通用插值接口..." << std::endl;
    
    // 创建测试元素类型
    ElementTypeStruct testElementType;
    testElementType.numberOfNodes = 4;
    testElementType.dimension = 2;
    testElementType.elementCode = 404; // 四边形单元
    
    // 设置双线性基函数
    BasisFunction basis1, basis2, basis3, basis4;
    
    // 基函数1: (1-u)(1-v)
    basis1.n = 4;
    basis1.p = {0, 1, 0, 1};
    basis1.q = {0, 0, 1, 1};
    basis1.r = {0, 0, 0, 0};
    basis1.coeff = {1.0, -1.0, -1.0, 1.0};
    
    // 基函数2: u(1-v)
    basis2.n = 4;
    basis2.p = {0, 1, 0, 1};
    basis2.q = {0, 0, 1, 1};
    basis2.r = {0, 0, 0, 0};
    basis2.coeff = {0.0, 1.0, 0.0, -1.0};
    
    // 基函数3: (1-u)v
    basis3.n = 4;
    basis3.p = {0, 1, 0, 1};
    basis3.q = {0, 0, 1, 1};
    basis3.r = {0, 0, 0, 0};
    basis3.coeff = {0.0, 0.0, 1.0, -1.0};
    
    // 基函数4: uv
    basis4.n = 4;
    basis4.p = {0, 1, 0, 1};
    basis4.q = {0, 0, 1, 1};
    basis4.r = {0, 0, 0, 0};
    basis4.coeff = {0.0, 0.0, 0.0, 1.0};
    
    testElementType.basisFunctions = {basis1, basis2, basis3, basis4};
    
    Element testElement;
    testElement.type = testElementType;
    
    std::vector<double> nodalValues = {1.0, 2.0, 3.0, 4.0};
    
    // 使用通用接口进行2D插值
    double result = Interpolation::InterpolateInElement(testElement, nodalValues, 0.5, 0.5, 0.0);
    
    ASSERT_NEAR(result, 2.5, 1e-12);
    
    std::cout << "通用插值接口测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试点元素检测功能
 */
bool TestPointInElement() {
    std::cout << "测试点元素检测功能..." << std::endl;
    
    // 创建测试元素类型
    ElementTypeStruct testElementType;
    testElementType.numberOfNodes = 4;
    testElementType.dimension = 2;
    testElementType.elementCode = 404; // 四边形单元
    
    // 设置双线性基函数
    BasisFunction basis1, basis2, basis3, basis4;
    
    // 基函数1: (1-u)(1-v)
    basis1.n = 4;
    basis1.p = {0, 1, 0, 1};
    basis1.q = {0, 0, 1, 1};
    basis1.r = {0, 0, 0, 0};
    basis1.coeff = {1.0, -1.0, -1.0, 1.0};
    
    // 基函数2: u(1-v)
    basis2.n = 4;
    basis2.p = {0, 1, 0, 1};
    basis2.q = {0, 0, 1, 1};
    basis2.r = {0, 0, 0, 0};
    basis2.coeff = {0.0, 1.0, 0.0, -1.0};
    
    // 基函数3: (1-u)v
    basis3.n = 4;
    basis3.p = {0, 1, 0, 1};
    basis3.q = {0, 0, 1, 1};
    basis3.r = {0, 0, 0, 0};
    basis3.coeff = {0.0, 0.0, 1.0, -1.0};
    
    // 基函数4: uv
    basis4.n = 4;
    basis4.p = {0, 1, 0, 1};
    basis4.q = {0, 0, 1, 1};
    basis4.r = {0, 0, 0, 0};
    basis4.coeff = {0.0, 0.0, 0.0, 1.0};
    
    testElementType.basisFunctions = {basis1, basis2, basis3, basis4};
    
    Element testElement;
    testElement.type = testElementType;
    
    // 创建测试节点坐标
    Nodes testNodes;
    testNodes.x = {0.0, 1.0, 0.0, 1.0};
    testNodes.y = {0.0, 0.0, 1.0, 1.0};
    testNodes.z = {0.0, 0.0, 0.0, 0.0};
    
    std::vector<double> point = {0.5, 0.5, 0.0}; // 元素内部的点
    std::vector<double> localCoords;
    
    bool isInElement = Interpolation::PointInElement(testElement, testNodes, point, localCoords);
    
    ASSERT_TRUE(isInElement);
    ASSERT_EQ(localCoords.size(), 3);
    
    // 测试元素外部的点
    std::vector<double> outsidePoint = {2.0, 2.0, 0.0};
    bool isOutside = Interpolation::PointInElement(testElement, testNodes, outsidePoint, localCoords);
    
    ASSERT_FALSE(isOutside);
    
    std::cout << "点元素检测测试通过" << std::endl;
    return true;
}

/**
 * @brief 测试边界情况
 */
bool TestEdgeCases() {
    std::cout << "测试边界情况..." << std::endl;
    
    // 创建测试元素类型
    ElementTypeStruct testElementType;
    testElementType.numberOfNodes = 4;
    testElementType.dimension = 2;
    testElementType.elementCode = 404; // 四边形单元
    
    // 设置双线性基函数
    BasisFunction basis1, basis2, basis3, basis4;
    
    // 基函数1: (1-u)(1-v)
    basis1.n = 4;
    basis1.p = {0, 1, 0, 1};
    basis1.q = {0, 0, 1, 1};
    basis1.r = {0, 0, 0, 0};
    basis1.coeff = {1.0, -1.0, -1.0, 1.0};
    
    // 基函数2: u(1-v)
    basis2.n = 4;
    basis2.p = {0, 1, 0, 1};
    basis2.q = {0, 0, 1, 1};
    basis2.r = {0, 0, 0, 0};
    basis2.coeff = {0.0, 1.0, 0.0, -1.0};
    
    // 基函数3: (1-u)v
    basis3.n = 4;
    basis3.p = {0, 1, 0, 1};
    basis3.q = {0, 0, 1, 1};
    basis3.r = {0, 0, 0, 0};
    basis3.coeff = {0.0, 0.0, 1.0, -1.0};
    
    // 基函数4: uv
    basis4.n = 4;
    basis4.p = {0, 1, 0, 1};
    basis4.q = {0, 0, 1, 1};
    basis4.r = {0, 0, 0, 0};
    basis4.coeff = {0.0, 0.0, 0.0, 1.0};
    
    testElementType.basisFunctions = {basis1, basis2, basis3, basis4};
    
    Element testElement;
    testElement.type = testElementType;
    
    // 测试空节点值
    std::vector<double> emptyValues;
    double result = Interpolation::InterpolateInElement2D(testElement, emptyValues, 0.5, 0.5);
    ASSERT_NEAR(result, 0.0, 1e-12);
    
    // 测试不匹配的节点值数量
    std::vector<double> wrongSizeValues = {1.0, 2.0}; // 只有2个值，需要4个
    result = Interpolation::InterpolateInElement2D(testElement, wrongSizeValues, 0.5, 0.5);
    ASSERT_NEAR(result, 0.0, 1e-12);
    
    std::cout << "边界情况测试通过" << std::endl;
    return true;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始Interpolation模块测试..." << std::endl;
    
    bool allTestsPassed = true;
    
    // 运行所有测试
    allTestsPassed &= Test1DInterpolation();
    allTestsPassed &= Test2DInterpolation();
    allTestsPassed &= TestGeneralInterpolation();
    allTestsPassed &= TestPointInElement();
    allTestsPassed &= TestEdgeCases();
    
    if (allTestsPassed) {
        std::cout << "所有Interpolation模块测试通过!" << std::endl;
        return 0;
    } else {
        std::cout << "部分Interpolation模块测试失败!" << std::endl;
        return 1;
    }
}