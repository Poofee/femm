/**
 * @file ElementDescription.h
 * @brief 元素描述模块
 * 
 * 对应Fortran模块：ElementDescription.F90
 * 提供有限元元素的类型定义、基函数描述等基础数据结构。
 */

#pragma once

#include "Types.h"
#include <vector>
#include <memory>

namespace ElmerCpp {

/**
 * @brief 基函数结构
 * 
 * 描述有限元基函数的系数和指数。
 */
struct BasisFunction {
    int n; ///< 项数
    std::vector<int> p; ///< u方向指数
    std::vector<int> q; ///< v方向指数  
    std::vector<int> r; ///< w方向指数
    std::vector<double> coeff; ///< 系数
    
    BasisFunction() : n(0) {}
    BasisFunction(int size) : n(size), p(size), q(size), r(size), coeff(size) {}
};

/**
 * @brief 元素类型结构
 * 
 * 描述有限元元素的类型信息，包括节点数、维度、基函数等。
 */
struct ElementType {
    int numberOfNodes; ///< 节点数量
    int dimension; ///< 维度
    int elementCode; ///< 元素代码
    std::vector<BasisFunction> basisFunctions; ///< 基函数数组
    
    ElementType() : numberOfNodes(0), dimension(0), elementCode(0) {}
};

/**
 * @brief 元素结构
 * 
 * 描述有限元元素的基本信息。
 */
struct Element {
    ElementType type; ///< 元素类型
    int index; ///< 元素索引
    
    Element() : index(0) {}
};

/**
 * @brief 节点坐标结构
 * 
 * 描述有限元节点的坐标信息。
 */
struct Nodes {
    std::vector<double> x; ///< x坐标数组
    std::vector<double> y; ///< y坐标数组
    std::vector<double> z; ///< z坐标数组
    
    Nodes() = default;
    Nodes(int size) : x(size), y(size), z(size) {}
    
    /**
     * @brief 获取节点坐标
     * @param index 节点索引
     * @return Node 节点对象
     */
    Node getNode(int index) const {
        if (index < 0 || index >= static_cast<int>(x.size())) {
            return Node();
        }
        return Node(x[index], y[index], z[index]);
    }
};

/**
 * @brief 网格结构
 * 
 * 描述有限元网格的基本信息。
 */
struct Mesh {
    int numberOfElements; ///< 元素数量
    int numberOfNodes; ///< 节点数量
    std::vector<Element> elements; ///< 元素数组
    std::vector<Node> nodes; ///< 节点对象数组
    
    Mesh() : numberOfElements(0), numberOfNodes(0) {}
};

} // namespace ElmerCpp