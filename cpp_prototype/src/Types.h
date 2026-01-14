/**
 * @file Types.h
 * @brief 基础类型定义
 * 
 * 定义ElmerCpp项目中使用的基础数据类型和结构。
 */

#pragma once

#include <vector>
#include <string>

namespace elmer {

/**
 * @brief 实数类型定义
 * 
 * 使用double作为实数类型，对应Fortran的REAL*8。
 */
using Real = double;

/**
 * @brief 整数类型定义
 * 
 * 使用int作为整数类型，对应Fortran的INTEGER。
 */
using Integer = int;

/**
 * @brief 节点结构
 * 
 * 描述有限元节点的基本信息。
 */
struct Node {
    double x; ///< x坐标
    double y; ///< y坐标
    double z; ///< z坐标
    
    Node() : x(0.0), y(0.0), z(0.0) {}
    Node(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

/**
 * @brief 矩阵结构
 * 
 * 描述矩阵的基本信息。
 */
struct Matrix {
    int rows; ///< 行数
    int cols; ///< 列数
    std::vector<double> values; ///< 矩阵值数组
    
    Matrix() : rows(0), cols(0) {}
    Matrix(int r, int c) : rows(r), cols(c), values(r * c, 0.0) {}
};

/**
 * @brief 向量结构
 * 
 * 描述向量的基本信息。
 */
struct Vector {
    int size; ///< 向量大小
    std::vector<double> values; ///< 向量值数组
    
    Vector() : size(0) {}
    Vector(int s) : size(s), values(s, 0.0) {}
};

/**
 * @brief 求解器结构
 * 
 * 描述求解器的基本信息。
 */
struct Solver {
    std::string name; ///< 求解器名称
    int maxIterations; ///< 最大迭代次数
    double tolerance; ///< 收敛容差
    
    Solver() : maxIterations(1000), tolerance(1.0e-6) {}
};

/**
 * @brief 模型结构
 * 
 * 描述有限元模型的基本信息。
 */
struct Model {
    std::string name; ///< 模型名称
    int numberOfVariables; ///< 变量数量
    
    Model() : numberOfVariables(0) {}
};

/**
 * @brief 变量结构
 * 
 * 描述变量的基本信息。
 */
struct Variable {
    std::string name; ///< 变量名称
    int dimension; ///< 变量维度
    std::vector<double> values; ///< 变量值数组
    
    Variable() : dimension(0) {}
};

} // namespace elmer