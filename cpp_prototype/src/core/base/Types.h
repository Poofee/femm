/**
 * @file Types.h
 * @brief 基础类型定义
 * 
 * 定义ElmerCpp项目中使用的基础数据类型和结构。
 */

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <cmath>

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

// Node structure is now defined in Mesh.h to maintain compatibility with Elmer Fortran code

// 前向声明，解决循环依赖问题
class Vector;

/**
 * @brief 矩阵接口类，用于线性代数操作
 */
class Matrix {
public:
    virtual ~Matrix() = default;
    
    virtual void Zero() = 0;
    virtual void SetElement(Integer i, Integer j, Real value) = 0;
    virtual Real GetElement(Integer i, Integer j) const = 0;
    virtual void AddToElement(Integer i, Integer j, Real value) = 0;
    
    // 矩阵大小信息
    virtual Integer GetNumRows() const = 0;
    virtual Integer GetNumCols() const = 0;
    
    // 矩阵-向量乘法
    virtual void Multiply(const Vector& x, Vector& y) const = 0;
    
    // 矩阵-向量乘法（带缩放）：y = alpha * A * x + beta * y
    virtual void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const = 0;
    
    // 工厂方法：创建不同类型的矩阵
    static std::unique_ptr<Matrix> CreateCRS(Integer nrows, Integer ncols);
    static std::unique_ptr<Matrix> CreateBand(Integer nrows, Integer bandwidth);
    static std::unique_ptr<Matrix> CreateDense(Integer nrows, Integer ncols);
};

/**
 * @brief 向量接口类，用于线性代数操作
 */
class Vector {
public:
    virtual ~Vector() = default;
    
    virtual Integer Size() const = 0;
    virtual Real& operator[](Integer i) = 0;
    virtual const Real& operator[](Integer i) const = 0;
    virtual void Zero() = 0;
    
    // 工厂方法：创建向量
    static std::unique_ptr<Vector> Create(Integer size);
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

/**
 * @brief 创建CRS矩阵的工厂函数
 */
std::unique_ptr<Matrix> CreateCRSMatrix(Integer nrows, Integer ncols);

} // namespace elmer