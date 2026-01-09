#ifndef ELMER_CPP_INTERPOLATION_H
#define ELMER_CPP_INTERPOLATION_H

#include <vector>
#include <memory>
#include <optional>
#include "ElementDescription.h"
#include "Types.h"

namespace ElmerCpp {

/**
 * @brief 插值模块，提供有限元插值功能
 * 
 * 该模块包含在元素内部进行变量插值的核心功能，
 * 支持1D、2D和3D元素的插值操作。
 */
namespace Interpolation {

/**
 * @brief 在1D元素内部进行插值
 * 
 * @param element 元素结构
 * @param nodalValues 节点值数组
 * @param u 局部坐标u
 * @return double 插值结果
 */
double InterpolateInElement1D(const Element& element, 
                             const std::vector<double>& nodalValues, 
                             double u);

/**
 * @brief 在2D元素内部进行插值
 * 
 * @param element 元素结构
 * @param nodalValues 节点值数组
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @return double 插值结果
 */
double InterpolateInElement2D(const Element& element,
                             const std::vector<double>& nodalValues,
                             double u, double v);

/**
 * @brief 在3D元素内部进行插值
 * 
 * @param element 元素结构
 * @param nodalValues 节点值数组
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @return double 插值结果
 */
double InterpolateInElement3D(const Element& element,
                             const std::vector<double>& nodalValues,
                             double u, double v, double w);

/**
 * @brief 在元素内部进行插值（通用接口）
 * 
 * @param element 元素结构
 * @param nodalValues 节点值数组
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param basisFunctions 可选的基函数值数组
 * @return double 插值结果
 */
double InterpolateInElement(const Element& element,
                           const std::vector<double>& nodalValues,
                           double u, double v, double w,
                           const std::vector<double>& basisFunctions = {});

/**
 * @brief 检查点是否在元素内部
 * 
 * @param element 元素结构
 * @param elementNodes 元素节点坐标
 * @param point 待检查的点坐标
 * @param localCoordinates 输出局部坐标
 * @param globalEps 全局坐标精度
 * @param localEps 局部坐标精度
 * @param numericEps 数值精度
 * @param globalDistance 输出全局距离
 * @param localDistance 输出局部距离
 * @return bool 点是否在元素内部
 */
bool PointInElement(const Element& element,
                   const Nodes& elementNodes,
                   const std::vector<double>& point,
                   std::vector<double>& localCoordinates,
                   double globalEps = 1.0e-4,
                   double localEps = 1.0e-10,
                   double numericEps = -1.0,
                   double* globalDistance = nullptr,
                   double* localDistance = nullptr);

/**
 * @brief 象限树节点结构
 */
struct Quadrant {
    std::vector<double> boundingBox; // 边界框 [minX, minY, minZ, maxX, maxY, maxZ]
    std::vector<std::shared_ptr<Quadrant>> childQuadrants; // 子象限
    std::vector<int> elementIndices; // 包含的元素索引
    
    Quadrant() = default;
    explicit Quadrant(const std::vector<double>& bbox) : boundingBox(bbox) {}
};

/**
 * @brief 查找点所在的叶子象限
 * 
 * @param point 点坐标
 * @param dimension 维度
 * @param rootQuadrant 根象限
 * @return std::shared_ptr<Quadrant> 叶子象限指针
 */
std::shared_ptr<Quadrant> FindLeafElements(const std::vector<double>& point,
                                          int dimension,
                                          const std::shared_ptr<Quadrant>& rootQuadrant);

/**
 * @brief 构建象限树
 * 
 * @param mesh 网格结构
 * @param boundingBox 边界框
 * @return std::shared_ptr<Quadrant> 根象限指针
 */
std::shared_ptr<Quadrant> BuildQuadrantTree(const Mesh& mesh,
                                           const std::vector<double>& boundingBox);

/**
 * @brief 变量到变量的插值（简化版本）
 * 
 * @param oldMesh 原网格
 * @param newMesh 新网格
 * @param variableName 变量名
 * @param variableDimensions 变量维度
 * @param globalEps 全局精度
 * @param localEps 局部精度
 * @param numericEps 数值精度
 * @return bool 插值是否成功
 */
bool InterpolateVarToVarReduced(const Mesh& oldMesh,
                               Mesh& newMesh,
                               const std::string& variableName,
                               const std::vector<int>& variableDimensions,
                               double globalEps = 1.0e-4,
                               double localEps = 1.0e-10,
                               double numericEps = -1.0);

} // namespace Interpolation
} // namespace ElmerCpp

#endif // ELMER_CPP_INTERPOLATION_H