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

/**
 * @brief 计算全局一阶导数（内部版本）
 * 
 * 给定元素结构，返回在局部坐标点(u,v,w)处，
 * 基于元素基函数计算的量f在全局坐标下的偏导数值。
 * 这是内部版本，通常不应直接由用户调用，而应通过包装例程GlobalFirstDerivatives调用。
 * 
 * @param elm 元素结构
 * @param nodes 元素节点坐标数组
 * @param df 节点值数组
 * @param gx 输出：f对x的偏导数
 * @param gy 输出：f对y的偏导数
 * @param gz 输出：f对z的偏导数
 * @param metric 逆变度量张量
 * @param dLBasisdx 基函数对局部坐标的偏导数
 */
void GlobalFirstDerivativesInternal(const Element& elm, const Nodes& nodes, 
                                   const std::vector<double>& df,
                                   double& gx, double& gy, double& gz,
                                   const std::vector<std::vector<double>>& metric,
                                   const std::vector<std::vector<double>>& dLBasisdx);

/**
 * @brief 计算全局一阶导数
 * 
 * 给定元素结构，返回在局部坐标点(u,v,w)处，
 * 基于元素基函数计算的量f在全局坐标下的偏导数值。
 * 
 * @param elm 元素结构
 * @param nodes 元素节点坐标数组
 * @param df 节点值数组
 * @param gx 输出：f对x的偏导数
 * @param gy 输出：f对y的偏导数
 * @param gz 输出：f对z的偏导数
 * @param metric 逆变度量张量
 * @param dLBasisdx 基函数对局部坐标的偏导数
 */
void GlobalFirstDerivatives(const Element& elm, const Nodes& nodes,
                           const std::vector<double>& df,
                           double& gx, double& gy, double& gz,
                           const std::vector<std::vector<double>>& metric,
                           const std::vector<std::vector<double>>& dLBasisdx);

/**
 * @brief 计算全局二阶导数
 * 
 * 在给定点(u,v,w)处计算全局坐标下的元素二阶偏导数矩阵。
 * 
 * @param elm 元素结构
 * @param nodes 元素节点坐标
 * @param values 输出：3x3二阶偏导数矩阵
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param metric 逆变度量张量
 * @param dBasisdx 基函数对局部坐标的一阶导数
 * @param ddLBasisddx 基函数对局部坐标的二阶导数
 * @param nd 节点数量
 */
void GlobalSecondDerivatives(const Element& elm, const Nodes& nodes,
                            std::vector<std::vector<std::vector<double>>>& values,
                            double u, double v, double w,
                            const std::vector<std::vector<double>>& metric,
                            const std::vector<std::vector<double>>& dBasisdx,
                            const std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                            int nd);

/**
 * @brief 获取坐标系维度
 * 
 * 返回当前坐标系的维度（1D、2D或3D）。
 * 
 * @return 坐标系维度（1、2或3）
 */
int CoordinateSystemDimension();

/**
 * @brief 棱柱单元P型基函数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有P型基函数值。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param phi 输出：基函数值数组（大小为6）
 */
void WedgeNodalPBasisAll(double u, double v, double w, std::vector<double>& phi);

/**
 * @brief 棱柱单元L型基函数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有L型基函数值。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param phi 输出：基函数值数组（大小为6）
 */
void WedgeNodalLBasisAll(double u, double v, double w, std::vector<double>& phi);

/**
 * @brief 棱柱单元P型基函数导数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有P型基函数导数。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param gradphi 输出：基函数导数矩阵（6x3）
 */
void dWedgeNodalPBasisAll(double u, double v, double w, std::vector<std::vector<double>>& gradphi);

/**
 * @brief 棱柱单元L型基函数导数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有L型基函数导数。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param gradphi 输出：基函数导数矩阵（6x3）
 */
void dWedgeNodalLBasisAll(double u, double v, double w, std::vector<std::vector<double>>& gradphi);

// =============================================================================
// P型元素支持函数
// =============================================================================

/**
 * @brief 检查元素是否为P型元素
 * 
 * 检查给定元素是否为P型（高阶）元素。
 * 
 * @param element 元素结构
 * @return true 如果是P型元素
 * @return false 如果不是P型元素
 */
bool isPElement(const Element& element);

/**
 * @brief 检查元素是否为活动P型元素
 * 
 * 检查给定元素在特定求解器中是否为活动的P型元素。
 * 
 * @param element 元素结构
 * @param solver 求解器指针（可选）
 * @return true 如果是活动的P型元素
 * @return false 如果不是活动的P型元素
 */
bool isActivePElement(const Element& element, void* solver = nullptr);

/**
 * @brief 获取参考P型元素节点
 * 
 * 获取P型元素的参考节点坐标。
 * 
 * @param element 元素类型结构
 * @param u 输出：u坐标数组
 * @param v 输出：v坐标数组
 * @param w 输出：w坐标数组
 */
void GetRefPElementNodes(const ElementType& element, std::vector<double>& u, 
                        std::vector<double>& v, std::vector<double>& w);

/**
 * @brief 线单元P型基函数（单个节点）
 * 
 * 计算线单元在给定局部坐标点u处的单个P型基函数值。
 * 
 * @param node 节点编号（1或2）
 * @param u 局部坐标u
 * @return 基函数值
 */
double LineNodalPBasis(int node, double u);

/**
 * @brief 线单元P型基函数（所有节点）
 * 
 * 计算线单元在给定局部坐标点u处的所有P型基函数值。
 * 
 * @param u 局部坐标u
 * @param phi 输出：基函数值数组（大小为2）
 */
void LineNodalPBasisAll(double u, std::vector<double>& phi);

/**
 * @brief 线单元P型基函数导数（单个节点）
 * 
 * 计算线单元在给定局部坐标点u处的单个P型基函数导数。
 * 
 * @param node 节点编号（1或2）
 * @param u 局部坐标u
 * @return 基函数导数值
 */
double dLineNodalPBasis(int node, double u);

/**
 * @brief 线单元P型基函数导数（所有节点）
 * 
 * 计算线单元在给定局部坐标点u处的所有P型基函数导数。
 * 
 * @param u 局部坐标u
 * @param gradphi 输出：基函数导数数组（大小为2）
 */
void dLineNodalPBasisAll(double u, std::vector<double>& gradphi);

/**
 * @brief 三角形单元P型基函数（所有节点）
 * 
 * 计算三角形单元在给定局部坐标点(u,v)处的所有P型基函数值。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param phi 输出：基函数值数组（大小为3）
 */
void TriangleNodalPBasisAll(double u, double v, std::vector<double>& phi);

/**
 * @brief 四边形单元P型基函数（所有节点）
 * 
 * 计算四边形单元在给定局部坐标点(u,v)处的所有P型基函数值。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param phi 输出：基函数值数组（大小为4）
 */
void QuadNodalPBasisAll(double u, double v, std::vector<double>& phi);

/**
 * @brief 四面体单元P型基函数（所有节点）
 * 
 * 计算四面体单元在给定局部坐标点(u,v,w)处的所有P型基函数值。
 * 
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param phi 输出：基函数值数组（大小为4）
 */
void TetraNodalPBasisAll(double u, double v, double w, std::vector<double>& phi);

/**
 * @brief 计算一维P型基函数
 * 
 * 基于Legendre多项式计算一维P型基函数。
 * 
 * @param basis 输出：基函数矩阵
 * @param n 基函数阶数
 */
void Compute1DPBasis(std::vector<std::vector<double>>& basis, int n);

// =============================================================================
// 边界处理函数
// =============================================================================

/**
 * @brief 获取三角形面的全局方向
 * 
 * 给定元素和其面映射（对于元素的某个三角形面），
 * 此例程返回三角形面的全局方向，以便函数在元素边界上连续。
 * 
 * @param element 元素结构
 * @param faceMap 元素三角形面映射（3个整数）
 * @param indexes 全局节点索引数组
 * @return 三角形面的全局方向（局部节点编号）
 */
std::vector<int> getTriangleFaceDirection(const Element& element, 
                                         const std::vector<int>& faceMap,
                                         const std::vector<int>& indexes);

/**
 * @brief 获取四边形面的全局方向
 * 
 * 给定元素和其面映射（对于元素的某个四边形面），
 * 此例程返回四边形面的全局方向，以便函数在元素边界上连续。
 * 
 * @param element 元素结构
 * @param faceMap 元素四边形面映射（4个整数）
 * @param indexes 全局节点索引数组
 * @return 四边形面的全局方向（局部节点编号）
 */
std::vector<int> getSquareFaceDirection(const Element& element,
                                       const std::vector<int>& faceMap,
                                       const std::vector<int>& indexes);

/**
 * @brief 检查楔形元素面的局部编号是否合法
 * 
 * 检查给定的四边形面的局部编号对于楔形元素是否合法。
 * 
 * @param ordering 楔形元素四边形面的局部编号
 * @return true 如果给定的编号对于楔形元素四边形面是合法的
 */
bool wedgeOrdering(const std::vector<int>& ordering);

} // namespace ElmerCpp