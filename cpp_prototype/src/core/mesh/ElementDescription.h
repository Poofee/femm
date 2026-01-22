/**
 * @file ElementDescription.h
 * @brief 元素描述模块
 * 
 * 对应Fortran模块：ElementDescription.F90
 * 提供有限元元素的类型定义、基函数描述等基础数据结构。
 */

#pragma once

#include "Types.h"
#include "BoundaryConditions.h"
#include "Mesh.h"
#include <vector>
#include <memory>

namespace elmer {

/**
 * @brief 边界类型枚举
 * 
 * 定义有限元边界的不同类型。
 */
enum class BoundaryType {
    UNKNOWN,    ///< 未知边界类型
    POINT,      ///< 点边界
    LINE,       ///< 线边界
    TRIANGLE,   ///< 三角形边界
    QUADRILATERAL ///< 四边形边界
};

/**
 * @brief 连续性类型枚举
 * 
 * 定义边界上的连续性要求。
 */
enum class ContinuityType {
    C0,         ///< C0连续性（函数连续）
    C1          ///< C1连续性（一阶导数连续）
};

/**
 * @brief 边界连续性结构
 * 
 * 描述边界上的连续性要求。
 */
struct BoundaryContinuity {
    BoundaryType type;       ///< 边界类型
    ContinuityType continuity; ///< 连续性要求
    
    BoundaryContinuity() : type(BoundaryType::POINT), continuity(ContinuityType::C0) {}
    BoundaryContinuity(BoundaryType t, ContinuityType c) : type(t), continuity(c) {}
};

/**
 * @brief 简单边界条件结构
 * 
 * 描述边界条件的简化信息，用于元素描述模块。
 */
struct SimpleBoundaryCondition {
    int boundaryIndex;          ///< 边界索引
    BoundaryConditionType type; ///< 边界条件类型
    double value;               ///< 边界值
    
    SimpleBoundaryCondition() : boundaryIndex(0), type(BoundaryConditionType::DIRICHLET), value(0.0) {}
    SimpleBoundaryCondition(int idx, BoundaryConditionType t, double v) : 
        boundaryIndex(idx), type(t), value(v) {}
};

// BasisFunction is now defined in Mesh.h to maintain compatibility with Elmer Fortran code

// ElementTypeStruct is now defined in Mesh.h to maintain compatibility with Elmer Fortran code





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
void GetRefPElementNodes(const ElementTypeStruct& element, std::vector<double>& u_coords, 
                        std::vector<double>& v_coords, std::vector<double>& w_coords);

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

/**
 * @brief 获取元素网格边信息
 * 
 * 获取元素的边度数和边方向信息，用于P型元素的边基函数计算。
 * 
 * @param mesh 网格结构
 * @param element 元素结构
 * @param edgeDegree 输出：每条边的度数
 * @param edgeDirection 输出：每条边的方向映射
 * @param edgeMaxDegree 输出：最大边度数
 */
void GetElementMeshEdgeInfo(const Mesh& mesh, const Element& element,
                           std::vector<int>& edgeDegree,
                           std::vector<std::vector<int>>& edgeDirection,
                           int& edgeMaxDegree);

/**
 * @brief 获取元素网格面信息
 * 
 * 获取元素的面度数和面方向信息，用于P型元素的面基函数计算。
 * 
 * @param mesh 网格结构
 * @param element 元素结构
 * @param faceDegree 输出：每个面的度数
 * @param faceDirection 输出：每个面的方向映射
 * @param faceMaxDegree 输出：最大面度数
 */
void GetElementMeshFaceInfo(const Mesh& mesh, const Element& element,
                           std::vector<int>& faceDegree,
                           std::vector<std::vector<int>>& faceDirection,
                           int& faceMaxDegree);

/**
 * @brief 面元素基函数排序
 * 
 * 产生关于面（向量）元素的基函数如何重新排序以符合全局排序约定的信息。
 * 
 * @param element 具有2D面的3D元素
 * @param fDofMap 输出：面基函数排列的面信息
 * @param faceIndex 可选：指定要检查的面索引
 * @param reverseSign 可选输出：关于符号反转的面信息
 */
void FaceElementBasisOrdering(const Element& element,
                             std::vector<std::vector<int>>& fDofMap,
                             int faceIndex = -1,
                             std::vector<bool>* reverseSign = nullptr);

/**
 * @brief 获取边基函数
 * 
 * 计算元素的边基函数及其旋转版本。
 * 
 * @param element 元素结构
 * @param wBasis 输出：边基函数值
 * @param rotWBasis 输出：旋转边基函数值
 * @param basis 输出：标量基函数值
 * @param dBasisdx 输出：基函数导数
 */
void GetEdgeBasis(const Element& element,
                 std::vector<double>& wBasis,
                 std::vector<double>& rotWBasis,
                 std::vector<double>& basis,
                 std::vector<std::vector<double>>& dBasisdx);

/**
 * @brief 检查面节点连续性
 * 
 * 检查两个相邻元素在共享面上的节点是否连续。
 * 
 * @param element1 第一个元素
 * @param element2 第二个元素
 * @param sharedFace 共享面索引
 * @param mesh 网格结构
 * @param tolerance 容差
 * @return true 如果节点连续
 * @return false 如果节点不连续
 */
bool CheckFaceNodeContinuity(const Element& element1,
                            const Element& element2,
                            int sharedFace,
                            const Mesh& mesh,
                            double tolerance = 1e-12);

/**
 * @brief 检查面几何连续性
 * 
 * 检查相邻元素在共享面上的几何连续性。
 * 
 * @param element1 第一个元素
 * @param element2 第二个元素
 * @param sharedFace 共享面索引
 * @param mesh 网格结构
 * @param nodes 节点坐标
 * @param tolerance 容差
 * @return true 如果几何连续
 * @return false 如果几何不连续
 */
bool CheckFaceGeometricContinuity(const Element& element1,
                                 const Element& element2,
                                 int sharedFace,
                                 const Mesh& mesh,
                                 const Nodes& nodes,
                                 double tolerance = 1e-12);

/**
 * @brief 获取元素边界连续性信息
 * 
 * 获取元素边界连续性信息，包括边界类型和连续性条件。
 * 
 * @param element 元素类型结构
 * @param continuity 输出：边界连续性信息数组
 */
void GetElementBoundaryContinuity(const ElementTypeStruct& element, std::vector<BoundaryContinuity>& continuity);

/**
 * @brief 获取边界节点
 * 
 * 获取元素边界上的节点索引。
 * 
 * @param element 元素类型结构
 * @param boundaryIndex 边界索引
 * @param boundaryNodes 输出：边界节点索引数组
 */
void GetBoundaryNodes(const ElementTypeStruct& element, int boundaryIndex, std::vector<int>& boundaryNodes);

/**
 * @brief 检查面函数连续性
 * 
 * 检查相邻元素在共享面上的函数连续性。
 * 
 * @param element1 第一个元素
 * @param element2 第二个元素
 * @param sharedFace 共享面索引
 * @param mesh 网格结构
 * @param nodes 节点坐标
 * @param field1 第一个元素的场值
 * @param field2 第二个元素的场值
 * @param tolerance 容差
 * @return true 如果函数连续
 * @return false 如果函数不连续
 */
bool CheckFaceFunctionContinuity(const Element& element1,
                                const Element& element2,
                                int sharedFace,
                                const Mesh& mesh,
                                const Nodes& nodes,
                                const std::vector<double>& field1,
                                const std::vector<double>& field2,
                                double tolerance = 1e-12);

/**
 * @brief 边界条件处理函数
 * 
 * 处理元素边界上的边界条件，包括Dirichlet和Neumann边界条件。
 * 
 * @param element 元素结构
 * @param boundaryType 边界类型（0: Dirichlet, 1: Neumann）
 * @param boundaryValue 边界值
 * @param boundaryNodes 边界节点索引
 * @param boundaryValues 输出：边界节点值
 */
void ProcessBoundaryConditions(const Element& element,
                              int boundaryType,
                              double boundaryValue,
                              const std::vector<int>& boundaryNodes,
                              std::vector<double>& boundaryValues);

/**
 * @brief 检查元素面连续性
 * 
 * 检查相邻元素在共享面上的连续性条件。
 * 
 * @param element1 第一个元素
 * @param element2 第二个元素
 * @param sharedFace 共享面索引
 * @param tolerance 容差
 * @return true 如果面连续
 * @return false 如果面不连续
 */
bool CheckFaceContinuity(const Element& element1,
                        const Element& element2,
                        int sharedFace,
                        double tolerance = 1e-12);

/**
 * @brief 获取面法向量
 * 
 * 计算元素面的法向量。
 * 
 * @param element 元素结构
 * @param faceIndex 面索引
 * @param nodes 节点坐标
 * @param normal 输出：法向量
 */
void GetFaceNormal(const Element& element,
                  int faceIndex,
                  const Nodes& nodes,
                  std::vector<double>& normal);

/**
 * @brief 计算面通量
 * 
 * 计算通过元素面的通量。
 * 
 * @param element 元素结构
 * @param faceIndex 面索引
 * @param nodes 节点坐标
 * @param field 场值
 * @param flux 输出：通量值
 */
void ComputeFaceFlux(const Element& element,
                    int faceIndex,
                    const Nodes& nodes,
                    const std::vector<double>& field,
                    double& flux);

/**
 * @brief 应用Dirichlet边界条件
 * 
 * 在系统矩阵和右端项中应用Dirichlet边界条件。
 * 
 * @param element 元素结构
 * @param boundaryNodes 边界节点索引
 * @param boundaryValues 边界值
 * @param stiffnessMatrix 刚度矩阵
 * @param forceVector 力向量
 */
void ApplyDirichletBoundaryConditions(const Element& element,
                                     const std::vector<int>& boundaryNodes,
                                     const std::vector<double>& boundaryValues,
                                     std::vector<std::vector<double>>& stiffnessMatrix,
                                     std::vector<double>& forceVector);

/**
 * @brief 应用Neumann边界条件
 * 
 * 在系统右端项中应用Neumann边界条件。
 * 
 * @param element 元素结构
 * @param faceIndex 面索引
 * @param boundaryValue 边界值
 * @param forceVector 力向量
 */
void ApplyNeumannBoundaryConditions(const Element& element,
                                   int faceIndex,
                                   double boundaryValue,
                                   std::vector<double>& forceVector);

/**
 * @brief 3D节点基函数计算
 * 
 * 对应Fortran函数：NodalBasisFunctions3D
 * 计算给定局部坐标点(u,v,w)处所有节点的基函数值。
 * 
 * @param element 元素结构
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param y 输出：基函数值数组
 */
void NodalBasisFunctions3D(const Element& element, double u, double v, double w, 
                          std::vector<double>& y);

/**
 * @brief 计算3D节点基函数导数
 * 
 * 计算给定局部坐标点(u,v,w)处所有节点的基函数导数。
 * 
 * @param element 元素结构
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param dphi_du 输出：基函数对u的导数
 * @param dphi_dv 输出：基函数对v的导数
 * @param dphi_dw 输出：基函数对w的导数
 */
void NodalBasisFunctions3DDerivatives(const Element& element, double u, double v, double w,
                                     std::vector<double>& dphi_du, std::vector<double>& dphi_dv, 
                                     std::vector<double>& dphi_dw);

/**
 * @brief 元素信息计算接口
 * 
 * 对应Fortran函数：ElementInfo
 * 计算给定局部坐标点(u,v,w)处的元素信息，包括基函数值、导数、雅可比矩阵等。
 * 
 * @param element 元素结构
 * @param nodes 节点坐标
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param detJ 输出：雅可比矩阵行列式的平方根
 * @param basis 输出：基函数值数组
 * @param dBasisdx 输出：基函数全局导数（可选）
 * @param ddBasisddx 输出：基函数全局二阶导数（可选）
 * @param secondDerivatives 输入：是否需要计算二阶导数（可选）
 * @param bubbles 输入：是否需要计算气泡基函数（可选）
 * @param basisDegree 输出：基函数度数（可选）
 * @param edgeBasis 输出：H(curl)相容边基函数（可选）
 * @param rotBasis 输出：边基函数的旋度（可选）
 * @param solver 输入：求解器指针（可选）
 * @return 是否成功计算
 */
bool ElementInfo(const Element& element, const Nodes& nodes, 
                double u, double v, double w,
                double& detJ, std::vector<double>& basis,
                std::vector<std::vector<double>>* dBasisdx = nullptr,
                std::vector<std::vector<std::vector<double>>>* ddBasisddx = nullptr,
                bool* secondDerivatives = nullptr,
                bool* bubbles = nullptr,
                std::vector<int>* basisDegree = nullptr,
                std::vector<std::vector<double>>* edgeBasis = nullptr,
                std::vector<std::vector<double>>* rotBasis = nullptr,
                void* solver = nullptr);

/**
 * @brief 处理线形P型元素
 */
bool processLinePElement(const Element& element, const Nodes& nodes,
                        double u, double v, double w, int bodyId,
                        std::vector<double>& basis, 
                        std::vector<std::vector<double>>& dLBasisdx,
                        std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                        bool compute2ndDerivatives,
                        std::vector<int>* basisDegree,
                        void* solver);

/**
 * @brief 处理三角形P型元素
 */
bool processTrianglePElement(const Element& element, const Nodes& nodes,
                           double u, double v, double w, int bodyId,
                           std::vector<double>& basis, 
                           std::vector<std::vector<double>>& dLBasisdx,
                           std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                           bool compute2ndDerivatives,
                           std::vector<int>* basisDegree,
                           void* solver);

/**
 * @brief 处理四边形P型元素
 */
bool processQuadPElement(const Element& element, const Nodes& nodes,
                        double u, double v, double w, int bodyId,
                        std::vector<double>& basis, 
                        std::vector<std::vector<double>>& dLBasisdx,
                        std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                        bool compute2ndDerivatives,
                        std::vector<int>* basisDegree,
                        void* solver);

/**
 * @brief 处理四面体P型元素
 */
bool processTetraPElement(const Element& element, const Nodes& nodes,
                         double u, double v, double w, int bodyId,
                         std::vector<double>& basis, 
                         std::vector<std::vector<double>>& dLBasisdx,
                         std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                         bool compute2ndDerivatives,
                         std::vector<int>* basisDegree,
                         void* solver);

/**
 * @brief 处理金字塔P型元素
 */
bool processPyramidPElement(const Element& element, const Nodes& nodes,
                           double u, double v, double w, int bodyId,
                           std::vector<double>& basis, 
                           std::vector<std::vector<double>>& dLBasisdx,
                           std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                           bool compute2ndDerivatives,
                           std::vector<int>* basisDegree,
                           void* solver);

/**
 * @brief 处理楔形P型元素
 */
bool processWedgePElement(const Element& element, const Nodes& nodes,
                         double u, double v, double w, int bodyId,
                         std::vector<double>& basis, 
                         std::vector<std::vector<double>>& dLBasisdx,
                         std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                         bool compute2ndDerivatives,
                         std::vector<int>* basisDegree,
                         void* solver);

/**
 * @brief 处理六面体P型元素
 */
bool processBrickPElement(const Element& element, const Nodes& nodes,
                         double u, double v, double w, int bodyId,
                         std::vector<double>& basis, 
                         std::vector<std::vector<double>>& dLBasisdx,
                         std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                         bool compute2ndDerivatives,
                         std::vector<int>* basisDegree,
                         void* solver);

/**
 * @brief 计算元素度量张量
 */
bool ElementMetric(int nBasis, const Element& element, const Nodes& nodes,
                  std::vector<std::vector<double>>& elmMetric, double& detJ,
                  const std::vector<std::vector<double>>& dLBasisdx,
                  std::vector<std::vector<double>>& ltoGMap);



/**
 * @brief 获取元素P值
 */
int getElementP(const Element& element, int bodyId, void* solver = nullptr);

/**
 * @brief 获取边自由度数量
 */
int getEdgeDOFs(const Element& element, int p);

/**
 * @brief 获取气泡自由度数量
 */
int getBubbleDOFs(const Element& element, int p);

/**
 * @brief 获取有效的气泡P值
 */
int getEffectiveBubbleP(const Element& element, int p, int bDofs);

/**
 * @brief 获取三角形边映射
 */
std::vector<int> getTriangleEdgeMap(int edgeIndex);

/**
 * @brief 线形气泡P型基函数
 */
double LineBubblePBasis(int i, double u, bool invert);

/**
 * @brief 线形气泡P型基函数导数
 */
double dLineBubblePBasis(int i, double u, bool invert);

/**
 * @brief 线形气泡P型基函数二阶导数
 */
double ddLineBubblePBasis(int i, double u, bool invert);

/**
 * @brief 三角形边P型基函数
 */
double TriangleEdgePBasis(int edge, int k, double u, double v, bool invert);

/**
 * @brief 三角形边P型基函数导数
 */
std::vector<double> dTriangleEdgePBasis(int edge, int k, double u, double v, bool invert);

/**
 * @brief 三角形边P型基函数二阶导数
 */
std::vector<std::vector<double>> ddTriangleEdgePBasis(int edge, int k, double u, double v, bool invert);

/**
 * @brief 三角形边界气泡P型基函数
 */
double TriangleEBubblePBasis(int i, int j, double u, double v, const std::vector<int>& direction);

/**
 * @brief 三角形边界气泡P型基函数导数
 */
std::vector<double> dTriangleEBubblePBasis(int i, int j, double u, double v, const std::vector<int>& direction);

/**
 * @brief 三角形边界气泡P型基函数二阶导数
 */
std::vector<std::vector<double>> ddTriangleEBubblePBasis(int i, int j, double u, double v, const std::vector<int>& direction);

/**
 * @brief 三角形气泡P型基函数
 */
double TriangleBubblePBasis(int i, int j, double u, double v);

/**
 * @brief 三角形气泡P型基函数导数
 */
std::vector<double> dTriangleBubblePBasis(int i, int j, double u, double v);

/**
 * @brief 三角形气泡P型基函数二阶导数
 */
std::vector<std::vector<double>> ddTriangleBubblePBasis(int i, int j, double u, double v);

/**
 * @brief 处理气泡基函数
 */
bool processBubbleBasis(const Element& element, const Nodes& nodes,
                       double u, double v, double w, double detJ,
                       std::vector<double>& basis, 
                       std::vector<std::vector<double>>& dBasisdx,
                       int cdim);

/**
 * @brief 边元素信息计算接口
 * 
 * 对应Fortran函数：EdgeElementInfo
 * 计算给定局部坐标点(u,v,w)处的边元素信息，包括H(curl)相容边基函数、旋度等。
 * 
 * @param element 元素结构
 * @param nodes 节点坐标
 * @param u 局部坐标u
 * @param v 局部坐标v
 * @param w 局部坐标w
 * @param F 输出：梯度矩阵F=Grad f（可选）
 * @param G 输出：梯度矩阵逆的转置（可选）
 * @param detF 输出：梯度矩阵行列式
 * @param basis 输出：H1相容基函数值
 * @param edgeBasis 输出：H(curl)相容边基函数
 * @param rotBasis 输出：边基函数的旋度（可选）
 * @param dBasisdx 输出：H1基函数的全局导数（可选）
 * @param secondFamily 输入：是否使用第二类Nedelec基函数（可选）
 * @param basisDegree 输入：近似度数（可选）
 * @param applyPiolaTransform 输入：是否应用Piola变换（可选）
 * @param readyEdgeBasis 输入：预计算的边基函数（可选）
 * @param readyRotBasis 输入：预计算的旋度基函数（可选）
 * @param tangentialTrMapping 输入：是否应用切向迹映射（可选）
 * @return 是否成功计算
 */
bool EdgeElementInfo(const Element& element, const Nodes& nodes,
                    double u, double v, double w,
                    std::vector<std::vector<double>>* F,
                    std::vector<std::vector<double>>* G,
                    double& detF,
                    std::vector<double>& basis,
                    std::vector<std::vector<double>>& edgeBasis,
                    std::vector<std::vector<double>>* rotBasis,
                    std::vector<std::vector<double>>* dBasisdx,
                    bool* secondFamily,
                    int* basisDegree,
                    bool* applyPiolaTransform,
                    std::vector<std::vector<double>>* readyEdgeBasis,
                    std::vector<std::vector<double>>* readyRotBasis,
                    bool* tangentialTrMapping);

} // namespace elmer