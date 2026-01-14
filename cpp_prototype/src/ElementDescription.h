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
 * @brief 边界条件类型枚举
 * 
 * 定义边界条件的类型。
 */
enum class BoundaryConditionType {
    DIRICHLET,  ///< Dirichlet边界条件
    NEUMANN,    ///< Neumann边界条件
    ROBIN       ///< Robin边界条件
};

/**
 * @brief 边界条件结构
 * 
 * 描述边界条件的信息。
 */
struct BoundaryCondition {
    int boundaryIndex;          ///< 边界索引
    BoundaryConditionType type; ///< 边界条件类型
    double value;               ///< 边界值
    
    BoundaryCondition() : boundaryIndex(0), type(BoundaryConditionType::DIRICHLET), value(0.0) {}
    BoundaryCondition(int idx, BoundaryConditionType t, double v) : 
        boundaryIndex(idx), type(t), value(v) {}
};

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
    int numberOfBoundaries; ///< 边界数量
    std::vector<BasisFunction> basisFunctions; ///< 基函数数组
    
    ElementType() : numberOfNodes(0), dimension(0), elementCode(0), numberOfBoundaries(0) {}
};

/**
 * @brief 元素结构
 * 
 * 描述有限元元素的基本信息。
 */
struct Element {
    ElementType type; ///< 元素类型
    int index; ///< 元素索引
    std::vector<int> nodeIndexes; ///< 元素节点索引
    std::vector<int> edgeIndexes; ///< 元素边索引
    std::vector<int> faceIndexes; ///< 元素面索引
    
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
 * @brief 边结构
 * 
 * 描述有限元网格中的边信息。
 */
struct Edge {
    int index; ///< 边索引
    std::vector<int> nodeIndexes; ///< 边上的节点索引
    int bDofs; ///< 边上的自由度数量
    int pDegree; ///< P型元素的度数
    
    Edge() : index(0), bDofs(0), pDegree(0) {}
    Edge(int idx, const std::vector<int>& nodes) : 
        index(idx), nodeIndexes(nodes), bDofs(0), pDegree(0) {}
};

/**
 * @brief 面结构
 * 
 * 描述有限元网格中的面信息。
 */
struct Face {
    int index; ///< 面索引
    std::vector<int> nodeIndexes; ///< 面上的节点索引
    int bDofs; ///< 面上的自由度数量
    int pDegree; ///< P型元素的度数
    
    Face() : index(0), bDofs(0), pDegree(0) {}
    Face(int idx, const std::vector<int>& nodes) : 
        index(idx), nodeIndexes(nodes), bDofs(0), pDegree(0) {}
};

/**
 * @brief 网格结构
 * 
 * 描述有限元网格的基本信息。
 */
struct Mesh {
    int numberOfElements; ///< 元素数量
    int numberOfNodes; ///< 节点数量
    int numberOfEdges; ///< 边数量
    int numberOfFaces; ///< 面数量
    int maxEdgeDofs; ///< 最大边自由度数量
    int maxFaceDofs; ///< 最大面自由度数量
    int minEdgeDofs; ///< 最小边自由度数量
    int minFaceDofs; ///< 最小面自由度数量
    std::vector<Element> elements; ///< 元素数组
    std::vector<Node> nodes; ///< 节点对象数组
    std::vector<Edge> edges; ///< 边数组
    std::vector<Face> faces; ///< 面数组
    
    Mesh() : numberOfElements(0), numberOfNodes(0), numberOfEdges(0), numberOfFaces(0),
             maxEdgeDofs(0), maxFaceDofs(0), minEdgeDofs(0), minFaceDofs(0) {}
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
void GetRefPElementNodes(const ElementType& element, std::vector<double>& u_coords, 
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

} // namespace ElmerCpp