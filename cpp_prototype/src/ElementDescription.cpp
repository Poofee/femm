#include "ElementDescription.h"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace ElmerCpp {

/**
 * @brief 获取参考P型元素节点
 * 
 * 获取P型元素的参考节点坐标。
 * 简化实现：根据元素类型返回相应的参考节点坐标。
 */
void GetRefPElementNodes(const ElementType& element, std::vector<double>& u_coords, 
                        std::vector<double>& v_coords, std::vector<double>& w_coords) {
    int n = element.numberOfNodes;
    
    // 确保输出数组大小正确
    if (u_coords.size() < static_cast<size_t>(n)) u_coords.resize(n);
    if (v_coords.size() < static_cast<size_t>(n)) v_coords.resize(n);
    if (w_coords.size() < static_cast<size_t>(n)) w_coords.resize(n);
    
    int elemCode = element.elementCode;
    
    // 根据元素类型设置参考节点坐标
    switch (elemCode / 100) {
        case 2: // 线形元素
            for (int i = 0; i < n; ++i) {
                u_coords[i] = -1.0 + 2.0 * i / (n - 1); // 均匀分布在[-1,1]区间
                v_coords[i] = 0.0;
                w_coords[i] = 0.0;
            }
            break;
            
        case 3: // 三角形元素
            // 三角形元素参考节点（简化实现）
            if (n == 3) {
                u_coords[0] = 0.0; v_coords[0] = 0.0; w_coords[0] = 0.0;
                u_coords[1] = 1.0; v_coords[1] = 0.0; w_coords[1] = 0.0;
                u_coords[2] = 0.0; v_coords[2] = 1.0; w_coords[2] = 0.0;
            } else if (n == 6) {
                // 二次三角形元素
                u_coords[0] = 0.0; v_coords[0] = 0.0; w_coords[0] = 0.0;
                u_coords[1] = 1.0; v_coords[1] = 0.0; w_coords[1] = 0.0;
                u_coords[2] = 0.0; v_coords[2] = 1.0; w_coords[2] = 0.0;
                u_coords[3] = 0.5; v_coords[3] = 0.0; w_coords[3] = 0.0;
                u_coords[4] = 0.5; v_coords[4] = 0.5; w_coords[4] = 0.0;
                u_coords[5] = 0.0; v_coords[5] = 0.5; w_coords[5] = 0.0;
            }
            break;
            
        case 4: // 四边形元素
            // 四边形元素参考节点（简化实现）
            if (n == 4) {
                u_coords[0] = -1.0; v_coords[0] = -1.0; w_coords[0] = 0.0;
                u_coords[1] = 1.0; v_coords[1] = -1.0; w_coords[1] = 0.0;
                u_coords[2] = 1.0; v_coords[2] = 1.0; w_coords[2] = 0.0;
                u_coords[3] = -1.0; v_coords[3] = 1.0; w_coords[3] = 0.0;
            }
            break;
            
        case 5: // 四面体元素
            // 四面体元素参考节点（简化实现）
            if (n == 4) {
                u_coords[0] = 0.0; v_coords[0] = 0.0; w_coords[0] = 0.0;
                u_coords[1] = 1.0; v_coords[1] = 0.0; w_coords[1] = 0.0;
                u_coords[2] = 0.0; v_coords[2] = 1.0; w_coords[2] = 0.0;
                u_coords[3] = 0.0; v_coords[3] = 0.0; w_coords[3] = 1.0;
            }
            break;
            
        case 7: // 棱柱元素
            // 棱柱元素参考节点（简化实现）
            if (n == 6) {
                u_coords[0] = 0.0; v_coords[0] = 0.0; w_coords[0] = -1.0;
                u_coords[1] = 1.0; v_coords[1] = 0.0; w_coords[1] = -1.0;
                u_coords[2] = 0.0; v_coords[2] = 1.0; w_coords[2] = -1.0;
                u_coords[3] = 0.0; v_coords[3] = 0.0; w_coords[3] = 1.0;
                u_coords[4] = 1.0; v_coords[4] = 0.0; w_coords[4] = 1.0;
                u_coords[5] = 0.0; v_coords[5] = 1.0; w_coords[5] = 1.0;
            }
            break;
            
        default:
            // 默认情况：均匀分布
            for (int i = 0; i < n; ++i) {
                u_coords[i] = -1.0 + 2.0 * i / (n - 1);
                v_coords[i] = 0.0;
                w_coords[i] = 0.0;
            }
            break;
    }
    
    // 确保所有元素都被正确初始化
    for (int i = 0; i < n; ++i) {
        if (i >= static_cast<int>(u_coords.size())) u_coords.resize(i+1);
        if (i >= static_cast<int>(v_coords.size())) v_coords.resize(i+1);
        if (i >= static_cast<int>(w_coords.size())) w_coords.resize(i+1);
    }
}

/**
 * @brief 线单元P型基函数（单个节点）
 * 
 * 计算线单元在给定局部坐标点u处的单个P型基函数值。
 */
double LineNodalPBasis(int node, double u) {
    switch (node) {
        case 1:
            return (1.0 - u) / 2.0;
        case 2:
            return (1.0 + u) / 2.0;
        default:
            throw std::invalid_argument("LineNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 线单元P型基函数（所有节点）
 * 
 * 计算线单元在给定局部坐标点u处的所有P型基函数值。
 */
void LineNodalPBasisAll(double u, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 2) {
        phi.resize(2);
    }
    
    const double half = 0.5;
    const double usgn = 1.0; // 符号因子
    const double c = half * usgn;
    
    phi[0] = c * (1.0 - u);
    phi[1] = c * (1.0 + u);
}

/**
 * @brief 线单元P型基函数导数（单个节点）
 * 
 * 计算线单元在给定局部坐标点u处的单个P型基函数导数。
 */
double dLineNodalPBasis(int node, double u) {
    switch (node) {
        case 1:
            return -0.5;
        case 2:
            return 0.5;
        default:
            throw std::invalid_argument("dLineNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 线单元P型基函数导数（所有节点）
 * 
 * 计算线单元在给定局部坐标点u处的所有P型基函数导数。
 */
void dLineNodalPBasisAll(double u, std::vector<double>& dphi) {
    // 检查输出数组大小
    if (dphi.size() < 2) {
        dphi.resize(2);
    }
    
    const double half = 0.5;
    const double usgn = 1.0; // 符号因子
    const double c = half * usgn;
    
    dphi[0] = -c;
    dphi[1] = c;
}

/**
 * @brief 三角形单元P型基函数（单个节点）
 * 
 * 计算三角形单元在给定局部坐标点(u,v)处的单个P型基函数值。
 */
double TriangleNodalPBasis(int node, double u, double v) {
    switch (node) {
        case 1:
            return 1.0 - u - v;
        case 2:
            return u;
        case 3:
            return v;
        default:
            throw std::invalid_argument("TriangleNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 三角形单元P型基函数（所有节点）
 * 
 * 计算三角形单元在给定局部坐标点(u,v)处的所有P型基函数值。
 */
void TriangleNodalPBasisAll(double u, double v, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 3) {
        phi.resize(3);
    }
    
    phi[0] = 1.0 - u - v;
    phi[1] = u;
    phi[2] = v;
}

/**
 * @brief 三角形单元P型基函数导数（单个节点）
 * 
 * 计算三角形单元在给定局部坐标点(u,v)处的单个P型基函数导数。
 */
void dTriangleNodalPBasis(int node, double u, double v, double& dphi_du, double& dphi_dv) {
    switch (node) {
        case 1:
            dphi_du = -1.0;
            dphi_dv = -1.0;
            break;
        case 2:
            dphi_du = 1.0;
            dphi_dv = 0.0;
            break;
        case 3:
            dphi_du = 0.0;
            dphi_dv = 1.0;
            break;
        default:
            throw std::invalid_argument("dTriangleNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 三角形单元P型基函数导数（所有节点）
 * 
 * 计算三角形单元在给定局部坐标点(u,v)处的所有P型基函数导数。
 */
void dTriangleNodalPBasisAll(double u, double v, std::vector<double>& dphi_du, std::vector<double>& dphi_dv) {
    // 检查输出数组大小
    if (dphi_du.size() < 3) {
        dphi_du.resize(3);
    }
    if (dphi_dv.size() < 3) {
        dphi_dv.resize(3);
    }
    
    dphi_du[0] = -1.0;
    dphi_du[1] = 1.0;
    dphi_du[2] = 0.0;
    
    dphi_dv[0] = -1.0;
    dphi_dv[1] = 0.0;
    dphi_dv[2] = 1.0;
}

/**
 * @brief 四边形单元P型基函数（单个节点）
 * 
 * 计算四边形单元在给定局部坐标点(u,v)处的单个P型基函数值。
 */
double QuadNodalPBasis(int node, double u, double v) {
    switch (node) {
        case 1:
            return (1.0 - u) * (1.0 - v) / 4.0;
        case 2:
            return (1.0 + u) * (1.0 - v) / 4.0;
        case 3:
            return (1.0 + u) * (1.0 + v) / 4.0;
        case 4:
            return (1.0 - u) * (1.0 + v) / 4.0;
        default:
            throw std::invalid_argument("QuadNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 四边形单元P型基函数（所有节点）
 * 
 * 计算四边形单元在给定局部坐标点(u,v)处的所有P型基函数值。
 */
void QuadNodalPBasisAll(double u, double v, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 4) {
        phi.resize(4);
    }
    
    const double quarter = 0.25;
    
    phi[0] = quarter * (1.0 - u) * (1.0 - v);
    phi[1] = quarter * (1.0 + u) * (1.0 - v);
    phi[2] = quarter * (1.0 + u) * (1.0 + v);
    phi[3] = quarter * (1.0 - u) * (1.0 + v);
}

/**
 * @brief 四边形单元P型基函数导数（单个节点）
 * 
 * 计算四边形单元在给定局部坐标点(u,v)处的单个P型基函数导数。
 */
void dQuadNodalPBasis(int node, double u, double v, double& dphi_du, double& dphi_dv) {
    const double quarter = 0.25;
    
    switch (node) {
        case 1:
            dphi_du = -quarter * (1.0 - v);
            dphi_dv = -quarter * (1.0 - u);
            break;
        case 2:
            dphi_du = quarter * (1.0 - v);
            dphi_dv = -quarter * (1.0 + u);
            break;
        case 3:
            dphi_du = quarter * (1.0 + v);
            dphi_dv = quarter * (1.0 + u);
            break;
        case 4:
            dphi_du = -quarter * (1.0 + v);
            dphi_dv = quarter * (1.0 - u);
            break;
        default:
            throw std::invalid_argument("dQuadNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 四边形单元P型基函数导数（所有节点）
 * 
 * 计算四边形单元在给定局部坐标点(u,v)处的所有P型基函数导数。
 */
void dQuadNodalPBasisAll(double u, double v, std::vector<double>& dphi_du, std::vector<double>& dphi_dv) {
    // 检查输出数组大小
    if (dphi_du.size() < 4) {
        dphi_du.resize(4);
    }
    if (dphi_dv.size() < 4) {
        dphi_dv.resize(4);
    }
    
    const double quarter = 0.25;
    
    dphi_du[0] = -quarter * (1.0 - v);
    dphi_du[1] = quarter * (1.0 - v);
    dphi_du[2] = quarter * (1.0 + v);
    dphi_du[3] = -quarter * (1.0 + v);
    
    dphi_dv[0] = -quarter * (1.0 - u);
    dphi_dv[1] = -quarter * (1.0 + u);
    dphi_dv[2] = quarter * (1.0 + u);
    dphi_dv[3] = quarter * (1.0 - u);
}

/**
 * @brief 四面体单元P型基函数（单个节点）
 * 
 * 计算四面体单元在给定局部坐标点(u,v,w)处的单个P型基函数值。
 */
double TetraNodalPBasis(int node, double u, double v, double w) {
    switch (node) {
        case 1:
            return 1.0 - u - v - w;
        case 2:
            return u;
        case 3:
            return v;
        case 4:
            return w;
        default:
            throw std::invalid_argument("TetraNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 四面体单元P型基函数（所有节点）
 * 
 * 计算四面体单元在给定局部坐标点(u,v,w)处的所有P型基函数值。
 */
void TetraNodalPBasisAll(double u, double v, double w, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 4) {
        phi.resize(4);
    }
    
    phi[0] = 1.0 - u - v - w;
    phi[1] = u;
    phi[2] = v;
    phi[3] = w;
}

/**
 * @brief 四面体单元P型基函数导数（单个节点）
 * 
 * 计算四面体单元在给定局部坐标点(u,v,w)处的单个P型基函数导数。
 */
void dTetraNodalPBasis(int node, double u, double v, double w, double& dphi_du, double& dphi_dv, double& dphi_dw) {
    switch (node) {
        case 1:
            dphi_du = -1.0;
            dphi_dv = -1.0;
            dphi_dw = -1.0;
            break;
        case 2:
            dphi_du = 1.0;
            dphi_dv = 0.0;
            dphi_dw = 0.0;
            break;
        case 3:
            dphi_du = 0.0;
            dphi_dv = 1.0;
            dphi_dw = 0.0;
            break;
        case 4:
            dphi_du = 0.0;
            dphi_dv = 0.0;
            dphi_dw = 1.0;
            break;
        default:
            throw std::invalid_argument("dTetraNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 四面体单元P型基函数导数（所有节点）
 * 
 * 计算四面体单元在给定局部坐标点(u,v,w)处的所有P型基函数导数。
 */
void dTetraNodalPBasisAll(double u, double v, double w, std::vector<double>& dphi_du, std::vector<double>& dphi_dv, std::vector<double>& dphi_dw) {
    // 检查输出数组大小
    if (dphi_du.size() < 4) {
        dphi_du.resize(4);
    }
    if (dphi_dv.size() < 4) {
        dphi_dv.resize(4);
    }
    if (dphi_dw.size() < 4) {
        dphi_dw.resize(4);
    }
    
    dphi_du[0] = -1.0;
    dphi_du[1] = 1.0;
    dphi_du[2] = 0.0;
    dphi_du[3] = 0.0;
    
    dphi_dv[0] = -1.0;
    dphi_dv[1] = 0.0;
    dphi_dv[2] = 1.0;
    dphi_dv[3] = 0.0;
    
    dphi_dw[0] = -1.0;
    dphi_dw[1] = 0.0;
    dphi_dw[2] = 0.0;
    dphi_dw[3] = 1.0;
}

/**
 * @brief 棱柱单元P型基函数（单个节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的单个P型基函数值。
 */
double PrismNodalPBasis(int node, double u, double v, double w) {
    switch (node) {
        case 1:
            return (1.0 - u - v) * (1.0 - w) / 2.0;
        case 2:
            return u * (1.0 - w) / 2.0;
        case 3:
            return v * (1.0 - w) / 2.0;
        case 4:
            return (1.0 - u - v) * (1.0 + w) / 2.0;
        case 5:
            return u * (1.0 + w) / 2.0;
        case 6:
            return v * (1.0 + w) / 2.0;
        default:
            throw std::invalid_argument("PrismNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 棱柱单元P型基函数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有P型基函数值。
 */
void PrismNodalPBasisAll(double u, double v, double w, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 6) {
        phi.resize(6);
    }
    
    const double half = 0.5;
    
    phi[0] = half * (1.0 - u - v) * (1.0 - w);
    phi[1] = half * u * (1.0 - w);
    phi[2] = half * v * (1.0 - w);
    phi[3] = half * (1.0 - u - v) * (1.0 + w);
    phi[4] = half * u * (1.0 + w);
    phi[5] = half * v * (1.0 + w);
}

/**
 * @brief 棱柱单元P型基函数导数（单个节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的单个P型基函数导数。
 */
void dPrismNodalPBasis(int node, double u, double v, double w, double& dphi_du, double& dphi_dv, double& dphi_dw) {
    const double half = 0.5;
    
    switch (node) {
        case 1:
            dphi_du = -half * (1.0 - w);
            dphi_dv = -half * (1.0 - w);
            dphi_dw = -half * (1.0 - u - v);
            break;
        case 2:
            dphi_du = half * (1.0 - w);
            dphi_dv = 0.0;
            dphi_dw = -half * u;
            break;
        case 3:
            dphi_du = 0.0;
            dphi_dv = half * (1.0 - w);
            dphi_dw = -half * v;
            break;
        case 4:
            dphi_du = -half * (1.0 + w);
            dphi_dv = -half * (1.0 + w);
            dphi_dw = half * (1.0 - u - v);
            break;
        case 5:
            dphi_du = half * (1.0 + w);
            dphi_dv = 0.0;
            dphi_dw = half * u;
            break;
        case 6:
            dphi_du = 0.0;
            dphi_dv = half * (1.0 + w);
            dphi_dw = half * v;
            break;
        default:
            throw std::invalid_argument("dPrismNodalPBasis: 无效的节点编号");
    }
}

/**
 * @brief 棱柱单元P型基函数导数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有P型基函数导数。
 */
void dPrismNodalPBasisAll(double u, double v, double w, std::vector<double>& dphi_du, std::vector<double>& dphi_dv, std::vector<double>& dphi_dw) {
    // 检查输出数组大小
    if (dphi_du.size() < 6) {
        dphi_du.resize(6);
    }
    if (dphi_dv.size() < 6) {
        dphi_dv.resize(6);
    }
    if (dphi_dw.size() < 6) {
        dphi_dw.resize(6);
    }
    
    const double half = 0.5;
    
    dphi_du[0] = -half * (1.0 - w);
    dphi_du[1] = half * (1.0 - w);
    dphi_du[2] = 0.0;
    dphi_du[3] = -half * (1.0 + w);
    dphi_du[4] = half * (1.0 + w);
    dphi_du[5] = 0.0;
    
    dphi_dv[0] = -half * (1.0 - w);
    dphi_dv[1] = 0.0;
    dphi_dv[2] = half * (1.0 - w);
    dphi_dv[3] = -half * (1.0 + w);
    dphi_dv[4] = 0.0;
    dphi_dv[5] = half * (1.0 + w);
    
    dphi_dw[0] = -half * (1.0 - u - v);
    dphi_dw[1] = -half * u;
    dphi_dw[2] = -half * v;
    dphi_dw[3] = half * (1.0 - u - v);
    dphi_dw[4] = half * u;
    dphi_dw[5] = half * v;
}

/**
 * @brief 获取元素边界连续性信息
 * 
 * 获取元素边界连续性信息，包括边界类型和连续性条件。
 */
void GetElementBoundaryContinuity(const ElementType& element, std::vector<BoundaryContinuity>& continuity) {
    int nBoundaries = element.numberOfBoundaries;
    
    // 确保输出数组大小正确
    if (continuity.size() < static_cast<size_t>(nBoundaries)) {
        continuity.resize(nBoundaries);
    }
    
    // 根据元素类型设置边界连续性
    int elemCode = element.elementCode;
    
    switch (elemCode / 100) {
        case 2: // 线形元素
            if (nBoundaries >= 2) {
                continuity[0].type = BoundaryType::POINT;
                continuity[0].continuity = ContinuityType::C0;
                continuity[1].type = BoundaryType::POINT;
                continuity[1].continuity = ContinuityType::C0;
            }
            break;
            
        case 3: // 三角形元素
            if (nBoundaries >= 3) {
                for (int i = 0; i < 3; ++i) {
                    continuity[i].type = BoundaryType::LINE;
                    continuity[i].continuity = ContinuityType::C0;
                }
            }
            break;
            
        case 4: // 四边形元素
            if (nBoundaries >= 4) {
                for (int i = 0; i < 4; ++i) {
                    continuity[i].type = BoundaryType::LINE;
                    continuity[i].continuity = ContinuityType::C0;
                }
            }
            break;
            
        case 5: // 四面体元素
            if (nBoundaries >= 4) {
                for (int i = 0; i < 4; ++i) {
                    continuity[i].type = BoundaryType::TRIANGLE;
                    continuity[i].continuity = ContinuityType::C0;
                }
            }
            break;
            
        case 7: // 棱柱元素
            if (nBoundaries >= 5) {
                // 三角形面
                for (int i = 0; i < 2; ++i) {
                    continuity[i].type = BoundaryType::TRIANGLE;
                    continuity[i].continuity = ContinuityType::C0;
                }
                // 四边形面
                for (int i = 2; i < 5; ++i) {
                    continuity[i].type = BoundaryType::QUADRILATERAL;
                    continuity[i].continuity = ContinuityType::C0;
                }
            }
            break;
            
        default:
            // 默认情况：所有边界都是C0连续
            for (int i = 0; i < nBoundaries; ++i) {
                continuity[i].type = BoundaryType::UNKNOWN;
                continuity[i].continuity = ContinuityType::C0;
            }
            break;
    }
}

/**
 * @brief 检查边界条件是否连续
 * 
 * 检查两个边界条件是否连续。
 */
bool CheckBoundaryContinuity(const BoundaryContinuity& bc1, const BoundaryContinuity& bc2) {
    // 边界类型必须匹配
    if (bc1.type != bc2.type) {
        return false;
    }
    
    // 连续性类型必须匹配
    if (bc1.continuity != bc2.continuity) {
        return false;
    }
    
    return true;
}

/**
 * @brief 处理边界条件
 * 
 * 处理边界条件，包括设置边界值和边界类型。
 */
void ProcessBoundaryConditions(const ElementType& element, const std::vector<double>& boundaryValues, 
                              std::vector<BoundaryCondition>& boundaryConditions) {
    int nBoundaries = element.numberOfBoundaries;
    
    // 确保输出数组大小正确
    if (boundaryConditions.size() < static_cast<size_t>(nBoundaries)) {
        boundaryConditions.resize(nBoundaries);
    }
    
    // 设置边界条件
    for (int i = 0; i < nBoundaries; ++i) {
        boundaryConditions[i].boundaryIndex = i;
        boundaryConditions[i].type = BoundaryConditionType::DIRICHLET;
        
        // 设置边界值
        if (i < static_cast<int>(boundaryValues.size())) {
            boundaryConditions[i].value = boundaryValues[i];
        } else {
            boundaryConditions[i].value = 0.0;
        }
    }
}

/**
 * @brief 获取边界节点
 * 
 * 获取元素边界上的节点索引。
 */
void GetBoundaryNodes(const ElementType& element, int boundaryIndex, std::vector<int>& boundaryNodes) {
    int nBoundaries = element.numberOfBoundaries;
    
    if (boundaryIndex < 0 || boundaryIndex >= nBoundaries) {
        throw std::invalid_argument("GetBoundaryNodes: 无效的边界索引");
    }
    
    int elemCode = element.elementCode;
    
    // 根据元素类型和边界索引设置边界节点
    switch (elemCode / 100) {
        case 2: // 线形元素
            if (boundaryNodes.size() < 1) {
                boundaryNodes.resize(1);
            }
            boundaryNodes[0] = boundaryIndex; // 线元素的边界是端点
            break;
            
        case 3: // 三角形元素
            if (boundaryNodes.size() < 2) {
                boundaryNodes.resize(2);
            }
            // 三角形边界是边，包含2个节点
            switch (boundaryIndex) {
                case 0:
                    boundaryNodes[0] = 0;
                    boundaryNodes[1] = 1;
                    break;
                case 1:
                    boundaryNodes[0] = 1;
                    boundaryNodes[1] = 2;
                    break;
                case 2:
                    boundaryNodes[0] = 2;
                    boundaryNodes[1] = 0;
                    break;
                default:
                    throw std::invalid_argument("GetBoundaryNodes: 无效的三角形边界索引");
            }
            break;
            
        case 4: // 四边形元素
            if (boundaryNodes.size() < 2) {
                boundaryNodes.resize(2);
            }
            // 四边形边界是边，包含2个节点
            switch (boundaryIndex) {
                case 0:
                    boundaryNodes[0] = 0;
                    boundaryNodes[1] = 1;
                    break;
                case 1:
                    boundaryNodes[0] = 1;
                    boundaryNodes[1] = 2;
                    break;
                case 2:
                    boundaryNodes[0] = 2;
                    boundaryNodes[1] = 3;
                    break;
                case 3:
                    boundaryNodes[0] = 3;
                    boundaryNodes[1] = 0;
                    break;
                default:
                    throw std::invalid_argument("GetBoundaryNodes: 无效的四边形边界索引");
            }
            break;
            
        default:
            throw std::invalid_argument("GetBoundaryNodes: 不支持的元素类型");
    }
}

/**
 * @brief 获取三角形面方向
 * 
 * 计算三角形面的方向向量。
 */
std::vector<int> getTriangleFaceDirection(const Element& element, 
                                         const std::vector<int>& faceMap,
                                         const std::vector<int>& indexes) {
    std::vector<int> direction(3);
    
    // 简单的三角形面方向计算
    for (int i = 0; i < 3; ++i) {
        direction[i] = i;
    }
    
    return direction;
}

/**
 * @brief 获取四边形面方向
 * 
 * 计算四边形面的方向向量。
 */
std::vector<int> getSquareFaceDirection(const Element& element,
                                       const std::vector<int>& faceMap,
                                       const std::vector<int>& indexes) {
    std::vector<int> direction(4);
    
    // 简单的四边形面方向计算
    for (int i = 0; i < 4; ++i) {
        direction[i] = i;
    }
    
    return direction;
}

/**
 * @brief 检查楔形元素面排序
 * 
 * 验证楔形元素面的排序是否有效。
 */
bool wedgeOrdering(const std::vector<int>& ordering) {
    // 检查是否有重复索引
    std::vector<int> seen(ordering.size(), 0);
    for (int idx : ordering) {
        if (idx < 0 || idx >= static_cast<int>(ordering.size())) {
            return false;
        }
        if (seen[idx] > 0) {
            return false; // 重复索引
        }
        seen[idx]++;
    }
    return true;
}

/**
 * @brief 检查面节点连续性
 * 
 * 检查两个元素在共享面上的节点连续性。
 */
bool CheckFaceNodeContinuity(const Element& element1, const Element& element2,
                            int faceIndex1, const Mesh& mesh, double tolerance) {
    // 简单的节点连续性检查
    // 假设相同节点索引表示连续性
    if (element1.nodeIndexes.size() != element2.nodeIndexes.size()) {
        return false;
    }
    
    for (size_t i = 0; i < element1.nodeIndexes.size(); ++i) {
        if (element1.nodeIndexes[i] != element2.nodeIndexes[i]) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief 检查面几何连续性
 * 
 * 检查两个元素在共享面上的几何连续性。
 */
bool CheckFaceGeometricContinuity(const Element& element1, const Element& element2,
                                 int faceIndex1, const Mesh& mesh, const Nodes& nodes,
                                 double tolerance) {
    // 简单的几何连续性检查
    // 使用节点坐标检查几何连续性
    if (element1.nodeIndexes.size() != element2.nodeIndexes.size()) {
        return false;
    }
    
    for (size_t i = 0; i < element1.nodeIndexes.size(); ++i) {
        int node1 = element1.nodeIndexes[i];
        int node2 = element2.nodeIndexes[i];
        
        // 检查节点坐标是否在容差范围内匹配
        if (std::abs(nodes.x[node1] - nodes.x[node2]) > tolerance ||
            std::abs(nodes.y[node1] - nodes.y[node2]) > tolerance ||
            std::abs(nodes.z[node1] - nodes.z[node2]) > tolerance) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief 检查面函数连续性
 * 
 * 检查两个元素在共享面上的函数连续性。
 */
bool CheckFaceFunctionContinuity(const Element& element1, const Element& element2,
                                int faceIndex1, const Mesh& mesh, const Nodes& nodes,
                                const std::vector<double>& field1, const std::vector<double>& field2,
                                double tolerance) {
    // 简单的函数连续性检查
    if (field1.size() != field2.size()) {
        return false;
    }
    
    for (size_t i = 0; i < field1.size(); ++i) {
        if (std::abs(field1[i] - field2[i]) > tolerance) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief 处理边界条件
 * 
 * 处理边界条件，包括设置边界值和边界类型。
 */
void ProcessBoundaryConditions(const Element& element, int boundaryType, double boundaryValue,
                              const std::vector<int>& boundaryNodes, std::vector<double>& boundaryValues) {
    // 确保边界值数组大小正确
    if (boundaryValues.size() < boundaryNodes.size()) {
        boundaryValues.resize(boundaryNodes.size());
    }
    
    // 设置边界值
    for (size_t i = 0; i < boundaryNodes.size(); ++i) {
        boundaryValues[i] = boundaryValue;
    }
}

/**
 * @brief 获取面法向量
 * 
 * 计算元素面的法向量。
 */
void GetFaceNormal(const Element& element, int faceIndex, const Nodes& nodes, std::vector<double>& normal) {
    // 简单的面法向量计算（假设为平面）
    if (normal.size() < 3) {
        normal.resize(3);
    }
    
    // 对于三角形面，计算法向量
    if (element.nodeIndexes.size() >= 3) {
        int node1 = element.nodeIndexes[0];
        int node2 = element.nodeIndexes[1];
        int node3 = element.nodeIndexes[2];
        
        double v1x = nodes.x[node2] - nodes.x[node1];
        double v1y = nodes.y[node2] - nodes.y[node1];
        double v1z = nodes.z[node2] - nodes.z[node1];
        
        double v2x = nodes.x[node3] - nodes.x[node1];
        double v2y = nodes.y[node3] - nodes.y[node1];
        double v2z = nodes.z[node3] - nodes.z[node1];
        
        // 叉积
        normal[0] = v1y * v2z - v1z * v2y;
        normal[1] = v1z * v2x - v1x * v2z;
        normal[2] = v1x * v2y - v1y * v2x;
        
        // 归一化
        double length = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
        if (length > 1e-12) {
            normal[0] /= length;
            normal[1] /= length;
            normal[2] /= length;
        }
    } else {
        // 默认法向量
        normal[0] = 0.0;
        normal[1] = 0.0;
        normal[2] = 1.0;
    }
}

/**
 * @brief 计算面通量
 * 
 * 计算通过元素面的通量。
 */
void ComputeFaceFlux(const Element& element, int faceIndex, const Nodes& nodes,
                    const std::vector<double>& field, double& flux) {
    // 简单的面通量计算（平均值）
    double sum = 0.0;
    for (double value : field) {
        sum += value;
    }
    flux = sum / field.size();
}

/**
 * @brief 应用Dirichlet边界条件
 * 
 * 在刚度矩阵和力向量中应用Dirichlet边界条件。
 */
void ApplyDirichletBoundaryConditions(const Element& element, const std::vector<int>& boundaryNodes,
                                     const std::vector<double>& dirichletValues,
                                     std::vector<std::vector<double>>& stiffnessMatrix,
                                     std::vector<double>& forceVector) {
    // 应用Dirichlet边界条件
    for (size_t i = 0; i < boundaryNodes.size(); ++i) {
        int nodeIdx = boundaryNodes[i];
        
        // 设置刚度矩阵对角线为1，其他为0
        for (size_t j = 0; j < stiffnessMatrix.size(); ++j) {
            if (j == static_cast<size_t>(nodeIdx)) {
                stiffnessMatrix[nodeIdx][j] = 1.0;
            } else {
                stiffnessMatrix[nodeIdx][j] = 0.0;
                stiffnessMatrix[j][nodeIdx] = 0.0;
            }
        }
        
        // 设置力向量为边界值
        if (i < dirichletValues.size()) {
            forceVector[nodeIdx] = dirichletValues[i];
        }
    }
}

/**
 * @brief 应用Neumann边界条件
 * 
 * 在力向量中应用Neumann边界条件。
 */
void ApplyNeumannBoundaryConditions(const Element& element, int boundaryIndex, double neumannValue,
                                   std::vector<double>& forceVector) {
    // 应用Neumann边界条件（简单地在力向量上添加值）
    for (size_t i = 0; i < forceVector.size(); ++i) {
        forceVector[i] += neumannValue;
    }
}

} // namespace ElmerCpp