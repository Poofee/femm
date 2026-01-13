/**
 * @file ElementDescription.cpp
 * @brief 元素描述模块实现
 * 
 * 对应Fortran模块：ElementDescription.F90
 * 提供有限元元素的类型定义、基函数描述等基础数据结构的实现。
 */

#include "ElementDescription.h"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace ElmerCpp {

/**
 * @brief 计算全局一阶导数（内部版本）
 * 
 * 给定元素结构，返回在局部坐标点(u,v,w)处，
 * 基于元素基函数计算的量f在全局坐标下的偏导数值。
 * 这是内部版本，通常不应直接由用户调用，而应通过包装例程GlobalFirstDerivatives调用。
 */
void GlobalFirstDerivativesInternal(const Element& elm, const Nodes& nodes, 
                                   const std::vector<double>& df,
                                   double& gx, double& gy, double& gz,
                                   const std::vector<std::vector<double>>& metric,
                                   const std::vector<std::vector<double>>& dLBasisdx) {
    
    int n = elm.type.numberOfNodes;
    int dim = elm.type.dimension;
    int cdim = CoordinateSystemDimension();

    // 检查输入参数的有效性
    if (n <= 0 || dim <= 0 || cdim <= 0) {
        throw std::invalid_argument("GlobalFirstDerivativesInternal: 无效的元素参数");
    }

    if (nodes.x.size() < static_cast<size_t>(n) || 
        nodes.y.size() < static_cast<size_t>(n) || 
        nodes.z.size() < static_cast<size_t>(n)) {
        throw std::invalid_argument("GlobalFirstDerivativesInternal: 节点坐标数组大小不足");
    }

    if (df.size() < static_cast<size_t>(dim)) {
        throw std::invalid_argument("GlobalFirstDerivativesInternal: 节点值数组大小不足");
    }

    // 全局坐标对局部坐标的偏导数矩阵
    std::vector<std::vector<double>> dx(3, std::vector<double>(3, 0.0));
    
    // 计算全局坐标对局部坐标的偏导数
    switch (cdim) {
        case 1:
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < n; ++j) {
                    dx[0][i] += nodes.x[j] * dLBasisdx[j][i];
                }
            }
            break;
            
        case 2:
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < n; ++j) {
                    dx[0][i] += nodes.x[j] * dLBasisdx[j][i];
                    dx[1][i] += nodes.y[j] * dLBasisdx[j][i];
                }
            }
            break;
            
        case 3:
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < n; ++j) {
                    dx[0][i] += nodes.x[j] * dLBasisdx[j][i];
                    dx[1][i] += nodes.y[j] * dLBasisdx[j][i];
                    dx[2][i] += nodes.z[j] * dLBasisdx[j][i];
                }
            }
            break;
            
        default:
            throw std::invalid_argument("GlobalFirstDerivativesInternal: 不支持的坐标系维度");
    }

    // 逆变分量：在元素坐标系下的偏导数
    std::vector<double> dfc(dim, 0.0);
    for (int i = 0; i < dim; ++i) {
        double s = 0.0;
        for (int j = 0; j < dim; ++j) {
            s += metric[i][j] * df[j];
        }
        dfc[i] = s;
    }

    // 变换到空间坐标
    gx = 0.0;
    gy = 0.0;
    gz = 0.0;
    
    switch (cdim) {
        case 1:
            for (int i = 0; i < dim; ++i) {
                gx += dx[0][i] * dfc[i];
            }
            break;
            
        case 2:
            for (int i = 0; i < dim; ++i) {
                gx += dx[0][i] * dfc[i];
                gy += dx[1][i] * dfc[i];
            }
            break;
            
        case 3:
            for (int i = 0; i < dim; ++i) {
                gx += dx[0][i] * dfc[i];
                gy += dx[1][i] * dfc[i];
                gz += dx[2][i] * dfc[i];
            }
            break;
    }
}

/**
 * @brief 计算全局一阶导数
 * 
 * 给定元素结构，返回在局部坐标点(u,v,w)处，
 * 基于元素基函数计算的量f在全局坐标下的偏导数值。
 */
void GlobalFirstDerivatives(const Element& elm, const Nodes& nodes,
                           const std::vector<double>& df,
                           double& gx, double& gy, double& gz,
                           const std::vector<std::vector<double>>& metric,
                           const std::vector<std::vector<double>>& dLBasisdx) {
    
    // 调用内部版本实现
    GlobalFirstDerivativesInternal(elm, nodes, df, gx, gy, gz, metric, dLBasisdx);
}

/**
 * @brief 计算全局二阶导数
 * 
 * 在给定点(u,v,w)处计算全局坐标下的元素二阶偏导数矩阵。
 */
void GlobalSecondDerivatives(const Element& elm, const Nodes& nodes,
                            std::vector<std::vector<std::vector<double>>>& values,
                            double u, double v, double w,
                            const std::vector<std::vector<double>>& metric,
                            const std::vector<std::vector<double>>& dBasisdx,
                            const std::vector<std::vector<std::vector<double>>>& ddLBasisddx,
                            int nd) {
    
    int n = elm.type.numberOfNodes;
    int dim = elm.type.dimension;
    int cdim = CoordinateSystemDimension();

    // 检查输入参数的有效性
    if (n <= 0 || dim <= 0 || cdim <= 0) {
        throw std::invalid_argument("GlobalSecondDerivatives: 无效的元素参数");
    }

    if (nodes.x.size() < static_cast<size_t>(n) || 
        nodes.y.size() < static_cast<size_t>(n) || 
        nodes.z.size() < static_cast<size_t>(n)) {
        throw std::invalid_argument("GlobalSecondDerivatives: 节点坐标数组大小不足");
    }

    // 初始化输出矩阵
    values.resize(3);
    for (auto& mat2d : values) {
        mat2d.resize(3);
        for (auto& vec : mat2d) {
            vec.resize(3, 0.0);
        }
    }

    // 全局坐标对局部坐标的一阶偏导数矩阵
    std::vector<std::vector<double>> dx(3, std::vector<double>(3, 0.0));
    
    // 全局坐标对局部坐标的二阶偏导数张量
    std::vector<std::vector<std::vector<double>>> ddx(3, 
        std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0)));

    // 计算一阶偏导数
    switch (cdim) {
        case 1:
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < nd; ++j) {
                    dx[0][i] += nodes.x[j] * dBasisdx[j][i];
                }
            }
            break;
            
        case 2:
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < nd; ++j) {
                    dx[0][i] += nodes.x[j] * dBasisdx[j][i];
                    dx[1][i] += nodes.y[j] * dBasisdx[j][i];
                }
            }
            break;
            
        case 3:
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < nd; ++j) {
                    dx[0][i] += nodes.x[j] * dBasisdx[j][i];
                    dx[1][i] += nodes.y[j] * dBasisdx[j][i];
                    dx[2][i] += nodes.z[j] * dBasisdx[j][i];
                }
            }
            break;
    }

    // 计算二阶偏导数
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < nd; ++k) {
                ddx[0][i][j] += ddLBasisddx[k][i][j] * nodes.x[k];
                ddx[1][i][j] += ddLBasisddx[k][i][j] * nodes.y[k];
                ddx[2][i][j] += ddLBasisddx[k][i][j] * nodes.z[k];
            }
        }
    }

    // Christoffel符号计算（元素坐标系下的第二类Christoffel符号）
    std::vector<std::vector<std::vector<double>>> C1(3, 
        std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0)));
    std::vector<std::vector<std::vector<double>>> C2(3, 
        std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0)));

    // 简化实现：对于线性元素，二阶导数为零
    if (elm.type.elementCode <= 202 || elm.type.elementCode == 303 || 
        elm.type.elementCode == 504) {
        // 线性元素，二阶导数为零，直接返回
        return;
    }

    // 计算Christoffel符号
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
                for (int l = 0; l < dim; ++l) {
                    C1[i][j][k] += metric[k][l] * ddx[l][i][j];
                }
            }
        }
    }

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
                for (int l = 0; l < dim; ++l) {
                    C2[i][j][k] += dx[k][l] * C1[i][j][l];
                }
            }
        }
    }

    // 计算二阶偏导数矩阵
    for (int q = 0; q < 3; ++q) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double s = 0.0;
                for (int k = 0; k < dim; ++k) {
                    s += C2[i][j][k] * dx[q][k];
                }
                values[q][i][j] = s;
            }
        }
    }
}

/**
 * @brief 获取坐标系维度
 * 
 * 返回当前坐标系的维度（1D、2D或3D）。
 * 这是一个简化实现，实际应用中可能需要根据具体配置返回正确的维度。
 * 
 * @return 坐标系维度（1、2或3）
 */
int CoordinateSystemDimension() {
    // 简化实现：默认返回3D
    // 在实际应用中，可能需要从配置文件或全局设置中获取
    return 3;
}

/**
 * @brief 棱柱单元P型基函数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有P型基函数值。
 * 棱柱单元有6个节点，基函数是三角形基函数和线形基函数的张量积。
 */
void WedgeNodalPBasisAll(double u, double v, double w, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 6) {
        phi.resize(6);
    }
    
    const double half = 0.5;
    const double c3 = 1.0 / std::sqrt(3.0);
    
    // 棱柱单元基函数：三角形基函数 × 线形基函数
    phi[0] = half * (1.0 - u - c3 * v) * (1.0 - w);
    phi[1] = half * (1.0 + u - c3 * v) * (1.0 - w);
    phi[2] = c3 * v * (1.0 - w);
    phi[3] = half * (1.0 - u - c3 * v) * (1.0 + w);
    phi[4] = half * (1.0 + u - c3 * v) * (1.0 + w);
    phi[5] = c3 * v * (1.0 + w);
}
    
    // 三角形基函数（P型）
    std::vector<double> tri(3);
    tri[0] = half * (1.0 - u - c3 * v);
    tri[1] = half * (1.0 + u - c3 * v);
    tri[2] = c3 * v;
    
    // 线形基函数
    std::vector<double> line(2);
    line[0] = half * (1.0 - w);
    line[1] = half * (1.0 + w);
    
    // 张量积：棱柱单元基函数 = 三角形基函数 × 线形基函数
    phi[0] = line[0] * tri[0];
    phi[1] = line[0] * tri[1];
    phi[2] = line[0] * tri[2];
    phi[3] = line[1] * tri[0];
    phi[4] = line[1] * tri[1];
    phi[5] = line[1] * tri[2];
}

/**
 * @brief 棱柱单元L型基函数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有L型基函数值。
 * L型基函数使用线性插值，适用于线性元素。
 */
void WedgeNodalLBasisAll(double u, double v, double w, std::vector<double>& phi) {
    // 检查输出数组大小
    if (phi.size() < 6) {
        phi.resize(6);
    }
    
    const double half = 0.5;
    
    // 三角形基函数（L型）
    std::vector<double> tri(3);
    tri[0] = 1.0 - u - v;
    tri[1] = u;
    tri[2] = v;
    
    // 线形基函数
    std::vector<double> line(2);
    line[0] = half * (1.0 - w);
    line[1] = half * (1.0 + w);
    
    // 张量积：棱柱单元基函数 = 三角形基函数 × 线形基函数
    phi[0] = line[0] * tri[0];
    phi[1] = line[0] * tri[1];
    phi[2] = line[0] * tri[2];
    phi[3] = line[1] * tri[0];
    phi[4] = line[1] * tri[1];
    phi[5] = line[1] * tri[2];
}

/**
 * @brief 棱柱单元P型基函数导数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有P型基函数导数。
 * 返回6x3矩阵，其中6个节点，3个方向导数（du, dv, dw）。
 */
void dWedgeNodalPBasisAll(double u, double v, double w, std::vector<std::vector<double>>& gradphi) {
    // 检查输出矩阵大小
    if (gradphi.size() < 6 || gradphi[0].size() < 3) {
        gradphi.resize(6, std::vector<double>(3, 0.0));
    }
    
    const double half = 0.5;
    const double c3 = 1.0 / std::sqrt(3.0);
    
    // 三角形基函数值
    std::vector<double> tri(3);
    tri[0] = half * (1.0 - u - c3 * v);
    tri[1] = half * (1.0 + u - c3 * v);
    tri[2] = c3 * v;
    
    // 线形基函数值
    std::vector<double> line(2);
    line[0] = half * (1.0 - w);
    line[1] = half * (1.0 + w);
    
    // 三角形基函数导数
    std::vector<std::vector<double>> gradtri(3, std::vector<double>(2, 0.0));
    gradtri[0][0] = -half;      // d(phi_tri1)/du
    gradtri[0][1] = -half * c3; // d(phi_tri1)/dv
    gradtri[1][0] = half;       // d(phi_tri2)/du
    gradtri[1][1] = -half * c3; // d(phi_tri2)/dv
    gradtri[2][0] = 0.0;        // d(phi_tri3)/du
    gradtri[2][1] = c3;         // d(phi_tri3)/dv
    
    // 线形基函数导数
    std::vector<double> gradline(2);
    gradline[0] = -half; // d(phi_line1)/dw
    gradline[1] = half;  // d(phi_line2)/dw
    
    // 计算棱柱单元基函数导数
    // 前三个节点（w = -1）
    for (int i = 0; i < 3; ++i) {
        // du导数：d(phi)/du = d(phi_tri)/du * phi_line
        gradphi[i][0] = gradtri[i][0] * line[0];
        // dv导数：d(phi)/dv = d(phi_tri)/dv * phi_line
        gradphi[i][1] = gradtri[i][1] * line[0];
        // dw导数：d(phi)/dw = phi_tri * d(phi_line)/dw
        gradphi[i][2] = tri[i] * gradline[0];
    }
    
    // 后三个节点（w = +1）
    for (int i = 0; i < 3; ++i) {
        // du导数
        gradphi[i+3][0] = gradtri[i][0] * line[1];
        // dv导数
        gradphi[i+3][1] = gradtri[i][1] * line[1];
        // dw导数
        gradphi[i+3][2] = tri[i] * gradline[1];
    }
}

/**
 * @brief 棱柱单元L型基函数导数（所有节点）
 * 
 * 计算棱柱单元在给定局部坐标点(u,v,w)处的所有L型基函数导数。
 * 返回6x3矩阵，其中6个节点，3个方向导数（du, dv, dw）。
 */
void dWedgeNodalLBasisAll(double u, double v, double w, std::vector<std::vector<double>>& gradphi) {
    // 检查输出矩阵大小
    if (gradphi.size() < 6 || gradphi[0].size() < 3) {
        gradphi.resize(6, std::vector<double>(3, 0.0));
    }
    
    const double half = 0.5;
    
    // 三角形基函数值
    std::vector<double> tri(3);
    tri[0] = 1.0 - u - v;
    tri[1] = u;
    tri[2] = v;
    
    // 线形基函数值
    std::vector<double> line(2);
    line[0] = half * (1.0 - w);
    line[1] = half * (1.0 + w);
    
    // 三角形基函数导数
    std::vector<std::vector<double>> gradtri(3, std::vector<double>(2, 0.0));
    gradtri[0][0] = -1.0; // d(phi_tri1)/du
    gradtri[0][1] = -1.0; // d(phi_tri1)/dv
    gradtri[1][0] = 1.0;  // d(phi_tri2)/du
    gradtri[1][1] = 0.0;  // d(phi_tri2)/dv
    gradtri[2][0] = 0.0;  // d(phi_tri3)/du
    gradtri[2][1] = 1.0;  // d(phi_tri3)/dv
    
    // 线形基函数导数
    std::vector<double> gradline(2);
    gradline[0] = -half; // d(phi_line1)/dw
    gradline[1] = half;  // d(phi_line2)/dw
    
    // 计算棱柱单元基函数导数
    // 前三个节点（w = -1）
    for (int i = 0; i < 3; ++i) {
        // du导数
        gradphi[i][0] = gradtri[i][0] * line[0];
        // dv导数
        gradphi[i][1] = gradtri[i][1] * line[0];
        // dw导数
        gradphi[i][2] = tri[i] * gradline[0];
    }
    
    // 后三个节点（w = +1）
    for (int i = 0; i < 3; ++i) {
        // du导数
        gradphi[i+3][0] = gradtri[i][0] * line[1];
        // dv导数
        gradphi[i+3][1] = gradtri[i][1] * line[1];
        // dw导数
    gradphi[i+3][2] = tri[i] * gradline[1];
  }
}

// =============================================================================
// P型元素支持函数实现
// =============================================================================

/**
 * @brief 检查元素是否为P型元素
 * 
 * 检查给定元素是否为P型（高阶）元素。
 * 简化实现：根据元素代码判断是否为P型元素。
 */
bool isPElement(const Element& element) {
    // 简化实现：根据元素代码判断
    // 在实际应用中，可能需要检查更复杂的条件
    int elemCode = element.type.elementCode;
    
    // P型元素通常具有特定的元素代码范围
    // 例如：2xx系列为线形P型元素，3xx系列为三角形P型元素等
    if (elemCode >= 200 && elemCode < 300) {
        return true; // 线形P型元素
    } else if (elemCode >= 300 && elemCode < 400) {
        return true; // 三角形P型元素
    } else if (elemCode >= 400 && elemCode < 500) {
        return true; // 四边形P型元素
    } else if (elemCode >= 500 && elemCode < 600) {
        return true; // 四面体P型元素
    } else if (elemCode >= 700 && elemCode < 800) {
        return true; // 棱柱P型元素
    }
    
    return false;
}

/**
 * @brief 检查元素是否为活动P型元素
 * 
 * 检查给定元素在特定求解器中是否为活动的P型元素。
 * 简化实现：暂时忽略求解器参数，仅检查是否为P型元素。
 */
bool isActivePElement(const Element& element, void* solver) {
    // 简化实现：暂时忽略求解器参数
    // 在实际应用中，可能需要检查求解器中的P型元素定义
    return isPElement(element);
}

/**
 * @brief 获取参考P型元素节点
 * 
 * 获取P型元素的参考节点坐标。
 * 简化实现：根据元素类型返回相应的参考节点坐标。
 */
void GetRefPElementNodes(const ElementType& element, std::vector<double>& u, 
                        std::vector<double>& v, std::vector<double>& w) {
    int n = element.numberOfNodes;
    
    // 确保输出数组大小正确
    if (u.size() < static_cast<size_t>(n)) u.resize(n);
    if (v.size() < static_cast<size_t>(n)) v.resize(n);
    if (w.size() < static_cast<size_t>(n)) w.resize(n);
    
    int elemCode = element.elementCode;
    
    // 根据元素类型设置参考节点坐标
    switch (elemCode / 100) {
        case 2: // 线形元素
            for (int i = 0; i < n; ++i) {
                u[i] = -1.0 + 2.0 * i / (n - 1); // 均匀分布在[-1,1]区间
                v[i] = 0.0;
                w[i] = 0.0;
            }
            break;
            
        case 3: // 三角形元素
            // 三角形元素参考节点（简化实现）
            if (n == 3) {
                u[0] = 0.0; v[0] = 0.0; w[0] = 0.0;
                u[1] = 1.0; v[1] = 0.0; w[1] = 0.0;
                u[2] = 0.0; v[2] = 1.0; w[2] = 0.0;
            } else if (n == 6) {
                // 二次三角形元素
                u[0] = 0.0; v[0] = 0.0; w[0] = 0.0;
                u[1] = 1.0; v[1] = 0.0; w[1] = 0.0;
                u[2] = 0.0; v[2] = 1.0; w[2] = 0.0;
                u[3] = 0.5; v[3] = 0.0; w[3] = 0.0;
                u[4] = 0.5; v[4] = 0.5; w[4] = 0.0;
                u[5] = 0.0; v[5] = 0.5; w[5] = 0.0;
            }
            break;
            
        case 4: // 四边形元素
            // 四边形元素参考节点（简化实现）
            if (n == 4) {
                u[0] = -1.0; v[0] = -1.0; w[0] = 0.0;
                u[1] = 1.0; v[1] = -1.0; w[1] = 0.0;
                u[2] = 1.0; v[2] = 1.0; w[2] = 0.0;
                u[3] = -1.0; v[3] = 1.0; w[3] = 0.0;
            }
            break;
            
        case 5: // 四面体元素
            // 四面体元素参考节点（简化实现）
            if (n == 4) {
                u[0] = 0.0; v[0] = 0.0; w[0] = 0.0;
                u[1] = 1.0; v[1] = 0.0; w[1] = 0.0;
                u[2] = 0.0; v[2] = 1.0; w[2] = 0.0;
                u[3] = 0.0; v[3] = 0.0; w[3] = 1.0;
            }
            break;
            
        case 7: // 棱柱元素
            // 棱柱元素参考节点（简化实现）
            if (n == 6) {
                u[0] = 0.0; v[0] = 0.0; w[0] = -1.0;
                u[1] = 1.0; v[1] = 0.0; w[1] = -1.0;
                u[2] = 0.0; v[2] = 1.0; w[2] = -1.0;
                u[3] = 0.0; v[3] = 0.0; w[3] = 1.0;
                u[4] = 1.0; v[4] = 0.0; w[4] = 1.0;
                u[5] = 0.0; v[5] = 1.0; w[5] = 1.0;
            }
            break;
            
        default:
            // 默认情况：均匀分布
            for (int i = 0; i < n; ++i) {
                u[i] = -1.0 + 2.0 * i / (n - 1);
                v[i] = 0.0;
                w[i] = 0.0;
            }
            break;
    }
    
    // 确保所有元素都被正确初始化
    for (int i = 0; i < n; ++i) {
        if (i >= static_cast<int>(u.size())) u.resize(i+1);
        if (i >= static_cast<int>(v.size())) v.resize(i+1);
        if (i >= static_cast<int>(w.size())) w.resize(i+1);
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
void dLineNodalPBasisAll(double u, std::vector<double>& gradphi) {
    // 检查输出数组大小
    if (gradphi.size() < 2) {
        gradphi.resize(2);
    }
    
    const double half = 0.5;
    const double usgn = 1.0; // 符号因子
    const double c = half * usgn;
    
    gradphi[0] = -c;
    gradphi[1] = c;
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
    
    const double half = 0.5;
    const double c3 = 1.0 / std::sqrt(3.0);
    
    phi[0] = half * (1.0 - u - c3 * v);
    phi[1] = half * (1.0 + u - c3 * v);
    phi[2] = c3 * v;
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
 * @brief 计算一维P型基函数
 * 
 * 基于Legendre多项式计算一维P型基函数。
 * 这是一个简化实现，实际应用中可能需要更复杂的Legendre多项式计算。
 */
void Compute1DPBasis(std::vector<std::vector<double>>& basis, int n) {
    // 检查输出矩阵大小
    if (basis.size() < static_cast<size_t>(n) || basis[0].size() < static_cast<size_t>(n)) {
        basis.resize(n, std::vector<double>(n, 0.0));
    }
    
    // 简化实现：使用线性基函数
    // 在实际应用中，这里应该使用Legendre多项式
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                basis[i][j] = 1.0;
            } else {
                basis[i][j] = 0.0;
            }
        }
    }
}

// =============================================================================
// 边界处理函数实现
// =============================================================================

/**
 * @brief 辅助函数：对整数数组进行排序
 * 
 * 使用冒泡排序算法对整数数组进行排序。
 * 
 * @param arr 要排序的数组
 */
void sort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (arr[j] > arr[j + 1]) {
                std::swap(arr[j], arr[j + 1]);
            }
        }
    }
}

/**
 * @brief 获取三角形面的全局方向
 * 
 * 给定元素和其面映射（对于元素的某个三角形面），
 * 此例程返回三角形面的全局方向，以便函数在元素边界上连续。
 */
std::vector<int> getTriangleFaceDirection(const Element& element, 
                                         const std::vector<int>& faceMap,
                                         const std::vector<int>& indexes) {
    // 检查输入参数
    if (faceMap.size() != 3) {
        throw std::invalid_argument("getTriangleFaceDirection: 面映射必须包含3个节点");
    }
    
    if (indexes.size() < static_cast<size_t>(element.type.numberOfNodes)) {
        throw std::invalid_argument("getTriangleFaceDirection: 索引数组大小不足");
    }
    
    // 将面的全局节点放入排序后的数组中
    std::vector<int> nodes(3);
    for (int i = 0; i < 3; ++i) {
        int nodeIndex = faceMap[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(indexes.size())) {
            throw std::invalid_argument("getTriangleFaceDirection: 无效的面映射节点索引");
        }
        nodes[i] = indexes[nodeIndex];
    }
    
    // 对节点进行排序
    sort(nodes);
    
    // 查找排序后节点的局部编号
    std::vector<int> globalDir(3, 0);
    for (int i = 0; i < element.type.numberOfNodes; ++i) {
        if (nodes[0] == indexes[i]) {
            globalDir[0] = i;
        } else if (nodes[1] == indexes[i]) {
            globalDir[1] = i;
        } else if (nodes[2] == indexes[i]) {
            globalDir[2] = i;
        }
    }
    
    return globalDir;
}

/**
 * @brief 获取四边形面的全局方向
 * 
 * 给定元素和其面映射（对于元素的某个四边形面），
 * 此例程返回四边形面的全局方向，以便函数在元素边界上连续。
 */
std::vector<int> getSquareFaceDirection(const Element& element,
                                       const std::vector<int>& faceMap,
                                       const std::vector<int>& indexes) {
    // 检查输入参数
    if (faceMap.size() != 4) {
        throw std::invalid_argument("getSquareFaceDirection: 面映射必须包含4个节点");
    }
    
    if (indexes.size() < static_cast<size_t>(element.type.numberOfNodes)) {
        throw std::invalid_argument("getSquareFaceDirection: 索引数组大小不足");
    }
    
    // 获取全局节点
    std::vector<int> nodes(4);
    for (int i = 0; i < 4; ++i) {
        int nodeIndex = faceMap[i];
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(indexes.size())) {
            throw std::invalid_argument("getSquareFaceDirection: 无效的面映射节点索引");
        }
        nodes[i] = indexes[nodeIndex];
    }
    
    // 查找最小全局节点
    int minGlobal = nodes[0];
    int A = 0;
    for (int i = 1; i < 4; ++i) {
        if (nodes[i] < minGlobal) {
            A = i;
            minGlobal = nodes[i];
        }
    }
    
    // 选择B节点作为最小节点旁边的次小节点
    int B = (A + 1) % 4;
    int C = (A + 3) % 4;
    int D = (A + 2) % 4;
    
    if (nodes[B] > nodes[C]) {
        std::swap(B, C);
    }
    
    // 查找节点A、B和C的局部编号
    std::vector<int> globalDir(4, 0);
    for (int i = 0; i < element.type.numberOfNodes; ++i) {
        if (nodes[A] == indexes[i]) {
            globalDir[0] = i;
        } else if (nodes[B] == indexes[i]) {
            globalDir[1] = i;
        } else if (nodes[C] == indexes[i]) {
            globalDir[3] = i;
        } else if (nodes[D] == indexes[i]) {
            globalDir[2] = i;
        }
    }
    
    return globalDir;
}

/**
 * @brief 检查楔形元素面的局部编号是否合法
 * 
 * 检查给定的四边形面的局部编号对于楔形元素是否合法。
 */
bool wedgeOrdering(const std::vector<int>& ordering) {
    if (ordering.size() != 4) {
        throw std::invalid_argument("wedgeOrdering: 排序数组必须包含4个元素");
    }
    
    // 检查前两个节点是否在同一个三角形面上
    // 对于楔形元素，合法的四边形面应该满足：
    // 前两个节点在同一个三角形面上（1-3或4-6）
    if ((ordering[0] >= 0 && ordering[0] <= 2 &&
         ordering[1] >= 0 && ordering[1] <= 2) ||
        (ordering[0] >= 3 && ordering[0] <= 5 &&
         ordering[1] >= 3 && ordering[1] <= 5)) {
        return true;
    }
    
    return false;
}

} // namespace ElmerCpp