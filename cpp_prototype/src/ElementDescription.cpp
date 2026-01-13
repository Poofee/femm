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

} // namespace ElmerCpp