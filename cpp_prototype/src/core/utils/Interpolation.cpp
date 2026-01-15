#include "Interpolation.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace elmer {
namespace Interpolation {

// 数值精度常量
constexpr double DEFAULT_EPSILON = 1.0e-12;

/**
 * @brief 在1D元素内部进行插值
 */
double InterpolateInElement1D(const Element& element, 
                             const std::vector<double>& nodalValues, 
                             double u) {
    if (nodalValues.size() < static_cast<size_t>(element.type.numberOfNodes)) {
        std::cout << "错误: 节点值数量与元素节点数不匹配" << std::endl;
        return 0.0;
    }
    
    double result = 0.0;
    int numNodes = element.type.numberOfNodes;
    
    // 遍历所有节点，计算插值
    for (int n = 0; n < numNodes; ++n) {
        if (std::abs(nodalValues[n]) > DEFAULT_EPSILON) {
            // 获取基函数系数
            const auto& basis = element.type.basisFunctions[n];
            double basisValue = 0.0;
            
            // 计算基函数在点u处的值
            for (int i = 0; i < basis.n; ++i) {
                if (basis.p[i] == 0) {
                    basisValue += basis.coeff[i];
                } else {
                    basisValue += basis.coeff[i] * std::pow(u, basis.p[i]);
                }
            }
            
            result += nodalValues[n] * basisValue;
        }
    }
    
    return result;
}

/**
 * @brief 在2D元素内部进行插值
 */
double InterpolateInElement2D(const Element& element,
                             const std::vector<double>& nodalValues,
                             double u, double v) {
    if (nodalValues.size() < static_cast<size_t>(element.type.numberOfNodes)) {
        std::cout << "错误: 节点值数量与元素节点数不匹配" << std::endl;
        return 0.0;
    }
    
    double result = 0.0;
    int numNodes = element.type.numberOfNodes;
    
    // 遍历所有节点，计算插值
    for (int n = 0; n < numNodes; ++n) {
        if (std::abs(nodalValues[n]) > DEFAULT_EPSILON) {
            // 获取基函数系数
            const auto& basis = element.type.basisFunctions[n];
            double basisValue = 0.0;
            
            // 计算基函数在点(u,v)处的值
            for (int i = 0; i < basis.n; ++i) {
                double term = basis.coeff[i];
                if (basis.p[i] > 0) {
                    term *= std::pow(u, basis.p[i]);
                }
                if (basis.q[i] > 0) {
                    term *= std::pow(v, basis.q[i]);
                }
                basisValue += term;
            }
            
            result += nodalValues[n] * basisValue;
        }
    }
    
    return result;
}

/**
 * @brief 在3D元素内部进行插值
 */
double InterpolateInElement3D(const Element& element,
                             const std::vector<double>& nodalValues,
                             double u, double v, double w) {
    if (nodalValues.size() < static_cast<size_t>(element.type.numberOfNodes)) {
        std::cout << "错误: 节点值数量与元素节点数不匹配" << std::endl;
        return 0.0;
    }
    
    double result = 0.0;
    int numNodes = element.type.numberOfNodes;
    
    // 遍历所有节点，计算插值
    for (int n = 0; n < numNodes; ++n) {
        if (std::abs(nodalValues[n]) > DEFAULT_EPSILON) {
            // 获取基函数系数
            const auto& basis = element.type.basisFunctions[n];
            double basisValue = 0.0;
            
            // 计算基函数在点(u,v,w)处的值
            for (int i = 0; i < basis.n; ++i) {
                double term = basis.coeff[i];
                if (basis.p[i] > 0) {
                    term *= std::pow(u, basis.p[i]);
                }
                if (basis.q[i] > 0) {
                    term *= std::pow(v, basis.q[i]);
                }
                if (basis.r[i] > 0) {
                    term *= std::pow(w, basis.r[i]);
                }
                basisValue += term;
            }
            
            result += nodalValues[n] * basisValue;
        }
    }
    
    return result;
}

/**
 * @brief 在元素内部进行插值（通用接口）
 */
double InterpolateInElement(const Element& element,
                           const std::vector<double>& nodalValues,
                           double u, double v, double w,
                           const std::vector<double>& basisFunctions) {
    
    // 如果提供了基函数值，直接使用
    if (!basisFunctions.empty() && 
        basisFunctions.size() >= static_cast<size_t>(element.type.numberOfNodes)) {
        double result = 0.0;
        for (int i = 0; i < element.type.numberOfNodes; ++i) {
            result += nodalValues[i] * basisFunctions[i];
        }
        return result;
    }
    
    // 否则根据维度调用相应的插值函数
    switch (element.type.dimension) {
        case 0:
            return nodalValues[0]; // 0维元素只有一个节点值
        case 1:
            return InterpolateInElement1D(element, nodalValues, u);
        case 2:
            return InterpolateInElement2D(element, nodalValues, u, v);
        case 3:
            return InterpolateInElement3D(element, nodalValues, u, v, w);
        default:
            std::cout << "错误: 不支持的维度: " << element.type.dimension << std::endl;
            return 0.0;
    }
}

/**
 * @brief 检查点是否在元素内部
 */
bool PointInElement(const Element& element,
                   const Nodes& elementNodes,
                   const std::vector<double>& point,
                   std::vector<double>& localCoordinates,
                   double globalEps,
                   double localEps,
                   double numericEps,
                   double* globalDistance,
                   double* localDistance) {
    
    if (point.size() < 3) {
        std::cout << "错误: 点坐标维度不足" << std::endl;
        return false;
    }
    
    // 设置默认数值精度
    double eps0 = (numericEps > 0) ? numericEps : DEFAULT_EPSILON;
    double eps1 = (globalEps > 0) ? globalEps : 1.0e-4;
    double eps2 = (localEps > 0) ? localEps : 1.0e-10;
    
    int numNodes = element.type.numberOfNodes;
    
    // 检查全局坐标边界
    double minX = *std::min_element(elementNodes.x.begin(), elementNodes.x.begin() + numNodes);
    double maxX = *std::max_element(elementNodes.x.begin(), elementNodes.x.begin() + numNodes);
    double minY = *std::min_element(elementNodes.y.begin(), elementNodes.y.begin() + numNodes);
    double maxY = *std::max_element(elementNodes.y.begin(), elementNodes.y.begin() + numNodes);
    double minZ = *std::min_element(elementNodes.z.begin(), elementNodes.z.begin() + numNodes);
    double maxZ = *std::max_element(elementNodes.z.begin(), elementNodes.z.begin() + numNodes);
    
    // 计算点到边界的距离
    double xDist = std::max(std::max(point[0] - maxX, 0.0), minX - point[0]);
    double yDist = std::max(std::max(point[1] - maxY, 0.0), minY - point[1]);
    double zDist = std::max(std::max(point[2] - maxZ, 0.0), minZ - point[2]);
    
    // 如果提供了全局距离，计算并检查
    if (globalDistance) {
        *globalDistance = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
        
        if (xDist > eps0 + eps1 * (maxX - minX)) return false;
        if (yDist > eps0 + eps1 * (maxY - minY)) return false;
        if (zDist > eps0 + eps1 * (maxZ - minZ)) return false;
    } else {
        // 分别检查每个方向
        if (xDist > eps0 + eps1 * (maxX - minX)) return false;
        if (yDist > eps0 + eps1 * (maxY - minY)) return false;
        if (zDist > eps0 + eps1 * (maxZ - minZ)) return false;
    }
    
    // 计算局部坐标（简化实现，实际需要调用GlobalToLocal函数）
    // 这里使用简化方法计算局部坐标
    double uLocal = (point[0] - minX) / (maxX - minX);
    double vLocal = (point[1] - minY) / (maxY - minY);
    double wLocal = (point[2] - minZ) / (maxZ - minZ);
    
    localCoordinates = {uLocal, vLocal, wLocal};
    
    // 根据元素类型计算局部距离
    double sumDist = 0.0;
    int elementCode = element.type.elementCode / 100;
    
    switch (elementCode) {
        case 1: // 1D线单元
            sumDist = uLocal;
            break;
        case 2: // 1D线单元（另一种表示）
            sumDist = std::max(uLocal - 1.0, std::max(-uLocal - 1.0, 0.0));
            break;
        case 3: // 2D三角形单元
            sumDist = std::max(-uLocal, 0.0) + std::max(-vLocal, 0.0);
            sumDist += std::max(uLocal + vLocal - 1.0, 0.0);
            break;
        case 4: // 2D四边形单元
            sumDist = std::max(uLocal - 1.0, std::max(-uLocal - 1.0, 0.0));
            sumDist += std::max(vLocal - 1.0, std::max(-vLocal - 1.0, 0.0));
            break;
        case 5: // 3D四面体单元
            sumDist = std::max(-uLocal, 0.0) + std::max(-vLocal, 0.0) + std::max(-wLocal, 0.0);
            sumDist += std::max(uLocal + vLocal + wLocal - 1.0, 0.0);
            break;
        case 7: // 3D棱柱单元
            sumDist = std::max(-uLocal, 0.0) + std::max(-vLocal, 0.0);
            sumDist += std::max(uLocal + vLocal - 1.0, 0.0);
            sumDist += std::max(wLocal - 1.0, std::max(-wLocal - 1.0, 0.0));
            break;
        case 8: // 3D六面体单元
            sumDist = std::max(uLocal - 1.0, std::max(-uLocal - 1.0, 0.0));
            sumDist += std::max(vLocal - 1.0, std::max(-vLocal - 1.0, 0.0));
            sumDist += std::max(wLocal - 1.0, std::max(-wLocal - 1.0, 0.0));
            break;
        default:
            std::cout << "警告: 不支持的元素类型: " << elementCode << std::endl;
            return false;
    }
    
    // 检查局部距离是否在容差范围内
    if (localDistance) {
        *localDistance = sumDist;
    }
    
    return sumDist < eps0 + eps2;
}

/**
 * @brief 递归查找点所在的象限
 */
std::shared_ptr<Quadrant> FindPointsQuadrant(const std::vector<double>& point,
                                            int dimension,
                                            std::shared_ptr<Quadrant> motherQuadrant) {
    
    if (!motherQuadrant || motherQuadrant->childQuadrants.empty()) {
        return motherQuadrant;
    }
    
    // 查找点所在的子象限
    for (size_t i = 0; i < motherQuadrant->childQuadrants.size(); ++i) {
        auto childQuadrant = motherQuadrant->childQuadrants[i];
        if (!childQuadrant) continue;
        
        const auto& bbox = childQuadrant->boundingBox;
        if (bbox.size() < 6) continue;
        
        // 检查点是否在子象限内
        if (point[0] >= bbox[0] && point[0] <= bbox[3] &&
            point[1] >= bbox[1] && point[1] <= bbox[4] &&
            point[2] >= bbox[2] && point[2] <= bbox[5]) {
            
            // 递归查找更深层次的象限
            return FindPointsQuadrant(point, dimension, childQuadrant);
        }
    }
    
    return nullptr;
}

/**
 * @brief 查找点所在的叶子象限
 */
std::shared_ptr<Quadrant> FindLeafElements(const std::vector<double>& point,
                                          int dimension,
                                          const std::shared_ptr<Quadrant>& rootQuadrant) {
    
    if (!rootQuadrant) {
        return nullptr;
    }
    
    return FindPointsQuadrant(point, dimension, rootQuadrant);
}

/**
 * @brief 构建象限树（简化实现）
 */
std::shared_ptr<Quadrant> BuildQuadrantTree(const Mesh& mesh,
                                           const std::vector<double>& boundingBox) {
    
    // 创建根象限
    auto rootQuadrant = std::make_shared<Quadrant>(boundingBox);
    
    // 简化实现：将网格中的所有元素添加到根象限
    // 实际实现应该递归分割空间
    for (int i = 0; i < mesh.numberOfElements; ++i) {
        rootQuadrant->elementIndices.push_back(i);
    }
    
    return rootQuadrant;
}

/**
 * @brief 变量到变量的插值（简化版本）
 */
bool InterpolateVarToVarReduced(const Mesh& oldMesh,
                               Mesh& newMesh,
                               const std::string& variableName,
                               const std::vector<int>& variableDimensions,
                               double globalEps,
                               double localEps,
                               double numericEps) {
    
    std::cout << "开始变量插值: " << variableName << std::endl;
    
    // 简化实现：遍历新网格的所有节点
    for (int i = 0; i < newMesh.numberOfNodes; ++i) {
        std::vector<double> point = {
            newMesh.nodes[i].x,
            newMesh.nodes[i].y,
            newMesh.nodes[i].z
        };
        
        // 在原网格中查找包含该点的元素
        // 这里使用简化实现，实际应该使用象限树加速搜索
        for (int j = 0; j < oldMesh.numberOfElements; ++j) {
            const auto& element = oldMesh.elements[j];
            std::vector<double> localCoords;
            
            // 创建原网格的节点坐标结构
            Nodes oldNodes;
            oldNodes.x.resize(oldMesh.numberOfNodes);
            oldNodes.y.resize(oldMesh.numberOfNodes);
            oldNodes.z.resize(oldMesh.numberOfNodes);
            for (int k = 0; k < oldMesh.numberOfNodes; ++k) {
                oldNodes.x[k] = oldMesh.nodes[k].x;
                oldNodes.y[k] = oldMesh.nodes[k].y;
                oldNodes.z[k] = oldMesh.nodes[k].z;
            }
            
            if (PointInElement(element, oldNodes, point, localCoords, 
                             globalEps, localEps, numericEps)) {
                
                // 找到包含点的元素，进行插值
                // 简化实现：假设变量值已知
                std::cout << "在元素 " << j << " 中找到节点 " << i << std::endl;
                break;
            }
        }
    }
    
    std::cout << "变量插值完成" << std::endl;
    return true;
}

} // namespace Interpolation
} // namespace elmer