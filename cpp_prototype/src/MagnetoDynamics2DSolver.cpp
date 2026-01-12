#include "MagnetoDynamics2DSolver.h"
#include "ElectromagneticMaterial.h"
#include "ShapeFunctions.h"
#include "GaussIntegration.h"
#include "ElementMatrix.h"
#include <iostream>
#include <cmath>

namespace elmer {

void MagnetoDynamics2DSolver::assembleElementContributions() {
    if (!mesh) {
        throw std::runtime_error("Mesh not set for element assembly");
    }
    
    const auto& elements = mesh->getElements();
    
    // 并行组装单元贡献（对应Fortran的OpenMP并行化）
    for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
        const auto& element = elements[elementIndex];
        
        // 跳过退化单元（如果启用）
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            // 在轴对称情况下检查退化单元
            if (element.isDegenerate()) {
                std::cout << "Skipping degenerate element: " << elementIndex << std::endl;
                continue;
            }
        }
        
        // 计算单元局部矩阵
        std::vector<std::vector<double>> localStiffness;
        std::vector<std::vector<double>> localMass;
        std::vector<std::vector<double>> localDamping;
        std::vector<double> localForce;
        
        computeLocalMatrix(element, localStiffness, localMass, localDamping, localForce);
        
        // 将局部矩阵组装到全局矩阵
        assembleLocalToGlobal(element, localStiffness, localMass, localDamping, localForce);
    }
}

void MagnetoDynamics2DSolver::computeLocalMatrix(const Element& element,
                                                std::vector<std::vector<double>>& stiffness,
                                                std::vector<std::vector<double>>& mass,
                                                std::vector<std::vector<double>>& damping,
                                                std::vector<double>& force) {
    int nNodes = element.getNumberOfNodes();
    int nDofs = nNodes; // 每个节点1个自由度：A_z
    
    // 初始化局部矩阵
    stiffness.resize(nDofs, std::vector<double>(nDofs, 0.0));
    mass.resize(nDofs, std::vector<double>(nDofs, 0.0));
    damping.resize(nDofs, std::vector<double>(nDofs, 0.0));
    force.resize(nDofs, 0.0);
    
    // 获取材料属性
    auto material = materialDB.getMaterial(element.getMaterialId());
    
    // 电导率
    double conductivity = material.getElectricConductivity();
    
    // 磁化强度
    std::array<double, 2> magnetization = {0.0, 0.0};
    if (material.hasProperty("Magnetization 1") && material.hasProperty("Magnetization 2")) {
        magnetization[0] = material.getProperty("Magnetization 1");
        magnetization[1] = material.getProperty("Magnetization 2");
    }
    
    // 电流密度（体载荷）
    double currentDensity = 0.0;
    if (element.hasBodyForce()) {
        auto bodyForce = element.getBodyForce();
        if (bodyForce.hasProperty("Current Density")) {
            currentDensity = bodyForce.getProperty("Current Density");
        }
    }
    
    // 洛伦兹速度（对流项）
    std::array<double, 2> lorentzVelocity = {0.0, 0.0};
    bool withVelocity = false;
    if (element.hasBodyForce()) {
        auto bodyForce = element.getBodyForce();
        if (bodyForce.hasProperty("Lorentz velocity X") && 
            bodyForce.hasProperty("Lorentz velocity Y")) {
            lorentzVelocity[0] = bodyForce.getProperty("Lorentz velocity X");
            lorentzVelocity[1] = bodyForce.getProperty("Lorentz velocity Y");
            withVelocity = true;
        }
    }
    
    // 角速度（旋转效应）
    double angularVelocity = 0.0;
    bool withAngularVelocity = false;
    if (element.hasBodyForce()) {
        auto bodyForce = element.getBodyForce();
        if (bodyForce.hasProperty("Angular velocity")) {
            angularVelocity = bodyForce.getProperty("Angular velocity");
            withAngularVelocity = true;
        }
    }
    
    // 获取高斯积分点
    auto integrationPoints = GaussIntegration::getPoints(element.getType(), 2); // 二阶积分
    
    // 数值积分循环
    for (const auto& ip : integrationPoints) {
        // 计算基函数及其导数
        std::vector<double> basis(nNodes);
        std::vector<std::array<double, 3>> dBasisdx(nNodes);
        double detJ;
        
        bool success = ShapeFunctions::evaluate(element, ip.xi, ip.eta, ip.zeta,
                                               basis, dBasisdx, detJ);
        if (!success) {
            continue; // 跳过无效积分点
        }
        
        // 坐标系修正（轴对称情况）
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            // 计算积分点的径向坐标
            double x = 0.0;
            for (int i = 0; i < nNodes; ++i) {
                x += basis[i] * element.getNode(i).getX();
            }
            detJ *= x; // 轴对称修正
        }
        
        double weight = ip.weight * detJ;
        
        // 计算磁阻率张量
        std::array<std::array<double, 2>, 2> reluctivityTensor;
        
        // 如果是非线性材料，需要计算当前磁场下的磁阻率
        if (material.isNonlinearMagnetic()) {
            // 需要当前磁矢势来计算磁通密度
            double magneticFluxDensity = 0.0; // 这里需要从当前解计算
            reluctivityTensor = computeReluctivity(material, magneticFluxDensity, element);
        } else {
            // 线性材料，使用恒定磁阻率
            double reluctivity = 1.0 / material.getMagneticPermeability();
            reluctivityTensor[0][0] = reluctivity;
            reluctivityTensor[0][1] = 0.0;
            reluctivityTensor[1][0] = 0.0;
            reluctivityTensor[1][1] = reluctivity;
        }
        
        // 计算B矩阵（磁通密度与基函数导数的关系）
        std::vector<std::array<double, 2>> Bt(nNodes);
        for (int i = 0; i < nNodes; ++i) {
            // 在2D情况下：B = curl A = (dA/dy, -dA/dx)
            Bt[i][0] = dBasisdx[i][1];  // dA/dy
            Bt[i][1] = -dBasisdx[i][0]; // -dA/dx
        }
        
        // 坐标系修正（轴对称情况）
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            double x = 0.0;
            for (int i = 0; i < nNodes; ++i) {
                x += basis[i] * element.getNode(i).getX();
            }
            
            for (int i = 0; i < nNodes; ++i) {
                Bt[i][0] = -Bt[i][0]; // 符号修正
                Bt[i][1] = -Bt[i][1];
                Bt[i][1] += basis[i] / x; // 轴对称附加项
            }
        }
        
        // 计算H矩阵（磁场强度与基函数的关系）
        std::vector<std::array<double, 2>> Ht(nNodes);
        for (int i = 0; i < nNodes; ++i) {
            Ht[i][0] = reluctivityTensor[0][0] * Bt[i][0] + reluctivityTensor[0][1] * Bt[i][1];
            Ht[i][1] = reluctivityTensor[1][0] * Bt[i][0] + reluctivityTensor[1][1] * Bt[i][1];
        }
        
        // 组装刚度矩阵贡献
        for (int i = 0; i < nNodes; ++i) {
            for (int j = 0; j < nNodes; ++j) {
                double contribution = 0.0;
                for (int k = 0; k < 2; ++k) {
                    contribution += Ht[i][k] * Bt[j][k];
                }
                stiffness[i][j] += weight * contribution;
            }
        }
        
        // 瞬态项（质量矩阵和阻尼矩阵）
        if (parameters.isTransient && conductivity != 0.0) {
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    // 质量矩阵（电导率相关）
                    mass[i][j] += weight * conductivity * basis[i] * basis[j];
                    
                    // 如果是电动力学模型，包含位移电流项
                    if (parameters.includeDisplacementCurrent) {
                        double permittivity = material.getElectricPermittivity();
                        damping[i][j] += weight * permittivity * basis[i] * basis[j];
                    }
                }
            }
        }
        
        // 对流项（洛伦兹速度或角速度）
        if (withVelocity || withAngularVelocity) {
            std::array<double, 2> velocity = {0.0, 0.0};
            
            if (withVelocity) {
                velocity = lorentzVelocity;
            } else if (withAngularVelocity) {
                // 计算积分点坐标
                double x = 0.0, y = 0.0;
                for (int i = 0; i < nNodes; ++i) {
                    x += basis[i] * element.getNode(i).getX();
                    y += basis[i] * element.getNode(i).getY();
                }
                
                // 角速度产生的速度场：v = ω × r
                velocity[0] = -angularVelocity * y;
                velocity[1] = angularVelocity * x;
            }
            
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    double convectionTerm = 0.0;
                    
                    if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
                        parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
                        convectionTerm = -velocity[1] * Bt[j][0] + velocity[0] * Bt[j][1];
                    } else {
                        convectionTerm = velocity[1] * Bt[j][0] - velocity[0] * Bt[j][1];
                    }
                    
                    stiffness[i][j] += weight * conductivity * basis[i] * convectionTerm;
                }
            }
        }
        
        // 右端向量（源项）
        double loadAtIP = currentDensity;
        
        for (int i = 0; i < nNodes; ++i) {
            if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
                parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
                force[i] += weight * (loadAtIP * basis[i] - 
                    magnetization[0] * dBasisdx[i][1] + 
                    magnetization[1] * (dBasisdx[i][0] + basis[i] / x));
            } else {
                force[i] += weight * (loadAtIP * basis[i] + 
                    magnetization[0] * dBasisdx[i][1] - 
                    magnetization[1] * dBasisdx[i][0]);
            }
        }
    }
}

std::array<std::array<double, 2>, 2> MagnetoDynamics2DSolver::computeReluctivity(
    const Material& material, double magneticFluxDensity, const Element& element) {
    
    std::array<std::array<double, 2>, 2> reluctivityTensor;
    
    // 如果是非线性磁性材料，使用B-H曲线
    if (material.isNonlinearMagnetic()) {
        // 从B-H曲线计算磁阻率
        double B = std::abs(magneticFluxDensity);
        double H = material.getMagneticFieldStrength(B);
        double mu = B / H; // 磁导率
        double nu = 1.0 / mu; // 磁阻率
        
        // 各向同性材料，假设为对角张量
        reluctivityTensor[0][0] = nu;
        reluctivityTensor[0][1] = 0.0;
        reluctivityTensor[1][0] = 0.0;
        reluctivityTensor[1][1] = nu;
        
        // 如果是牛顿-拉夫逊迭代，还需要计算导数
        if (parameters.useNewtonRaphson) {
            // 这里需要计算dnu/dB，用于雅可比矩阵
            // 简化处理：使用有限差分近似
            double dB = 1.0e-6;
            double H_plus = material.getMagneticFieldStrength(B + dB);
            double mu_plus = (B + dB) / H_plus;
            double nu_plus = 1.0 / mu_plus;
            
            double dnu_dB = (nu_plus - nu) / dB;
            // 雅可比矩阵贡献将在后续实现中处理
        }
    } else {
        // 线性材料
        double mu = material.getMagneticPermeability();
        double nu = 1.0 / mu;
        
        reluctivityTensor[0][0] = nu;
        reluctivityTensor[0][1] = 0.0;
        reluctivityTensor[1][0] = 0.0;
        reluctivityTensor[1][1] = nu;
    }
    
    return reluctivityTensor;
}

void MagnetoDynamics2DSolver::assembleLocalToGlobal(const Element& element,
                                                   const std::vector<std::vector<double>>& localStiffness,
                                                   const std::vector<std::vector<double>>& localMass,
                                                   const std::vector<std::vector<double>>& localDamping,
                                                   const std::vector<double>& localForce) {
    
    int nNodes = element.getNumberOfNodes();
    
    for (int i = 0; i < nNodes; ++i) {
        int globalI = element.getNode(i).getGlobalIndex();
        
        // 组装右端向量
        rhsVector->SetElement(globalI, rhsVector->GetElement(globalI) + localForce[i]);
        
        for (int j = 0; j < nNodes; ++j) {
            int globalJ = element.getNode(j).getGlobalIndex();
            
            // 组装刚度矩阵
            stiffnessMatrix->AddToElement(globalI, globalJ, localStiffness[i][j]);
            
            // 组装质量矩阵和阻尼矩阵（如果存在）
            if (massMatrix) {
                massMatrix->AddToElement(globalI, globalJ, localMass[i][j]);
            }
            if (dampingMatrix) {
                dampingMatrix->AddToElement(globalI, globalJ, localDamping[i][j]);
            }
        }
    }
}

void MagnetoDynamics2DSolver::assembleBoundaryContributions() {
    if (!mesh) {
        throw std::runtime_error("Mesh not set for boundary assembly");
    }
    
    const auto& boundaryElements = mesh->getBoundaryElements();
    
    for (const auto& element : boundaryElements) {
        auto bc = bcManager.getBoundaryConditionForElement(element);
        
        if (bc) {
            // 处理不同类型的边界条件
            switch (bc->getType()) {
                case BoundaryConditionType::INFINITY_BC:
                    applyInfinityBoundaryCondition(element, stiffnessMatrix, rhsVector);
                    break;
                case BoundaryConditionType::AIR_GAP_BC:
                    applyAirGapBoundaryCondition(element, stiffnessMatrix, rhsVector);
                    break;
                default:
                    // 其他边界条件由基类处理
                    break;
            }
        }
    }
}

void MagnetoDynamics2DSolver::applyInfinityBoundaryCondition(const Element& element,
                                                            std::shared_ptr<Matrix> stiffness,
                                                            std::shared_ptr<Vector> force) {
    // 无限远边界条件的实现
    // 对应Fortran的LocalMatrixInfinityBC子程序
    
    int nNodes = element.getNumberOfNodes();
    
    // 获取高斯积分点
    auto integrationPoints = GaussIntegration::getPoints(element.getType(), 2);
    
    for (const auto& ip : integrationPoints) {
        // 计算基函数及其导数
        std::vector<double> basis(nNodes);
        std::vector<std::array<double, 3>> dBasisdx(nNodes);
        double detJ;
        
        bool success = ShapeFunctions::evaluate(element, ip.xi, ip.eta, ip.zeta,
                                               basis, dBasisdx, detJ);
        if (!success) continue;
        
        double weight = ip.weight * detJ;
        
        // 无限远边界条件的典型实现：添加一个惩罚项
        // 这里使用简化的实现，实际应根据具体物理模型调整
        double penaltyParameter = 1.0e6; // 惩罚参数
        
        for (int i = 0; i < nNodes; ++i) {
            int globalI = element.getNode(i).getGlobalIndex();
            
            for (int j = 0; j < nNodes; ++j) {
                int globalJ = element.getNode(j).getGlobalIndex();
                
                stiffness->AddToElement(globalI, globalJ, 
                    weight * penaltyParameter * basis[i] * basis[j]);
            }
        }
    }
}

void MagnetoDynamics2DSolver::applyAirGapBoundaryCondition(const Element& element,
                                                          std::shared_ptr<Matrix> stiffness,
                                                          std::shared_ptr<Vector> force) {
    // 气隙边界条件的实现
    // 对应Fortran的LocalMatrixAirGapBC子程序
    
    int nNodes = element.getNumberOfNodes();
    
    // 获取高斯积分点
    auto integrationPoints = GaussIntegration::getPoints(element.getType(), 2);
    
    for (const auto& ip : integrationPoints) {
        // 计算基函数及其导数
        std::vector<double> basis(nNodes);
        std::vector<std::array<double, 3>> dBasisdx(nNodes);
        double detJ;
        
        bool success = ShapeFunctions::evaluate(element, ip.xi, ip.eta, ip.zeta,
                                               basis, dBasisdx, detJ);
        if (!success) continue;
        
        double weight = ip.weight * detJ;
        
        // 气隙边界条件的典型实现：考虑气隙的磁导
        // 这里使用简化的实现，实际应根据具体物理模型调整
        double airGapPermeability = 4.0 * M_PI * 1e-7; // 真空磁导率
        double airGapLength = 0.001; // 假设气隙长度1mm
        
        double couplingCoefficient = airGapPermeability / airGapLength;
        
        for (int i = 0; i < nNodes; ++i) {
            int globalI = element.getNode(i).getGlobalIndex();
            
            for (int j = 0; j < nNodes; ++j) {
                int globalJ = element.getNode(j).getGlobalIndex();
                
                stiffness->AddToElement(globalI, globalJ, 
                    weight * couplingCoefficient * basis[i] * basis[j]);
            }
        }
    }
}

void MagnetoDynamics2DSolver::applyBoundaryConditions() {
    // 应用边界条件
    bcManager.applyBoundaryConditions(*stiffnessMatrix, *rhsVector);
    
    // 如果是瞬态分析，也需要对质量矩阵应用边界条件
    if (massMatrix) {
        bcManager.applyBoundaryConditions(*massMatrix);
    }
    if (dampingMatrix) {
        bcManager.applyBoundaryConditions(*dampingMatrix);
    }
}

std::vector<double> MagnetoDynamics2DSolver::solveLinearSystem() {
    // 使用迭代求解器求解线性系统
    IterativeSolver solver;
    
    // 设置求解器参数
    solver.setTolerance(parameters.tolerance);
    solver.setMaxIterations(parameters.maxIterations);
    
    // 求解线性系统
    auto solution = solver.solve(*stiffnessMatrix, *rhsVector);
    
    return solution;
}

bool MagnetoDynamics2DSolver::checkConvergence(const MagnetoDynamics2DResults& results) {
    // 简化的收敛检查
    // 实际实现应该计算残差并与容差比较
    
    if (results.nonlinearIterations >= parameters.maxNonlinearIterations) {
        std::cout << "Maximum nonlinear iterations reached" << std::endl;
        return true;
    }
    
    // 这里应该实现更复杂的收敛检查逻辑
    // 包括残差计算、相对变化检查等
    
    return results.nonlinearIterations > 1; // 简化实现
}

void MagnetoDynamics2DSolver::calculateDerivedFields(MagnetoDynamics2DResults& results) {
    if (!mesh || results.vectorPotential.empty()) {
        return;
    }
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    results.magneticFluxDensity.resize(nNodes);
    results.magneticFieldStrength.resize(nNodes);
    results.currentDensity.resize(nNodes, 0.0);
    
    // 计算每个节点的磁场量
    for (size_t nodeIndex = 0; nodeIndex < nNodes; ++nodeIndex) {
        // 这里需要实现磁场量的计算
        // 实际实现应该使用基函数导数来计算磁场
        
        // 简化实现：使用相邻单元的平均值
        results.magneticFluxDensity[nodeIndex] = {0.0, 0.0};
        results.magneticFieldStrength[nodeIndex] = {0.0, 0.0};
    }
}

void MagnetoDynamics2DSolver::calculateLumpedParameters(MagnetoDynamics2DResults& results) {
    // 计算集总参数（转矩、磁能、电感等）
    // 对应Fortran的CalculateLumpedTransient子程序
    
    // 简化实现
    results.torque = 0.0;
    results.magneticEnergy = 0.0;
    results.inductance = 0.0;
    
    // 实际实现需要积分计算这些参数
}

void MagnetoDynamics2DSolver::reassembleNonlinearSystem() {
    // 重新组装非线性系统（用于牛顿-拉夫逊迭代）
    
    // 重置系统矩阵
    stiffnessMatrix->Zero();
    rhsVector->Zero();
    
    if (massMatrix) massMatrix->Zero();
    if (dampingMatrix) dampingMatrix->Zero();
    
    // 重新组装
    assembleElementContributions();
    assembleBoundaryContributions();
    applyBoundaryConditions();
}

void MagnetoDynamics2DSolver::initializeBasisFunctionCache() {
    // 初始化基函数缓存（性能优化）
    // 对应Fortran的TabulateBasisFunctions子程序
    
    if (!mesh) {
        return;
    }
    
    const auto& elements = mesh->getElements();
    
    // 计算总积分点数
    size_t totalIntegrationPoints = 0;
    for (const auto& element : elements) {
        auto integrationPoints = GaussIntegration::getPoints(element.getType(), 2);
        totalIntegrationPoints += integrationPoints.size();
    }
    
    basisCache.resize(totalIntegrationPoints);
    
    size_t cacheIndex = 0;
    for (const auto& element : elements) {
        auto integrationPoints = GaussIntegration::getPoints(element.getType(), 2);
        
        for (const auto& ip : integrationPoints) {
            int nNodes = element.getNumberOfNodes();
            
            std::vector<double> basis(nNodes);
            std::vector<std::array<double, 3>> dBasisdx(nNodes);
            double detJ;
            
            bool success = ShapeFunctions::evaluate(element, ip.xi, ip.eta, ip.zeta,
                                                   basis, dBasisdx, detJ);
            
            if (success) {
                basisCache[cacheIndex].basis = basis;
                basisCache[cacheIndex].dBasisdx = dBasisdx;
                basisCache[cacheIndex].weight = ip.weight * detJ;
                
                // 坐标系修正
                if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
                    parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
                    double x = 0.0;
                    for (int i = 0; i < nNodes; ++i) {
                        x += basis[i] * element.getNode(i).getX();
                    }
                    basisCache[cacheIndex].weight *= x;
                }
            }
            
            ++cacheIndex;
        }
    }
    
    useBasisFunctionsCache = true;
}

} // namespace elmer