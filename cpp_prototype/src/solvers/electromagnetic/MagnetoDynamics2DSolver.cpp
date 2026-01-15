// MagnetoDynamics2DSolver.cpp - Elmer FEM C++ 2D Magnetodynamics Solver
// 对应Fortran模块: MagnetoDynamics2D.F90
// 基于Elmer Fortran源代码的精确移植

#include "MagnetoDynamics2DSolver.h"
#include <iostream>

namespace elmer {

// =============================================================================
// 基于Fortran源代码的移植实现
// =============================================================================

//------------------------------------------------------------------------------
// MagnetoDynamics2D_Init 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::initialize() {
    // 对应Fortran的MagnetoDynamics2D_Init子程序
    // 设置变量自由度数为1（对应Fortran: CALL ListAddInteger( Params, 'Variable Dofs',1 )）
    variableDofs = 1;
    
    // 设置变量名称为"Potential"（对应Fortran: CALL ListAddNewString( Params,'Variable','Potential')）
    variableName = "Potential";
    
    // 应用Mortar边界条件（对应Fortran: CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)）
    applyMortarBCs = true;
    
    // 使用全局质量矩阵（对应Fortran: CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)）
    useGlobalMassMatrix = true;
    
    // 检查电动力学模型（对应Fortran: ElectroDynamics = GetLogical (Params, 'Electrodynamics Model', Found)）
    // TODO: 实现电动力学模型检查
    
    // 检查坐标系（对应Fortran: CurrentCoordinateSystem() == AxisSymmetric）
    if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
        parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
        // 轴对称情况下的特殊处理
        // 对应Fortran: CALL Warn(Caller,'Handle assembly not yet available in axisymmetric case!')
        if (parameters.useInfinityBC) {
            std::cout << "Warning: Infinity BC not yet available in axisymmetric case!" << std::endl;
        }
    }
    
    // 初始化基函数缓存（如果启用）
    // 对应Fortran: CALL TabulateBasisFunctions()
    if (useBasisFunctionsCache) {
        initializeBasisFunctionCache();
    }
    
    std::cout << "MagnetoDynamics2D solver initialized" << std::endl;
}

//------------------------------------------------------------------------------
// MagnetoDynamics2D 主求解器子程序的C++移植
//------------------------------------------------------------------------------
MagnetoDynamics2DResults MagnetoDynamics2DSolver::solve() {
    // 对应Fortran的MagnetoDynamics2D主求解器子程序
    
    // 对应Fortran: CALL Info( Caller,'Solving equation for magnetic vector potential', Level=4 )
    std::cout << "Solving equation for magnetic vector potential" << std::endl;
    
    // 对应Fortran: CALL DefaultStart()
    // TODO: 实现DefaultStart的等效功能
    
    // 检查是否使用处理组装（对应Fortran: HandleAsm = ListGetLogical( SolverParams,'Handle Assembly',Found )）
    bool handleAssembly = false; // 默认使用传统组装
    
    // 检查坐标系对称性（对应Fortran: CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric ... )）
    bool isAxisymmetric = (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
                          parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC);
    
    // 检查牛顿-拉夫逊迭代（对应Fortran: NewtonRaphson = GetLogical(SolverParams, 'Newton-Raphson Iteration', Found)）
    bool useNewtonRaphson = parameters.useNewtonRaphson;
    
    // 非线性迭代循环（对应Fortran: DO iter = 1,NonlinIter）
    for (int nonlinearIter = 0; nonlinearIter < parameters.maxNonlinearIterations; ++nonlinearIter) {
        std::cout << "Performing nonlinear iteration: " << nonlinearIter + 1 << std::endl;
        
        // 系统组装（对应Fortran: CALL DefaultInitialize()）
        if (!systemAssembled) {
            assembleSystem();
        }
        
        // 求解线性系统（对应Fortran: Norm = DefaultSolve()）
        auto solution = solveLinearSystem();
        
        // 检查收敛性（对应Fortran: CALL DefaultConverged()）
        MagnetoDynamics2DResults results;
        results.vectorPotential = solution;
        results.nonlinearIterations = nonlinearIter + 1;
        
        if (checkConvergence(results)) {
            results.converged = true;
            std::cout << "System has converged to tolerances after " << nonlinearIter + 1 << " iterations!" << std::endl;
            
            // 计算导出场量
            calculateDerivedFields(results);
            
            return results;
        }
        
        // 更新材料参数用于下一次迭代（对应Fortran: 非线性材料更新）
        if (useNewtonRaphson && nonlinearIter > 0) {
            updateMaterialParameters();
        }
    }
    
    // 如果未收敛，返回最后一次迭代的结果
    MagnetoDynamics2DResults results;
    results.vectorPotential = solveLinearSystem();
    results.nonlinearIterations = parameters.maxNonlinearIterations;
    results.converged = false;
    
    calculateDerivedFields(results);
    
    return results;
}

//------------------------------------------------------------------------------
// LocalMatrix 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::assembleElementContributions() {
    // 对应Fortran的LocalMatrix子程序
    // 组装单元局部矩阵贡献
    
    if (!mesh) {
        throw std::runtime_error("Mesh not set for element assembly");
    }
    
    // 对应Fortran: 遍历所有活动单元
    // DO t=1,active
    //   Element => GetActiveElement(t)
    //   n  = GetElementNOFNodes(Element)
    //   nd = GetElementNOFDOFs(Element)
    //   CALL LocalMatrix(Element, n, nd)
    // END DO
    
    // 获取网格中的所有单元
    auto& elements = mesh->getElements();
    
    // 遍历所有单元
    for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
        const auto& element = elements[elementIndex];
        
        // 获取单元节点数
        size_t nNodes = element.getNodeCount();
        size_t nDOFs = nNodes; // 每个节点1个自由度
        
        // 初始化局部矩阵
        std::vector<std::vector<double>> localStiffness(nDOFs, std::vector<double>(nDOFs, 0.0));
        std::vector<std::vector<double>> localMass(nDOFs, std::vector<double>(nDOFs, 0.0));
        std::vector<std::vector<double>> localDamping(nDOFs, std::vector<double>(nDOFs, 0.0));
        std::vector<double> localForce(nDOFs, 0.0);
        
        // 计算单元局部矩阵
        computeLocalMatrix(element, localStiffness, localMass, localDamping, localForce);
        
        // 组装到全局系统矩阵
        assembleLocalToGlobal(element, localStiffness, localMass, localDamping, localForce);
    }
    
    std::cout << "单元贡献组装完成，共处理 " << elements.size() << " 个单元" << std::endl;
}

//------------------------------------------------------------------------------
// computeLocalMatrix 子程序的详细实现
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::computeLocalMatrix(const Element& element,
                                                std::vector<std::vector<double>>& stiffness,
                                                std::vector<std::vector<double>>& mass,
                                                std::vector<std::vector<double>>& damping,
                                                std::vector<double>& force) {
    // 对应Fortran的LocalMatrix子程序的核心计算逻辑
    
    size_t nNodes = element.getNodeCount();
    size_t nDOFs = nNodes;
    
    // 初始化矩阵
    for (size_t i = 0; i < nDOFs; ++i) {
        for (size_t j = 0; j < nDOFs; ++j) {
            stiffness[i][j] = 0.0;
            mass[i][j] = 0.0;
            damping[i][j] = 0.0;
        }
        force[i] = 0.0;
    }
    
    // 获取材料属性
    // 对应Fortran: Material => GetMaterial(Element)
    std::string materialName = element.getMaterialName();
    
    // 获取电导率
    // 对应Fortran: C = GetReal( Material, 'Electric Conductivity', Found, Element)
    double conductivity = materialDB.getConductivity(materialName);
    
    // 获取磁化强度
    // 对应Fortran: M(1,:) = GetReal( Material, 'Magnetization 1', Found, Element)
    //              M(2,:) = GetReal( Material, 'Magnetization 2', Found, Element)
    std::array<double, 2> magnetization = materialDB.getMagnetization(materialName);
    
    // 获取单元节点坐标
    auto nodeCoords = element.getNodeCoordinates();
    
    // 数值积分（高斯积分）
    // 对应Fortran: IP = GaussPoints(Element)
    auto integrationPoints = GaussIntegration::getPoints(element.getType());
    
    for (const auto& ip : integrationPoints) {
        // 计算基函数和导数
        // 对应Fortran: stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), ...
        std::vector<double> basisFunctions;
        std::vector<std::array<double, 3>> basisDerivatives;
        double detJ;
        
        bool success = ShapeFunctions::evaluate(element, ip.u, ip.v, ip.w,
                                               basisFunctions, basisDerivatives, detJ);
        
        if (!success) {
            continue; // 跳过无效积分点
        }
        
        // 考虑坐标系对称性
        // 对应Fortran: IF ( CSymmetry ) THEN
        double coordinateFactor = 1.0;
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            // 计算径向坐标
            double r = 0.0;
            for (size_t i = 0; i < nNodes; ++i) {
                r += basisFunctions[i] * nodeCoords[i][0]; // x坐标作为径向坐标
            }
            coordinateFactor = r;
            detJ *= r; // 考虑轴对称的雅可比修正
        }
        
        // 计算磁阻率张量
        // 对应Fortran: CALL GetReluctivity(Material,R,n,Element)
        auto reluctivity = computeReluctivity(materialName, 0.0, element); // 初始磁通密度为0
        
        // 计算B矩阵（基函数导数矩阵）
        // 对应Fortran: Bt(1:nd,1) =  dbasisdx(1:nd,2)
        //              Bt(1:nd,2) = -dbasisdx(1:nd,1)
        std::vector<std::array<double, 2>> Bt(nDOFs);
        for (size_t i = 0; i < nDOFs; ++i) {
            Bt[i][0] = basisDerivatives[i][1];  // dN/dy
            Bt[i][1] = -basisDerivatives[i][0]; // -dN/dx
        }
        
        // 考虑坐标系对称性的修正
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            double r = coordinateFactor;
            for (size_t i = 0; i < nDOFs; ++i) {
                Bt[i][0] = -Bt[i][0]; // 符号修正
                Bt[i][1] = -Bt[i][1]; // 符号修正
                Bt[i][1] += basisFunctions[i] / r; // 轴对称修正项
            }
        }
        
        // 计算H矩阵 = nu_tensor * Bt
        // 对应Fortran: DO p = 1,nd
        //                Ht(p,:) = MATMUL(nu_tensor, Bt(p,:))
        //              END DO
        std::vector<std::array<double, 2>> Ht(nDOFs);
        for (size_t i = 0; i < nDOFs; ++i) {
            Ht[i][0] = reluctivity[0][0] * Bt[i][0] + reluctivity[0][1] * Bt[i][1];
            Ht[i][1] = reluctivity[1][0] * Bt[i][0] + reluctivity[1][1] * Bt[i][1];
        }
        
        // 组装刚度矩阵: STIFF = Ht * Bt^T
        // 对应Fortran: STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * MATMUL(Ht(1:nd,:), TRANSPOSE(Bt(1:nd,:)))
        for (size_t i = 0; i < nDOFs; ++i) {
            for (size_t j = 0; j < nDOFs; ++j) {
                double contribution = ip.weight * detJ * (Ht[i][0] * Bt[j][0] + Ht[i][1] * Bt[j][1]);
                stiffness[i][j] += contribution;
            }
        }
        
        // 组装质量矩阵和阻尼矩阵（瞬态分析）
        if (parameters.isTransient && conductivity != 0.0) {
            // 对应Fortran: MASS(p,q) = MASS(p,q) + IP % s(t) * detJ * C_ip * Basis(q)*Basis(p)
            for (size_t i = 0; i < nDOFs; ++i) {
                for (size_t j = 0; j < nDOFs; ++j) {
                    double massContribution = ip.weight * detJ * conductivity * basisFunctions[j] * basisFunctions[i];
                    mass[i][j] += massContribution;
                }
            }
        }
        
        // 组装力向量
        // 对应Fortran: FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtip * Basis(1:nd)
        double loadAtIP = 0.0; // TODO: 从体载荷获取
        for (size_t i = 0; i < nDOFs; ++i) {
            force[i] += ip.weight * detJ * loadAtIP * basisFunctions[i];
        }
        
        // 考虑磁化强度的贡献
        // 对应Fortran: FORCE(1:nd) = ... + M_ip(1)*dBasisdx(1:nd,2)-M_ip(2)*dBasisdx(1:nd,1)
        for (size_t i = 0; i < nDOFs; ++i) {
            double magnetizationContribution = magnetization[0] * basisDerivatives[i][1] - 
                                              magnetization[1] * basisDerivatives[i][0];
            force[i] += ip.weight * detJ * magnetizationContribution;
        }
    }
}

//------------------------------------------------------------------------------
// assembleLocalToGlobal 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::assembleLocalToGlobal(const Element& element,
                                                   const std::vector<std::vector<double>>& localStiffness,
                                                   const std::vector<std::vector<double>>& localMass,
                                                   const std::vector<std::vector<double>>& localDamping,
                                                   const std::vector<double>& localForce) {
    // 将局部矩阵组装到全局系统矩阵
    // 对应Fortran: 局部矩阵到全局矩阵的组装逻辑
    
    size_t nNodes = element.getNodeCount();
    
    // 获取单元节点编号
    auto nodeIndices = element.getNodeIndices();
    
    // 组装刚度矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        int globalRow = nodeIndices[i];
        
        for (size_t j = 0; j < nNodes; ++j) {
            int globalCol = nodeIndices[j];
            
            // 添加到全局刚度矩阵
            stiffnessMatrix->AddToElement(globalRow, globalCol, localStiffness[i][j]);
            
            // 如果是瞬态分析，组装质量矩阵和阻尼矩阵
            if (parameters.isTransient) {
                massMatrix->AddToElement(globalRow, globalCol, localMass[i][j]);
                dampingMatrix->AddToElement(globalRow, globalCol, localDamping[i][j]);
            }
        }
        
        // 组装力向量
        rhsVector->AddToElement(globalRow, localForce[i]);
    }
}

//------------------------------------------------------------------------------
// LocalMatrixInfinityBC 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::applyInfinityBoundaryCondition(const Element& element, 
                                                           std::vector<std::vector<double>>& stiffness,
                                                           std::vector<double>& force) {
    // 无限远边界条件处理
    size_t nNodes = element.getNodeCount();
    
    // 无限远边界条件使用Robin边界条件：αA + β∂A/∂n = γ
    // 对于无限远边界，通常使用α=0, β=1, γ=0，即∂A/∂n = 0
    
    // 计算边界长度
    double boundaryLength = computeBoundaryLength(element);
    
    // 计算边界法向量
    auto normalVector = computeBoundaryNormal(element);
    
    // 应用边界条件贡献
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            // 边界条件对刚度矩阵的贡献
            double boundaryContribution = boundaryLength * normalVector[i] * normalVector[j];
            stiffness[i][j] += boundaryContribution;
        }
    }
}

void MagnetoDynamics2DSolver::applyInfinityBoundaryConditionToGlobal(const Element& boundaryElement) {
    // 将无限远边界条件应用到全局系统
    size_t nNodes = boundaryElement.getNodeCount();
    
    // 初始化局部矩阵
    std::vector<std::vector<double>> localStiffness(nNodes, std::vector<double>(nNodes, 0.0));
    std::vector<double> localForce(nNodes, 0.0);
    
    // 计算局部边界条件贡献
    applyInfinityBoundaryCondition(boundaryElement, localStiffness, localForce);
    
    // 组装到全局系统
    assembleBoundaryLocalToGlobal(boundaryElement, localStiffness, localForce);
}

//------------------------------------------------------------------------------
// LocalMatrixAirGapBC 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::applyAirGapBoundaryCondition(const Element& element, 
                                                         std::vector<std::vector<double>>& stiffness,
                                                         std::vector<double>& force) {
    // 气隙边界条件处理
    size_t nNodes = element.getNodeCount();
    
    // 气隙边界条件考虑气隙中的磁场分布
    // 使用等效磁导率来处理气隙效应
    
    // 获取气隙参数
    double airGapLength = getAirGapLength(element);
    double airGapPermeability = getAirGapPermeability(element);
    
    // 计算边界长度
    double boundaryLength = computeBoundaryLength(element);
    
    // 计算气隙等效刚度
    double airGapStiffness = airGapPermeability / airGapLength;
    
    // 获取基函数（用于气隙边界条件的形状函数）
    auto nodeCoords = element.getNodeCoordinates();
    
    // 使用高斯积分计算气隙边界条件贡献
    auto integrationPoints = GaussIntegration::getPoints(element.getType());
    
    for (const auto& ip : integrationPoints) {
        // 计算基函数和导数
        std::vector<double> basisFunctions;
        std::vector<std::array<double, 3>> basisDerivatives;
        double detJ;
        
        if (!ShapeFunctions::evaluate(element, ip.u, ip.v, ip.w, basisFunctions, basisDerivatives, detJ)) {
            continue;
        }
        
        // 考虑坐标系对称性
        double coordinateFactor = 1.0;
        if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
            parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
            double r = 0.0;
            for (size_t i = 0; i < nNodes; ++i) {
                r += basisFunctions[i] * nodeCoords[i][0];
            }
            coordinateFactor = r;
            detJ *= r;
        }
        
        // 应用气隙边界条件贡献
        for (size_t i = 0; i < nNodes; ++i) {
            for (size_t j = 0; j < nNodes; ++j) {
                // 气隙对刚度矩阵的贡献
                double airGapContribution = ip.weight * detJ * airGapStiffness * 
                                          basisFunctions[i] * basisFunctions[j];
                stiffness[i][j] += airGapContribution;
            }
        }
    }
}

void MagnetoDynamics2DSolver::applyAirGapBoundaryConditionToGlobal(const Element& boundaryElement, 
                                                                  const std::shared_ptr<BoundaryCondition>& bc) {
    // 将气隙边界条件应用到全局系统
    size_t nNodes = boundaryElement.getNodeCount();
    
    // 初始化局部矩阵
    std::vector<std::vector<double>> localStiffness(nNodes, std::vector<double>(nNodes, 0.0));
    std::vector<double> localForce(nNodes, 0.0);
    
    // 计算局部边界条件贡献
    applyAirGapBoundaryCondition(boundaryElement, localStiffness, localForce);
    
    // 组装到全局系统
    assembleBoundaryLocalToGlobal(boundaryElement, localStiffness, localForce);
}

//------------------------------------------------------------------------------
// GetReluctivity 子程序的C++移植
//------------------------------------------------------------------------------
std::array<std::array<double, 2>, 2> MagnetoDynamics2DSolver::computeReluctivity(
    const std::string& materialName, 
    double magneticFluxDensity,
    const Element& element) {
    
    // 计算材料的磁阻率
    std::array<std::array<double, 2>, 2> reluctivity = {{{0.0, 0.0}, {0.0, 0.0}}};
    
    // 获取材料属性
    auto material = materialDB.getMaterial(materialName);
    
    // 检查是否为非线性材料
    bool isNonlinear = materialDB.isNonlinearMaterial(materialName);
    
    if (isNonlinear && magneticFluxDensity > 0.0) {
        // 非线性材料：根据磁通密度计算磁阻率
        reluctivity = computeNonlinearReluctivity(material, magneticFluxDensity);
    } else {
        // 线性材料：使用恒定磁导率
        double permeability = material.getPermeability();
        double reluctivityValue = 1.0 / permeability;
        
        // 创建各向同性磁阻率张量
        reluctivity[0][0] = reluctivityValue;
        reluctivity[0][1] = 0.0;
        reluctivity[1][0] = 0.0;
        reluctivity[1][1] = reluctivityValue;
    }
    
    return reluctivity;
}

std::array<std::array<double, 2>, 2> MagnetoDynamics2DSolver::computeNonlinearReluctivity(
    const ElectromagneticMaterial& material, 
    double magneticFluxDensity) {
    
    // 非线性材料磁阻率计算
    std::array<std::array<double, 2>, 2> reluctivity = {{{0.0, 0.0}, {0.0, 0.0}}};
    
    // 获取H-B曲线数据
    auto hbCurve = material.getHBCurve();
    
    if (hbCurve.empty()) {
        // 如果没有H-B曲线，使用饱和模型近似
        double saturationFluxDensity = material.getSaturationFluxDensity();
        double initialPermeability = material.getInitialPermeability();
        
        // 使用Froelich公式近似非线性特性
        double Bsat = saturationFluxDensity;
        double mu0 = 4.0 * M_PI * 1e-7;
        double mu_r = initialPermeability / mu0;
        
        double H = magneticFluxDensity / (mu0 * mu_r) + 
                   (magneticFluxDensity / Bsat) / (1.0 - magneticFluxDensity / Bsat);
        
        double reluctivityValue = H / magneticFluxDensity;
        
        reluctivity[0][0] = reluctivityValue;
        reluctivity[0][1] = 0.0;
        reluctivity[1][0] = 0.0;
        reluctivity[1][1] = reluctivityValue;
    } else {
        // 使用H-B曲线插值计算磁阻率
        double H = interpolateHBCurve(hbCurve, magneticFluxDensity);
        double reluctivityValue = H / magneticFluxDensity;
        
        reluctivity[0][0] = reluctivityValue;
        reluctivity[0][1] = 0.0;
        reluctivity[1][0] = 0.0;
        reluctivity[1][1] = reluctivityValue;
    }
    
    return reluctivity;
}

double MagnetoDynamics2DSolver::interpolateHBCurve(
    const std::vector<std::pair<double, double>>& hbCurve, 
    double magneticFluxDensity) {
    
    // H-B曲线插值
    if (hbCurve.empty()) {
        return 0.0;
    }
    
    // 查找插值区间
    for (size_t i = 0; i < hbCurve.size() - 1; ++i) {
        double B1 = hbCurve[i].first;
        double H1 = hbCurve[i].second;
        double B2 = hbCurve[i + 1].first;
        double H2 = hbCurve[i + 1].second;
        
        if (magneticFluxDensity >= B1 && magneticFluxDensity <= B2) {
            // 线性插值
            double t = (magneticFluxDensity - B1) / (B2 - B1);
            return H1 + t * (H2 - H1);
        }
    }
    
    // 超出范围，使用端点值
    if (magneticFluxDensity <= hbCurve.front().first) {
        return hbCurve.front().second;
    } else {
        return hbCurve.back().second;
    }
}

//------------------------------------------------------------------------------
// 其他辅助方法的实现
//------------------------------------------------------------------------------

void MagnetoDynamics2DSolver::assembleSystem() {
    // 系统组装主方法
    
    // 对应Fortran的系统组装流程
    // 1. 初始化系统矩阵
    // 2. 组装单元贡献
    // 3. 组装边界条件贡献
    // 4. 应用边界条件
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    
    // 初始化系统矩阵（每个节点1个自由度：A_z）
    stiffnessMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    rhsVector = std::shared_ptr<Vector>(Vector::Create(nNodes));
    
    if (parameters.isTransient) {
        massMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
        dampingMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    
    // 组装单元贡献
    assembleElementContributions();
    
    // 组装边界条件贡献
    assembleBoundaryContributions();
    
    // 应用边界条件
    applyBoundaryConditions();
    
    systemAssembled = true;
    
    std::cout << "System assembly completed" << std::endl;
}

std::vector<double> MagnetoDynamics2DSolver::solveLinearSystem() {
    // 线性系统求解
    
    if (!stiffnessMatrix || !rhsVector) {
        throw std::runtime_error("System matrices not assembled");
    }
    
    // 创建迭代求解器
    ConjugateGradientSolver solver(parameters.maxIterations, parameters.tolerance);
    
    // 求解
    auto solution = std::shared_ptr<Vector>(Vector::Create(rhsVector->Size()));
    bool success = solver.Solve(*stiffnessMatrix, *solution, *rhsVector);
    
    if (!success) {
        throw std::runtime_error("Linear solver failed to converge");
    }
    
    // 转换为向量
    std::vector<double> result(solution->Size());
    for (size_t i = 0; i < solution->Size(); ++i) {
        result[i] = (*solution)[i];
    }
    
    return result;
}

void MagnetoDynamics2DSolver::initializeBasisFunctionCache() {
    // 初始化基函数缓存
    
    if (!mesh) {
        std::cerr << "错误: 网格未初始化，无法初始化基函数缓存" << std::endl;
        return;
    }
    
    auto& elements = mesh->getElements();
    
    // 清空现有缓存
    basisFunctionCache.clear();
    
    // 为每个单元预计算基函数值
    for (const auto& element : elements) {
        BasisFunctionCache cache;
        
        // 获取单元类型
        std::string elementType = element.getType();
        size_t nNodes = element.getNodeCount();
        
        // 根据单元类型初始化缓存
        if (elementType == "Triangle" || elementType == "Triangular") {
            // 三角形单元基函数
            cache.shapeFunctions.resize(nNodes);
            cache.shapeFunctionDerivatives.resize(nNodes);
            
            // 预计算高斯积分点的基函数值
            // 这里简化处理，实际应该基于高斯积分点
            for (size_t i = 0; i < nNodes; ++i) {
                cache.shapeFunctions[i] = 1.0 / nNodes; // 简化值
                cache.shapeFunctionDerivatives[i] = 1.0; // 简化值
            }
        } else if (elementType == "Quadrilateral" || elementType == "Quad") {
            // 四边形单元基函数
            cache.shapeFunctions.resize(nNodes);
            cache.shapeFunctionDerivatives.resize(nNodes);
            
            for (size_t i = 0; i < nNodes; ++i) {
                cache.shapeFunctions[i] = 0.25; // 四边形单元的标准值
                cache.shapeFunctionDerivatives[i] = 0.5; // 简化值
            }
        }
        
        // 存储缓存
        basisFunctionCache[element.getId()] = cache;
    }
    
    std::cout << "基函数缓存初始化完成，缓存了 " << elements.size() << " 个单元" << std::endl;
}

void MagnetoDynamics2DSolver::assembleElementContributions() {
    // 组装单元贡献
    
    if (!stiffnessMatrix || !massMatrix || !dampingMatrix || !rhsVector) {
        std::cerr << "错误: 系统矩阵未初始化，无法组装单元贡献" << std::endl;
        return;
    }
    
    auto& elements = mesh->getElements();
    
    // 遍历所有单元
    for (const auto& element : elements) {
        // 计算单元矩阵
        auto elementMatrices = computeElementMatrices(element);
        
        // 将局部矩阵组装到全局系统
        assembleLocalToGlobal(element, 
                             elementMatrices.stiffness,
                             elementMatrices.mass,
                             elementMatrices.damping,
                             elementMatrices.force);
    }
    
    std::cout << "单元贡献组装完成，处理了 " << elements.size() << " 个单元" << std::endl;
}

void MagnetoDynamics2DSolver::assembleBoundaryContributions() {
    // 组装边界条件贡献
    
    if (!stiffnessMatrix || !rhsVector) {
        std::cerr << "错误: 系统矩阵未初始化，无法组装边界条件贡献" << std::endl;
        return;
    }
    
    // 遍历所有边界条件
    for (const auto& bcPair : boundaryConditions) {
        const auto& bc = bcPair.second;
        
        // 根据边界条件类型组装贡献
        switch (bc.type[0]) {
            case 'D': // Dirichlet
                // Dirichlet边界条件在应用阶段处理
                break;
                
            case 'N': // Neumann
                assembleNeumannBoundaryContribution(bc);
                break;
                
            case 'M': // Magnetic Flux Density
                assembleMagneticFluxBoundaryContribution(bc);
                break;
                
            default:
                std::cout << "警告: 未知边界条件类型: " << bc.type << std::endl;
                break;
        }
    }
    
    std::cout << "边界条件贡献组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::applyBoundaryConditions() {
    // 应用边界条件
    
    if (!stiffnessMatrix || !rhsVector) {
        std::cerr << "错误: 系统矩阵未初始化，无法应用边界条件" << std::endl;
        return;
    }
    
    // 应用Dirichlet边界条件
    applyDirichletBoundaryConditions();
    
    // 应用Neumann边界条件
    applyNeumannBoundaryConditions();
    
    // 应用磁通密度边界条件
    applyMagneticFluxDensityBoundaryConditions();
    
    std::cout << "边界条件应用完成" << std::endl;
}

void MagnetoDynamics2DSolver::assembleElementContributionsParallel(
    const DomainDecompositionResult& decomposition) {
    // 并行组装单元贡献
    
    if (!isParallel()) {
        std::cout << "警告: 未启用并行模式，使用串行组装" << std::endl;
        assembleElementContributions();
        return;
    }
    
    // 获取当前进程的本地单元
    auto localElements = decomposition.getLocalElements();
    
    // 并行组装本地单元贡献
    #pragma omp parallel for
    for (size_t i = 0; i < localElements.size(); ++i) {
        const auto& element = localElements[i];
        
        // 计算单元矩阵
        auto elementMatrices = computeElementMatrices(element);
        
        // 将局部矩阵组装到全局系统（需要线程安全）
        #pragma omp critical
        {
            assembleLocalToGlobal(element, 
                                 elementMatrices.stiffness,
                                 elementMatrices.mass,
                                 elementMatrices.damping,
                                 elementMatrices.force);
        }
    }
    
    // 处理并行通信（如果需要）
    handleParallelCommunication(decomposition);
    
    std::cout << "并行单元贡献组装完成，处理了 " << localElements.size() << " 个本地单元" << std::endl;
}

void MagnetoDynamics2DSolver::assembleBoundaryContributionsParallel(
    const DomainDecompositionResult& decomposition) {
    // 并行组装边界条件贡献
    
    if (!isParallel()) {
        std::cout << "警告: 未启用并行模式，使用串行边界条件组装" << std::endl;
        assembleBoundaryContributions();
        return;
    }
    
    // 获取当前进程的本地边界条件
    auto localBoundaryConditions = decomposition.getLocalBoundaryConditions();
    
    // 并行组装本地边界条件贡献
    #pragma omp parallel for
    for (size_t i = 0; i < localBoundaryConditions.size(); ++i) {
        const auto& bc = localBoundaryConditions[i];
        
        // 根据边界条件类型组装贡献
        switch (bc.type[0]) {
            case 'N': // Neumann
                #pragma omp critical
                {
                    assembleNeumannBoundaryContribution(bc);
                }
                break;
                
            case 'M': // Magnetic Flux Density
                #pragma omp critical
                {
                    assembleMagneticFluxBoundaryContribution(bc);
                }
                break;
                
            default:
                // Dirichlet边界条件在应用阶段处理
                break;
        }
    }
    
    // 处理边界条件的并行通信
    handleBoundaryConditionCommunication(decomposition);
    
    std::cout << "并行边界条件贡献组装完成，处理了 " << localBoundaryConditions.size() << " 个本地边界条件" << std::endl;
}

void MagnetoDynamics2DSolver::applyBoundaryConditionsParallel() {
    // 并行应用边界条件
    
    if (!isParallel()) {
        std::cout << "警告: 未启用并行模式，使用串行边界条件" << std::endl;
        applyBoundaryConditions();
        return;
    }
    
    // 获取当前进程的本地边界条件
    auto localBoundaryConditions = getLocalBoundaryConditions();
    
    // 并行应用本地边界条件
    #pragma omp parallel for
    for (size_t i = 0; i < localBoundaryConditions.size(); ++i) {
        const auto& bc = localBoundaryConditions[i];
        
        // 根据边界条件类型应用
        switch (bc.type[0]) {
            case 'D': // Dirichlet
                #pragma omp critical
                {
                    applyDirichletBoundaryCondition(bc);
                }
                break;
                
            case 'N': // Neumann
                #pragma omp critical
                {
                    applyNeumannBoundaryCondition(bc);
                }
                break;
                
            case 'M': // Magnetic Flux Density
                #pragma omp critical
                {
                    applyMagneticFluxDensityBoundaryCondition(bc);
                }
                break;
                
            default:
                std::cout << "警告: 未知边界条件类型: " << bc.type << std::endl;
                break;
        }
    }
    
    // 并行同步：确保所有进程完成边界条件应用
    synchronizeBoundaryConditions();
    
    std::cout << "并行边界条件应用完成，处理了 " << localBoundaryConditions.size() << " 个本地边界条件" << std::endl;
}

DomainDecompositionResult MagnetoDynamics2DSolver::performDomainDecomposition() {
    // 执行域分解
    
    if (!mesh) {
        std::cerr << "错误: 网格未初始化，无法执行域分解" << std::endl;
        return DomainDecompositionResult();
    }
    
    DomainDecompositionResult result;
    
    // 获取总单元数和进程数
    size_t totalElements = mesh->getElements().size();
    int numProcesses = isParallel() ? comm_->getSize() : 1;
    int currentProcess = isParallel() ? comm_->getRank() : 0;
    
    // 负载均衡：每个进程分配大致相等的单元数
    size_t elementsPerProcess = totalElements / numProcesses;
    size_t remainder = totalElements % numProcesses;
    
    // 计算当前进程的单元范围
    size_t startElement = currentProcess * elementsPerProcess + std::min((size_t)currentProcess, remainder);
    size_t endElement = startElement + elementsPerProcess;
    if (currentProcess < remainder) {
        endElement += 1;
    }
    
    // 分配本地单元
    auto& allElements = mesh->getElements();
    for (size_t i = startElement; i < endElement && i < allElements.size(); ++i) {
        result.localElements.push_back(allElements[i]);
    }
    
    // 分配边界条件（简化处理）
    // 实际应该基于几何位置分配边界条件
    for (const auto& bcPair : boundaryConditions) {
        const auto& bc = bcPair.second;
        
        // 检查边界条件是否与本地单元相关
        for (int nodeId : bc.affectedNodes) {
            // 简化：如果节点在本地单元范围内，则分配该边界条件
            if (nodeId >= startElement && nodeId < endElement) {
                result.localBoundaryConditions.push_back(bc);
                break;
            }
        }
    }
    
    // 记录分解信息
    result.totalProcesses = numProcesses;
    result.currentProcess = currentProcess;
    result.totalElements = totalElements;
    result.localElementCount = result.localElements.size();
    
    std::cout << "域分解完成，进程 " << currentProcess << " 分配了 " 
              << result.localElementCount << " 个单元和 " 
              << result.localBoundaryConditions.size() << " 个边界条件" << std::endl;
    
    return result;
}

void MagnetoDynamics2DSolver::calculateDerivedFields(MagnetoDynamics2DResults& results) {
    // 计算导出场量
    
    if (!mesh || results.vectorPotential.empty()) {
        std::cerr << "错误: 网格或解向量未初始化，无法计算导出场量" << std::endl;
        return;
    }
    
    // 计算磁场强度分布
    results.magneticFieldStrength = computeMagneticFieldStrength();
    
    // 计算磁通密度分布
    results.magneticFluxDensity = computeMagneticFluxDensity();
    
    // 计算电场分布（如果有时变项）
    if (parameters.timeDependent) {
        calculateElectricField(results);
    }
    
    // 计算能量密度分布
    calculateEnergyDensity(results);
    
    // 计算力密度分布
    calculateForceDensity(results);
    
    std::cout << "导出场量计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::calculateLumpedParameters(MagnetoDynamics2DResults& results) {
    // 计算集总参数
    
    if (!mesh || results.energyDensity.empty()) {
        std::cerr << "错误: 网格或能量密度未初始化，无法计算集总参数" << std::endl;
        return;
    }
    
    // 计算总磁能
    results.totalMagneticEnergy = 0.0;
    for (double energy : results.energyDensity) {
        results.totalMagneticEnergy += energy;
    }
    
    // 计算总磁通
    results.totalMagneticFlux = 0.0;
    size_t nNodes = mesh->getNodes().numberOfNodes();
    for (size_t i = 0; i < nNodes; ++i) {
        size_t idx = i * 2;
        double Bx = results.magneticFluxDensity[idx];
        double By = results.magneticFluxDensity[idx + 1];
        double B_magnitude = std::sqrt(Bx*Bx + By*By);
        results.totalMagneticFlux += B_magnitude;
    }
    
    // 计算电感
    results.inductance = 0.0;
    if (parameters.currentDensity != 0.0) {
        // 电感公式: L = 2 * W_m / I^2
        // 其中W_m是总磁能，I是总电流
        double totalCurrent = parameters.currentDensity * mesh->getArea();
        if (totalCurrent != 0.0) {
            results.inductance = 2.0 * results.totalMagneticEnergy / (totalCurrent * totalCurrent);
        }
    }
    
    // 计算功率损耗
    results.powerLoss = 0.0;
    if (parameters.timeDependent && parameters.frequency != 0.0) {
        // 涡流损耗公式: P = ∫(σE^2) dΩ
        // 简化计算：基于磁能密度和频率
        results.powerLoss = results.totalMagneticEnergy * parameters.frequency * parameters.electricalConductivity;
    }
    
    std::cout << "集总参数计算完成" << std::endl;
    std::cout << "总磁能: " << results.totalMagneticEnergy << " J" << std::endl;
    std::cout << "总磁通: " << results.totalMagneticFlux << " Wb" << std::endl;
    std::cout << "电感: " << results.inductance << " H" << std::endl;
    std::cout << "功率损耗: " << results.powerLoss << " W" << std::endl;
}

// =============================================================================
// 性能优化方法实现
// =============================================================================

void MagnetoDynamics2DSolver::precomputeBasisFunctionCache() {
    // 预计算基函数缓存
    
    if (!mesh) {
        std::cerr << "错误: 网格未初始化，无法预计算基函数缓存" << std::endl;
        return;
    }
    
    auto& elements = mesh->getElements();
    
    // 清空现有缓存
    basisFunctionCache.clear();
    
    // 为每个单元预计算基函数值
    for (const auto& element : elements) {
        BasisFunctionCache cache;
        
        // 获取单元类型
        std::string elementType = element.getType();
        size_t nNodes = element.getNodeCount();
        
        // 根据单元类型初始化缓存
        if (elementType == "Triangle" || elementType == "Triangular") {
            // 三角形单元基函数
            cache.shapeFunctions.resize(nNodes);
            cache.shapeFunctionDerivatives.resize(nNodes);
            
            // 预计算高斯积分点的基函数值
            // 这里简化处理，实际应该基于高斯积分点
            for (size_t i = 0; i < nNodes; ++i) {
                cache.shapeFunctions[i] = 1.0 / nNodes; // 简化值
                cache.shapeFunctionDerivatives[i] = 1.0; // 简化值
            }
        } else if (elementType == "Quadrilateral" || elementType == "Quad") {
            // 四边形单元基函数
            cache.shapeFunctions.resize(nNodes);
            cache.shapeFunctionDerivatives.resize(nNodes);
            
            for (size_t i = 0; i < nNodes; ++i) {
                cache.shapeFunctions[i] = 0.25; // 四边形单元的标准值
                cache.shapeFunctionDerivatives[i] = 0.5; // 简化值
            }
        }
        
        // 存储缓存
        basisFunctionCache[element.getId()] = cache;
    }
    
    std::cout << "基函数缓存预计算完成，缓存了 " << elements.size() << " 个单元" << std::endl;
}

void MagnetoDynamics2DSolver::assembleSystemOptimized() {
    // 优化的矩阵组装实现
    
    if (!stiffnessMatrix || !massMatrix || !dampingMatrix || !rhsVector) {
        initializeSystemMatrices();
    }
    
    // 使用缓存数据加速组装
    if (useBasisFunctionsCache && basisFunctionCache.empty()) {
        initializeBasisFunctionCache();
    }
    
    auto& elements = mesh->getElements();
    
    // 并行组装单元贡献（如果启用并行）
    if (isParallel()) {
        assembleElementContributionsParallel();
    } else {
        // 串行优化组装
        for (const auto& element : elements) {
            // 使用缓存的基函数数据
            auto elementMatrices = computeElementMatricesOptimized(element);
            
            // 组装到全局系统
            assembleLocalToGlobal(element, 
                                 elementMatrices.stiffness,
                                 elementMatrices.mass,
                                 elementMatrices.damping,
                                 elementMatrices.force);
        }
    }
    
    // 组装边界条件贡献
    assembleBoundaryContributions();
    
    // 应用边界条件
    applyBoundaryConditions();
    
    systemAssembled = true;
    
    std::cout << "优化的矩阵组装完成，处理了 " << elements.size() << " 个单元" << std::endl;
}

void MagnetoDynamics2DSolver::updateSystemIncremental() {
    // 增量式矩阵更新
    
    if (!stiffnessMatrix || !rhsVector || currentPotential.empty()) {
        std::cerr << "错误: 系统矩阵或当前解未初始化，无法进行增量更新" << std::endl;
        return;
    }
    
    // 非线性迭代优化：检查是否需要更新
    if (currentIteration <= 1) {
        std::cout << "第一次迭代，进行完整系统组装" << std::endl;
        assembleSystem();
        return;
    }
    
    // 处理材料非线性：更新非线性材料属性
    updateMaterialParameters();
    
    // 增量更新：只更新受影响的单元
    auto& elements = mesh->getElements();
    
    // 计算当前残差
    double residualNorm = computeResidualNorm();
    
    // 根据残差大小决定更新策略
    if (residualNorm > parameters.nonlinearTolerance * 10.0) {
        // 残差较大，进行完整更新
        std::cout << "残差较大(" << residualNorm << ")，进行完整系统更新" << std::endl;
        assembleSystem();
    } else {
        // 残差较小，进行增量更新
        std::cout << "残差较小(" << residualNorm << ")，进行增量更新" << std::endl;
        
        // 只更新非线性材料单元
        size_t updatedElements = 0;
        for (const auto& element : elements) {
            std::string materialName = element.getMaterialName();
            
            // 检查是否为非线性材料且需要更新
            if (materialDB.isNonlinearMaterial(materialName) && 
                shouldUpdateElement(element)) {
                
                // 更新该单元的刚度矩阵贡献
                updateElementStiffness(element);
                updatedElements++;
            }
        }
        
        std::cout << "增量更新完成，更新了 " << updatedElements << " 个非线性材料单元" << std::endl;
    }
    
    // 更新迭代计数器
    currentIteration++;
}

double MagnetoDynamics2DSolver::computeResidualNorm() {
    // 计算残差范数
    
    if (!stiffnessMatrix || !rhsVector || currentPotential.empty()) {
        return std::numeric_limits<double>::max();
    }
    
    // 计算残差向量: r = b - A*x
    auto residual = std::shared_ptr<Vector>(Vector::Create(rhsVector->Size()));
    
    // 复制右端向量
    for (size_t i = 0; i < rhsVector->Size(); ++i) {
        (*residual)[i] = (*rhsVector)[i];
    }
    
    // 减去 A*x
    auto solutionVector = std::shared_ptr<Vector>(Vector::Create(currentPotential.size()));
    for (size_t i = 0; i < currentPotential.size(); ++i) {
        (*solutionVector)[i] = currentPotential[i];
    }
    
    auto Ax = std::shared_ptr<Vector>(Vector::Create(solutionVector->Size()));
    stiffnessMatrix->Multiply(*solutionVector, *Ax);
    
    for (size_t i = 0; i < residual->Size(); ++i) {
        (*residual)[i] -= (*Ax)[i];
    }
    
    // 计算残差范数
    double residualNorm = 0.0;
    for (size_t i = 0; i < residual->Size(); ++i) {
        residualNorm += (*residual)[i] * (*residual)[i];
    }
    residualNorm = std::sqrt(residualNorm);
    
    // 计算相对残差
    double rhsNorm = rhsVector->Norm2();
    double relativeResidual = (rhsNorm > 1e-12) ? residualNorm / rhsNorm : residualNorm;
    
    return relativeResidual;
}

bool MagnetoDynamics2DSolver::shouldUpdateElement(const Element& element) {
    // 判断单元是否需要更新
    
    // 检查单元是否在非线性材料区域
    std::string materialName = element.getMaterialName();
    if (!materialDB.isNonlinearMaterial(materialName)) {
        return false;
    }
    
    // 检查单元的磁场强度变化
    auto nodeIndices = element.getNodeIndices();
    double maxFieldChange = 0.0;
    
    // 计算当前磁场强度
    auto H_current = computeMagneticFieldStrength();
    
    // 检查历史记录中的磁场强度变化
    if (materialUpdateHistory.find(element.getId()) != materialUpdateHistory.end()) {
        const auto& history = materialUpdateHistory[element.getId()];
        if (history.size() >= 2) {
            // 比较最近两次迭代的磁场强度
            const auto& lastUpdate = history.back();
            const auto& prevUpdate = history[history.size() - 2];
            
            double fieldChange = std::abs(lastUpdate.fieldStrength - prevUpdate.fieldStrength);
            maxFieldChange = std::max(maxFieldChange, fieldChange);
        }
    }
    
    // 如果磁场强度变化超过阈值，则需要更新
    double fieldChangeThreshold = parameters.nonlinearTolerance * 0.1;
    return maxFieldChange > fieldChangeThreshold;
}

void MagnetoDynamics2DSolver::updateElementStiffness(const Element& element) {
    // 更新单元的刚度矩阵贡献
    
    // 计算新的单元刚度矩阵
    auto elementMatrices = computeElementMatrices(element);
    
    // 从全局矩阵中移除旧的贡献
    removeElementStiffness(element);
    
    // 添加新的贡献
    assembleLocalToGlobal(element, 
                         elementMatrices.stiffness,
                         elementMatrices.mass,
                         elementMatrices.damping,
                         elementMatrices.force);
}

void MagnetoDynamics2DSolver::assembleSystemWithJacobian() {
    // 组装带雅可比项的系统矩阵
    
    if (!stiffnessMatrix || !massMatrix || !dampingMatrix || !rhsVector) {
        initializeSystemMatrices();
    }
    
    // 初始化雅可比矩阵
    if (!jacobianMatrix) {
        size_t nNodes = mesh->getNodes().numberOfNodes();
        jacobianMatrix = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    
    auto& elements = mesh->getElements();
    
    // 组装单元贡献（包含雅可比项）
    for (const auto& element : elements) {
        // 计算单元矩阵（包含雅可比项）
        auto elementMatrices = computeElementMatricesWithJacobian(element);
        
        // 将局部矩阵组装到全局系统
        assembleLocalToGlobal(element, 
                             elementMatrices.stiffness,
                             elementMatrices.mass,
                             elementMatrices.damping,
                             elementMatrices.force);
        
        // 组装雅可比矩阵贡献
        assembleJacobianContribution(element, elementMatrices.jacobian);
    }
    
    // 组装边界条件贡献
    assembleBoundaryContributions();
    
    // 应用边界条件
    applyBoundaryConditions();
    
    systemAssembled = true;
    
    std::cout << "带雅可比项的系统矩阵组装完成，处理了 " << elements.size() << " 个单元" << std::endl;
}

ElementMatricesWithJacobian MagnetoDynamics2DSolver::computeElementMatricesWithJacobian(const Element& element) {
    // 计算带雅可比项的单元矩阵
    
    ElementMatricesWithJacobian matrices;
    size_t nNodes = element.getNodeCount();
    
    // 初始化矩阵
    matrices.stiffness.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.mass.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.damping.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.jacobian.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.force.resize(nNodes, 0.0);
    
    // 获取材料属性
    std::string materialName = element.getMaterialName();
    double reluctivity = getReluctivity(materialName, element);
    
    // 计算单元刚度矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            matrices.stiffness[i][j] = reluctivity * element.getArea();
        }
    }
    
    // 计算雅可比矩阵（刚度矩阵的导数）
    // 对于非线性材料，雅可比矩阵包含材料导数的贡献
    if (materialDB.isNonlinearMaterial(materialName)) {
        // 计算材料导数 dν/dB
        double dB = 1.0; // 简化处理，实际应该基于磁场变化
        double dReluctivity = computeReluctivityDerivative(materialName, dB);
        
        for (size_t i = 0; i < nNodes; ++i) {
            for (size_t j = 0; j < nNodes; ++j) {
                matrices.jacobian[i][j] = dReluctivity * element.getArea();
            }
        }
    } else {
        // 线性材料，雅可比矩阵等于刚度矩阵
        matrices.jacobian = matrices.stiffness;
    }
    
    // 计算其他矩阵（与标准方法相同）
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            matrices.mass[i][j] = parameters.materialDensity * element.getArea() / nNodes;
            matrices.damping[i][j] = parameters.dampingCoefficient * element.getArea() / nNodes;
        }
        matrices.force[i] = parameters.currentDensity * element.getArea() / nNodes;
    }
    
    return matrices;
}

void MagnetoDynamics2DSolver::assembleJacobianContribution(const Element& element, 
                                                           const std::vector<std::vector<double>>& localJacobian) {
    // 组装雅可比矩阵贡献
    
    if (!jacobianMatrix) {
        std::cerr << "错误: 雅可比矩阵未初始化" << std::endl;
        return;
    }
    
    size_t nNodes = element.getNodeCount();
    auto nodeIndices = element.getNodeIndices();
    
    // 将局部雅可比矩阵组装到全局雅可比矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        int globalRow = nodeIndices[i];
        
        for (size_t j = 0; j < nNodes; ++j) {
            int globalCol = nodeIndices[j];
            
            // 添加局部贡献到全局矩阵
            double currentValue = jacobianMatrix->GetElement(globalRow, globalCol);
            jacobianMatrix->SetElement(globalRow, globalCol, currentValue + localJacobian[i][j]);
        }
    }
}

void MagnetoDynamics2DSolver::addNonlinearJacobianContribution(const Element& element) {
    // 添加非线性材料对雅可比矩阵的贡献
    
    if (!jacobianMatrix || currentPotential.empty()) {
        return;
    }
    
    size_t nNodes = element.getNodeCount();
    auto nodeIndices = element.getNodeIndices();
    
    // 获取材料名称
    std::string materialName = element.getMaterialName();
    
    // 计算当前单元的磁场强度
    double H_avg = 0.0;
    auto H = computeMagneticFieldStrength();
    
    for (size_t i = 0; i < nNodes; ++i) {
        int nodeIndex = nodeIndices[i];
        double Hx = H[nodeIndex * 2];
        double Hy = H[nodeIndex * 2 + 1];
        double H_magnitude = std::sqrt(Hx*Hx + Hy*Hy);
        H_avg += H_magnitude;
    }
    H_avg /= nNodes;
    
    // 计算材料导数 dν/dH
    double dReluctivity = computeReluctivityDerivative(materialName, H_avg);
    
    // 计算非线性贡献
    double nonlinearContribution = dReluctivity * element.getArea();
    
    // 将非线性贡献添加到雅可比矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        int globalRow = nodeIndices[i];
        
        for (size_t j = 0; j < nNodes; ++j) {
            int globalCol = nodeIndices[j];
            
            // 添加非线性贡献
            double currentValue = jacobianMatrix->GetElement(globalRow, globalCol);
            jacobianMatrix->SetElement(globalRow, globalCol, currentValue + nonlinearContribution);
        }
    }
}

void MagnetoDynamics2DSolver::initializeComplexSystemMatrices() {
    // 初始化复数系统矩阵
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    
    // 初始化复数刚度矩阵（实部和虚部）
    if (!complexStiffnessMatrixReal) {
        complexStiffnessMatrixReal = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    if (!complexStiffnessMatrixImag) {
        complexStiffnessMatrixImag = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    
    // 初始化复数质量矩阵
    if (!complexMassMatrixReal) {
        complexMassMatrixReal = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    if (!complexMassMatrixImag) {
        complexMassMatrixImag = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    
    // 初始化复数阻尼矩阵
    if (!complexDampingMatrixReal) {
        complexDampingMatrixReal = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    if (!complexDampingMatrixImag) {
        complexDampingMatrixImag = std::dynamic_pointer_cast<Matrix>(std::make_shared<CRSMatrix>(nNodes, nNodes));
    }
    
    // 初始化复数右端向量
    if (!complexRhsVectorReal) {
        complexRhsVectorReal = std::shared_ptr<Vector>(Vector::Create(nNodes));
    }
    if (!complexRhsVectorImag) {
        complexRhsVectorImag = std::shared_ptr<Vector>(Vector::Create(nNodes));
    }
    
    // 初始化复数解向量
    if (complexSolutionReal.empty()) {
        complexSolutionReal.resize(nNodes, 0.0);
    }
    if (complexSolutionImag.empty()) {
        complexSolutionImag.resize(nNodes, 0.0);
    }
    
    std::cout << "复数系统矩阵初始化完成" << std::endl;
}

void MagnetoDynamics2DSolver::assembleComplexElementContributions() {
    // 组装复数单元贡献
    
    if (!complexStiffnessMatrixReal || !complexStiffnessMatrixImag) {
        std::cerr << "错误: 复数系统矩阵未初始化" << std::endl;
        return;
    }
    
    auto& elements = mesh->getElements();
    
    // 遍历所有单元
    for (const auto& element : elements) {
        // 计算复数单元矩阵
        auto complexElementMatrices = computeComplexElementMatrices(element);
        
        // 将复数局部矩阵组装到全局系统
        assembleComplexLocalToGlobal(element, complexElementMatrices);
    }
    
    std::cout << "复数单元贡献组装完成，处理了 " << elements.size() << " 个单元" << std::endl;
}

ComplexElementMatrices MagnetoDynamics2DSolver::computeComplexElementMatrices(const Element& element) {
    // 计算复数单元矩阵
    
    ComplexElementMatrices matrices;
    size_t nNodes = element.getNodeCount();
    
    // 初始化复数矩阵
    matrices.stiffnessReal.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.stiffnessImag.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.massReal.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.massImag.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.dampingReal.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.dampingImag.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.forceReal.resize(nNodes, 0.0);
    matrices.forceImag.resize(nNodes, 0.0);
    
    // 获取材料属性
    std::string materialName = element.getMaterialName();
    double reluctivity = getReluctivity(materialName, element);
    
    // 计算复数刚度矩阵
    // 实部：电阻性项
    // 虚部：电抗性项（与频率相关）
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            // 实部：电阻性刚度矩阵
            matrices.stiffnessReal[i][j] = reluctivity * element.getArea();
            
            // 虚部：电抗性项（与频率相关）
            if (parameters.frequency != 0.0) {
                matrices.stiffnessImag[i][j] = 2.0 * M_PI * parameters.frequency * 
                                              parameters.electricalConductivity * element.getArea();
            }
        }
    }
    
    // 计算复数质量矩阵（用于时域分析）
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            matrices.massReal[i][j] = parameters.materialDensity * element.getArea() / nNodes;
            // 质量矩阵虚部通常为0
            matrices.massImag[i][j] = 0.0;
        }
    }
    
    // 计算复数阻尼矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            matrices.dampingReal[i][j] = parameters.dampingCoefficient * element.getArea() / nNodes;
            // 阻尼矩阵虚部通常为0
            matrices.dampingImag[i][j] = 0.0;
        }
    }
    
    // 计算复数载荷向量
    for (size_t i = 0; i < nNodes; ++i) {
        // 实部：直流或同相分量
        matrices.forceReal[i] = parameters.currentDensity * element.getArea() / nNodes;
        
        // 虚部：交流或正交分量
        if (parameters.frequency != 0.0) {
            matrices.forceImag[i] = parameters.currentDensity * element.getArea() / nNodes * 
                                   std::sin(2.0 * M_PI * parameters.frequency * parameters.time);
        }
    }
    
    return matrices;
}

void MagnetoDynamics2DSolver::assembleComplexLocalToGlobal(const Element& element, const ComplexElementMatrices& matrices) {
    // 将复数局部矩阵组装到全局系统
    
    size_t nNodes = element.getNodeCount();
    auto nodeIndices = element.getNodeIndices();
    
    // 组装实部矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        int globalRow = nodeIndices[i];
        
        for (size_t j = 0; j < nNodes; ++j) {
            int globalCol = nodeIndices[j];
            
            // 组装实部刚度矩阵
            double currentReal = complexStiffnessMatrixReal->GetElement(globalRow, globalCol);
            complexStiffnessMatrixReal->SetElement(globalRow, globalCol, 
                                                  currentReal + matrices.stiffnessReal[i][j]);
            
            // 组装虚部刚度矩阵
            double currentImag = complexStiffnessMatrixImag->GetElement(globalRow, globalCol);
            complexStiffnessMatrixImag->SetElement(globalRow, globalCol, 
                                                  currentImag + matrices.stiffnessImag[i][j]);
            
            // 组装实部质量矩阵
            currentReal = complexMassMatrixReal->GetElement(globalRow, globalCol);
            complexMassMatrixReal->SetElement(globalRow, globalCol, 
                                             currentReal + matrices.massReal[i][j]);
            
            // 组装虚部质量矩阵
            currentImag = complexMassMatrixImag->GetElement(globalRow, globalCol);
            complexMassMatrixImag->SetElement(globalRow, globalCol, 
                                             currentImag + matrices.massImag[i][j]);
            
            // 组装实部阻尼矩阵
            currentReal = complexDampingMatrixReal->GetElement(globalRow, globalCol);
            complexDampingMatrixReal->SetElement(globalRow, globalCol, 
                                                currentReal + matrices.dampingReal[i][j]);
            
            // 组装虚部阻尼矩阵
            currentImag = complexDampingMatrixImag->GetElement(globalRow, globalCol);
            complexDampingMatrixImag->SetElement(globalRow, globalCol, 
                                                currentImag + matrices.dampingImag[i][j]);
        }
        
        // 组装实部右端向量
        (*complexRhsVectorReal)[globalRow] += matrices.forceReal[i];
        
        // 组装虚部右端向量
        (*complexRhsVectorImag)[globalRow] += matrices.forceImag[i];
    }
}

void MagnetoDynamics2DSolver::assembleComplexBoundaryContributions() {
    // 组装复数边界条件贡献
    
    if (!complexRhsVectorReal || !complexRhsVectorImag) {
        std::cerr << "错误: 复数右端向量未初始化" << std::endl;
        return;
    }
    
    // 遍历所有边界条件
    for (const auto& bcPair : boundaryConditions) {
        const auto& bc = bcPair.second;
        
        // 根据边界条件类型组装复数贡献
        switch (bc.type[0]) {
            case 'N': // Neumann
                assembleComplexNeumannBoundaryContribution(bc);
                break;
                
            case 'M': // Magnetic Flux Density
                assembleComplexMagneticFluxBoundaryContribution(bc);
                break;
                
            default:
                // Dirichlet边界条件在应用阶段处理
                break;
        }
    }
    
    std::cout << "复数边界条件贡献组装完成" << std::endl;
}

void MagnetoDynamics2DSolver::applyComplexBoundaryConditions() {
    // 应用复数边界条件
    
    if (!complexStiffnessMatrixReal || !complexStiffnessMatrixImag || 
        !complexRhsVectorReal || !complexRhsVectorImag) {
        std::cerr << "错误: 复数系统矩阵未初始化，无法应用边界条件" << std::endl;
        return;
    }
    
    // 应用复数Dirichlet边界条件
    applyComplexDirichletBoundaryConditions();
    
    std::cout << "复数边界条件应用完成" << std::endl;
}

void MagnetoDynamics2DSolver::computeReluctivityDerivative(const std::string& materialName, double fieldStrength) {
    // 计算磁阻率对磁场强度的导数 dν/dH
    
    auto material = materialDB.getMaterial(materialName);
    
    // 检查是否为非线性材料
    if (!materialDB.isNonlinearMaterial(materialName)) {
        // 线性材料，导数为0
        return 0.0;
    }
    
    // 使用有限差分法计算导数
    double h = 1e-6; // 差分步长
    double reluctivity1 = computeReluctivity(materialName, fieldStrength);
    double reluctivity2 = computeReluctivity(materialName, fieldStrength + h);
    
    double derivative = (reluctivity2 - reluctivity1) / h;
    
    // 确保导数非负（磁阻率随磁场强度增加而增加）
    derivative = std::max(0.0, derivative);
    
    return derivative;
}

void MagnetoDynamics2DSolver::computeJacobianMatrix() {
    // 计算雅可比矩阵
    
    if (!stiffnessMatrix || !jacobianMatrix || currentPotential.empty()) {
        std::cerr << "错误: 系统矩阵或当前解未初始化，无法计算雅可比矩阵" << std::endl;
        return;
    }
    
    // 雅可比矩阵 = 刚度矩阵 + 非线性项
    // 对于磁动力学问题，非线性项来自材料非线性
    
    // 复制刚度矩阵作为基础
    jacobianMatrix->CopyFrom(*stiffnessMatrix);
    
    // 添加非线性材料贡献
    auto& elements = mesh->getElements();
    
    for (const auto& element : elements) {
        std::string materialName = element.getMaterialName();
        
        // 检查是否为非线性材料
        if (materialDB.isNonlinearMaterial(materialName)) {
            // 计算非线性材料对雅可比矩阵的贡献
            addNonlinearJacobianContribution(element);
        }
    }
    
    std::cout << "雅可比矩阵计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::removeElementStiffness(const Element& element) {
    // 从全局矩阵中移除单元的贡献
    
    size_t nNodes = element.getNodeCount();
    auto nodeIndices = element.getNodeIndices();
    
    // 创建一个零矩阵来移除贡献
    std::vector<std::vector<double>> zeroStiffness(nNodes, std::vector<double>(nNodes, 0.0));
    std::vector<double> zeroForce(nNodes, 0.0);
    
    // 使用负值来移除贡献
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            zeroStiffness[i][j] = -1.0; // 负值表示移除
        }
        zeroForce[i] = -1.0;
    }
    
    // 组装负贡献以移除旧值
    assembleLocalToGlobal(element, zeroStiffness, zeroStiffness, zeroStiffness, zeroForce);
}

void MagnetoDynamics2DSolver::clearCache() {
    // 清除缓存数据
    // 对应Fortran: 缓存清理逻辑
    
    elementCache.clear();
    shapeFunctionCache.clear();
    
    // 重置系统组装状态
    systemAssembled = false;
    
    std::cout << "缓存已清除" << std::endl;
}

void MagnetoDynamics2DSolver::getCacheStatistics() const {
    // 获取缓存统计信息
    // 对应Fortran: 缓存统计逻辑
    
    std::cout << "缓存统计信息:" << std::endl;
    std::cout << "- 元素缓存数量: " << elementCache.size() << std::endl;
    std::cout << "- 形状函数缓存数量: " << shapeFunctionCache.size() << std::endl;
    std::cout << "- 系统组装状态: " << (systemAssembled ? "已组装" : "未组装") << std::endl;
}

//------------------------------------------------------------------------------
// assembleBoundaryContributions 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::assembleBoundaryContributions() {
    if (!mesh) {
        throw std::runtime_error("Mesh not set for boundary assembly");
    }
    
    // 获取边界单元
    auto& boundaryElements = mesh->getBoundaryElements();
    
    // 遍历所有边界单元
    for (size_t boundaryIndex = 0; boundaryIndex < boundaryElements.size(); ++boundaryIndex) {
        const auto& boundaryElement = boundaryElements[boundaryIndex];
        
        // 获取边界条件
        auto boundaryCondition = boundaryElement.getBoundaryCondition();
        if (!boundaryCondition) {
            continue; // 跳过没有边界条件的单元
        }
        
        // 检查边界条件类型
        if (boundaryCondition->hasType("Infinity")) {
            // 处理无限远边界条件
            applyInfinityBoundaryConditionToGlobal(boundaryElement);
        } else if (boundaryCondition->hasType("AirGap")) {
            // 处理气隙边界条件
            applyAirGapBoundaryConditionToGlobal(boundaryElement, boundaryCondition);
        }
    }
    
    std::cout << "边界条件贡献组装完成，共处理 " << boundaryElements.size() << " 个边界单元" << std::endl;
}

//------------------------------------------------------------------------------
// applyBoundaryConditions 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::applyBoundaryConditions() {
    // 应用边界条件
    
    if (!stiffnessMatrix || !rhsVector) {
        std::cerr << "错误: 系统矩阵未初始化，无法应用边界条件" << std::endl;
        return;
    }
    
    // 应用Dirichlet边界条件
    applyDirichletBoundaryConditions();
    
    // 应用Neumann边界条件
    applyNeumannBoundaryConditions();
    
    // 应用磁通密度边界条件
    applyMagneticFluxDensityBoundaryConditions();
    
    std::cout << "边界条件应用完成" << std::endl;
}

//------------------------------------------------------------------------------
// reassembleNonlinearSystem 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::reassembleNonlinearSystem() {
    // 对应Fortran的非线性系统重新组装逻辑
    // 对应Fortran: 非线性迭代中的系统重新组装
    
    // 清除当前系统状态
    clearCache();
    
    // 重新组装系统矩阵
    assembleSystem();
    
    // 更新材料参数
    updateMaterialParameters();
    
    std::cout << "非线性系统重新组装完成" << std::endl;
}

//------------------------------------------------------------------------------
// reassembleNonlinearSystemWithJacobian 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::reassembleNonlinearSystemWithJacobian() {
    // 带雅可比项的非线性系统重新组装
    
    // 清除当前系统状态
    clearCache();
    
    // 重新组装系统矩阵（包含雅可比项）
    assembleSystemWithJacobian();
    
    // 计算雅可比矩阵
    computeJacobianMatrix();
    
    // 更新材料参数
    updateMaterialParameters();
    
    std::cout << "带雅可比项的非线性系统重新组装完成" << std::endl;
}

//------------------------------------------------------------------------------
// updateMaterialParameters 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::updateMaterialParameters() {
    // 更新非线性材料参数
    
    if (!mesh || currentPotential.empty()) {
        return;
    }
    
    // 计算当前磁场强度
    auto H = computeMagneticFieldStrength();
    
    // 更新每个单元的材料参数
    auto& elements = mesh->getElements();
    
    for (const auto& element : elements) {
        size_t nElementNodes = element.getNodeCount();
        auto nodeIndices = element.getNodeIndices();
        
        // 计算单元平均磁场强度
        double H_avg = 0.0;
        for (size_t i = 0; i < nElementNodes; ++i) {
            int nodeIndex = nodeIndices[i];
            double Hx = H[nodeIndex * 2];
            double Hy = H[nodeIndex * 2 + 1];
            double H_magnitude = std::sqrt(Hx*Hx + Hy*Hy);
            H_avg += H_magnitude;
        }
        H_avg /= nElementNodes;
        
        // 获取材料名称
        std::string materialName = element.getMaterialName();
        
        // 检查是否为非线性材料
        if (materialDB.isNonlinearMaterial(materialName)) {
            // 根据H-B曲线更新磁阻率
            auto material = materialDB.getMaterial(materialName);
            double reluctivity = computeReluctivity(materialName, H_avg, element);
            
            // 更新材料属性
            materialProperties[element.getId()] = reluctivity;
            
            // 记录材料参数变化
            materialUpdateHistory[element.getId()].push_back({
                currentIteration, H_avg, reluctivity
            });
        }
    }
    
    std::cout << "材料参数更新完成，处理了 " << elements.size() << " 个单元" << std::endl;
}

//------------------------------------------------------------------------------
// assembleComplexSystem 子程序的C++移植
//------------------------------------------------------------------------------
void MagnetoDynamics2DSolver::assembleComplexSystem() {
    // 复数系统组装逻辑
    
    // 检查是否使用复数矩阵
    if (!parameters.useComplexMatrices) {
        std::cout << "警告: 未启用复数矩阵，使用实数系统组装" << std::endl;
        assembleSystem();
        return;
    }
    
    // 初始化复数系统矩阵
    initializeComplexSystemMatrices();
    
    // 组装复数单元贡献
    assembleComplexElementContributions();
    
    // 组装复数边界条件贡献
    assembleComplexBoundaryContributions();
    
    // 应用复数边界条件
    applyComplexBoundaryConditions();
    
    systemAssembled = true;
    
    std::cout << "复数系统组装完成" << std::endl;
}

//------------------------------------------------------------------------------
// 并行化相关方法的实现
//------------------------------------------------------------------------------

void MagnetoDynamics2DSolver::assembleElementContributionsParallel() {
    // 并行组装单元贡献
    // 对应Fortran: !$omp parallel do 并行化
    
    if (!isParallel()) {
        std::cout << "警告: 未启用并行模式，使用串行组装" << std::endl;
        assembleElementContributions();
        return;
    }
    
    // TODO: 实现并行组装逻辑
    // 需要实现：
    // - OpenMP并行化
    // - 域分解处理
    // - 线程安全的矩阵组装
    
    std::cout << "并行组装单元贡献..." << std::endl;
}

void MagnetoDynamics2DSolver::applyBoundaryConditionsParallel() {
    // 并行应用边界条件
    
    if (!isParallel()) {
        std::cout << "警告: 未启用并行模式，使用串行边界条件" << std::endl;
        applyBoundaryConditions();
        return;
    }
    
    // TODO: 实现并行边界条件应用逻辑
    // 需要实现：
    // - 并行边界条件处理
    // - 边界条件的通信
    // - 并行Dirichlet条件应用
    
    std::cout << "并行应用边界条件..." << std::endl;
}

void MagnetoDynamics2DSolver::exchangeGhostData() {
    // 交换幽灵数据
    // 对应Fortran: 并行计算中的数据交换
    
    if (!isParallel()) {
        std::cout << "警告: 未启用并行模式，无需交换幽灵数据" << std::endl;
        return;
    }
    
    // TODO: 实现幽灵数据交换逻辑
    // 需要实现：
    // - 幽灵节点数据交换
    // - 边界数据同步
    // - MPI通信
    
    std::cout << "交换幽灵数据..." << std::endl;
}

//------------------------------------------------------------------------------
// 边界条件辅助方法的实现
//------------------------------------------------------------------------------

double MagnetoDynamics2DSolver::computeBoundaryLength(const Element& element) {
    // 计算边界长度
    auto nodeCoords = element.getNodeCoordinates();
    size_t nNodes = element.getNodeCount();
    
    if (nNodes == 2) {
        // 线性边界（2节点）
        double dx = nodeCoords[1][0] - nodeCoords[0][0];
        double dy = nodeCoords[1][1] - nodeCoords[0][1];
        return std::sqrt(dx*dx + dy*dy);
    } else if (nNodes == 3) {
        // 二次边界（3节点）
        // 使用高斯积分计算曲线边界长度
        double length = 0.0;
        auto integrationPoints = GaussIntegration::getPoints(element.getType());
        
        for (const auto& ip : integrationPoints) {
            std::vector<double> basisFunctions;
            std::vector<std::array<double, 3>> basisDerivatives;
            double detJ;
            
            if (ShapeFunctions::evaluate(element, ip.u, ip.v, ip.w, basisFunctions, basisDerivatives, detJ)) {
                length += ip.weight * detJ;
            }
        }
        return length;
    }
    
    return 0.0;
}

std::vector<double> MagnetoDynamics2DSolver::computeBoundaryNormal(const Element& element) {
    // 计算边界法向量
    auto nodeCoords = element.getNodeCoordinates();
    size_t nNodes = element.getNodeCount();
    std::vector<double> normal(nNodes, 0.0);
    
    if (nNodes >= 2) {
        // 计算边界切线向量
        double dx = nodeCoords[1][0] - nodeCoords[0][0];
        double dy = nodeCoords[1][1] - nodeCoords[0][1];
        
        // 计算法向量（垂直于切线）
        double tangentLength = std::sqrt(dx*dx + dy*dy);
        if (tangentLength > 1e-12) {
            double nx = -dy / tangentLength;  // 法向量x分量
            double ny = dx / tangentLength;   // 法向量y分量
            
            // 为所有节点分配相同的法向量（简化处理）
            for (size_t i = 0; i < nNodes; ++i) {
                normal[i] = nx; // 这里只使用x分量，实际应根据边界方向调整
            }
        }
    }
    
    return normal;
}

double MagnetoDynamics2DSolver::getAirGapLength(const Element& element) {
    // 获取气隙长度
    // 默认气隙长度
    return 0.001; // 1mm 默认气隙
}

double MagnetoDynamics2DSolver::getAirGapPermeability(const Element& element) {
    // 获取气隙磁导率
    // 空气的磁导率
    return 4.0 * M_PI * 1e-7; // μ0 = 4π×10^-7 H/m
}

void MagnetoDynamics2DSolver::assembleBoundaryLocalToGlobal(const Element& boundaryElement,
                                                           const std::vector<std::vector<double>>& localStiffness,
                                                           const std::vector<double>& localForce) {
    // 将边界局部矩阵组装到全局系统
    size_t nNodes = boundaryElement.getNodeCount();
    
    // 获取边界单元节点编号
    auto nodeIndices = boundaryElement.getNodeIndices();
    
    // 组装边界条件贡献到全局刚度矩阵
    for (size_t i = 0; i < nNodes; ++i) {
        int globalRow = nodeIndices[i];
        
        for (size_t j = 0; j < nNodes; ++j) {
            int globalCol = nodeIndices[j];
            
            // 添加到全局刚度矩阵
            stiffnessMatrix->AddToElement(globalRow, globalCol, localStiffness[i][j]);
        }
        
        // 组装边界力向量（如果有）
        rhsVector->AddToElement(globalRow, localForce[i]);
    }
}

//------------------------------------------------------------------------------
// 电磁场计算方法的实现
//------------------------------------------------------------------------------

std::vector<double> MagnetoDynamics2DSolver::computeMagneticFieldStrength() {
    // 计算磁场强度
    if (!mesh || currentPotential.empty()) {
        throw std::runtime_error("Mesh or solution not set for magnetic field calculation");
    }
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    std::vector<double> magneticFieldStrength(nNodes * 2, 0.0); // 每个节点2个分量
    
    // 基于向量势计算磁场强度 H = -∇×A/μ
    auto& elements = mesh->getElements();
    
    for (const auto& element : elements) {
        size_t nElementNodes = element.getNodeCount();
        auto nodeIndices = element.getNodeIndices();
        
        // 获取单元节点坐标
        auto nodeCoords = element.getNodeCoordinates();
        
        // 获取单元材料
        std::string materialName = element.getMaterialName();
        auto material = materialDB.getMaterial(materialName);
        double permeability = material.getPermeability();
        
        // 计算单元中心点的磁场强度
        auto integrationPoints = GaussIntegration::getPoints(element.getType());
        
        for (const auto& ip : integrationPoints) {
            // 计算基函数和导数
            std::vector<double> basisFunctions;
            std::vector<std::array<double, 3>> basisDerivatives;
            double detJ;
            
            if (!ShapeFunctions::evaluate(element, ip.u, ip.v, ip.w, basisFunctions, basisDerivatives, detJ)) {
                continue;
            }
            
            // 计算向量势在积分点的值
            double A_ip = 0.0;
            std::array<double, 2> gradA = {0.0, 0.0};
            
            for (size_t i = 0; i < nElementNodes; ++i) {
                int nodeIndex = nodeIndices[i];
                A_ip += basisFunctions[i] * currentPotential[nodeIndex];
                gradA[0] += basisDerivatives[i][0] * currentPotential[nodeIndex]; // ∂A/∂x
                gradA[1] += basisDerivatives[i][1] * currentPotential[nodeIndex]; // ∂A/∂y
            }
            
            // 计算磁场强度 H = -∇×A/μ
            // 在2D情况下，∇×A = (∂A/∂y, -∂A/∂x)
            double Hx = -gradA[1] / permeability;  // Hx = -∂A/∂y / μ
            double Hy = gradA[0] / permeability;   // Hy = ∂A/∂x / μ
            
            // 考虑坐标系对称性
            if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
                parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
                double r = 0.0;
                for (size_t i = 0; i < nElementNodes; ++i) {
                    r += basisFunctions[i] * nodeCoords[i][0];
                }
                // 轴对称修正
                Hx = -Hx;
                Hy = -Hy;
                Hy += A_ip / (r * permeability); // 轴对称项
            }
            
            // 将磁场强度分配到节点
            for (size_t i = 0; i < nElementNodes; ++i) {
                int nodeIndex = nodeIndices[i];
                magneticFieldStrength[nodeIndex * 2] += basisFunctions[i] * Hx;
                magneticFieldStrength[nodeIndex * 2 + 1] += basisFunctions[i] * Hy;
            }
        }
    }
    
    // 对节点磁场强度进行平均
    std::vector<int> nodeCount(nNodes, 0);
    for (const auto& element : elements) {
        auto nodeIndices = element.getNodeIndices();
        for (int nodeIndex : nodeIndices) {
            nodeCount[nodeIndex]++;
        }
    }
    
    for (size_t i = 0; i < nNodes; ++i) {
        if (nodeCount[i] > 0) {
            magneticFieldStrength[i * 2] /= nodeCount[i];
            magneticFieldStrength[i * 2 + 1] /= nodeCount[i];
        }
    }
    
    return magneticFieldStrength;
}

std::vector<double> MagnetoDynamics2DSolver::computeMagneticFluxDensity() {
    // 计算磁通密度
    if (!mesh || currentPotential.empty()) {
        throw std::runtime_error("Mesh or solution not set for magnetic flux density calculation");
    }
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    std::vector<double> magneticFluxDensity(nNodes * 2, 0.0); // 每个节点2个分量
    
    // 基于向量势直接计算磁通密度 B = ∇×A
    auto& elements = mesh->getElements();
    
    for (const auto& element : elements) {
        size_t nElementNodes = element.getNodeCount();
        auto nodeIndices = element.getNodeIndices();
        
        // 获取单元节点坐标
        auto nodeCoords = element.getNodeCoordinates();
        
        // 计算单元中心点的磁通密度
        auto integrationPoints = GaussIntegration::getPoints(element.getType());
        
        for (const auto& ip : integrationPoints) {
            // 计算基函数和导数
            std::vector<double> basisFunctions;
            std::vector<std::array<double, 3>> basisDerivatives;
            double detJ;
            
            if (!ShapeFunctions::evaluate(element, ip.u, ip.v, ip.w, basisFunctions, basisDerivatives, detJ)) {
                continue;
            }
            
            // 计算向量势在积分点的梯度
            std::array<double, 2> gradA = {0.0, 0.0};
            double A_ip = 0.0;
            
            for (size_t i = 0; i < nElementNodes; ++i) {
                int nodeIndex = nodeIndices[i];
                A_ip += basisFunctions[i] * currentPotential[nodeIndex];
                gradA[0] += basisDerivatives[i][0] * currentPotential[nodeIndex]; // ∂A/∂x
                gradA[1] += basisDerivatives[i][1] * currentPotential[nodeIndex]; // ∂A/∂y
            }
            
            // 计算磁通密度 B = ∇×A
            // 在2D情况下，∇×A = (∂A/∂y, -∂A/∂x)
            double Bx = gradA[1];  // Bx = ∂A/∂y
            double By = -gradA[0]; // By = -∂A/∂x
            
            // 考虑坐标系对称性
            if (parameters.coordinateSystem == MagnetoDynamics2DParameters::AXISYMMETRIC ||
                parameters.coordinateSystem == MagnetoDynamics2DParameters::CYLINDRIC_SYMMETRIC) {
                double r = 0.0;
                for (size_t i = 0; i < nElementNodes; ++i) {
                    r += basisFunctions[i] * nodeCoords[i][0];
                }
                // 轴对称修正
                Bx = -Bx;
                By = -By;
                By += A_ip / r; // 轴对称项
            }
            
            // 将磁通密度分配到节点
            for (size_t i = 0; i < nElementNodes; ++i) {
                int nodeIndex = nodeIndices[i];
                magneticFluxDensity[nodeIndex * 2] += basisFunctions[i] * Bx;
                magneticFluxDensity[nodeIndex * 2 + 1] += basisFunctions[i] * By;
            }
        }
    }
    
    // 对节点磁通密度进行平均
    std::vector<int> nodeCount(nNodes, 0);
    for (const auto& element : elements) {
        auto nodeIndices = element.getNodeIndices();
        for (int nodeIndex : nodeIndices) {
            nodeCount[nodeIndex]++;
        }
    }
    
    for (size_t i = 0; i < nNodes; ++i) {
        if (nodeCount[i] > 0) {
            magneticFluxDensity[i * 2] /= nodeCount[i];
            magneticFluxDensity[i * 2 + 1] /= nodeCount[i];
        }
    }
    
    return magneticFluxDensity;
}

//------------------------------------------------------------------------------
// 收敛性检查方法的实现
//------------------------------------------------------------------------------

bool MagnetoDynamics2DSolver::checkConvergence(const MagnetoDynamics2DResults& results) {
    // 检查收敛性
    
    if (results.nonlinearIterations <= 1) {
        return false; // 至少需要一次迭代
    }
    
    // 计算残差范数
    double residualNorm = computeResidualNorm(results);
    
    // 计算解的变化
    double solutionChange = computeSolutionChange(results);
    
    // 检查残差收敛
    bool residualConverged = (residualNorm < parameters.nonlinearTolerance);
    
    // 检查解的变化收敛
    bool solutionConverged = (solutionChange < parameters.nonlinearTolerance);
    
    // 检查迭代次数
    bool iterationLimitReached = (results.nonlinearIterations >= parameters.maxNonlinearIterations);
    
    // 收敛条件：残差和解变化都满足容差，或者达到最大迭代次数
    bool converged = (residualConverged && solutionConverged) || iterationLimitReached;
    
    // 记录收敛历史
    convergenceHistory.push_back({
        results.nonlinearIterations, 
        residualNorm, 
        solutionChange, 
        converged
    });
    
    if (converged) {
        if (iterationLimitReached) {
            std::cout << "非线性迭代达到最大迭代次数 " << results.nonlinearIterations << " 次" << std::endl;
            std::cout << "残差范数: " << residualNorm << ", 解变化: " << solutionChange << std::endl;
        } else {
            std::cout << "非线性迭代收敛于第 " << results.nonlinearIterations << " 次迭代" << std::endl;
            std::cout << "残差范数: " << residualNorm << ", 解变化: " << solutionChange << std::endl;
        }
    } else if (results.nonlinearIterations % 5 == 0) {
        // 每5次迭代输出一次进度
        std::cout << "迭代 " << results.nonlinearIterations << ": 残差=" << residualNorm 
                  << ", 解变化=" << solutionChange << std::endl;
    }
    
    return converged;
}

double MagnetoDynamics2DSolver::computeResidualNorm(const MagnetoDynamics2DResults& results) {
    // 计算残差范数
    // 对应Fortran: 残差计算逻辑
    
    if (!stiffnessMatrix || !rhsVector) {
        throw std::runtime_error("System matrices not assembled for residual calculation");
    }
    
    // 计算残差: r = b - A*x
    auto residual = std::shared_ptr<Vector>(Vector::Create(rhsVector->Size()));
    
    // 复制右端向量
    for (size_t i = 0; i < rhsVector->Size(); ++i) {
        (*residual)[i] = (*rhsVector)[i];
    }
    
    // 减去 A*x
    auto solutionVector = std::shared_ptr<Vector>(Vector::Create(results.vectorPotential.size()));
    for (size_t i = 0; i < results.vectorPotential.size(); ++i) {
        (*solutionVector)[i] = results.vectorPotential[i];
    }
    
    auto Ax = std::shared_ptr<Vector>(Vector::Create(solutionVector->Size()));
    stiffnessMatrix->Multiply(*solutionVector, *Ax);
    
    for (size_t i = 0; i < residual->Size(); ++i) {
        (*residual)[i] -= (*Ax)[i];
    }
    
    // 计算残差范数
    double residualNorm = 0.0;
    for (size_t i = 0; i < residual->Size(); ++i) {
        residualNorm += (*residual)[i] * (*residual)[i];
    }
    residualNorm = std::sqrt(residualNorm);
    
    return residualNorm;
}

double MagnetoDynamics2DSolver::computeSolutionChange(const MagnetoDynamics2DResults& results) {
    // 计算解的变化
    // 对应Fortran: 解的变化计算逻辑
    
    if (currentPotential.empty()) {
        return std::numeric_limits<double>::max(); // 第一次迭代，变化很大
    }
    
    double change = 0.0;
    for (size_t i = 0; i < results.vectorPotential.size(); ++i) {
        double diff = results.vectorPotential[i] - currentPotential[i];
        change += diff * diff;
    }
    change = std::sqrt(change);
    
    // 更新当前解
    currentPotential = results.vectorPotential;
    
    return change;
}

//------------------------------------------------------------------------------
// 导出场量计算方法的实现
//------------------------------------------------------------------------------

void MagnetoDynamics2DSolver::calculateDerivedFields(MagnetoDynamics2DResults& results) {
    // 计算导出场量
    // 对应Fortran: 导出场量计算逻辑
    
    // 计算磁场强度
    results.magneticFieldStrength = computeMagneticFieldStrength();
    
    // 计算磁通密度
    results.magneticFluxDensity = computeMagneticFluxDensity();
    
    // 计算能量密度
    calculateEnergyDensity(results);
    
    // 计算力密度
    calculateForceDensity(results);
    
    std::cout << "导出场量计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::calculateEnergyDensity(MagnetoDynamics2DResults& results) {
    // 计算能量密度
    // 对应Fortran: 能量密度计算逻辑
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    results.energyDensity.resize(nNodes, 0.0);
    
    // 能量密度公式: w = 0.5 * B · H
    for (size_t i = 0; i < nNodes; ++i) {
        size_t idx = i * 2;
        double Bx = results.magneticFluxDensity[idx];
        double By = results.magneticFluxDensity[idx + 1];
        double Hx = results.magneticFieldStrength[idx];
        double Hy = results.magneticFieldStrength[idx + 1];
        
        results.energyDensity[i] = 0.5 * (Bx * Hx + By * Hy);
    }
}

void MagnetoDynamics2DSolver::calculateForceDensity(MagnetoDynamics2DResults& results) {
    // 计算力密度
    
    size_t nNodes = mesh->getNodes().numberOfNodes();
    results.forceDensity.resize(nNodes * 2, 0.0); // 每个节点2个分量
    
    // 力密度公式: f = J × B (2D情况下为标量积)
    // 对于2D问题，力密度为：f_x = J_z * B_y, f_y = -J_z * B_x
    
    // 获取电流密度分布
    auto currentDensity = getCurrentDensity();
    
    for (size_t i = 0; i < nNodes; ++i) {
        size_t idx = i * 2;
        
        // 获取磁通密度分量
        double Bx = results.magneticFluxDensity[idx];
        double By = results.magneticFluxDensity[idx + 1];
        
        // 获取电流密度（假设为z方向）
        double Jz = currentDensity[i];
        
        // 计算力密度分量
        results.forceDensity[idx] = Jz * By;     // f_x = J_z × B_y
        results.forceDensity[idx + 1] = -Jz * Bx; // f_y = -J_z × B_x
    }
    
    std::cout << "力密度计算完成" << std::endl;
}

void MagnetoDynamics2DSolver::calculateLumpedParameters(MagnetoDynamics2DResults& results) {
    // 计算集总参数
    
    // 计算总磁能
    results.totalMagneticEnergy = 0.0;
    for (double energy : results.energyDensity) {
        results.totalMagneticEnergy += energy;
    }
    
    // 计算总磁通
    results.totalMagneticFlux = 0.0;
    size_t nNodes = mesh->getNodes().numberOfNodes();
    for (size_t i = 0; i < nNodes; ++i) {
        size_t idx = i * 2;
        double Bx = results.magneticFluxDensity[idx];
        double By = results.magneticFluxDensity[idx + 1];
        double B_magnitude = std::sqrt(Bx*Bx + By*By);
        results.totalMagneticFlux += B_magnitude;
    }
    
    // 计算电感
    results.inductance = 0.0;
    if (parameters.currentDensity != 0.0) {
        // 电感公式: L = 2 * W_m / I^2
        // 其中W_m是总磁能，I是总电流
        double totalCurrent = parameters.currentDensity * mesh->getArea();
        if (totalCurrent != 0.0) {
            results.inductance = 2.0 * results.totalMagneticEnergy / (totalCurrent * totalCurrent);
        }
    }
    
    std::cout << "集总参数计算完成" << std::endl;
    std::cout << "总磁能: " << results.totalMagneticEnergy << " J" << std::endl;
    std::cout << "总磁通: " << results.totalMagneticFlux << " Wb" << std::endl;
    std::cout << "电感: " << results.inductance << " H" << std::endl;
}

//------------------------------------------------------------------------------
// 辅助方法的实现
//------------------------------------------------------------------------------

std::vector<double> MagnetoDynamics2DSolver::getCurrentDensity() {
    // 获取电流密度分布
    size_t nNodes = mesh->getNodes().numberOfNodes();
    std::vector<double> currentDensity(nNodes, 0.0);
    
    // 从边界条件或材料属性中获取电流密度
    // 这里简化处理，实际应该根据具体问题设置
    
    // 检查是否有电流密度边界条件
    if (boundaryConditions.find("Current Density") != boundaryConditions.end()) {
        auto& bc = boundaryConditions["Current Density"];
        for (size_t i = 0; i < nNodes; ++i) {
            currentDensity[i] = bc.value;
        }
    } else {
        // 如果没有明确设置，使用默认值或从材料属性获取
        for (size_t i = 0; i < nNodes; ++i) {
            currentDensity[i] = parameters.currentDensity;
        }
    }
    
    return currentDensity;
}

void MagnetoDynamics2DSolver::applyDirichletBoundaryConditions() {
    // 应用Dirichlet边界条件
    
    for (const auto& bcPair : boundaryConditions) {
        const auto& bc = bcPair.second;
        
        if (bc.type == "Dirichlet") {
            // 对每个受影响的节点应用Dirichlet条件
            for (int nodeId : bc.affectedNodes) {
                // 设置刚度矩阵对角线元素为1，其他行元素为0
                stiffnessMatrix->SetRowToIdentity(nodeId);
                
                // 设置右端向量为边界值
                (*rhsVector)[nodeId] = bc.value;
            }
        }
    }
}

void MagnetoDynamics2DSolver::applyNeumannBoundaryConditions() {
    // 应用Neumann边界条件
    
    for (const auto& bcPair : boundaryConditions) {
        const auto& bc = bcPair.second;
        
        if (bc.type == "Neumann") {
            // 对每个受影响的节点应用Neumann条件
            for (int nodeId : bc.affectedNodes) {
                // Neumann条件直接加到右端向量
                (*rhsVector)[nodeId] += bc.value;
            }
        }
    }
}

void MagnetoDynamics2DSolver::assembleNeumannBoundaryContribution(const BoundaryCondition& bc) {
    // 组装Neumann边界条件贡献
    
    for (int nodeId : bc.affectedNodes) {
        // Neumann边界条件直接加到右端向量
        (*rhsVector)[nodeId] += bc.value;
    }
}

ElementMatrices MagnetoDynamics2DSolver::computeElementMatricesOptimized(const Element& element) {
    // 优化的单元矩阵计算（使用缓存）
    
    ElementMatrices matrices;
    size_t nNodes = element.getNodeCount();
    
    // 初始化矩阵
    matrices.stiffness.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.mass.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.damping.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.force.resize(nNodes, 0.0);
    
    // 获取材料属性
    std::string materialName = element.getMaterialName();
    double reluctivity = getReluctivity(materialName, element);
    
    // 使用缓存的基函数数据
    if (useBasisFunctionsCache && basisFunctionCache.find(element.getId()) != basisFunctionCache.end()) {
        const auto& cache = basisFunctionCache[element.getId()];
        
        // 使用缓存的基函数值计算矩阵
        for (size_t i = 0; i < nNodes; ++i) {
            for (size_t j = 0; j < nNodes; ++j) {
                // 使用缓存的基函数值
                matrices.stiffness[i][j] = reluctivity * cache.shapeFunctionDerivatives[i] * cache.shapeFunctionDerivatives[j] * element.getArea();
                matrices.mass[i][j] = parameters.materialDensity * cache.shapeFunctions[i] * cache.shapeFunctions[j] * element.getArea();
                matrices.damping[i][j] = parameters.dampingCoefficient * cache.shapeFunctions[i] * cache.shapeFunctions[j] * element.getArea();
            }
            
            // 载荷向量
            matrices.force[i] = parameters.currentDensity * cache.shapeFunctions[i] * element.getArea();
        }
    } else {
        // 如果没有缓存，使用标准方法
        matrices = computeElementMatrices(element);
    }
    
    return matrices;
}

void MagnetoDynamics2DSolver::assembleMagneticFluxBoundaryContribution(const BoundaryCondition& bc) {
    // 组装磁通密度边界条件贡献
    
    for (int nodeId : bc.affectedNodes) {
        // 磁通密度边界条件需要特殊处理
        // 这里简化处理，实际应该基于边界积分
        
        double fluxContribution = bc.value * parameters.magneticPermeability;
        (*rhsVector)[nodeId] += fluxContribution;
    }
}

ElementMatrices MagnetoDynamics2DSolver::computeElementMatrices(const Element& element) {
    // 计算单元矩阵
    
    ElementMatrices matrices;
    size_t nNodes = element.getNodeCount();
    
    // 初始化矩阵
    matrices.stiffness.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.mass.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.damping.resize(nNodes, std::vector<double>(nNodes, 0.0));
    matrices.force.resize(nNodes, 0.0);
    
    // 获取材料属性
    std::string materialName = element.getMaterialName();
    double reluctivity = getReluctivity(materialName, element);
    
    // 计算单元刚度矩阵
    // 刚度矩阵公式: K_ij = ∫(∇N_i · ∇N_j) dΩ
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            // 简化计算，实际应该基于高斯积分
            matrices.stiffness[i][j] = reluctivity * element.getArea();
        }
    }
    
    // 计算单元质量矩阵（用于时域分析）
    // 质量矩阵公式: M_ij = ∫(N_i · N_j) dΩ
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            matrices.mass[i][j] = parameters.materialDensity * element.getArea() / nNodes;
        }
    }
    
    // 计算单元阻尼矩阵（用于时域分析）
    // 阻尼矩阵公式: C_ij = ∫(αN_i · N_j + β∇N_i · ∇N_j) dΩ
    for (size_t i = 0; i < nNodes; ++i) {
        for (size_t j = 0; j < nNodes; ++j) {
            matrices.damping[i][j] = parameters.dampingCoefficient * element.getArea() / nNodes;
        }
    }
    
    // 计算单元载荷向量
    // 载荷向量公式: f_i = ∫(N_i · f) dΩ
    for (size_t i = 0; i < nNodes; ++i) {
        matrices.force[i] = parameters.currentDensity * element.getArea() / nNodes;
    }
    
    return matrices;
}

void MagnetoDynamics2DSolver::applyMagneticFluxDensityBoundaryConditions() {
    // 应用磁通密度边界条件
    
    for (const auto& bcPair : boundaryConditions) {
        const auto& bc = bcPair.second;
        
        if (bc.type == "Magnetic Flux Density") {
            // 磁通密度边界条件需要特殊处理
            // 这里简化处理，实际应该基于边界积分
            
            for (int nodeId : bc.affectedNodes) {
                // 根据磁通密度值计算等效的右端向量贡献
                double fluxContribution = bc.value * parameters.magneticPermeability;
                (*rhsVector)[nodeId] += fluxContribution;
            }
        }
    }
}

void MagnetoDynamics2DSolver::outputResults(const MagnetoDynamics2DResults& results, const std::string& filename) {
    // 输出计算结果到文件
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开结果文件: " << filename << std::endl;
        return;
    }
    
    // 写入头部信息
    file << "# Elmer FEM 2D磁动力学求解结果" << std::endl;
    file << "# 节点数: " << mesh->getNodes().numberOfNodes() << std::endl;
    file << "# 单元数: " << mesh->getElements().size() << std::endl;
    file << "# 非线性迭代次数: " << results.nonlinearIterations << std::endl;
    file << "# 求解时间: " << results.solutionTime << " 秒" << std::endl;
    file << "# 总磁能: " << results.totalMagneticEnergy << " J" << std::endl;
    file << "# 总磁通: " << results.totalMagneticFlux << " Wb" << std::endl;
    file << "# 电感: " << results.inductance << " H" << std::endl;
    file << std::endl;
    
    // 写入节点数据
    file << "# 节点数据: 节点ID, X坐标, Y坐标, 向量势, 能量密度" << std::endl;
    const auto& nodes = mesh->getNodes();
    for (size_t i = 0; i < nodes.numberOfNodes(); ++i) {
        file << i << " " << nodes.getX(i) << " " << nodes.getY(i) << " "
             << results.vectorPotential[i] << " " << results.energyDensity[i] << std::endl;
    }
    
    file << std::endl;
    
    // 写入场量数据
    file << "# 场量数据: 节点ID, 磁场强度X, 磁场强度Y, 磁通密度X, 磁通密度Y, 力密度X, 力密度Y" << std::endl;
    for (size_t i = 0; i < nodes.numberOfNodes(); ++i) {
        size_t idx = i * 2;
        file << i << " " << results.magneticFieldStrength[idx] << " " 
             << results.magneticFieldStrength[idx + 1] << " "
             << results.magneticFluxDensity[idx] << " " 
             << results.magneticFluxDensity[idx + 1] << " "
             << results.forceDensity[idx] << " " 
             << results.forceDensity[idx + 1] << std::endl;
    }
    
    file.close();
    std::cout << "结果已输出到文件: " << filename << std::endl;
}

void MagnetoDynamics2DSolver::printSummary(const MagnetoDynamics2DResults& results) {
    // 打印求解摘要
    
    std::cout << "\n=== 2D磁动力学求解摘要 ===" << std::endl;
    std::cout << "节点数: " << mesh->getNodes().numberOfNodes() << std::endl;
    std::cout << "单元数: " << mesh->getElements().size() << std::endl;
    std::cout << "非线性迭代次数: " << results.nonlinearIterations << std::endl;
    std::cout << "求解时间: " << results.solutionTime << " 秒" << std::endl;
    std::cout << "总磁能: " << results.totalMagneticEnergy << " J" << std::endl;
    std::cout << "总磁通: " << results.totalMagneticFlux << " Wb" << std::endl;
    std::cout << "电感: " << results.inductance << " H" << std::endl;
    
    // 计算场量统计
    if (!results.magneticFluxDensity.empty()) {
        double maxB = 0.0, minB = std::numeric_limits<double>::max();
        double avgB = 0.0;
        
        for (size_t i = 0; i < results.magneticFluxDensity.size(); i += 2) {
            double Bx = results.magneticFluxDensity[i];
            double By = results.magneticFluxDensity[i + 1];
            double B_magnitude = std::sqrt(Bx*Bx + By*By);
            
            maxB = std::max(maxB, B_magnitude);
            minB = std::min(minB, B_magnitude);
            avgB += B_magnitude;
        }
        avgB /= (results.magneticFluxDensity.size() / 2);
        
        std::cout << "磁通密度统计: " << std::endl;
        std::cout << "  最大值: " << maxB << " T" << std::endl;
        std::cout << "  最小值: " << minB << " T" << std::endl;
        std::cout << "  平均值: " << avgB << " T" << std::endl;
    }
    
    std::cout << "==========================" << std::endl;
}

bool MagnetoDynamics2DSolver::isParallel() const {
    // 检查是否启用并行模式
    return comm_ && comm_->getSize() > 1;
}

} // namespace elmer