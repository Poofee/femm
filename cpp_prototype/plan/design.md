# Elmer FEM C++移植项目架构设计文档

## 系统架构概述

本项目采用分层架构设计，将Elmer FEM的功能模块化重构为现代C++架构。整体架构分为核心层、求解器层和应用层，支持渐进式迁移策略。

## 核心模块划分

### 1. 核心基础库 (Core Foundation)

#### 1.1 线性代数模块 (LinearAlgebra)
```cpp
namespace elmer {
    class Matrix;           // 抽象矩阵接口
    class Vector;           // 向量类
    class CRSMatrix;        // 压缩行存储矩阵
    class DenseMatrix;      // 稠密矩阵
    class LinearSolver;     // 线性求解器基类
}
```

**功能职责**:
- 提供统一的矩阵/向量操作接口
- 实现稀疏和稠密矩阵存储
- 支持基本的线性代数运算
- 为求解器提供数值计算基础

#### 1.2 网格处理模块 (Mesh)
```cpp
namespace elmer {
    class Mesh;             // 网格数据结构
    class Element;          // 单元基类
    class Node;            // 节点类
    class MeshIO;          // 网格I/O接口
    class MeshUtils;       // 网格工具函数
}
```

**功能职责**:
- 管理网格拓扑和几何信息
- 支持多种单元类型
- 提供网格质量检查工具
- 实现网格文件读写功能

#### 1.3 有限元核心模块 (FEMCore)
```cpp
namespace elmer {
    class ShapeFunction;    // 形函数
    class IntegrationRule; // 积分规则
    class ElementMatrix;   // 单元矩阵组装
    class BoundaryCondition; // 边界条件
}
```

**功能职责**:
- 实现有限元方法的核心算法
- 提供形函数和数值积分
- 管理单元矩阵组装过程
- 处理各类边界条件

### 2. 求解器层 (Solver Layer)

#### 2.1 迭代求解器模块 (IterativeSolvers)
```cpp
namespace elmer {
    class CGSolver;         // 共轭梯度法
    class BiCGStabSolver;   // 稳定双共轭梯度法
    class GMRESSolver;      // 广义最小残差法
    class Preconditioner;   // 预处理器基类
}
```

**功能职责**:
- 实现主流迭代求解算法
- 提供预处理技术支持
- 管理求解器参数和收敛控制

#### 2.2 直接求解器模块 (DirectSolvers)
```cpp
namespace elmer {
    class LUSolver;         // LU分解求解器
    class CholeskySolver;   // Cholesky分解求解器
    class SparseLUSolver;   // 稀疏LU求解器
}
```

**功能职责**:
- 实现直接求解方法
- 集成第三方求解库（如SuperLU、MUMPS）
- 提供小规模问题的精确解

### 3. 物理场求解器模块 (PhysicsSolvers)

#### 3.1 电磁场求解器 (Electromagnetics)
```cpp
namespace elmer {
    class MagnetoDynamicsSolver;  // 磁动力学求解器
    class ElectrostaticsSolver;   // 静电求解器
    class WaveEquationSolver;      // 波动方程求解器
}
```

#### 3.2 力学求解器 (Mechanics)
```cpp
namespace elmer {
    class ElasticitySolver;        // 弹性力学求解器
    class PlasticitySolver;        // 塑性力学求解器
    class ContactSolver;           // 接触问题求解器
}
```

#### 3.3 流体求解器 (Fluids)
```cpp
namespace elmer {
    class NavierStokesSolver;      // 纳维-斯托克斯求解器
    class HeatTransferSolver;      // 热传导求解器
    class MultiphaseSolver;        // 多相流求解器
}
```

## 数据结构设计

### 4.1 核心数据结构

#### 4.1.1 稀疏矩阵存储 (CRS格式)
```cpp
class CRSMatrix : public Matrix {
private:
    std::vector<double> values;      // 非零元素值
    std::vector<int> col_indices;    // 列索引
    std::vector<int> row_pointers;   // 行指针
    int rows, cols;                  // 矩阵维度
    
public:
    // 矩阵-向量乘法
    void Multiply(const Vector& x, Vector& y) const override;
    // 设置矩阵元素
    void SetElement(int row, int col, double value) override;
    // 获取矩阵元素
    double GetElement(int row, int col) const override;
};
```

#### 4.1.2 网格数据结构
```cpp
class Mesh {
private:
    std::vector<Node> nodes;                    // 节点列表
    std::vector<std::unique_ptr<Element>> elements; // 单元列表
    std::vector<BoundaryCondition> bcs;         // 边界条件
    
public:
    // 网格操作接口
    int GetNodeCount() const;
    int GetElementCount() const;
    const Element& GetElement(int index) const;
    void AddElement(std::unique_ptr<Element> element);
};
```

### 4.2 算法设计模式

#### 4.2.1 策略模式 - 求解器选择
```cpp
class SolverContext {
private:
    std::unique_ptr<LinearSolver> solver_;
    
public:
    void SetSolver(std::unique_ptr<LinearSolver> solver) {
        solver_ = std::move(solver);
    }
    
    bool Solve(const Matrix& A, Vector& x, const Vector& b) {
        return solver_->Solve(A, x, b);
    }
};
```

#### 4.2.2 工厂模式 - 对象创建
```cpp
class SolverFactory {
public:
    static std::unique_ptr<LinearSolver> CreateSolver(SolverType type) {
        switch (type) {
            case SolverType::CG: return std::make_unique<CGSolver>();
            case SolverType::BiCGStab: return std::make_unique<BiCGStabSolver>();
            case SolverType::GMRES: return std::make_unique<GMRESSolver>();
            default: throw std::invalid_argument("Unknown solver type");
        }
    }
};
```

## 接口设计

### 5.1 C++抽象接口

#### 5.1.1 矩阵接口
```cpp
class Matrix {
public:
    virtual ~Matrix() = default;
    virtual int Rows() const = 0;
    virtual int Cols() const = 0;
    virtual void Multiply(const Vector& x, Vector& y) const = 0;
    virtual void SetElement(int row, int col, double value) = 0;
    virtual double GetElement(int row, int col) const = 0;
};
```

#### 5.1.2 求解器接口
```cpp
class LinearSolver {
public:
    virtual ~LinearSolver() = default;
    virtual bool Solve(const Matrix& A, Vector& x, const Vector& b) = 0;
    virtual void SetTolerance(double tol) = 0;
    virtual void SetMaxIterations(int max_iter) = 0;
    virtual int GetIterationCount() const = 0;
    virtual double GetResidual() const = 0;
};
```

### 5.2 Fortran-C++互操作接口

#### 5.2.1 类型映射
```cpp
extern "C" {
    // Fortran矩阵到C++矩阵的转换
    void fortran_to_cpp_matrix(double* values, int* col_indices, 
                              int* row_pointers, int n, int m,
                              CRSMatrix* cpp_matrix);
    
    // C++矩阵到Fortran矩阵的转换
    void cpp_to_fortran_matrix(const CRSMatrix& cpp_matrix,
                              double* values, int* col_indices,
                              int* row_pointers);
}
```

## 构建和依赖管理

### 6.1 CMake构建系统
```cmake
# 核心模块
add_library(elmer_core
    src/LinearAlgebra/Matrix.cpp
    src/LinearAlgebra/Vector.cpp
    src/Mesh/Mesh.cpp
    src/Mesh/Element.cpp
    src/FEMCore/ShapeFunction.cpp
)

# 求解器模块
add_library(elmer_solvers
    src/Solvers/IterativeSolvers.cpp
    src/Solvers/DirectSolvers.cpp
)

# 物理场求解器
add_library(elmer_physics
    src/Physics/Electromagnetics.cpp
    src/Physics/Mechanics.cpp
    src/Physics/Fluids.cpp
)
```

### 6.2 外部依赖管理
```cmake
# Eigen库依赖
find_package(Eigen3 REQUIRED)
target_link_libraries(elmer_core Eigen3::Eigen)

# MPI支持（可选）
if(USE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(elmer_solvers MPI::MPI_CXX)
endif()

# VTK可视化支持（可选）
if(USE_VTK)
    find_package(VTK REQUIRED)
    target_link_libraries(elmer_physics VTK::CommonCore)
endif()
```

## 性能优化策略

### 7.1 内存优化
- 使用移动语义避免不必要的拷贝
- 实现自定义内存分配器
- 优化稀疏矩阵存储格式

### 7.2 计算优化
- 利用SIMD指令集优化关键循环
- 实现多线程并行计算
- 使用表达式模板优化矩阵运算

### 7.3 缓存优化
- 优化数据访问模式提高缓存命中率
- 使用分块算法减少缓存失效
- 实现数据预取机制

## 测试策略

### 8.1 单元测试
```cpp
TEST(MatrixTest, MatrixVectorMultiplication) {
    auto matrix = CreateTestMatrix(3);
    auto x = CreateTestVector(3);
    auto y = Vector::Create(3);
    
    matrix->Multiply(*x, *y);
    
    // 验证计算结果
    EXPECT_NEAR((*y)[0], 2.0, 1e-10);
    EXPECT_NEAR((*y)[1], 4.0, 1e-10);
    EXPECT_NEAR((*y)[2], 10.0, 1e-10);
}
```

### 8.2 集成测试
- 验证Fortran-C++互操作性
- 测试完整求解流程
- 性能基准测试

### 8.3 回归测试
- 与原始Fortran版本结果对比
- 确保数值精度一致性
- 验证边界条件处理正确性

---

*文档版本: 1.0*  
*最后更新: 2026-01-08*  
*负责人: 项目架构师*