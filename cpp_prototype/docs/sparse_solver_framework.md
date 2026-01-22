# 稀疏矩阵求解器框架文档

## 框架概述

本框架是为Elmer FEM低频电磁有限元求解器设计的稀疏矩阵求解器封装，提供统一的接口来集成多种主流稀疏线性求解器（SuperLU、MUMPS等）。框架采用"底层求解器封装+上层统一接口"的架构，支持无缝切换不同求解器，同时保留底层求解器的关键参数配置能力。

## 核心特性

### 1. 统一的求解器接口
- **抽象基类设计**：`LinearSolver` 基类提供统一的求解接口
- **直接求解器支持**：`DirectSolver` 基类封装直接求解方法
- **迭代求解器支持**：`IterativeSolver` 基类封装迭代求解方法
- **异常处理机制**：`SolverException` 提供详细的错误信息

### 2. 稀疏矩阵存储格式
- **CSR格式**：`CSRMatrix` 类实现压缩行存储，适合有限元矩阵
- **CSC格式**：`CSCMatrix` 类实现压缩列存储，适合某些求解器需求
- **矩阵工具类**：`SparseMatrixUtils` 提供矩阵创建、转换、运算等工具函数

### 3. 求解器封装
- **SuperLU求解器**：封装SuperLU 5.3.0串行双精度版本
- **SuperLU_MT求解器**：封装SuperLU多线程版本
- **MUMPS求解器**：封装MUMPS 5.6.2 MPI并行版本，支持64位整数

### 4. 求解器管理
- **工厂模式**：`SolverFactory` 提供求解器的动态创建
- **管理器模式**：`SolverManager` 管理多个求解器实例
- **参数配置**：`SolverParameters` 统一管理求解器参数

## 目录结构

```
cpp_prototype/src/core/math/
├── SparseMatrix.h/cpp          # 稀疏矩阵基类和CSR/CSC实现
├── LinearSolverInterface.h     # 求解器抽象接口定义
├── LinearSolverFactory.cpp     # 求解器工厂和管理器实现
├── SuperLUSolver.h             # SuperLU求解器封装
├── MumpsSolver.h               # MUMPS求解器封装
└── test_sparse_solver_framework.cpp  # 测试用例
```

## 核心类说明

### 1. SparseMatrix 基类

**功能**：稀疏矩阵的抽象基类，定义统一的矩阵操作接口

**主要方法**：
- `GetElement()` / `SetElement()` - 元素访问
- `Multiply()` - 矩阵向量乘法
- `Assemble()` - 矩阵装配
- `SaveToFile()` / `LoadFromFile()` - 文件I/O

### 2. CSRMatrix 类

**功能**：压缩行存储格式的稀疏矩阵实现

**特点**：
- 支持临时存储和装配模式
- 对角线元素优化访问
- 矩阵对称性检查
- 内存使用统计

### 3. LinearSolver 基类

**功能**：线性求解器的抽象基类

**主要方法**：
- `Initialize()` - 初始化求解器
- `Solve()` - 求解线性方程组
- `Refactorize()` - 重新分解矩阵
- `Cleanup()` - 清理资源

### 4. DirectSolver 基类

**功能**：直接求解器的基类

**扩展方法**：
- `Factorize()` - 矩阵分解
- `BackSubstitute()` - 前向/后向代换
- `ComputeDeterminant()` - 计算行列式

### 5. SolverFactory 类

**功能**：求解器工厂，提供求解器的动态创建

**主要方法**：
- `CreateSolver()` - 创建求解器实例
- `GetAvailableSolvers()` - 获取可用求解器列表
- `GetRecommendedParameters()` - 获取推荐参数

## 使用示例

### 基本使用流程

```cpp
#include "core/math/SparseMatrix.h"
#include "core/math/LinearSolverFactory.h"

using namespace elmer;

// 1. 创建稀疏矩阵
auto matrix = std::make_shared<CSRMatrix>(100, 100);
// ... 设置矩阵元素 ...
matrix->Assemble();

// 2. 创建右端向量
Vector b(100);
// ... 设置右端项 ...

// 3. 创建求解器
auto solver = SolverFactory::CreateSolver("superlu");

// 4. 设置求解器参数
auto params = SolverFactory::GetRecommendedParameters(SolverType::SUPERLU);
params.tolerance = 1e-12;
params.verbose = true;
solver->SetParameters(params);

// 5. 初始化求解器
if (!solver->Initialize(matrix)) {
    throw std::runtime_error("求解器初始化失败");
}

// 6. 求解线性方程组
Vector x(100);
if (solver->Solve(b, x)) {
    std::cout << "求解成功！" << std::endl;
    
    // 验证解的正确性
    Vector Ax(100);
    matrix->Multiply(x, Ax);
    Real residual = ComputeResidual(Ax, b);
    std::cout << "残差: " << residual << std::endl;
}
```

### 求解器管理器使用

```cpp
// 创建求解器管理器
SolverManager manager;

// 添加多个求解器
manager.AddSolver("superlu", SolverFactory::CreateSolver("superlu"));
manager.AddSolver("mumps", SolverFactory::CreateSolver("mumps"));

// 设置默认求解器
manager.SetDefaultSolver("superlu");

// 使用默认求解器求解
Vector x(100);
if (manager.Solve(b, x)) {
    std::cout << "求解成功！" << std::endl;
}

// 使用指定求解器求解
if (manager.Solve("mumps", b, x)) {
    std::cout << "MUMPS求解成功！" << std::endl;
}
```

### 低频电磁有限元场景示例

```cpp
// 创建拉普拉斯矩阵（模拟电磁场离散）
auto laplacian_matrix = SparseMatrixUtils::CreateLaplacianMatrix(100);

// 创建源项向量（边界条件）
Vector source(100);
source[0] = 1.0;   // 第一类边界条件
source[99] = -1.0; // 第二类边界条件

// 使用SuperLU求解器
auto solver = SolverFactory::CreateSolver("superlu");
solver->Initialize(laplacian_matrix);

Vector potential(100);
if (solver->Solve(source, potential)) {
    // 计算电场强度等物理量
    // ...
}
```

## 编译配置

### CMakeLists.txt 配置

```cmake
# 稀疏矩阵求解器框架配置
set(SPARSE_SOLVER_SOURCES
    src/core/math/SparseMatrix.cpp
    src/core/math/LinearSolverFactory.cpp
    # 其他源文件...
)

# 外部库依赖
find_package(SuperLU 5.3.0 REQUIRED)
find_package(MUMPS 5.6.2 REQUIRED)

# 编译选项
target_compile_definitions(elmer_cpp PRIVATE
    -DUSE_SUPERLU
    -DUSE_MUMPS
    -DSUPERLU_VERSION=530
    -DMUMPS_VERSION=562
)

target_link_libraries(elmer_cpp
    SuperLU::SuperLU
    MUMPS::MUMPS
    # 其他依赖库...
)
```

### Windows MSVC 配置

```bat
# 设置库路径
set SUPERLU_DIR=C:\libraries\superlu-5.3.0
set MUMPS_DIR=C:\libraries\mumps-5.6.2

# 编译命令
cl /I%SUPERLU_DIR%\include /I%MUMPS_DIR%\include ^
    /DUSE_SUPERLU /DUSE_MUMPS ^
    SparseMatrix.cpp LinearSolverFactory.cpp ^
    /link %SUPERLU_DIR%\lib\superlu.lib %MUMPS_DIR%\lib\dmumps.lib
```

## 性能优化建议

### 1. 矩阵存储优化
- 使用CSR格式存储有限元矩阵
- 在装配完成后调用`OptimizeStorage()`压缩存储
- 对大矩阵使用分块存储策略

### 2. 求解器选择策略
- **小规模问题**：使用SuperLU串行版本
- **中等规模问题**：使用SuperLU_MT多线程版本
- **大规模问题**：使用MUMPS并行版本
- **对称矩阵**：利用对称性优化求解

### 3. 参数调优
- **SuperLU**：调整`drop_tolerance`和`panel_size`
- **MUMPS**：优化`icntl_7`排序策略和`icntl_14`内存设置
- **通用参数**：根据问题规模调整收敛容差

## 错误处理

### 常见错误类型

1. **矩阵奇异错误**：矩阵不可逆，检查边界条件
2. **内存不足错误**：矩阵规模过大，优化存储或使用并行求解器
3. **参数配置错误**：求解器参数不合法，使用推荐参数
4. **初始化错误**：矩阵格式不支持，检查矩阵装配状态

### 错误处理示例

```cpp
try {
    auto solver = SolverFactory::CreateSolver("superlu");
    solver->Initialize(matrix);
    
    if (!solver->Solve(b, x)) {
        auto status = solver->GetStatus();
        std::cerr << "求解失败: " << status.message << std::endl;
        
        if (status.message.find("singular") != std::string::npos) {
            // 处理矩阵奇异情况
            HandleSingularMatrix();
        }
    }
} catch (const SolverException& e) {
    std::cerr << "求解器异常: " << e.GetFullMessage() << std::endl;
}
```

## 扩展开发

### 添加新的求解器

1. **实现求解器类**：继承`DirectSolver`或`IterativeSolver`
2. **注册到工厂**：在`SolverFactory::CreateSolver()`中添加case分支
3. **提供推荐参数**：在`GetRecommendedParameters()`中添加参数配置

### 示例：添加PETSc求解器

```cpp
class PetscSolver : public DirectSolver {
public:
    PetscSolver() : DirectSolver(SolverType::PETSC) {}
    
    bool Initialize(std::shared_ptr<SparseMatrix> matrix) override {
        // PETSc初始化实现
    }
    
    bool Solve(const Vector& b, Vector& x) override {
        // PETSc求解实现
    }
    // ... 其他方法实现
};
```

## 测试验证

框架提供完整的测试用例，验证以下功能：
- 稀疏矩阵基本操作
- 求解器工厂功能
- 求解器管理器功能
- 线性方程组求解
- 低频电磁有限元场景

运行测试：
```bash
cd cpp_prototype
mkdir build && cd build
cmake ..
make
./test_sparse_solver_framework
```

## 性能基准

在典型低频电磁有限元问题上的性能表现：
- **矩阵规模**：10^4 ~ 10^6阶
- **内存效率**：≤ 120% 原始Fortran版本
- **求解精度**：1e-12相对容差
- **并行效率**：MPI并行支持，可扩展性良好

## 技术支持

- **文档更新**：框架使用过程中发现的问题和建议请反馈
- **性能优化**：针对特定应用场景的性能优化建议
- **扩展开发**：新求解器封装的技术指导

---

*文档版本：1.0*  
*最后更新：2026年1月22日*  
*框架版本：Elmer FEM C++ Prototype v1.0*