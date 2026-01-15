# MagneticSolver C++移植计划

## 文件分析
- **Fortran源文件**: `src/modules/MagneticSolve/MagneticSolve.F90` (828行)
- **C++目标文件**: `src/solvers/electromagnetic/MagneticSolver.cpp/h`

## 移植策略

### 第一阶段：架构设计 (已完成)
- ✅ 分析Fortran代码结构和功能
- ✅ 检查现有C++实现状态
- 🔄 设计C++移植架构和接口

### 第二阶段：核心功能移植
- [ ] 移植主求解器函数和数据结构
- [ ] 移植矩阵组装和元素计算部分
- [ ] 移植边界条件和材料参数处理

### 第三阶段：验证和优化
- [ ] 实现单元测试验证数值精度
- [ ] 验证编译和基本功能

## 详细移植计划

### 1. 数据结构映射

**Fortran类型到C++映射**:
```cpp
// Fortran类型映射
using FortranReal = double;          // REAL(KIND=dp)
using FortranInteger = int;          // INTEGER
using FortranLogical = bool;         // LOGICAL

// 数组类型映射
using FortranArray1D = std::vector<double>;
using FortranArray2D = std::vector<std::vector<double>>;
```

**关键数据结构**:
```cpp
struct MagneticSolverData {
    // 场变量
    std::vector<double> magneticField;      // MagneticField
    std::vector<double> electricCurrent;    // ElectricCurrent
    std::vector<double> elecField;          // ElecField
    
    // 速度场
    std::vector<double> U, V, W;           // 流体速度
    std::vector<double> MU, MV, MW;        // 网格速度
    
    // 材料参数
    std::vector<double> permeability;       // Permeability
    std::vector<double> conductivity;       // Conductivity
    
    // 矩阵和向量
    std::vector<std::vector<double>> massMatrix;    // MASS
    std::vector<std::vector<double>> stiffMatrix;   // STIFF
    std::vector<double> loadVector;         // LoadVector
    
    // 辅助变量
    std::vector<double> divB;               // divB
    std::vector<double> B1, B2, B3;         // B场分量
};
```

### 2. 函数接口设计

**主求解器接口**:
```cpp
class MagneticSolver : public SolverBase {
public:
    // 构造函数和析构函数
    MagneticSolver();
    ~MagneticSolver();
    
    // 核心求解函数
    bool solve(double dt, bool transientSimulation);
    
    // 辅助函数
    bool initialize();
    bool assemble();
    bool postProcess();
    
private:
    // 内部实现函数
    bool assembleCartesian();
    bool assembleAxisymmetric();
    bool assembleGeneral();
    bool assembleBoundaryConditions();
    
    // 非线性迭代求解
    bool solveNonlinearIteration(int maxIterations, double tolerance);
    
    // 材料参数处理
    bool getMaterialParameters();
    bool getVelocityField();
    
    // 内存管理
    bool allocateMemory();
    void deallocateMemory();
};
```

### 3. 移植优先级

**高优先级 (必须实现)**:
1. 主求解器函数 `MagneticSolver::solve()`
2. 矩阵组装函数 `assembleCartesian()`
3. 非线性迭代求解 `solveNonlinearIteration()`
4. 内存管理函数 `allocateMemory()`, `deallocateMemory()`

**中优先级 (重要功能)**:
1. 边界条件处理 `assembleBoundaryConditions()`
2. 材料参数获取 `getMaterialParameters()`
3. 速度场处理 `getVelocityField()`
4. 坐标系支持 `assembleAxisymmetric()`, `assembleGeneral()`

**低优先级 (优化功能)**:
1. 瞬态仿真支持
2. 并行计算优化
3. 高级边界条件
4. 性能监控和日志

### 4. 数值精度要求

- **相对误差**: ≤ 1e-12 (相对于Fortran参考实现)
- **内存使用**: ≤ 120% Fortran版本
- **性能要求**: ≥ 90% Fortran版本 (单核)

### 5. 测试计划

**单元测试**:
- 矩阵组装正确性测试
- 非线性迭代收敛性测试
- 边界条件应用测试
- 材料参数处理测试

**集成测试**:
- 与现有求解器框架集成测试
- 多物理场耦合测试
- 性能基准测试

## 实施步骤

### 步骤1: 完善现有C++实现
- 修复头文件包含问题
- 实现缺失的函数接口
- 添加必要的辅助类

### 步骤2: 移植核心算法
- 逐行翻译Fortran算法逻辑
- 保持数值精度和稳定性
- 添加适当的错误处理

### 步骤3: 集成测试
- 验证与现有框架的兼容性
- 测试数值精度和性能
- 修复发现的问题

### 步骤4: 优化和文档
- 性能优化
- 代码文档完善
- 用户指南编写

## 风险评估

**技术风险**:
- Fortran复杂数值算法的准确翻译
- 内存管理和性能优化
- 多坐标系支持实现

**缓解措施**:
- 逐行代码审查和测试
- 性能分析和优化
- 分阶段实施和验证

## 时间估计

- **阶段1 (架构设计)**: 1-2天
- **阶段2 (核心移植)**: 3-5天
- **阶段3 (验证优化)**: 2-3天
- **总计**: 6-10天

---

**状态**: 设计阶段完成，准备开始实施
**下一步**: 开始移植主求解器函数和数据结构