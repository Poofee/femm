# Elmer FEM C++ 移植项目开发规范

## 项目概述

本项目是Elmer FEM有限元软件的C++17移植版本，旨在将Fortran代码库转换为现代、安全的C++代码，同时保持数值精度和性能。

## 开发环境配置

### 编译环境
- **编译器**: MSVC (Windows) / GCC (Linux) / Clang (macOS)
- **标准**: C++17
- **构建系统**: CMake 3.15+

### 依赖库
- **Eigen**: 线性代数计算
- **spdlog**: 日志系统
- **Google Test**: 单元测试框架
- **MPI**: 并行计算支持（可选）

## 项目编译

### 编译命令
```bash
# 创建构建目录
mkdir build && cd build

# 配置CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# 编译项目
cmake --build . --config Release

# 运行测试
ctest --output-on-failure
```

### 编译选项
- `-DUSE_MPI=ON`: 启用MPI并行支持
- `-DBUILD_TESTS=ON`: 构建测试套件
- `-DCMAKE_BUILD_TYPE=Debug/Release`: 构建类型

## 代码结构框架

```
cpp_prototype/
├── src/
│   ├── core/                    # 核心功能模块
│   │   ├── base/               # 基础类（求解器基类等）
│   │   ├── mesh/               # 网格处理
│   │   ├── utils/              # 工具函数
│   │   ├── parallel/           # 并行计算
│   │   └── logging/            # 日志系统
│   ├── solvers/                # 求解器实现
│   │   ├── ElmerSolver.cpp     # 主求解程序
│   │   └── electromagnetic/    # 电磁求解器
│   ├── boundary/               # 边界条件
│   └── materials/              # 材料模型
├── tests/                      # 单元测试
├── docs/                       # 文档
└── CMakeLists.txt              # 构建配置
```

## 代码规范

### 1. 头文件包含规范
只包含文件名，不要包含路径。路径已经配置好了。

**正确示例**:
```cpp
#include "ElmerSolver.h"
#include "LoggerFactory.h"
#include "Mesh.h"
```

**错误示例**:
```cpp
#include "../solvers/ElmerSolver.h"           // 禁止相对路径
去掉 ../solvers/
#include "d:/project/src/solvers/ElmerSolver.h" // 禁止绝对路径
去掉 d:/project/src/solvers/
```

**规范要求**:
- 只包含头文件（.h文件）
- 使用项目根目录的相对路径
- 禁止使用`../`等相对路径
- 禁止使用绝对路径
- 确保CMakeLists.txt中正确配置include路径

### 2. 变量使用规范

**正确示例**:
```cpp
// 声明后使用
const auto& nodesContainer = mesh_->getNodes();
const auto& nodes = nodesContainer.getNodes();
for (const auto& node : nodes) {
    // 使用node
}
```

**错误示例**:
```cpp
// 未声明直接使用
for (const auto& node : nodes) {  // nodes未声明
    // ...
}
```

**规范要求**:
- 所有变量在使用前必须声明
- 检查变量作用域和生命周期
- 避免重复定义变量

### 3. 类成员函数声明与实现对应规范

**正确示例**:
```cpp
// 头文件声明
class ElmerSolver {
public:
    virtual bool initialize();
    virtual bool loadModel();
};

// 源文件实现
bool ElmerSolver::initialize() {
    // 实现代码
}

bool ElmerSolver::loadModel() {
    // 实现代码
}
```

**规范要求**:
- 头文件中的每个public/protected成员函数都必须在源文件中实现
- 检查函数签名（返回类型、参数列表）完全一致
- 虚函数必须正确实现
- 使用`override`关键字明确重写

### 4. 日志接口使用规范

#### 日志级别定义
- `ELMER_TRACE`: 跟踪级别，最详细的日志
- `ELMER_DEBUG`: 调试级别，用于调试信息  
- `ELMER_INFO`: 信息级别，常规信息
- `ELMER_WARN`: 警告级别，潜在问题
- `ELMER_ERROR`: 错误级别，错误信息
- `ELMER_CRITICAL`: 严重级别，严重错误

#### 使用示例
```cpp
// 替换std::cout为日志接口
ELMER_INFO("开始初始化ElmerSolver...");

// 替换std::cerr为日志接口
ELMER_ERROR("错误: 网格目录不存在: {}", meshDir);

// 格式化输出
ELMER_DEBUG("当前残差: {}", residual);
ELMER_INFO("时间步 {}/{}, 当前时间: {}s", step, totalSteps, currentTime);

// 警告信息
ELMER_WARN("警告: 达到最大迭代次数，可能未完全收敛");
```

#### 日志接口替换规则
- `std::cout << "信息" << std::endl;` → `ELMER_INFO("信息");`
- `std::cerr << "错误" << std::endl;` → `ELMER_ERROR("错误");`
- 包含变量的输出使用格式化字符串`{}`

### 5. 命名规范

#### 类名
- 使用PascalCase: `ElmerSolver`, `MagneticSolver`

#### 函数名
- 使用camelCase: `initialize()`, `solveSteadyState()`

#### 变量名
- 使用snake_case: `mesh_`, `solver_manager_`
- 成员变量以`_`结尾: `initialized_`, `mesh_loaded_`

#### 常量
- 使用UPPER_CASE: `MAX_ITERATIONS`, `TOLERANCE`

### 6. 内存管理规范

#### 智能指针使用
```cpp
// 推荐使用智能指针
std::shared_ptr<Mesh> mesh_ = std::make_shared<Mesh>();
std::unique_ptr<SolverManager> solver_manager_ = std::make_unique<SolverManager>();

// 避免裸指针
Mesh* mesh_;  // 不推荐
```

#### RAII原则
- 资源在构造函数中获取，在析构函数中释放
- 使用智能指针管理资源生命周期

### 7. 错误处理规范

#### 返回值检查
```cpp
if (!mesh_->loadFromFile(filepath)) {
    ELMER_ERROR("网格文件加载失败: {}", filepath);
    return false;
}
```

#### 异常处理
```cpp
try {
    // 可能抛出异常的代码
} catch (const std::exception& e) {
    ELMER_ERROR("异常发生: {}", e.what());
    return false;
}
```

## 开发流程

### 1. 代码修改流程
1. 从Fortran源代码选择要移植的子程序
2. 编写或更新对应的单元测试
3. 实现C++版本，确保数值精度在1e-12以内
4. 运行测试套件验证正确性
5. 记录性能指标和修改日志

### 2. 提交规范
- 每次提交对应一个逻辑完整的移植单元
- 提交信息遵循Conventional Commits规范
- 提交前确保编译通过且测试通过

### 3. 测试要求
- 每个求解器必须有对应的单元测试
- 测试覆盖率应达到80%以上
- 数值结果应与Fortran版本在1e-12容差内一致

## 常见问题解决

### 编译错误处理
1. **未声明的标识符**: 检查头文件包含和变量声明
2. **函数未实现**: 检查头文件声明和源文件实现是否对应
3. **类型不匹配**: 检查函数签名和参数类型
4. **链接错误**: 检查库依赖和链接顺序

### 调试技巧
1. 使用`ELMER_DEBUG`输出调试信息
2. 启用Debug构建进行调试
3. 使用单元测试隔离问题
4. 检查数值精度差异

## 性能要求

- **单核性能**: ≥ 90% 原始Fortran版本
- **内存占用**: ≤ 120% 原始Fortran版本
- **并行效率**: MPI和OpenMP并行化支持

## 文档维护

### 必须维护的文档
- `docs/development_guide.md`: 本开发规范文档
- `docs/fortran_cpp_faq.md`: Fortran到C++常见模式转换FAQ
- `docs/design_decisions.md`: 重要设计决策记录
- `daily_log.md`: 每日开发日志
- `project_status.md`: 项目状态跟踪

### 文档更新要求
- 每次架构变更后更新设计文档
- 每日工作结束后更新开发日志
- 发现新模式时更新FAQ

## 联系方式

- 项目负责人: [负责人姓名]
- 技术讨论: [讨论组/邮件列表]
- 问题报告: [Issue跟踪系统]

---

*最后更新: 2026-01-16*
*版本: 1.0*