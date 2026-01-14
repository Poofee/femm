# Elmer FEM C++ 项目编译指南

## 项目概述
Elmer FEM C++ 移植项目是将Fortran代码库转换为现代C++17的项目，包含电磁学求解器、有限元分析等功能。

## 编译环境要求

### 必需软件
- **Visual Studio 2022** (Enterprise/Community版)
- **CMake 3.20+**
- **Windows 10/11**

### 可选依赖
- **Eigen 5.0.1** (已包含在项目中)
- **MPI** (并行计算支持)
- **HDF5** (数据存储支持)

## 编译流程

### 步骤1: 环境准备

#### 方法A: 使用环境设置脚本 (推荐)
```powershell
# 进入项目目录
cd d:\downloads\elmerfem-devel\fem\cpp_prototype

# 运行环境设置脚本
.\setup_vs_env.bat
```

**注意**: 脚本会设置Visual Studio编译器环境变量，完成后按任意键继续。

#### 方法B: 手动设置环境
```powershell
# 手动调用Visual Studio环境设置
call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
```

### 步骤2: 配置CMake项目

```powershell
# 确保在项目根目录
cd d:\downloads\elmerfem-devel\fem\cpp_prototype

# 创建构建目录 (如果不存在)
if (!(Test-Path "build")) { mkdir build }

# 配置CMake项目
cd build
cmake ..  # 注意：必须使用 .. 指向父目录的CMakeLists.txt
```

**重要提醒**: 
- ❌ **错误**: 在根目录运行 `cmake .` - 这会将生成文件放在根目录
- ✅ **正确**: 在build目录运行 `cmake ..` - 生成文件隔离在build目录中

**CMake最佳实践**:
- 所有CMake生成文件都应在build目录中
- 根目录只包含源代码和配置文件
- 支持多配置构建（Debug/Release等）

### 步骤3: 编译项目

#### 编译整个项目
```powershell
# 在build目录中执行
cmake --build . --config Release
```

#### 编译特定目标
```powershell
# 编译主库
cmake --build . --config Release --target ElmerCppMatrixAssembly

# 编译测试程序
cmake --build . --config Release --target test_magnetic_solve
cmake --build . --config Release --target test_linear_algebra
```

#### 编译调试版本
```powershell
cmake --build . --config Debug
```

### 步骤4: 运行测试

```powershell
# 运行简单测试
.\test_simple.exe

# 运行磁动力学求解器测试
.\test_magnetic_solve.exe

# 运行线性代数测试
.\test_linear_algebra.exe
```

## 常见编译问题及解决方案

### 问题1: 编译器未识别 (cl命令不存在)
**症状**: `cl : 无法将"cl"项识别为 cmdlet、函数、脚本文件或可运行程序的名称`

**解决方案**:
```powershell
# 确保已运行环境设置脚本
.\setup_vs_env.bat

# 验证编译器路径
where cl
# 应该输出: C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.39.33519\bin\Hostx64\x64\cl.exe
```

### 问题2: CMake配置失败
**症状**: `CMake Error: The source directory does not appear to contain CMakeLists.txt`

**解决方案**:
```powershell
# 确保在正确的目录
cd d:\downloads\elmerfem-devel\fem\cpp_prototype

# 重新配置
mkdir build
cd build
cmake ..
```

### 问题3: 中文注释编码警告 (C4819)
**症状**: `warning C4819: 该文件包含不能在当前代码页(936)中表示的字符`

**解决方案**: 项目已配置使用/utf-8选项，无需额外处理。

### 问题4: 链接错误
**症状**: `LNK2019: 无法解析的外部符号`

**解决方案**:
```powershell
# 清理构建目录并重新配置
cd build
rm -Recurse *
cmake ..
cmake --build . --config Release
```

## 项目结构说明

### 主要目录
- `src/` - 源代码目录
- `test/` - 测试代码目录
- `build/` - 构建输出目录
- `docs/` - 文档目录
- `eigen-5.0.1/` - Eigen数学库

### 主要编译目标
- `ElmerCppMatrixAssembly` - 主库，包含矩阵组装、求解器等核心功能
- `test_*` - 各种测试程序
- `simple_test` - 简单测试程序

## 开发工作流程

### 新增功能开发
1. 在`src/`目录添加新源文件
2. 在`CMakeLists.txt`中更新库文件列表
3. 在`test/`目录添加对应的测试
4. 重新编译并运行测试

### 调试模式编译
```powershell
# 配置调试版本
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . --config Debug
```

### 性能优化编译
```powershell
# 配置发布版本并启用优化
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

## 注意事项

1. **环境变量**: 每次打开新的PowerShell会话都需要重新运行`setup_vs_env.bat`
2. **构建目录**: 建议始终在`build`目录中进行构建操作
3. **编码问题**: 项目使用UTF-8编码，确保编辑器设置正确
4. **依赖管理**: Eigen库已包含在项目中，无需额外安装

## 快速参考命令

```powershell
# 完整编译流程 (从零开始)
cd d:\downloads\elmerfem-devel\fem\cpp_prototype
.\setup_vs_env.bat
mkdir build
cd build
cmake ..
cmake --build . --config Release

# 快速重新编译
cd build
cmake --build . --config Release --clean-first

# 仅编译主库
cd build
cmake --build . --config Release --target ElmerCppMatrixAssembly
```

## 故障排除

如果遇到编译问题，请按以下步骤排查：

1. 检查环境变量是否正确设置
2. 清理构建目录并重新配置
3. 检查CMakeLists.txt文件是否完整
4. 查看具体的编译错误信息
5. 参考`docs/compilation_issues.md`中的解决方案

---

**文档版本**: 1.0  
**最后更新**: 2026-01-14  
**维护者**: Elmer FEM C++移植团队