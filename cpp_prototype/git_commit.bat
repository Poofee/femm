@echo off
chcp 65001 >nul
echo Committing code changes...

REM 1. 添加所有新文件和修改的文件
echo Adding new files and modifications...
git add .

REM 2. 删除旧位置的文件（这些文件已经移动到新目录）
echo Removing old file locations...
git rm src/BoundaryConditions.cpp src/BoundaryConditions.h src/CRSMatrix.cpp src/CRSMatrix.h src/DefUtils.cpp src/DefUtils.h src/DistributedLinearAlgebra.cpp src/DistributedLinearAlgebra.h src/DomainDecomposition.cpp src/DomainDecomposition.h src/ElectromagneticMaterial.h src/ElementDescription.cpp src/ElementDescription.h src/ElementUtils.cpp src/ElementUtils.h src/ElmerSolver.cpp src/ElmerSolver.h src/GaussIntegration.cpp src/GaussIntegration.h src/HeatSolver.cpp src/HeatSolver.h src/HeatTransferSolver.cpp src/HeatTransferSolver.h src/InputFileParser.cpp src/InputFileParser.h src/Integration.cpp src/Integration.h src/Interpolation.cpp src/Interpolation.h src/IterativeSolver.cpp src/IterativeSolver.h src/LinearAlgebra.cpp src/LinearAlgebra.h src/Lists.cpp src/Lists.h src/LoggerFactory.cpp src/LoggerFactory.h src/LoggerInterface.h src/MPICommunicator.h src/MPIConfig.cpp src/MPIConfig.h src/MPIUtils.h src/MagneticSolver.cpp src/MagneticSolver.h src/MagnetoDynamics2DSolver.cpp src/MagnetoDynamics2DSolver.h src/MagnetodynamicsSolver.h src/Matrix.cpp src/Mesh.h src/MeshIO.cpp src/MeshIO.h src/MeshUtils.cpp src/MeshUtils.h src/NonlinearSolver.h src/ParallelLinearSolver.cpp src/ParallelLinearSolver.h src/ParallelMatrixAssembly.cpp src/ParallelMatrixAssembly.h src/ShapeFunctions.cpp src/ShapeFunctions.h src/SolverBase.h src/SolverRegistry.h src/SpdlogAdapter.h src/Types.h

REM 3. 提交更改
echo Committing changes...
git commit -m "refactor: 重构代码文件组织结构

- 将代码文件按照功能模块重新组织到新的目录结构中
- 创建 core/ 目录包含基础类、数学计算、网格处理、工具类和日志系统
- 创建 solvers/ 目录包含求解器模块，按物理域分类
- 创建 parallel/ 目录包含并行计算模块
- 创建 io/ 目录包含输入输出模块
- 创建 boundary/ 目录包含边界条件模块
- 更新 CMakeLists.txt 以反映新的文件结构
- 修复头文件包含路径问题
- 修复 HeatTransferSolver.cpp 中的类成员引用错误
- 验证编译通过，确保模块化结构正确"

REM 4. 显示提交状态
echo Commit completed. Current status:
git status --short

echo Code committed successfully!
pause