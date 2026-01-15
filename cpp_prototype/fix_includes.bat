@echo off
chcp 65001 >nul
echo Fixing header include paths...

cd src

REM 修复并行计算模块的头文件包含
REM MPI模块
sed -i "s/#include \"MPIConfig.h\"/#include \"..\\..\\parallel\\mpi\\MPIConfig.h\"/g" parallel/mpi/DistributedLinearAlgebra.h
sed -i "s/#include \"LinearAlgebra.h\"/#include \"..\\..\\core\\math\\LinearAlgebra.h\"/g" parallel/mpi/DistributedLinearAlgebra.h

sed -i "s/#include \"DistributedLinearAlgebra.h\"/#include \"..\\..\\parallel\\mpi\\DistributedLinearAlgebra.h\"/g" parallel/openmp/ParallelMatrixAssembly.h

sed -i "s/#include \"MPIConfig.h\"/#include \"..\\..\\parallel\\mpi\\MPIConfig.h\"/g" parallel/domain/DomainDecomposition.h
sed -i "s/#include \"DistributedLinearAlgebra.h\"/#include \"..\\..\\parallel\\mpi\\DistributedLinearAlgebra.h\"/g" parallel/domain/ParallelLinearSolver.h

REM 修复核心模块的头文件包含
sed -i "s/#include \"Types.h\"/#include \"..\\core\\base\\Types.h\"/g" MatrixAssembly.h
sed -i "s/#include \"Types.h\"/#include \"..\\core\\base\\Types.h\"/g" ElementMatrix.h

REM 修复边界条件模块的头文件包含
sed -i "s/#include \"LinearAlgebra.h\"/#include \"..\\core\\math\\LinearAlgebra.h\"/g" boundary/BoundaryConditions.h

REM 修复求解器模块的头文件包含
sed -i "s/#include \"SolverBase.h\"/#include \"..\\core\\base\\SolverBase.h\"/g" solvers/ElmerSolver.h
sed -i "s/#include \"SolverRegistry.h\"/#include \"..\\core\\base\\SolverRegistry.h\"/g" solvers/ElmerSolver.h
sed -i "s/#include \"Types.h\"/#include \"..\\core\\base\\Types.h\"/g" solvers/ElmerSolver.h

REM 修复电磁求解器的头文件包含
sed -i "s/#include \"SolverBase.h\"/#include \"..\\..\\core\\base\\SolverBase.h\"/g" solvers/electromagnetic/MagneticSolver.h

REM 修复热求解器的头文件包含
sed -i "s/#include \"SolverBase.h\"/#include \"..\\..\\core\\base\\SolverBase.h\"/g" solvers/thermal/HeatSolver.h

REM 修复结构求解器的头文件包含
sed -i "s/#include \"SolverBase.h\"/#include \"..\\core\\base\\SolverBase.h\"/g" StructuralSolver.h

REM 修复多物理场求解器的头文件包含
sed -i "s/#include \"SolverBase.h\"/#include \"..\\core\\base\\SolverBase.h\"/g" MultiphysicsSolver.h

echo Header include paths fixed!
cd ..
pause