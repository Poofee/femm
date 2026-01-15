@echo off
chcp 65001 >nul
echo Fixing compilation issues...

REM 1. 修复头文件包含路径
cd src

echo Fixing header include paths...

REM 修复IterativeSolver.h
powershell -Command "(Get-Content 'core\math\IterativeSolver.h') -replace '#include \"Types.h\"', '#include \"../base/Types.h\"' | Set-Content 'core\math\IterativeSolver.h'"

REM 修复ShapeFunctions.h
powershell -Command "(Get-Content 'core\utils\ShapeFunctions.h') -replace '#include \"Types.h\"', '#include \"../base/Types.h\"' | Set-Content 'core\utils\ShapeFunctions.h'"

REM 修复Mesh.h
powershell -Command "(Get-Content 'core\mesh\Mesh.h') -replace '#include \"Types.h\"', '#include \"../base/Types.h\"' | Set-Content 'core\mesh\Mesh.h'"

REM 修复MatrixAssembly.h
powershell -Command "(Get-Content 'MatrixAssembly.h') -replace '#include \"../core/base/Types.h\"', '#include \"core/base/Types.h\"' | Set-Content 'MatrixAssembly.h'"

REM 修复ElementMatrix.h
powershell -Command "(Get-Content 'ElementMatrix.h') -replace '#include \"../core/base/Types.h\"', '#include \"core/base/Types.h\"' | Set-Content 'ElementMatrix.h'"

REM 修复BoundaryConditions.h
powershell -Command "(Get-Content 'boundary\BoundaryConditions.h') -replace '#include \"Mesh.h\"', '#include \"../core/mesh/Mesh.h\"' | Set-Content 'boundary\BoundaryConditions.h'"

REM 2. 暂时禁用MPI相关模块
cd ..

echo Temporarily disabling MPI modules...

REM 创建一个简化的CMakeLists.txt，不包含MPI模块
copy CMakeLists.txt CMakeLists_backup.txt >nul

echo Creating simplified CMakeLists.txt without MPI modules...
(
echo cmake_minimum_required(VERSION 3.10)
echo project(ElmerCppMatrixAssembly)
echo.
echo # Set C++ standard
echo set(CMAKE_CXX_STANDARD 17)
echo set(CMAKE_CXX_STANDARD_REQUIRED ON)
echo.
echo # Set compiler flags
echo if(MSVC)
echo     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /utf-8")
echo else()
echo     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
echo endif()
echo.
echo # Find OpenMP for parallel computing
echo find_package(OpenMP)
echo if(OpenMP_CXX_FOUND)
echo     message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
echo     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
echo else()
echo     message(WARNING "OpenMP not found - parallel features will be disabled")
echo endif()
echo.
echo # Include directories
echo include_directories(${CMAKE_SOURCE_DIR}/src)
echo include_directories(${CMAKE_SOURCE_DIR}/spdlog/include)
echo.
echo # Create library with simplified structure (no MPI modules)
echo add_library(ElmerCppMatrixAssembly
echo     # Core infrastructure modules
echo     src/core/base/SolverBase.h
echo     src/core/base/SolverRegistry.h
echo     src/core/base/Types.h
echo     
echo     # Math computation
echo     src/core/math/LinearAlgebra.h
echo     src/core/math/LinearAlgebra.cpp
echo     src/core/math/CRSMatrix.h
echo     src/core/math/CRSMatrix.cpp
echo     src/core/math/Matrix.cpp
echo     src/core/math/IterativeSolver.h
echo     src/core/math/IterativeSolver.cpp
echo     
echo     # Mesh processing
echo     src/core/mesh/Mesh.h
echo     src/core/mesh/MeshIO.h
echo     src/core/mesh/MeshIO.cpp
echo     src/core/mesh/MeshUtils.h
echo     src/core/mesh/MeshUtils.cpp
echo     src/core/mesh/ElementDescription.h
echo     src/core/mesh/ElementDescription.cpp
echo     
echo     # Utility classes
echo     src/core/utils/DefUtils.h
echo     src/core/utils/DefUtils.cpp
echo     src/core/utils/Lists.h
echo     src/core/utils/Lists.cpp
echo     src/core/utils/Integration.h
echo     src/core/utils/Integration.cpp
echo     src/core/utils/Interpolation.h
echo     src/core/utils/Interpolation.cpp
echo     src/core/utils/ShapeFunctions.h
echo     src/core/utils/ShapeFunctions.cpp
echo     src/core/utils/ElementUtils.h
echo     src/core/utils/ElementUtils.cpp
echo     src/core/utils/GaussIntegration.h
echo     src/core/utils/GaussIntegration.cpp
echo     
echo     # Logging system
echo     src/core/logging/LoggerInterface.h
echo     src/core/logging/SpdlogAdapter.h
echo     src/core/logging/LoggerFactory.h
echo     src/core/logging/LoggerFactory.cpp
echo     
echo     # Solver modules
echo     src/solvers/ElmerSolver.h
echo     src/solvers/ElmerSolver.cpp
echo     src/solvers/base/NonlinearSolver.h
echo     
echo     # Electromagnetic solvers
echo     src/solvers/electromagnetic/MagneticSolver.h
echo     src/solvers/electromagnetic/MagneticSolver.cpp
echo     src/solvers/electromagnetic/MagnetoDynamics2DSolver.h
echo     src/solvers/electromagnetic/MagnetoDynamics2DSolver.cpp
echo     src/solvers/electromagnetic/MagnetodynamicsSolver.h
echo     src/solvers/electromagnetic/ElectromagneticMaterial.h
echo     
echo     # Thermal solvers
echo     src/solvers/thermal/HeatSolver.h
echo     src/solvers/thermal/HeatSolver.cpp
echo     src/solvers/thermal/HeatTransferSolver.h
echo     src/solvers/thermal/HeatTransferSolver.cpp
echo     
echo     # I/O modules
echo     src/io/input/InputFileParser.h
echo     src/io/input/InputFileParser.cpp
echo     
echo     # Material database (temporary location)
echo     src/Material.h
echo     src/MaterialDatabase.h
echo     
echo     # Matrix assembly (temporary location)
echo     src/MatrixAssembly.h
echo     src/MatrixAssembly.cpp
echo     
echo     # Element matrix (temporary location)
echo     src/ElementMatrix.h
echo     src/ElementMatrix.cpp
echo     
echo     # Boundary conditions
echo     src/boundary/BoundaryConditions.h
echo     src/boundary/BoundaryConditions.cpp
echo )
echo.
echo # Create simple test executable
echo add_executable(test_simple
echo     test/test_simple.cpp
echo )
echo.
echo target_link_libraries(test_simple ElmerCppMatrixAssembly)
) > CMakeLists_simple.txt

copy CMakeLists_simple.txt CMakeLists.txt >nul

echo Compilation issues fixed!
echo.
echo Now run clean_and_build.bat to test compilation...
pause