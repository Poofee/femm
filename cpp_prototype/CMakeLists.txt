cmake_minimum_required(VERSION 3.10)
project(ElmerCppMatrixAssembly)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /utf-8")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

# Find OpenMP for parallel computing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP not found - parallel features will be disabled")
endif()

# Create library with simplified structure (no MPI modules)
add_library(ElmerCppMatrixAssembly
    # Core infrastructure modules
    src/core/base/SolverBase.h
    src/core/base/SolverRegistry.h
    src/core/base/Types.h
    
    # Math computation
    src/core/math/LinearAlgebra.h
    src/core/math/LinearAlgebra.cpp
    src/core/math/CRSMatrix.h
    src/core/math/CRSMatrix.cpp
    src/core/math/Matrix.cpp
    src/core/math/IterativeSolver.h
    src/core/math/IterativeSolver.cpp
    
    # Mesh processing
    src/core/mesh/Mesh.h
    src/core/mesh/MeshIO.h
    src/core/mesh/MeshIO.cpp
    src/core/mesh/MeshUtils.h
    src/core/mesh/MeshUtils.cpp
    src/core/mesh/ElementDescription.h
    src/core/mesh/ElementDescription.cpp
    
    # Utility classes
    src/core/utils/DefUtils.h
    src/core/utils/DefUtils.cpp
    src/core/utils/Lists.h
    src/core/utils/Lists.cpp
    src/core/utils/Integration.h
    src/core/utils/Integration.cpp
    src/core/utils/Interpolation.h
    src/core/utils/Interpolation.cpp
    src/core/utils/ShapeFunctions.h
    src/core/utils/ShapeFunctions.cpp
    src/core/utils/ElementUtils.h
    src/core/utils/ElementUtils.cpp
    src/core/utils/GaussIntegration.h
    src/core/utils/GaussIntegration.cpp
    
    # Logging system
    src/core/logging/LoggerInterface.h
    src/core/logging/SpdlogAdapter.h
    src/core/logging/LoggerFactory.h
    src/core/logging/LoggerFactory.cpp
    
    # Solver modules
    src/solvers/ElmerSolver.h
    src/solvers/ElmerSolver.cpp
    src/solvers/base/NonlinearSolver.h
    
    # Electromagnetic solvers
    src/solvers/electromagnetic/MagneticSolver.h
    src/solvers/electromagnetic/MagneticSolver.cpp
    src/solvers/electromagnetic/ElectromagneticMaterial.h
    # TODO: 重新移植WhitneyAVSolver
    # src/solvers/electromagnetic/MagnetoDynamics2DSolver.h
    src/solvers/electromagnetic/MagnetodynamicsSolver.h
    
    # Thermal solvers
    src/solvers/thermal/HeatSolver.h
    src/solvers/thermal/HeatSolver.cpp
    
    # I/O modules
    src/io/input/InputFileParser.h
    src/io/input/InputFileParser.cpp
    
    # Material database (temporary location)
    src/Material.h
    src/MaterialDatabase.h
    
    # Matrix assembly (temporary location)
    src/MatrixAssembly.h
    src/MatrixAssembly.cpp
    
    # Element matrix (temporary location)
    src/ElementMatrix.h
    src/ElementMatrix.cpp
    
    # Boundary conditions
    src/boundary/BoundaryConditions.h
    src/boundary/BoundaryConditions.cpp
)

# 设置头文件搜索路径 - 这是关键修复
target_include_directories(ElmerCppMatrixAssembly PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/core/base
    ${CMAKE_SOURCE_DIR}/src/core/math
    ${CMAKE_SOURCE_DIR}/src/core/mesh
    ${CMAKE_SOURCE_DIR}/src/core/utils
    ${CMAKE_SOURCE_DIR}/src/core/logging
    ${CMAKE_SOURCE_DIR}/src/solvers
    ${CMAKE_SOURCE_DIR}/src/solvers/base
    ${CMAKE_SOURCE_DIR}/src/solvers/electromagnetic
    ${CMAKE_SOURCE_DIR}/src/solvers/thermal
    ${CMAKE_SOURCE_DIR}/src/io/input
    ${CMAKE_SOURCE_DIR}/src/boundary
    ${CMAKE_SOURCE_DIR}/spdlog/include
)

# 添加googletest子目录
add_subdirectory(googletest)

# Create simple test executable
add_executable(test_simple
    test/test_simple.cpp
)

target_link_libraries(test_simple ElmerCppMatrixAssembly)

# 设置测试程序的头文件搜索路径
target_include_directories(test_simple PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/spdlog/include
)