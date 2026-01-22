cmake_minimum_required(VERSION 3.10)
project(Elmer)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
if(MSVC)
    # /W4: 启用所有警告
    # /utf-8: 使用UTF-8编码
    # 禁用常见无害警告:
    # /wd4100: 未引用形参警告 (C4100)
    # /wd4244: 类型转换警告 (C4244)
    # /wd4267: size_t转换警告 (C4267)
    # /wd4101: 未引用局部变量警告 (C4101)
    # /wd4189: 局部变量已初始化但不引用警告 (C4189)
    # 添加宏定义:
    # _CRT_SECURE_NO_WARNINGS: 禁用CRT安全警告
    # _SCL_SECURE_NO_WARNINGS: 禁用STL安全警告
    # NOMINMAX: 禁用min/max宏定义冲突
    # WIN32_LEAN_AND_MEAN: 减少Windows头文件包含
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /utf-8 /wd4100 /wd4244 /wd4267 /wd4101 /wd4189")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS -DNOMINMAX -DWIN32_LEAN_AND_MEAN)
    
    # 设置运行时库为静态链接 (/MT) - 确保所有库使用相同的运行时库
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    
    # 强制使用静态运行时库
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    else()
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    endif()
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

# MPI support removed - serial computation only

# Create library (MPI support removed)
add_library(Elmer
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
    
    # Sparse matrix solver framework
    src/core/math/SparseMatrix.h
    src/core/math/SparseMatrix.cpp
    src/core/math/LinearSolverInterface.h
    src/core/math/LinearSolverFactory.cpp
    src/core/math/SuperLUSolver.h
    src/core/math/SuperLUSolver.cpp
    src/core/math/SuperLUMTSolver.h
    src/core/math/SuperLUMTSolver.cpp
    src/core/math/MumpsSolver.h
    src/core/math/MumpsSolver.cpp
    
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
    # 低频电磁求解器（整合2D和3D）
    src/solvers/electromagnetic/MagnetoDynamicsSolverBase.h
    src/solvers/electromagnetic/MagnetoDynamicsSolverBase.cpp
    src/solvers/electromagnetic/MagnetoDynamics2DSolver.h
    src/solvers/electromagnetic/MagnetoDynamics2DSolver.cpp
    src/solvers/electromagnetic/MagnetoDynamics3DSolver.h
    src/solvers/electromagnetic/MagnetoDynamics3DSolver.cpp
    src/solvers/electromagnetic/MagnetoDynamicsSolverFactory.h
    src/solvers/electromagnetic/MagnetoDynamicsSolverFactory.cpp
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

# 设置头文件搜索路径
target_include_directories(Elmer PUBLIC
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

# MPI支持已移除

# 添加googletest子目录
add_subdirectory(googletest)

# Create simple test executable
add_executable(test_simple
    test/test_simple.cpp
)

target_link_libraries(test_simple Elmer)

# 设置测试程序的头文件搜索路径
target_include_directories(test_simple PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/spdlog/include
)

# Create element description test executable
add_executable(test_element_description
    test/test_element_description.cpp
)

target_link_libraries(test_element_description Elmer gtest gtest_main)

# 设置测试程序的头文件搜索路径
target_include_directories(test_element_description PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/spdlog/include
)

# Create sparse solver framework test executable
add_executable(test_sparse_solver_framework
    test/test_sparse_solver_framework.cpp
)

target_link_libraries(test_sparse_solver_framework Elmer gtest gtest_main)

# 设置测试程序的头文件搜索路径
target_include_directories(test_sparse_solver_framework PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/spdlog/include
)

# Create simple sparse solver framework test executable
add_executable(test_sparse_solver_framework_simple
    test/test_sparse_solver_framework_simple.cpp
)

target_link_libraries(test_sparse_solver_framework_simple Elmer gtest gtest_main)

# 设置测试程序的头文件搜索路径
target_include_directories(test_sparse_solver_framework_simple PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/spdlog/include
)