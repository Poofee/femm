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

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# Create library
add_library(ElmerCppMatrixAssembly
    src/ElmerCpp.h
    src/ShapeFunctions.h
    src/ShapeFunctions.cpp
    src/MatrixAssembly.h
    src/MatrixAssembly.cpp
    src/Matrix.cpp
    src/CRSMatrix.h
    src/CRSMatrix.cpp
    src/IterativeSolver.h
    src/IterativeSolver.cpp
    src/BoundaryConditions.h
    src/BoundaryConditions.cpp
    src/LinearAlgebra.h
    src/LinearAlgebra.cpp
    src/Mesh.h
    src/MeshUtils.h
    src/MeshUtils.cpp
    src/DefUtils.h
    src/DefUtils.cpp
    src/Integration.h
    src/Integration.cpp
    src/ElementUtils.h
    src/ElementUtils.cpp
    src/Lists.h
    src/Lists.cpp
    src/Interpolation.h
    src/Interpolation.cpp
    src/MagnetoDynamics2DSolver.h
    src/MagnetoDynamics2DSolver.cpp
    src/MagnetodynamicsSolver.h
    src/ElectromagneticMaterial.h
    src/MaterialDatabase.h
)

# Create test executables
add_executable(test_mesh
    test/test_mesh.cpp
)

add_executable(test_integration
    test/test_integration.cpp
)

add_executable(test_shapefunctions
    test/test_shapefunctions.cpp
)

add_executable(test_matrix_assembly
    test/test_matrix_assembly.cpp
)

add_executable(test_linear_algebra
    test/test_linear_algebra.cpp
)

add_executable(test_def_utils
    test/test_def_utils.cpp
)

add_executable(test_simple
    test/test_simple.cpp
)

add_executable(simple_test
    test/simple_test.cpp
)

add_executable(test_element_utils
    test/test_element_utils.cpp
)

add_executable(test_crs_matrix
    test/test_crs_matrix.cpp
)

add_executable(test_lists
    test/test_lists.cpp
)

add_executable(test_interpolation
    test/test_interpolation.cpp
)

add_executable(test_boundary_conditions_priority
    test/test_boundary_conditions_priority.cpp
)

add_executable(test_magnetodynamics2d_solver
    test/test_magnetodynamics2d_solver.cpp
)

add_executable(test_nonlinear_materials
    test/test_nonlinear_materials.cpp
)

add_executable(benchmark_magnetodynamics2d
    test/benchmark_magnetodynamics2d.cpp
)

target_include_directories(test_mesh PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_integration PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_shapefunctions PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_crs_matrix PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_lists PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_interpolation PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_boundary_conditions_priority PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_magnetodynamics2d_solver PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(test_nonlinear_materials PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(benchmark_magnetodynamics2d PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Link test executables with library
target_link_libraries(test_mesh ElmerCppMatrixAssembly)
target_link_libraries(test_integration ElmerCppMatrixAssembly)
target_link_libraries(test_shapefunctions ElmerCppMatrixAssembly)
target_link_libraries(test_matrix_assembly ElmerCppMatrixAssembly)
target_link_libraries(test_linear_algebra ElmerCppMatrixAssembly)
target_link_libraries(test_def_utils ElmerCppMatrixAssembly)
target_link_libraries(test_simple ElmerCppMatrixAssembly)
target_link_libraries(simple_test ElmerCppMatrixAssembly)
target_link_libraries(test_element_utils ElmerCppMatrixAssembly)
target_link_libraries(test_crs_matrix ElmerCppMatrixAssembly)
target_link_libraries(test_lists ElmerCppMatrixAssembly)
target_link_libraries(test_interpolation ElmerCppMatrixAssembly)
target_link_libraries(test_boundary_conditions_priority ElmerCppMatrixAssembly)
target_link_libraries(test_magnetodynamics2d_solver ElmerCppMatrixAssembly)
target_link_libraries(test_nonlinear_materials ElmerCppMatrixAssembly)
target_link_libraries(benchmark_magnetodynamics2d ElmerCppMatrixAssembly)

target_link_libraries(test_matrix_assembly ElmerCppMatrixAssembly)
target_link_libraries(test_linear_algebra ElmerCppMatrixAssembly)
target_link_libraries(test_def_utils ElmerCppMatrixAssembly)
target_link_libraries(test_simple ElmerCppMatrixAssembly)
target_link_libraries(simple_test ElmerCppMatrixAssembly)
target_link_libraries(test_integration ElmerCppMatrixAssembly)
target_link_libraries(test_element_utils ElmerCppMatrixAssembly)
target_link_libraries(test_crs_matrix ElmerCppMatrixAssembly)
target_link_libraries(test_lists ElmerCppMatrixAssembly)
target_link_libraries(test_interpolation ElmerCppMatrixAssembly)
target_link_libraries(test_boundary_conditions_priority ElmerCppMatrixAssembly)
target_link_libraries(test_magnetodynamics2d_solver ElmerCppMatrixAssembly)

# Enable testing if CTest is available
if(BUILD_TESTING)
    enable_testing()
    add_test(NAME ShapeFunctionsTest COMMAND test_shape_functions)
    add_test(NAME InterpolationTest COMMAND test_interpolation)
endif()

# Installation (optional)
install(TARGETS ElmerCppMatrixAssembly
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES src/ElmerCpp.h src/ShapeFunctions.h src/MatrixAssembly.h
    DESTINATION include/elmer
)