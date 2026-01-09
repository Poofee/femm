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
)

# Create test executables
add_executable(test_shape_functions
    test/test_shape_functions.cpp
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

add_executable(test_integration
    test/test_integration.cpp
)

add_executable(test_element_utils
    test/test_element_utils.cpp
)

add_executable(test_crs_matrix
    test/test_crs_matrix.cpp
)

# Set include directories for test executables
target_include_directories(test_crs_matrix PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Link test executables with library
target_link_libraries(test_shape_functions ElmerCppMatrixAssembly)
target_link_libraries(test_matrix_assembly ElmerCppMatrixAssembly)
target_link_libraries(test_linear_algebra ElmerCppMatrixAssembly)
target_link_libraries(test_def_utils ElmerCppMatrixAssembly)
target_link_libraries(test_simple ElmerCppMatrixAssembly)
target_link_libraries(simple_test ElmerCppMatrixAssembly)
target_link_libraries(test_integration ElmerCppMatrixAssembly)
target_link_libraries(test_element_utils ElmerCppMatrixAssembly)
target_link_libraries(test_crs_matrix ElmerCppMatrixAssembly)

# Enable testing if CTest is available
if(BUILD_TESTING)
    enable_testing()
    add_test(NAME ShapeFunctionsTest COMMAND test_shape_functions)
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