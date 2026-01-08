# Elmer FEM C++移植版CMake配置文件
cmake_minimum_required(VERSION 3.15)
project(ElmerCpp VERSION 0.1.0 LANGUAGES C CXX Fortran)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 编译器选项
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Od /Zi")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
endif()

# 查找依赖包
find_package(MPI REQUIRED)

# 查找Eigen库（用于线性代数计算）
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 查找Boost库（可选，用于高级功能）
find_package(Boost COMPONENTS system filesystem program_options)

# 查找VTK库（可选，用于可视化）
find_package(VTK)

# 包含目录
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Eigen3_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_PATH}
)

# 添加Elmer FEM Fortran库路径
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/modules
)

# 定义目标
add_library(ElmerCppCore
    src/ElmerCpp.cpp
    src/Matrix.cpp
    src/Vector.cpp
    src/Mesh.cpp
    src/Solver.cpp
    src/Model.cpp
)

# 链接库
target_link_libraries(ElmerCppCore
    PUBLIC
        Eigen3::Eigen
        MPI::MPI_CXX
)

if(Boost_FOUND)
    target_link_libraries(ElmerCppCore PUBLIC Boost::boost)
endif()

if(VTK_FOUND)
    target_link_libraries(ElmerCppCore PUBLIC ${VTK_LIBRARIES})
endif()

# Fortran互操作性库
add_library(ElmerFortranInterface
    fortran/ElmerFortranWrapper.cpp
)

# 链接Fortran库（需要先编译Elmer FEM）
target_link_libraries(ElmerFortranInterface
    # 这里需要链接编译后的Elmer FEM Fortran库
    # -lelmerfem 或类似
)

# 主可执行文件
add_executable(elmer_cpp
    src/main.cpp
    src/ElmerSolver.cpp
    src/Matrix.cpp
    src/IterativeSolver.cpp
    src/LinearAlgebra.cpp
)

# 测试可执行文件
add_executable(test_linear_algebra
    test/test_linear_algebra.cpp
    src/LinearAlgebra.cpp
)

target_link_libraries(elmer_cpp
    ElmerCppCore
    ElmerFortranInterface
    ${EIGEN3_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

target_link_libraries(test_linear_algebra
    ${EIGEN3_LIBRARIES}
)

# 启用测试
enable_testing()
add_test(NAME LinearAlgebraTest COMMAND test_linear_algebra)

# 测试目标
if(BUILD_TESTING)
    enable_testing()
    
    add_executable(test_matrix
        tests/test_matrix.cpp
    )
    
    target_link_libraries(test_matrix
        ElmerCppCore
    )
    
    add_test(NAME MatrixTest COMMAND test_matrix)
endif()

# 安装配置
install(TARGETS elmer_cpp
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES include/ElmerCpp.h
    DESTINATION include/elmer
)

# 包配置
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ElmerCppConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT ElmerCppTargets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/ElmerCppTargets.cmake
)

configure_file(cmake/ElmerCppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ElmerCppConfig.cmake
    @ONLY
)