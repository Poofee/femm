cmake_minimum_required(VERSION 3.10)
project(ElmerCppShapeFunctions)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# Create library
add_library(ElmerCppShapeFunctions
    src/ElmerCpp.h
    src/ShapeFunctions.h
    src/ShapeFunctions.cpp
)

# Create test executable
add_executable(test_shape_functions
    test/test_shape_functions.cpp
)

# Link test executable with library
target_link_libraries(test_shape_functions ElmerCppShapeFunctions)

# Enable testing if CTest is available
if(BUILD_TESTING)
    enable_testing()
    add_test(NAME ShapeFunctionsTest COMMAND test_shape_functions)
endif()

# Installation (optional)
install(TARGETS ElmerCppShapeFunctions
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES src/ElmerCpp.h src/ShapeFunctions.h
    DESTINATION include/elmer
)