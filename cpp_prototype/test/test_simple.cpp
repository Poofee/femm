#include <iostream>
#include "Mesh.h"
#include "MeshUtils.h"

int main() {
    std::cout << "Testing basic C++ compilation..." << std::endl;
    
    // Test basic functionality
    auto mesh = elmer::MeshUtils::allocateMesh("test");
    std::cout << "Mesh created successfully!" << std::endl;
    
    return 0;
}