// test_matrix_assembly.cpp - Matrix Assembly Module Test Program

#include "../src/MatrixAssembly.h"
#include "../src/Types.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

using namespace elmer;

// Test basic matrix operations
void testBasicMatrixOperations() {
    std::cout << "Testing basic matrix operations..." << std::endl;
    
    // Create a 10x10 CRS matrix
    std::cout << "Creating CRS matrix..." << std::endl;
    auto matrix = Matrix::CreateCRS(10, 10);
    std::cout << "CRS matrix created successfully" << std::endl;
    
    auto sharedMatrix = std::shared_ptr<elmer::Matrix>(matrix.release());
    std::cout << "Shared pointer conversion successful" << std::endl;
    
    // Test setting and getting matrix elements
    MatrixAssembly::SetMatrixElement(sharedMatrix, 0, 0, 1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 1, 1, 2.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 2, 2, 3.0);
    
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 0) == 1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 1, 1) == 2.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 2, 2) == 3.0);
    
    // Test modifying matrix elements
    Real oldValue = MatrixAssembly::ChangeMatrixElement(sharedMatrix, 0, 0, 5.0);
    assert(oldValue == 1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 0) == 5.0);
    
    // Test adding to matrix elements
    MatrixAssembly::AddToMatrixElement(sharedMatrix, 0, 0, 1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 0) == 6.0);
    
    std::cout << "Basic matrix operations test passed!" << std::endl;
}

// Test matrix-vector multiplication
void testMatrixVectorMultiplication() {
    std::cout << "Testing matrix-vector multiplication..." << std::endl;
    
    // Create a 3x3 dense matrix
    auto matrix = Matrix::CreateDense(3, 3);
    auto sharedMatrix = std::shared_ptr<elmer::Matrix>(matrix.release());
    
    // Set matrix elements (identity matrix)
    MatrixAssembly::SetMatrixElement(sharedMatrix, 0, 0, 1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 1, 1, 1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 2, 2, 1.0);
    
    // Create test vectors
    std::vector<Real> x = {1.0, 2.0, 3.0};
    std::vector<Real> y(3, 0.0);
    
    // Perform matrix-vector multiplication
    auto xVector = Vector::Create(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        (*xVector)[i] = x[i];
    }
    auto yVector = Vector::Create(3);
    auto sharedXVector = std::shared_ptr<elmer::Vector>(xVector.release());
    auto sharedYVector = std::shared_ptr<elmer::Vector>(yVector.release());
    MatrixAssembly::MatrixVectorMultiply(sharedMatrix, sharedXVector, sharedYVector);
    
    // Convert result back to vector
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] = (*sharedYVector)[i];
    }
    
    // Verify results
    assert(y[0] == 1.0);
    assert(y[1] == 2.0);
    assert(y[2] == 3.0);
    
    // Test non-identity matrix
    MatrixAssembly::SetMatrixElement(sharedMatrix, 0, 1, 1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 1, 2, 1.0);
    
    y = {0.0, 0.0, 0.0};
    auto xVector2 = Vector::Create(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        (*xVector2)[i] = x[i];
    }
    auto yVector2 = Vector::Create(3);
    auto sharedXVector2 = std::shared_ptr<elmer::Vector>(xVector2.release());
    auto sharedYVector2 = std::shared_ptr<elmer::Vector>(yVector2.release());
    MatrixAssembly::MatrixVectorMultiply(sharedMatrix, sharedXVector2, sharedYVector2);
    
    // Convert result back to vector
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] = (*sharedYVector2)[i];
    }
    
    assert(y[0] == 1.0 + 2.0);  // 1*1 + 1*2
    assert(y[1] == 2.0 + 3.0);  // 1*2 + 1*3
    assert(y[2] == 3.0);        // 1*3
    
    std::cout << "Matrix-vector multiplication test passed!" << std::endl;
}

// Test element matrix assembly
void testElementMatrixAssembly() {
    std::cout << "Testing element matrix assembly..." << std::endl;
    
    // Create a 6x6 global matrix
    auto globalMatrix = Matrix::CreateCRS(6, 6);
    auto sharedGlobalMatrix = std::shared_ptr<elmer::Matrix>(globalMatrix.release());
    
    // Create a 2x2 element stiffness matrix
    std::vector<std::vector<Real>> elementMatrix = {
        {2.0, -1.0},
        {-1.0, 2.0}
    };
    
    // DOF indices (assemble to global matrix positions 1 and 3)
    std::vector<Integer> dofIndices = {1, 3};
    
    // Assemble element matrix
    MatrixAssembly::AssembleElementMatrix(sharedGlobalMatrix, dofIndices, elementMatrix);
    
    // Verify assembly results
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 1, 1) == 2.0);
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 1, 3) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 3, 1) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 3, 3) == 2.0);
    
    // Test multiple assembly (superposition)
    MatrixAssembly::AssembleElementMatrix(sharedGlobalMatrix, dofIndices, elementMatrix);
    
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 1, 1) == 4.0);
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 1, 3) == -2.0);
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 3, 1) == -2.0);
    assert(MatrixAssembly::GetMatrixElement(sharedGlobalMatrix, 3, 3) == 4.0);
    
    std::cout << "Element matrix assembly test passed!" << std::endl;
}

// Test element vector assembly
void testElementVectorAssembly() {
    std::cout << "Testing element vector assembly..." << std::endl;
    
    // Create a 6-dimensional global vector
    auto globalVector = Vector::Create(6);
    globalVector->Zero();
    auto sharedGlobalVector = std::shared_ptr<elmer::Vector>(globalVector.release());
    
    // Create a 2-dimensional element load vector
    std::vector<Real> elementVector = {1.5, 2.5};
    
    // DOF indices
    std::vector<Integer> dofIndices = {2, 4};
    
    // Assemble element vector
    MatrixAssembly::AssembleElementVector(sharedGlobalVector, dofIndices, elementVector);
    
    // Verify assembly results
    assert((*sharedGlobalVector)[2] == 1.5);
    assert((*sharedGlobalVector)[4] == 2.5);
    
    // Test multiple assembly (superposition)
    MatrixAssembly::AssembleElementVector(sharedGlobalVector, dofIndices, elementVector);
    
    assert((*sharedGlobalVector)[2] == 3.0);
    assert((*sharedGlobalVector)[4] == 5.0);
    
    std::cout << "Element vector assembly test passed!" << std::endl;
}

// Test Dirichlet boundary conditions application
void testDirichletBoundaryConditions() {
    std::cout << "Testing Dirichlet boundary conditions..." << std::endl;
    
    // Create a 3x3 matrix
    auto matrix = Matrix::CreateDense(3, 3);
    auto sharedMatrix = std::shared_ptr<elmer::Matrix>(matrix.release());
    
    // Set matrix elements (symmetric positive definite matrix)
    MatrixAssembly::SetMatrixElement(sharedMatrix, 0, 0, 4.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 0, 1, -1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 0, 2, -1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 1, 0, -1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 1, 1, 4.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 1, 2, -1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 2, 0, -1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 2, 1, -1.0);
    MatrixAssembly::SetMatrixElement(sharedMatrix, 2, 2, 4.0);
    
    // Create right-hand side vector
    std::vector<Real> rhs = {1.0, 2.0, 3.0};
    
    // Apply Dirichlet boundary condition (DOF 0 fixed to 5.0)
    MatrixAssembly::ApplyDirichletBC(sharedMatrix, rhs, 0, 5.0);
    
    // Verify boundary condition application
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 0) == 1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 1) == 0.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 2) == 0.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 1, 0) == 0.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 2, 0) == 0.0);
    assert(rhs[0] == 5.0);
    
    // Verify other rows remain unchanged
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 1, 1) == 4.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 1, 2) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 2, 1) == -1.0);
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 2, 2) == 4.0);
    
    std::cout << "Dirichlet boundary conditions test passed!" << std::endl;
}

// Test matrix properties check
void testMatrixProperties() {
    std::cout << "Testing matrix properties check..." << std::endl;
    
    // Create a symmetric matrix
    auto symmetricMatrix = Matrix::CreateDense(3, 3);
    auto sharedSymmetricMatrix = std::shared_ptr<elmer::Matrix>(symmetricMatrix.release());
    MatrixAssembly::SetMatrixElement(sharedSymmetricMatrix, 0, 0, 4.0);
    MatrixAssembly::SetMatrixElement(sharedSymmetricMatrix, 0, 1, 1.0);
    MatrixAssembly::SetMatrixElement(sharedSymmetricMatrix, 1, 0, 1.0);
    MatrixAssembly::SetMatrixElement(sharedSymmetricMatrix, 1, 1, 4.0);
    MatrixAssembly::SetMatrixElement(sharedSymmetricMatrix, 2, 2, 4.0);
    
    // Test symmetry check
    assert(MatrixAssembly::IsSymmetric(sharedSymmetricMatrix));
    
    // Create a non-symmetric matrix
    auto nonSymmetricMatrix = Matrix::CreateDense(2, 2);
    auto sharedNonSymmetricMatrix = std::shared_ptr<elmer::Matrix>(nonSymmetricMatrix.release());
    MatrixAssembly::SetMatrixElement(sharedNonSymmetricMatrix, 0, 1, 1.0);
    MatrixAssembly::SetMatrixElement(sharedNonSymmetricMatrix, 1, 0, 2.0);
    
    assert(!MatrixAssembly::IsSymmetric(sharedNonSymmetricMatrix));
    
    // Test positive definiteness check (simplified implementation)
    auto positiveDefiniteMatrix = Matrix::CreateDense(2, 2);
    auto sharedPositiveDefiniteMatrix = std::shared_ptr<elmer::Matrix>(positiveDefiniteMatrix.release());
    MatrixAssembly::SetMatrixElement(sharedPositiveDefiniteMatrix, 0, 0, 4.0);
    MatrixAssembly::SetMatrixElement(sharedPositiveDefiniteMatrix, 1, 1, 4.0);
    
    assert(MatrixAssembly::IsPositiveDefinite(sharedPositiveDefiniteMatrix));
    
    std::cout << "Matrix properties check test passed!" << std::endl;
}

// Test complex matrix operations
void testComplexMatrixOperations() {
    std::cout << "Testing complex matrix operations..." << std::endl;
    
    // Create a 4x4 matrix (for storing 2x2 complex matrix)
    auto matrix = Matrix::CreateDense(4, 4);
    auto sharedMatrix = std::shared_ptr<elmer::Matrix>(matrix.release());
    
    // Add complex matrix element (real part=2.0, imaginary part=1.0)
    MatrixAssembly::AddToComplexMatrixElement(sharedMatrix, 0, 0, 2.0, 1.0);
    
    // Verify complex matrix storage format
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 0) == 2.0);  // Real part
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 0, 1) == -1.0); // -Imaginary part
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 1, 0) == 1.0);  // Imaginary part
    assert(MatrixAssembly::GetMatrixElement(sharedMatrix, 1, 1) == 2.0);  // Real part
    
    std::cout << "Complex matrix operations test passed!" << std::endl;
}

// Main test function
int main() {
    std::cout << "Starting matrix assembly module tests..." << std::endl;
    std::cout << "=======================================" << std::endl;
    
    try {
        testBasicMatrixOperations();
        testMatrixVectorMultiplication();
        testElementMatrixAssembly();
        testElementVectorAssembly();
        testDirichletBoundaryConditions();
        testMatrixProperties();
        testComplexMatrixOperations();
        
        std::cout << "=======================================" << std::endl;
        std::cout << "All matrix assembly module tests passed!" << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}