// main.cpp - Elmer FEM C++ Prototype Main Program
// Demonstrates the integration of matrix, solver, and mesh functionality

#include "src/Types.h"
#include "src/Matrix.h"
#include "src/IterativeSolver.h"
#include "src/Mesh.h"
#include "src/MeshIO.h"
#include "src/MeshUtils.h"
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

using namespace elmer;

// Function to create a test matrix (Poisson equation on a simple grid)
std::unique_ptr<Matrix> CreateTestMatrix(Integer size) {
    auto matrix = Matrix::CreateCRS(size, size);
    
    // Create a simple 1D Poisson matrix: -u'' = f
    for (Integer i = 0; i < size; ++i) {
        // Diagonal element
        matrix->SetElement(i, i, 2.0);
        
        // Off-diagonal elements
        if (i > 0) {
            matrix->SetElement(i, i - 1, -1.0);
        }
        if (i < size - 1) {
            matrix->SetElement(i, i + 1, -1.0);
        }
    }
    
    return matrix;
}

// Function to create a test right-hand side vector
std::unique_ptr<Vector> CreateTestRHS(Integer size) {
    auto rhs = Vector::Create(size);
    
    // Simple test: f(x) = 1 for all x
    for (Integer i = 0; i < size; ++i) {
        (*rhs)[i] = 1.0;
    }
    
    return rhs;
}

// Test matrix operations
void TestMatrixOperations() {
    std::cout << "=== Testing Matrix Operations ===" << std::endl;
    
    // Create a 3x3 test matrix
    auto matrix = Matrix::CreateCRS(3, 3);
    
    // Set matrix elements
    matrix->SetElement(0, 0, 4.0);
    matrix->SetElement(0, 1, -1.0);
    matrix->SetElement(1, 0, -1.0);
    matrix->SetElement(1, 1, 4.0);
    matrix->SetElement(1, 2, -1.0);
    matrix->SetElement(2, 1, -1.0);
    matrix->SetElement(2, 2, 4.0);
    
    // Create test vectors
    auto x = Vector::Create(3);
    auto y = Vector::Create(3);
    
    // Set x vector
    (*x)[0] = 1.0;
    (*x)[1] = 2.0;
    (*x)[2] = 3.0;
    
    // Test matrix-vector multiplication
    matrix->Multiply(*x, *y);
    
    std::cout << "Matrix-vector multiplication result:" << std::endl;
    for (Integer i = 0; i < 3; ++i) {
        std::cout << "y[" << i << "] = " << (*y)[i] << std::endl;
    }
    
    std::cout << std::endl;
}

// Test iterative solvers
void TestIterativeSolvers() {
    std::cout << "=== Testing Iterative Solvers ===" << std::endl;
    
    Integer problem_size = 100;
    
    // Create test problem
    auto A = CreateTestMatrix(problem_size);
    auto b = CreateTestRHS(problem_size);
    auto x = Vector::Create(problem_size);
    x->Zero();  // Initial guess
    
    // Test different solvers
    std::vector<std::string> solver_names = {"CG", "BiCGStab", "GMRES"};
    
    for (const auto& name : solver_names) {
        std::cout << "Testing " << name << " solver:" << std::endl;
        
        auto solver = SolverFactory::CreateSolver(name);
        solver->SetMaxIterations(1000);
        solver->SetTolerance(1e-8);
        
        // Reset initial guess
        x->Zero();
        
        bool success = solver->Solve(*A, *x, *b);
        
        std::cout << "  Converged: " << (success ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << solver->GetIterations() << std::endl;
        std::cout << "  Final residual: " << solver->GetResidual() << std::endl;
        
        // Verify solution by computing residual
        auto r = Vector::Create(problem_size);
        A->Multiply(*x, *r);
        
        Real residual = 0.0;
        for (Integer i = 0; i < problem_size; ++i) {
            Real diff = (*b)[i] - (*r)[i];
            residual += diff * diff;
        }
        residual = std::sqrt(residual);
        
        std::cout << "  Computed residual: " << residual << std::endl;
        std::cout << std::endl;
    }
}

// Test mesh functionality
void TestMeshFunctionality() {
    std::cout << "=== Testing Mesh Functionality ===" << std::endl;
    
    // Create a simple 2D mesh
    Mesh mesh("TestMesh");
    
    // Add nodes
    mesh.getNodes().addNode(0.0, 0.0, 0.0);
    mesh.getNodes().addNode(1.0, 0.0, 0.0);
    mesh.getNodes().addNode(0.0, 1.0, 0.0);
    mesh.getNodes().addNode(1.0, 1.0, 0.0);
    
    // Add bulk elements (triangles)
    Element tri1(ElementType::TETRAHEDRON, 1);
    tri1.addNodeIndex(0);
    tri1.addNodeIndex(1);
    tri1.addNodeIndex(2);
    mesh.addBulkElement(tri1);
    
    Element tri2(ElementType::TETRAHEDRON, 1);
    tri2.addNodeIndex(1);
    tri2.addNodeIndex(3);
    tri2.addNodeIndex(2);
    mesh.addBulkElement(tri2);
    
    // Add boundary elements
    Element boundary1(ElementType::LINEAR, 1);
    boundary1.setBoundaryId(1);
    boundary1.addNodeIndex(0);
    boundary1.addNodeIndex(1);
    mesh.addBoundaryElement(boundary1);
    
    // Test mesh properties
    std::cout << "Mesh name: " << mesh.getName() << std::endl;
    std::cout << "Number of nodes: " << mesh.numberOfNodes() << std::endl;
    std::cout << "Number of bulk elements: " << mesh.numberOfBulkElements() << std::endl;
    std::cout << "Number of boundary elements: " << mesh.numberOfBoundaryElements() << std::endl;
    std::cout << "Total elements: " << mesh.totalElements() << std::endl;
    
    // Validate mesh
    std::cout << "Mesh validation: " << (mesh.validate() ? "Passed" : "Failed") << std::endl;
    
    std::cout << std::endl;
}

// Test Eigen integration
void TestEigenIntegration() {
    std::cout << "=== Testing Eigen Integration ===" << std::endl;
    
    // Create a small test matrix
    auto matrix = CreateTestMatrix(5);
    
    // Convert to Eigen format
    try {
        Eigen::SparseMatrix<Real> eigen_matrix = ConvertToEigen(*matrix);
        
        std::cout << "Eigen matrix conversion successful" << std::endl;
        std::cout << "Matrix size: " << eigen_matrix.rows() << " x " << eigen_matrix.cols() << std::endl;
        std::cout << "Non-zero elements: " << eigen_matrix.nonZeros() << std::endl;
        
        // Test Eigen solver
        Eigen::VectorXd b(5);
        b << 1.0, 1.0, 1.0, 1.0, 1.0;
        
        Eigen::ConjugateGradient<Eigen::SparseMatrix<Real>> cg;
        cg.compute(eigen_matrix);
        Eigen::VectorXd x = cg.solve(b);
        
        std::cout << "Eigen CG solver iterations: " << cg.iterations() << std::endl;
        std::cout << "Eigen CG solver error: " << cg.error() << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "Eigen integration test failed: " << e.what() << std::endl;
    }
    
    std::cout << std::endl;
}

int main() {
    std::cout << "Elmer FEM C++ Prototype Test Program" << std::endl;
    std::cout << "====================================" << std::endl;
    std::cout << std::endl;
    
    try {
        TestMatrixOperations();
        TestIterativeSolvers();
        TestMeshFunctionality();
        TestEigenIntegration();
        
        std::cout << "All tests completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}