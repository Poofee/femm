# Elmer FEM C++ Shape Functions System

This project implements a C++ shape function system for finite element analysis, based on Elmer FEM's Fortran implementation. The system provides shape functions and their derivatives for various element types used in finite element analysis.

## Features

- **Multiple Element Types**: Support for 1D, 2D, and 3D elements
- **Linear and Quadratic Elements**: Both linear and quadratic shape functions
- **Isoparametric Mapping**: Coordinate transformation between natural and physical coordinates
- **Jacobian Computation**: Automatic computation of Jacobian matrices and determinants
- **Global Derivatives**: Transformation of derivatives from natural to global coordinates

## Supported Element Types

- **1D Elements**: Linear bar elements (2 nodes) and quadratic bar elements (3 nodes)
- **2D Elements**: 
  - Quadrilateral elements (4 nodes linear, 8 nodes quadratic)
  - Triangular elements (3 nodes linear, 6 nodes quadratic)
- **3D Elements**:
  - Hexahedral elements (8 nodes linear, 20 nodes quadratic)
  - Tetrahedral elements (4 nodes linear, 10 nodes quadratic)

## Project Structure

```
cpp_prototype/
├── src/
│   ├── ElmerCpp.h          # Basic type definitions
│   ├── ShapeFunctions.h    # Shape function class interface
│   └── ShapeFunctions.cpp  # Shape function implementation
├── test/
│   └── test_shape_functions.cpp  # Test program
├── CMakeLists.txt          # Build configuration
└── README.md              # This documentation
```

## Usage

### Basic Usage

```cpp
#include "src/ShapeFunctions.h"

using namespace elmer;

// Create element nodes
std::vector<Node> nodes = {
    {0.0, 0.0, 0.0},  // Node 0
    {2.0, 0.0, 0.0},  // Node 1
    {2.0, 1.0, 0.0},  // Node 2
    {0.0, 1.0, 0.0}   // Node 3
};

// Compute shape functions at natural coordinates
auto result = ShapeFunctions::computeShapeFunctions(
    ElementType::QUADRILATERAL, nodes, 0.0, 0.0, 0.0);

// Access results
std::cout << "Shape function values: ";
for (double val : result.values) {
    std::cout << val << " ";
}
std::cout << "\nJacobian determinant: " << result.detJ << std::endl;
```

### Individual Shape Function Evaluation

```cpp
// Evaluate linear quadrilateral shape functions
auto N = ShapeFunctions::linearQuadrilateral(0.5, 0.5, 4);

// Evaluate derivatives
auto [dNdxi, dNdeta] = ShapeFunctions::linearQuadrilateralDerivatives(0.5, 0.5, 4);
```

## Building the Project

### Prerequisites

- CMake 3.10 or higher
- C++17 compatible compiler

### Build Steps

```bash
# Create build directory
mkdir build
cd build

# Configure with CMake
cmake ..

# Build the project
cmake --build .

# Run tests
./Debug/test_shape_functions.exe
```

## API Reference

### ShapeFunctions Class

#### Main Methods

- `computeShapeFunctions()`: Compute complete shape function results including global derivatives
- `linear1D()`, `linearQuadrilateral()`, `linearHexahedron()`, etc.: Evaluate shape function values
- `linear1DDerivatives()`, `linearQuadrilateralDerivatives()`, etc.: Evaluate shape function derivatives

#### Utility Methods

- `computeJacobianMatrix()`: Compute Jacobian matrix for coordinate transformation
- `computeInverseJacobian()`: Compute inverse of Jacobian matrix
- `transformDerivatives()`: Transform derivatives from natural to global coordinates
- `computeJacobianDeterminant()`: Compute determinant of Jacobian matrix
- `getNumberOfNodes()`: Get number of nodes for given element type and order
- `isSupported()`: Check if element type is supported
- `getDimension()`: Get dimension of element type

### ShapeResult Structure

Contains all computed shape function information:

- `values`: Shape function values at evaluation point
- `dNdxi`, `dNdeta`, `dNdzeta`: Derivatives in natural coordinates
- `dNdx`, `dNdy`, `dNdz`: Derivatives in global coordinates
- `detJ`: Jacobian determinant

## Mathematical Background

### Shape Functions

Shape functions are used in finite element analysis to interpolate field variables within elements. They satisfy the partition of unity property:

\[ \sum_{i=1}^{n} N_i(\xi, \eta, \zeta) = 1 \]

### Isoparametric Mapping

The mapping between natural coordinates (ξ, η, ζ) and physical coordinates (x, y, z) is given by:

\[ x = \sum_{i=1}^{n} N_i(\xi, \eta, \zeta) x_i \]

### Jacobian Matrix

The Jacobian matrix relates derivatives in natural and physical coordinates:

\[ J = \begin{bmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} & \frac{\partial z}{\partial \xi} \\
\frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta} & \frac{\partial z}{\partial \eta} \\
\frac{\partial x}{\partial \zeta} & \frac{\partial y}{\partial \zeta} & \frac{\partial z}{\partial \zeta}
\end{bmatrix} \]

### Derivative Transformation

Derivatives in global coordinates are computed using the inverse Jacobian:

\[ \frac{\partial N_i}{\partial x} = J^{-1} \frac{\partial N_i}{\partial \xi} \]

## Testing

The test program verifies:

1. **Partition of Unity**: Sum of shape functions equals 1.0
2. **Jacobian Determinant**: Positive values for valid elements
3. **Coordinate Mapping**: Correct transformation between coordinate systems
4. **Derivative Computation**: Accurate computation of shape function derivatives

## Implementation Notes

- The implementation follows Elmer FEM's mathematical formulation
- All computations use double precision floating-point arithmetic
- Error handling includes checks for invalid Jacobian determinants
- The code is designed to be easily extensible for additional element types

## Future Enhancements

- Support for higher-order elements (cubic, etc.)
- Additional element types (wedge, pyramid, etc.)
- Integration with numerical quadrature
- Parallel computation support
- GPU acceleration

## License

This project is part of the Elmer FEM ecosystem and follows the same licensing terms.

## Contributing

Contributions are welcome! Please ensure that:

1. Code follows the existing style and conventions
2. All tests pass
3. New features include appropriate test cases
4. Documentation is updated accordingly