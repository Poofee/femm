#pragma once

#include <memory>
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>

namespace elmer {

/**
 * @brief Node class representing a point in 3D space
 */
class Node {
public:
    double x, y, z;
    
    Node(double x_val = 0.0, double y_val = 0.0, double z_val = 0.0)
        : x(x_val), y(y_val), z(z_val) {}
    
    // Copy constructor
    Node(const Node& other) = default;
    
    // Assignment operator
    Node& operator=(const Node& other) = default;
    
    // Comparison operators
    bool operator==(const Node& other) const {
        return std::abs(x - other.x) < 1e-15 && 
               std::abs(y - other.y) < 1e-15 && 
               std::abs(z - other.z) < 1e-15;
    }
    
    bool operator!=(const Node& other) const {
        return !(*this == other);
    }
    
    // Distance calculation
    double distance(const Node& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    // Vector operations
    Node operator+(const Node& other) const {
        return Node(x + other.x, y + other.y, z + other.z);
    }
    
    Node operator-(const Node& other) const {
        return Node(x - other.x, y - other.y, z - other.z);
    }
    
    Node operator*(double scalar) const {
        return Node(x * scalar, y * scalar, z * scalar);
    }
    
    Node operator/(double scalar) const {
        if (std::abs(scalar) < 1e-15) {
            throw std::runtime_error("Division by zero");
        }
        return Node(x / scalar, y / scalar, z / scalar);
    }
    
    // Dot product
    double dot(const Node& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    // Cross product
    Node cross(const Node& other) const {
        return Node(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
    
    // Magnitude
    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    
    // Normalize
    Node normalize() const {
        double mag = magnitude();
        if (mag < 1e-15) {
            throw std::runtime_error("Cannot normalize zero vector");
        }
        return *this / mag;
    }
};

/**
 * @brief Element type enumeration
 */
enum class ElementType {
    LINEAR = 0,
    QUADRATIC = 1,
    TETRAHEDRON = 2,
    HEXAHEDRON = 3,
    TRIANGLE = 4,
    QUADRILATERAL = 5
};

} // namespace elmer