#include "MeshUtils.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace elmer {

// Helper functions for volume calculations
static double calculateTetrahedronVolume(const Node& a, const Node& b, const Node& c, const Node& d) {
    // Calculate volume using scalar triple product
    double v1x = b.x - a.x;
    double v1y = b.y - a.y;
    double v1z = b.z - a.z;
    
    double v2x = c.x - a.x;
    double v2y = c.y - a.y;
    double v2z = c.z - a.z;
    
    double v3x = d.x - a.x;
    double v3y = d.y - a.y;
    double v3z = d.z - a.z;
    
    // Calculate cross product of v2 and v3
    double crossX = v2y * v3z - v2z * v3y;
    double crossY = v2z * v3x - v2x * v3z;
    double crossZ = v2x * v3y - v2y * v3x;
    
    // Calculate dot product of v1 with cross product
    double volume = std::abs(v1x * crossX + v1y * crossY + v1z * crossZ) / 6.0;
    
    return volume;
}

// Helper function for hexahedron volume calculation
static double calculateHexahedronVolume(const std::vector<Node>& nodes, const std::vector<size_t>& indices) {
    if (indices.size() != 8) {
        return 0.0;
    }
    
    // Divide hexahedron into 5 tetrahedrons and sum their volumes
    double volume = 0.0;
    
    // Tetrahedron 1: 0,1,3,4
    volume += calculateTetrahedronVolume(
        nodes[indices[0]], nodes[indices[1]], nodes[indices[3]], nodes[indices[4]]
    );
    
    // Tetrahedron 2: 1,2,3,6
    volume += calculateTetrahedronVolume(
        nodes[indices[1]], nodes[indices[2]], nodes[indices[3]], nodes[indices[6]]
    );
    
    // Tetrahedron 3: 1,3,4,6
    volume += calculateTetrahedronVolume(
        nodes[indices[1]], nodes[indices[3]], nodes[indices[4]], nodes[indices[6]]
    );
    
    // Tetrahedron 4: 1,4,5,6
    volume += calculateTetrahedronVolume(
        nodes[indices[1]], nodes[indices[4]], nodes[indices[5]], nodes[indices[6]]
    );
    
    // Tetrahedron 5: 3,4,6,7
    volume += calculateTetrahedronVolume(
        nodes[indices[3]], nodes[indices[4]], nodes[indices[6]], nodes[indices[7]]
    );
    
    return volume;
}

} // namespace elmer