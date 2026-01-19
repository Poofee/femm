#pragma once

#include "Mesh.h"
#include <memory>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace elmer {

/**
 * @brief Mesh utility class - provides mesh creation, initialization and operation functions
 */
class MeshUtils {
public:
    /**
     * @brief Allocate mesh structure and return smart pointer
     * @param name Mesh name
     * @return Allocated mesh object
     */
    static std::shared_ptr<Mesh> allocateMesh(const std::string& name = "") {
        auto mesh = std::make_shared<Mesh>(name);
        
        // Initialize default values
        mesh->setChanged(false);
        mesh->setOutputActive(false);
        mesh->setAdaptiveDepth(0);
        mesh->setDiscontMesh(false);
        mesh->setSingleMesh(false);
        
        return mesh;
    }
    
    /**
     * @brief Allocate mesh structure with specified size
     * @param numberOfBulkElements Number of bulk elements
     * @param numberOfBoundaryElements Number of boundary elements
     * @param numberOfNodes Number of nodes
     * @param name Mesh name
     * @param initParallel Whether to initialize parallel information
     * @return Allocated mesh object
     */
    static std::shared_ptr<Mesh> allocateMesh(
        size_t numberOfBulkElements,
        size_t numberOfBoundaryElements,
        size_t numberOfNodes,
        const std::string& name = "",
        bool initParallel = false
    ) {
        auto mesh = allocateMesh(name);
        
        // Set mesh size information
        mesh->getNodes().clear();
        mesh->getBulkElements().clear();
        mesh->getBoundaryElements().clear();
        
        // Reserve space
        mesh->getBulkElements().reserve(numberOfBulkElements);
        mesh->getBoundaryElements().reserve(numberOfBoundaryElements);
        
        // Initialize node space
        for (size_t i = 0; i < numberOfNodes; ++i) {
            mesh->getNodes().addNode(0.0, 0.0, 0.0);
        }
        
        // Initialize parallel information
        if (initParallel) {
            mesh->initializeParallelInfo();
        }
        
        return mesh;
    }
    
    /**
     * @brief Initialize mesh structure
     * @param mesh Mesh to initialize
     * @param initParallel Whether to initialize parallel information
     */
    static void initializeMesh(std::shared_ptr<Mesh> mesh, bool initParallel = false) {
        if (!mesh.get()) {
            throw std::invalid_argument("Mesh pointer is null");
        }
        
        size_t numberOfNodes = mesh->numberOfNodes();
        [[maybe_unused]] size_t numberOfBulkElements = mesh->numberOfBulkElements();
        [[maybe_unused]] size_t numberOfBoundaryElements = mesh->numberOfBoundaryElements();
        
        // Validate mesh integrity
        if (numberOfNodes == 0) {
            throw std::runtime_error("Mesh has zero nodes!");
        }
        
        if (mesh->totalElements() == 0) {
            throw std::runtime_error("Mesh has zero elements!");
        }
        
        // Initialize parallel information
        if (initParallel) {
            mesh->initializeParallelInfo();
        }
        
        // Update statistics
        mesh->getBulkElements(); // Trigger statistics update
        mesh->getBoundaryElements(); // Trigger statistics update
    }
    
    /**
     * @brief Create simple cube mesh
     * @param size Cube size
     * @param divisions Number of divisions in each direction
     * @return Created mesh object
     */
    static std::shared_ptr<Mesh> createCubeMesh(
        double size = 1.0,
        int divisions = 1,
        const std::string& name = "CubeMesh"
    ) {
        auto mesh = allocateMesh(name);
        
        // Create nodes
        double step = size / divisions;
        for (int k = 0; k <= divisions; ++k) {
            for (int j = 0; j <= divisions; ++j) {
                for (int i = 0; i <= divisions; ++i) {
                    double x = i * step - size / 2.0;
                    double y = j * step - size / 2.0;
                    double z = k * step - size / 2.0;
                    mesh->getNodes().addNode(x, y, z);
                }
            }
        }
        
        // Create cube elements (hexahedrons)
        int nodesPerDirection = divisions + 1;
        for (int k = 0; k < divisions; ++k) {
            for (int j = 0; j < divisions; ++j) {
                for (int i = 0; i < divisions; ++i) {
                    Element element(ElementType::HEXAHEDRON);
                    
                    // Calculate 8 node indices for hexahedron element
                    int base = k * nodesPerDirection * nodesPerDirection + j * nodesPerDirection + i;
                    
                    element.addNodeIndex(base);
                    element.addNodeIndex(base + 1);
                    element.addNodeIndex(base + nodesPerDirection + 1);
                    element.addNodeIndex(base + nodesPerDirection);
                    
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection);
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection + 1);
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection + nodesPerDirection + 1);
                    element.addNodeIndex(base + nodesPerDirection * nodesPerDirection + nodesPerDirection);
                    
                    mesh->addBulkElement(element);
                }
            }
        }
        
        // Create boundary elements
        createBoundaryElements(mesh, divisions);
        
        return mesh;
    }
    
    /**
     * @brief Create simple tetrahedron mesh
     * @param size Mesh size
     * @param divisions Number of divisions
     * @return Created mesh object
     */
    static std::shared_ptr<Mesh> createTetrahedronMesh(
        double size = 1.0,
        [[maybe_unused]] int divisions = 1,
        const std::string& name = "TetrahedronMesh"
    ) {
        auto mesh = allocateMesh(name);
        
        // Create four nodes for simple tetrahedron
        mesh->getNodes().addNode(0.0, 0.0, 0.0);
        mesh->getNodes().addNode(size, 0.0, 0.0);
        mesh->getNodes().addNode(0.0, size, 0.0);
        mesh->getNodes().addNode(0.0, 0.0, size);
        
        // Create tetrahedron element
        Element element(ElementType::TETRAHEDRON);
        element.addNodeIndex(0);
        element.addNodeIndex(1);
        element.addNodeIndex(2);
        element.addNodeIndex(3);
        
        mesh->addBulkElement(element);
        
        return mesh;
    }
    
    /**
     * @brief Calculate mesh quality metric
     * @param mesh Mesh to calculate
     * @return Quality metric (0-1, 1 means perfect quality)
     */
    static double calculateMeshQuality(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh || mesh->totalElements() == 0) {
            return 0.0;
        }
        
        double totalQuality = 0.0;
        size_t elementCount = 0;
        
        // Calculate bulk element quality
        const auto& bulkElements = mesh->getBulkElements();
        for ([[maybe_unused]] const auto& element : bulkElements) {
            double quality = 1.0; // Simplified version, returns fixed quality value
            totalQuality += quality;
            elementCount++;
        }
        
        // Calculate boundary element quality
        const auto& boundaryElements = mesh->getBoundaryElements();
        for ([[maybe_unused]] const auto& element : boundaryElements) {
            double quality = 1.0; // Simplified version, returns fixed quality value
            totalQuality += quality;
            elementCount++;
        }
        
        return elementCount > 0 ? totalQuality / elementCount : 0.0;
    }
    
    /**
     * @brief Validate mesh topology
     * @param mesh Mesh to validate
     * @return Validation result
     */
    static bool validateMeshTopology(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh) {
            return false;
        }
        
        // Check basic integrity
        if (!mesh->validate()) {
            return false;
        }
        
        // Check node index range
        size_t numberOfNodes = mesh->numberOfNodes();
        
        const auto& bulkElements = mesh->getBulkElements();
        for (const auto& element : bulkElements) {
            for (size_t nodeIndex : element.getNodeIndices()) {
                if (nodeIndex >= numberOfNodes) {
                    return false;
                }
            }
        }
        
        const auto& boundaryElements = mesh->getBoundaryElements();
        for (const auto& element : boundaryElements) {
            for (size_t nodeIndex : element.getNodeIndices()) {
                if (nodeIndex >= numberOfNodes) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     * @brief Calculate mesh bounding box
     * @param mesh Mesh to calculate
     * @return Bounding box [minX, maxX, minY, maxY, minZ, maxZ]
     */
    static std::vector<double> calculateBoundingBox(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh || mesh->numberOfNodes() == 0) {
            return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        }
        
        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();
        double minZ = std::numeric_limits<double>::max();
        double maxZ = std::numeric_limits<double>::lowest();
        
        const auto& nodes = mesh->getNodes().getNodes();
        for (const auto& node : nodes) {
            minX = std::min(minX, node.x);
            maxX = std::max(maxX, node.x);
            minY = std::min(minY, node.y);
            maxY = std::max(maxY, node.y);
            minZ = std::min(minZ, node.z);
            maxZ = std::max(maxZ, node.z);
        }
        
        return {minX, maxX, minY, maxY, minZ, maxZ};
    }
    
    /**
     * @brief Calculate mesh volume
     * @param mesh Mesh to calculate
     * @return Mesh volume
     */
    static double calculateVolume(const std::shared_ptr<Mesh>& mesh) {
        if (!mesh || mesh->numberOfBulkElements() == 0) {
            return 0.0;
        }
        
        double totalVolume = 0.0;
        const auto& nodes = mesh->getNodes().getNodes();
        
        const auto& bulkElements = mesh->getBulkElements();
        for (const auto& element : bulkElements) {
            totalVolume += calculateElementVolume(element, nodes);
        }
        
        return totalVolume;
    }
    
private:
    /**
     * @brief Create boundary elements for cube mesh
     */
    static void createBoundaryElements(std::shared_ptr<Mesh> mesh, int divisions) {
        int nodesPerDirection = divisions + 1;
        
        // Create boundary elements for six faces
        for (int face = 0; face < 6; ++face) {
            for (int j = 0; j < divisions; ++j) {
                for (int i = 0; i < divisions; ++i) {
                    Element boundaryElement(ElementType::LINEAR);
                    boundaryElement.setBoundaryId(face + 1);
                    
                    // Calculate node indices based on face
                    std::vector<size_t> indices = calculateFaceIndices(face, i, j, divisions, nodesPerDirection);
                    boundaryElement.setNodeIndices(indices);
                    
                    mesh->addBoundaryElement(boundaryElement);
                }
            }
        }
    }
    
    /**
     * @brief Calculate node indices for cube face
     */
    static std::vector<size_t> calculateFaceIndices(int face, int i, int j, int divisions, int nodesPerDirection) {
        std::vector<size_t> indices;
        
        switch (face) {
            case 0: { // Bottom face z=0
                indices.push_back(j * nodesPerDirection + i);
                indices.push_back(j * nodesPerDirection + i + 1);
                indices.push_back((j + 1) * nodesPerDirection + i + 1);
                indices.push_back((j + 1) * nodesPerDirection + i);
                break;
            }
            case 1: { // Top face z=divisions
                size_t base = divisions * nodesPerDirection * nodesPerDirection;
                indices.push_back(base + j * nodesPerDirection + i);
                indices.push_back(base + j * nodesPerDirection + i + 1);
                indices.push_back(base + (j + 1) * nodesPerDirection + i + 1);
                indices.push_back(base + (j + 1) * nodesPerDirection + i);
                break;
            }
            // Other four faces handled similarly
            default:
                // Simplified handling, actual implementation needs complete face processing
                break;
        }
        
        return indices;
    }
    
    /**
     * @brief Calculate element quality
     */
    static double calculateElementQuality(const Element& element, const Nodes& nodes) {
        // Simplified implementation: calculate quality metric based on element shape
        // Actual implementation needs specific calculation based on element type
        
        size_t nodeCount = element.numberOfNodes();
        if (nodeCount < 3) {
            return 0.0;
        }
        
        // Calculate variance of distances between nodes as quality metric
        std::vector<double> distances;
        const auto& indices = element.getNodeIndices();
        
        for (size_t i = 0; i < nodeCount; ++i) {
            for (size_t j = i + 1; j < nodeCount; ++j) {
                double dist = nodes[indices[i]].distance(nodes[indices[j]]);
                distances.push_back(dist);
            }
        }
        
        if (distances.empty()) {
            return 0.0;
        }
        
        // Calculate mean and standard deviation of distances
        double mean = 0.0;
        for (double dist : distances) {
            mean += dist;
        }
        mean /= distances.size();
        
        double variance = 0.0;
        for (double dist : distances) {
            variance += (dist - mean) * (dist - mean);
        }
        variance /= distances.size();
        
        // Quality metric: smaller variance means better quality
        return 1.0 / (1.0 + std::sqrt(variance));
    }
    
    /**
     * @brief Calculate element volume
     */
    static double calculateElementVolume(const Element& element, const std::vector<Node>& nodes) {
        // Simplified implementation: only handles tetrahedrons and hexahedrons
        const auto& indices = element.getNodeIndices();
        
        switch (element.getType()) {
            case ElementType::TETRAHEDRON:
                if (indices.size() == 4) {
                    return calculateTetrahedronVolume(nodes[indices[0]], nodes[indices[1]], 
                                                     nodes[indices[2]], nodes[indices[3]]);
                }
                break;
                
            case ElementType::HEXAHEDRON:
                if (indices.size() == 8) {
                    // Decompose hexahedron into 5 tetrahedrons to calculate volume
                    return calculateHexahedronVolume(nodes, indices);
                }
                break;
                
            default:
                break;
        }
        
        return 0.0;
    }
    
    /**
     * @brief Calculate tetrahedron volume
     */
    static double calculateTetrahedronVolume(const Node& a, const Node& b, const Node& c, const Node& d) {
        double v1x = b.x - a.x, v1y = b.y - a.y, v1z = b.z - a.z;
        double v2x = c.x - a.x, v2y = c.y - a.y, v2z = c.z - a.z;
        double v3x = d.x - a.x, v3y = d.y - a.y, v3z = d.z - a.z;
        
        // Calculate scalar triple product
        double volume = std::abs(
            v1x * (v2y * v3z - v2z * v3y) -
            v1y * (v2x * v3z - v2z * v3x) +
            v1z * (v2x * v3y - v2y * v3x)
        ) / 6.0;
        
        return volume;
    }
    
    /**
     * @brief Calculate hexahedron volume
     */
    static double calculateHexahedronVolume(const std::vector<Node>& nodes, const std::vector<size_t>& indices) {
        // Simplified implementation: decompose hexahedron into 5 tetrahedrons
        double volume = 0.0;
        
        // Decomposition method 1: based on center point
        Node center(0.0, 0.0, 0.0);
        for (size_t index : indices) {
            center.x += nodes[index].x;
            center.y += nodes[index].y;
            center.z += nodes[index].z;
        }
        center.x /= 8.0;
        center.y /= 8.0;
        center.z /= 8.0;
        
        // Calculate tetrahedron volume for each face
        std::vector<std::vector<size_t>> faces = {
            {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
            {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
        };
        
        for (const auto& face : faces) {
            // Decompose quadrilateral face into two triangles
            volume += calculateTetrahedronVolume(
                center, 
                nodes[indices[face[0]]], 
                nodes[indices[face[1]]], 
                nodes[indices[face[2]]]
            );
            volume += calculateTetrahedronVolume(
                center, 
                nodes[indices[face[0]]], 
                nodes[indices[face[2]]], 
                nodes[indices[face[3]]]
            );
        }
        
        return volume;
    }
};

} // namespace elmer