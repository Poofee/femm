#pragma once

#include <vector>
#include <memory>
#include <string>
#include <cstddef>
#include <cmath>
#include <iostream>

namespace elmer {

/**
 * @brief Node data structure containing 3D coordinates
 */
struct Node {
    double x, y, z;
    
    Node() : x(0.0), y(0.0), z(0.0) {}
    Node(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    // Calculate distance between two points
    double distance(const Node& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

/**
 * @brief Node collection management class
 */
class Nodes {
private:
    std::vector<Node> nodes_;
    
public:
    Nodes() = default;
    
    // Add node
    void addNode(const Node& node) {
        nodes_.push_back(node);
    }
    
    void addNode(double x, double y, double z) {
        nodes_.emplace_back(x, y, z);
    }
    
    // Get number of nodes
    size_t numberOfNodes() const {
        return nodes_.size();
    }
    
    // Get node reference
    Node& operator[](size_t index) {
        return nodes_[index];
    }
    
    const Node& operator[](size_t index) const {
        return nodes_[index];
    }
    
    // Get all nodes
    const std::vector<Node>& getNodes() const {
        return nodes_;
    }
    
    // Clear nodes
    void clear() {
        nodes_.clear();
    }
};

/**
 * @brief Element type definition
 */
enum class ElementType {
    LINEAR,     // Linear element
    QUADRATIC,  // Quadratic element
    CUBIC,      // Cubic element
    TETRAHEDRON, // Tetrahedron
    HEXAHEDRON,  // Hexahedron
    PRISM,       // Prism
    PYRAMID      // Pyramid
};

/**
 * @brief Element data structure
 */
class Element {
private:
    std::vector<size_t> nodeIndices_;  // Node indices
    ElementType type_;                 // Element type
    int bodyId_;                       // Body ID
    int boundaryId_;                   // Boundary ID (used for boundary elements)
    
public:
    Element() : type_(ElementType::LINEAR), bodyId_(0), boundaryId_(0) {}
    
    Element(ElementType type, int bodyId = 0) 
        : type_(type), bodyId_(bodyId), boundaryId_(0) {}
    
    // Add node index
    void addNodeIndex(size_t nodeIndex) {
        nodeIndices_.push_back(nodeIndex);
    }
    
    // Set node indices
    void setNodeIndices(const std::vector<size_t>& indices) {
        nodeIndices_ = indices;
    }
    
    // Get node indices
    const std::vector<size_t>& getNodeIndices() const {
        return nodeIndices_;
    }
    
    // Get element type
    ElementType getType() const {
        return type_;
    }
    
    // Set element type
    void setType(ElementType type) {
        type_ = type;
    }
    
    // Get body ID
    int getBodyId() const {
        return bodyId_;
    }
    
    // Set body ID
    void setBodyId(int bodyId) {
        bodyId_ = bodyId;
    }
    
    // Get boundary ID
    int getBoundaryId() const {
        return boundaryId_;
    }
    
    // Set boundary ID
    void setBoundaryId(int boundaryId) {
        boundaryId_ = boundaryId;
    }
    
    // Get number of nodes
    size_t numberOfNodes() const {
        return nodeIndices_.size();
    }
    
    // Check if it's a boundary element
    bool isBoundaryElement() const {
        return boundaryId_ > 0;
    }
    
    // Check if it's a bulk element
    bool isBulkElement() const {
        return boundaryId_ == 0;
    }
    
    // Calculate element centroid (requires node collection as parameter)
    Node calculateCentroid(const Nodes& nodes) const {
        Node centroid;
        size_t nodeCount = nodeIndices_.size();
        
        if (nodeCount == 0) {
            return centroid; // Return origin if no nodes
        }
        
        for (size_t index : nodeIndices_) {
            const Node& node = nodes[index];
            centroid.x += node.x;
            centroid.y += node.y;
            centroid.z += node.z;
        }
        
        centroid.x /= nodeCount;
        centroid.y /= nodeCount;
        centroid.z /= nodeCount;
        
        return centroid;
    }
    
    // Calculate element bounding box
    void getBoundingBox(const Nodes& nodes, Node& minCorner, Node& maxCorner) const {
        if (nodeIndices_.empty()) {
            minCorner = Node();
            maxCorner = Node();
            return;
        }
        
        minCorner = nodes[nodeIndices_[0]];
        maxCorner = nodes[nodeIndices_[0]];
        
        for (size_t i = 1; i < nodeIndices_.size(); ++i) {
            const Node& node = nodes[nodeIndices_[i]];
            
            if (node.x < minCorner.x) minCorner.x = node.x;
            if (node.y < minCorner.y) minCorner.y = node.y;
            if (node.z < minCorner.z) minCorner.z = node.z;
            
            if (node.x > maxCorner.x) maxCorner.x = node.x;
            if (node.y > maxCorner.y) maxCorner.y = node.y;
            if (node.z > maxCorner.z) maxCorner.z = node.z;
        }
    }
    
    // Calculate element characteristic length (maximum edge length)
    double calculateCharacteristicLength(const Nodes& nodes) const {
        if (nodeIndices_.size() < 2) {
            return 0.0;
        }
        
        double maxDistance = 0.0;
        
        // Calculate distances between all node pairs
        for (size_t i = 0; i < nodeIndices_.size(); ++i) {
            for (size_t j = i + 1; j < nodeIndices_.size(); ++j) {
                const Node& node1 = nodes[nodeIndices_[i]];
                const Node& node2 = nodes[nodeIndices_[j]];
                double distance = node1.distance(node2);
                
                if (distance > maxDistance) {
                    maxDistance = distance;
                }
            }
        }
        
        return maxDistance;
    }
    
    // Get element type string representation
    std::string getTypeString() const {
        switch (type_) {
            case ElementType::LINEAR: return "LINEAR";
            case ElementType::QUADRATIC: return "QUADRATIC";
            case ElementType::CUBIC: return "CUBIC";
            case ElementType::TETRAHEDRON: return "TETRAHEDRON";
            case ElementType::HEXAHEDRON: return "HEXAHEDRON";
            case ElementType::PRISM: return "PRISM";
            case ElementType::PYRAMID: return "PYRAMID";
            default: return "UNKNOWN";
        }
    }
    
    // Get element information string
    std::string getInfoString(const Nodes& nodes) const {
        std::string info = "Element Type: " + getTypeString();
        info += ", Nodes: " + std::to_string(nodeIndices_.size());
        info += ", Body ID: " + std::to_string(bodyId_);
        
        if (isBoundaryElement()) {
            info += ", Boundary ID: " + std::to_string(boundaryId_);
        }
        
        if (!nodeIndices_.empty()) {
            Node centroid = calculateCentroid(nodes);
            info += ", Centroid: [" + std::to_string(centroid.x) + ", " + 
                    std::to_string(centroid.y) + ", " + std::to_string(centroid.z) + "]";
        }
        
        return info;
    }
};

/**
 * @brief Parallel information structure
 */
struct ParallelInfo {
    std::vector<int> globalDOFs;           // Global degrees of freedom
    std::vector<int> gInterface;           // Global interface
    std::vector<std::vector<int>> neighbourList; // Neighbor list
    int numberOfIfDOFs = 0;                // Interface degrees of freedom count
    
    ParallelInfo() = default;
    
    // Initialize parallel information
    void initialize(size_t numberOfNodes) {
        globalDOFs.resize(numberOfNodes, 0);
        gInterface.resize(numberOfNodes, 0);
        neighbourList.resize(numberOfNodes);
    }
    
    // Clear parallel information
    void clear() {
        globalDOFs.clear();
        gInterface.clear();
        neighbourList.clear();
        numberOfIfDOFs = 0;
    }
};

/**
 * @brief 网格数据结构 - 核心类
 */
class Mesh {
private:
    std::string name_;                    // Mesh name
    Nodes nodes_;                         // Node collection
    std::vector<Element> bulkElements_;   // Bulk elements
    std::vector<Element> boundaryElements_; // Boundary elements
    ParallelInfo parallelInfo_;           // Parallel information
    
    // Mesh attributes
    bool changed_ = false;                // Whether mesh has changed
    bool outputActive_ = false;           // Whether output is active
    int adaptiveDepth_ = 0;               // Adaptive depth
    bool discontMesh_ = false;            // Discontinuous mesh
    bool singleMesh_ = false;             // Single mesh
    
    // Statistical information
    size_t maxElementDOFs_ = 0;           // Maximum element DOFs
    size_t maxElementNodes_ = 0;          // Maximum element nodes
    size_t minEdgeDOFs_ = 1000;           // Minimum edge DOFs
    size_t minFaceDOFs_ = 1000;           // Minimum face DOFs
    size_t maxEdgeDOFs_ = 0;              // Maximum edge DOFs
    size_t maxFaceDOFs_ = 0;              // Maximum face DOFs
    size_t maxBDOFs_ = 0;                 // Maximum boundary DOFs
    
public:
    Mesh() = default;
    explicit Mesh(const std::string& name) : name_(name) {}
    
    // 网格名称操作
    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }
    
    // 节点操作
    Nodes& getNodes() { return nodes_; }
    const Nodes& getNodes() const { return nodes_; }
    size_t numberOfNodes() const { return nodes_.numberOfNodes(); }
    
    // 体单元操作
    void addBulkElement(const Element& element) {
        bulkElements_.push_back(element);
        updateStatistics();
    }
    
    std::vector<Element>& getBulkElements() { return bulkElements_; }
    const std::vector<Element>& getBulkElements() const { return bulkElements_; }
    size_t numberOfBulkElements() const { return bulkElements_.size(); }
    
    // 边界单元操作
    void addBoundaryElement(const Element& element) {
        boundaryElements_.push_back(element);
        updateStatistics();
    }
    
    std::vector<Element>& getBoundaryElements() { return boundaryElements_; }
    const std::vector<Element>& getBoundaryElements() const { return boundaryElements_; }
    size_t numberOfBoundaryElements() const { return boundaryElements_.size(); }
    
    // 总单元数
    size_t totalElements() const {
        return bulkElements_.size() + boundaryElements_.size();
    }
    
    // 并行信息操作
    ParallelInfo& getParallelInfo() { return parallelInfo_; }
    const ParallelInfo& getParallelInfo() const { return parallelInfo_; }
    
    // 网格属性操作
    bool isChanged() const { return changed_; }
    void setChanged(bool changed) { changed_ = changed; }
    
    bool isOutputActive() const { return outputActive_; }
    void setOutputActive(bool active) { outputActive_ = active; }
    
    int getAdaptiveDepth() const { return adaptiveDepth_; }
    void setAdaptiveDepth(int depth) { adaptiveDepth_ = depth; }
    
    bool isDiscontMesh() const { return discontMesh_; }
    void setDiscontMesh(bool discont) { discontMesh_ = discont; }
    
    bool isSingleMesh() const { return singleMesh_; }
    void setSingleMesh(bool single) { singleMesh_ = single; }
    
    // 统计信息操作
    size_t getMaxElementDOFs() const { return maxElementDOFs_; }
    size_t getMaxElementNodes() const { return maxElementNodes_; }
    size_t getMinEdgeDOFs() const { return minEdgeDOFs_; }
    size_t getMinFaceDOFs() const { return minFaceDOFs_; }
    size_t getMaxEdgeDOFs() const { return maxEdgeDOFs_; }
    size_t getMaxFaceDOFs() const { return maxFaceDOFs_; }
    size_t getMaxBDOFs() const { return maxBDOFs_; }
    
    // 清空网格
    void clear() {
        nodes_.clear();
        bulkElements_.clear();
        boundaryElements_.clear();
        parallelInfo_.clear();
        resetStatistics();
    }
    
    // 初始化并行信息
    void initializeParallelInfo() {
        parallelInfo_.initialize(nodes_.numberOfNodes());
    }
    
    // Validate mesh integrity
    bool validate() const {
        if (nodes_.numberOfNodes() == 0) {
            return false;
        }
        
        if (totalElements() == 0) {
            return false;
        }
        
        // Validate validity of all element node indices
        for (const auto& element : bulkElements_) {
            for (size_t nodeIndex : element.getNodeIndices()) {
                if (nodeIndex >= nodes_.numberOfNodes()) {
                    return false;
                }
            }
        }
        
        for (const auto& element : boundaryElements_) {
            for (size_t nodeIndex : element.getNodeIndices()) {
                if (nodeIndex >= nodes_.numberOfNodes()) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    // Calculate element centroid
    Node calculateElementCentroid(const Element& element) const {
        Node centroid;
        size_t nodeCount = element.numberOfNodes();
        
        if (nodeCount == 0) {
            return centroid; // Return origin if no nodes
        }
        
        for (size_t nodeIndex : element.getNodeIndices()) {
            const Node& node = nodes_[nodeIndex];
            centroid.x += node.x;
            centroid.y += node.y;
            centroid.z += node.z;
        }
        
        centroid.x /= nodeCount;
        centroid.y /= nodeCount;
        centroid.z /= nodeCount;
        
        return centroid;
    }
    
    // Calculate element volume (simplified version, suitable for tetrahedrons and hexahedrons)
    double calculateElementVolume(const Element& element) const {
        size_t nodeCount = element.numberOfNodes();
        
        if (nodeCount < 3) {
            return 0.0; // At least 3 nodes are needed to form a volume
        }
        
        // Simplified calculation: for linear elements, use average distance from centroid to nodes
        Node centroid = calculateElementCentroid(element);
        double totalDistance = 0.0;
        
        for (size_t nodeIndex : element.getNodeIndices()) {
            const Node& node = nodes_[nodeIndex];
            totalDistance += centroid.distance(node);
        }
        
        // Approximate volume calculation (simplified)
        double avgDistance = totalDistance / nodeCount;
        
        // Return different volume approximations based on element type
        switch (element.getType()) {
            case ElementType::TETRAHEDRON:
                return (avgDistance * avgDistance * avgDistance) / 6.0; // V = a³/6
            case ElementType::HEXAHEDRON:
                return avgDistance * avgDistance * avgDistance; // V = a³
            case ElementType::PRISM:
                return (avgDistance * avgDistance * avgDistance) / 2.0; // V = a³/2
            case ElementType::PYRAMID:
                return (avgDistance * avgDistance * avgDistance) / 3.0; // V = a³/3
            default:
                return avgDistance * avgDistance; // For 2D elements, return area approximation
        }
    }
    
    // Find element containing given point
    Element* findElementContainingPoint(const Node& point, double tolerance = 1e-6) {
        for (auto& element : bulkElements_) {
            Node centroid = calculateElementCentroid(element);
            double distance = centroid.distance(point);
            
            // Simplified judgment: if point is near element centroid, consider it contained
            if (distance < tolerance) {
                return &element;
            }
        }
        
        return nullptr; // Element containing the point not found
    }
    
    // Get mesh bounding box
    void getBoundingBox(Node& minCorner, Node& maxCorner) const {
        if (nodes_.numberOfNodes() == 0) {
            minCorner = Node();
            maxCorner = Node();
            return;
        }
        
        minCorner = nodes_[0];
        maxCorner = nodes_[0];
        
        for (size_t i = 1; i < nodes_.numberOfNodes(); ++i) {
            const Node& node = nodes_[i];
            
            if (node.x < minCorner.x) minCorner.x = node.x;
            if (node.y < minCorner.y) minCorner.y = node.y;
            if (node.z < minCorner.z) minCorner.z = node.z;
            
            if (node.x > maxCorner.x) maxCorner.x = node.x;
            if (node.y > maxCorner.y) maxCorner.y = node.y;
            if (node.z > maxCorner.z) maxCorner.z = node.z;
        }
    }
    
    // Calculate total mesh volume
    double calculateTotalVolume() const {
        double totalVolume = 0.0;
        
        for (const auto& element : bulkElements_) {
            totalVolume += calculateElementVolume(element);
        }
        
        return totalVolume;
    }
    
    // Get mesh statistics
    void getMeshStatistics(std::ostream& os) const {
        os << "=== Mesh Statistics ===" << std::endl;
        os << "Name: " << name_ << std::endl;
        os << "Number of nodes: " << numberOfNodes() << std::endl;
        os << "Number of bulk elements: " << numberOfBulkElements() << std::endl;
        os << "Number of boundary elements: " << numberOfBoundaryElements() << std::endl;
        os << "Total elements: " << totalElements() << std::endl;
        os << "Max element nodes: " << maxElementNodes_ << std::endl;
        os << "Total volume: " << calculateTotalVolume() << std::endl;
        
        // Calculate bounding box
        Node minCorner, maxCorner;
        getBoundingBox(minCorner, maxCorner);
        os << "Bounding box: [" << minCorner.x << ", " << minCorner.y << ", " << minCorner.z << "] - ["
           << maxCorner.x << ", " << maxCorner.y << ", " << maxCorner.z << "]" << std::endl;
    }
    
private:
    // Update statistics
    void updateStatistics() {
        maxElementNodes_ = 0;
        
        // Update bulk element statistics
        for (const auto& element : bulkElements_) {
            size_t nodeCount = element.numberOfNodes();
            if (nodeCount > maxElementNodes_) {
                maxElementNodes_ = nodeCount;
            }
        }
        
        // Update boundary element statistics
        for (const auto& element : boundaryElements_) {
            size_t nodeCount = element.numberOfNodes();
            if (nodeCount > maxElementNodes_) {
                maxElementNodes_ = nodeCount;
            }
        }
    }
    
    // Reset statistics
    void resetStatistics() {
        maxElementDOFs_ = 0;
        maxElementNodes_ = 0;
        minEdgeDOFs_ = 1000;
        minFaceDOFs_ = 1000;
        maxEdgeDOFs_ = 0;
        maxFaceDOFs_ = 0;
        maxBDOFs_ = 0;
    }
};

} // namespace elmer