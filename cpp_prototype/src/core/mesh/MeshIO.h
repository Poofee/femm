#pragma once

#include "Mesh.h"
#include "MeshUtils.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

namespace elmer {

/**
 * @brief 网格文件读写类 - 支持Elmer网格格式
 */
class MeshIO {
public:
    /**
     * @brief 从Elmer网格文件加载网格
     * @param filename 网格文件名
     * @return 加载的网格对象
     */
    static std::shared_ptr<Mesh> loadMesh(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("无法打开网格文件: " + filename);
        }
        
        std::string line;
        
        // 读取文件头
        if (!std::getline(file, line)) {
            throw std::runtime_error("网格文件为空或格式错误");
        }
        
        // 检查文件格式
        if (line.find("ElmerGrid") == std::string::npos) {
            throw std::runtime_error("不支持的网格文件格式");
        }
        
        auto mesh = MeshUtils::allocateMesh(filename);
        
        // 读取节点和单元信息
        readMeshData(file, mesh);
        
        file.close();
        return mesh;
    }
    
    /**
     * @brief 保存网格到Elmer网格文件
     * @param mesh 要保存的网格
     * @param filename 输出文件名
     */
    static void saveMesh(const std::shared_ptr<Mesh>& mesh, const std::string& filename) {
        if (!mesh.get()) {
            throw std::invalid_argument("网格指针为空");
        }
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("无法创建输出文件: " + filename);
        }
        
        // 写入文件头
        file << "! ElmerGrid file created by ElmerCpp\n";
        file << "! Mesh name: " << mesh->getName() << "\n";
        file << "! Number of nodes: " << mesh->numberOfNodes() << "\n";
        file << "! Number of bulk elements: " << mesh->numberOfBulkElements() << "\n";
        file << "! Number of boundary elements: " << mesh->numberOfBoundaryElements() << "\n\n";
        
        // 写入节点数据
        writeNodes(file, mesh);
        
        // 写入体单元数据
        writeBulkElements(file, mesh);
        
        // 写入边界单元数据
        writeBoundaryElements(file, mesh);
        
        file.close();
    }
    
    /**
     * @brief 从Gmsh .msh文件加载网格
     * @param filename Gmsh网格文件名
     * @return 加载的网格对象
     */
    static std::shared_ptr<Mesh> loadGmshMesh(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("无法打开Gmsh文件: " + filename);
        }
        
        auto mesh = MeshUtils::allocateMesh(filename);
        std::string line;
        
        // 读取Gmsh文件格式
        while (std::getline(file, line)) {
            if (line == "$Nodes") {
                readGmshNodes(file, mesh);
            } else if (line == "$Elements") {
                readGmshElements(file, mesh);
            }
        }
        
        file.close();
        return mesh;
    }
    
    /**
     * @brief 导出网格到VTK格式
     * @param mesh 要导出的网格
     * @param filename VTK文件名
     */
    static void exportToVTK(const std::shared_ptr<Mesh>& mesh, const std::string& filename) {
        if (!mesh.get()) {
            throw std::invalid_argument("网格指针为空");
        }
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("无法创建VTK文件: " + filename);
        }
        
        // 写入VTK文件头
        file << "# vtk DataFile Version 3.0\n";
        file << "ElmerCpp Mesh Export\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n\n";
        
        // 写入节点
        file << "POINTS " << mesh->numberOfNodes() << " double\n";
        const auto& nodes = mesh->getNodes().getNodes();
        for (const auto& node : nodes) {
            file << node.x << " " << node.y << " " << node.z << "\n";
        }
        file << "\n";
        
        // 写入体单元
        writeVTKCells(file, mesh);
        
        file.close();
    }
    
    /**
     * @brief 创建简单测试网格
     * @return 测试网格对象
     */
    static std::shared_ptr<Mesh> createTestMesh() {
        auto mesh = MeshUtils::allocateMesh("TestMesh");
        
        // 创建4个节点（四面体）
        mesh->getNodes().addNode(0.0, 0.0, 0.0);
        mesh->getNodes().addNode(1.0, 0.0, 0.0);
        mesh->getNodes().addNode(0.0, 1.0, 0.0);
        mesh->getNodes().addNode(0.0, 0.0, 1.0);
        
        // 创建四面体单元
        Element tetra(ElementType::TETRAHEDRON);
        tetra.addNodeIndex(0);
        tetra.addNodeIndex(1);
        tetra.addNodeIndex(2);
        tetra.addNodeIndex(3);
        tetra.setBodyId(1);
        
        mesh->addBulkElement(tetra);
        
        // 创建边界单元（三角形面）
        std::vector<std::vector<size_t>> boundaryFaces = {
            {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}
        };
        
        for (size_t i = 0; i < boundaryFaces.size(); ++i) {
            Element boundary(ElementType::LINEAR);
            boundary.setNodeIndices(boundaryFaces[i]);
            boundary.setBoundaryId(static_cast<int>(i + 1));
            mesh->addBoundaryElement(boundary);
        }
        
        return mesh;
    }
    
private:
    /**
     * @brief 读取Elmer网格数据
     */
    static void readMeshData(std::ifstream& file, std::shared_ptr<Mesh> mesh) {
        std::string line;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!') {
                continue; // 跳过注释和空行
            }
            
            std::istringstream iss(line);
            std::string keyword;
            iss >> keyword;
            
            if (keyword == "Nodes") {
                readNodes(file, mesh, iss);
            } else if (keyword == "Elements") {
                readElements(file, mesh, iss);
            } else if (keyword == "BoundaryElements") {
                readBoundaryElements(file, mesh, iss);
            }
        }
    }
    
    /**
     * @brief 读取节点数据
     */
    static void readNodes(std::ifstream& file, std::shared_ptr<Mesh> mesh, std::istringstream& iss) {
        size_t numberOfNodes;
        iss >> numberOfNodes;
        
        std::string line;
        for (size_t i = 0; i < numberOfNodes; ++i) {
            if (!std::getline(file, line)) {
                throw std::runtime_error("节点数据不完整");
            }
            
            std::istringstream nodeStream(line);
            double x, y, z;
            nodeStream >> x >> y >> z;
            
            mesh->getNodes().addNode(x, y, z);
        }
    }
    
    /**
     * @brief 读取体单元数据
     */
    static void readElements(std::ifstream& file, std::shared_ptr<Mesh> mesh, std::istringstream& iss) {
        size_t numberOfElements;
        iss >> numberOfElements;
        
        std::string line;
        for (size_t i = 0; i < numberOfElements; ++i) {
            if (!std::getline(file, line)) {
                throw std::runtime_error("体单元数据不完整");
            }
            
            std::istringstream elemStream(line);
            size_t elementId, bodyId;
            elemStream >> elementId >> bodyId;
            
            Element element(ElementType::HEXAHEDRON); // 默认类型
            element.setBodyId(static_cast<int>(bodyId));
            
            size_t nodeIndex;
            while (elemStream >> nodeIndex) {
                element.addNodeIndex(nodeIndex - 1); // Elmer使用1-based索引
            }
            
            mesh->addBulkElement(element);
        }
    }
    
    /**
     * @brief 读取边界单元数据
     */
    static void readBoundaryElements(std::ifstream& file, std::shared_ptr<Mesh> mesh, std::istringstream& iss) {
        size_t numberOfBoundaryElements;
        iss >> numberOfBoundaryElements;
        
        std::string line;
        for (size_t i = 0; i < numberOfBoundaryElements; ++i) {
            if (!std::getline(file, line)) {
                throw std::runtime_error("边界单元数据不完整");
            }
            
            std::istringstream elemStream(line);
            size_t elementId, boundaryId;
            elemStream >> elementId >> boundaryId;
            
            Element boundary(ElementType::LINEAR);
            boundary.setBoundaryId(static_cast<int>(boundaryId));
            
            size_t nodeIndex;
            while (elemStream >> nodeIndex) {
                boundary.addNodeIndex(nodeIndex - 1); // Elmer使用1-based索引
            }
            
            mesh->addBoundaryElement(boundary);
        }
    }
    
    /**
     * @brief 写入节点数据
     */
    static void writeNodes(std::ofstream& file, const std::shared_ptr<Mesh>& mesh) {
        file << "Nodes " << mesh->numberOfNodes() << "\n";
        
        const auto& nodes = mesh->getNodes().getNodes();
        for (size_t i = 0; i < nodes.size(); ++i) {
            file << i + 1 << " " << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z << "\n";
        }
        file << "\n";
    }
    
    /**
     * @brief 写入体单元数据
     */
    static void writeBulkElements(std::ofstream& file, const std::shared_ptr<Mesh>& mesh) {
        file << "Elements " << mesh->numberOfBulkElements() << "\n";
        
        const auto& elements = mesh->getBulkElements();
        for (size_t i = 0; i < elements.size(); ++i) {
            file << i + 1 << " " << elements[i].getBodyId();
            
            const auto& nodeIndices = elements[i].getNodeIndices();
            for (size_t nodeIndex : nodeIndices) {
                file << " " << nodeIndex + 1; // 转换为1-based索引
            }
            file << "\n";
        }
        file << "\n";
    }
    
    /**
     * @brief 写入边界单元数据
     */
    static void writeBoundaryElements(std::ofstream& file, const std::shared_ptr<Mesh>& mesh) {
        file << "BoundaryElements " << mesh->numberOfBoundaryElements() << "\n";
        
        const auto& elements = mesh->getBoundaryElements();
        for (size_t i = 0; i < elements.size(); ++i) {
            file << i + 1 << " " << elements[i].getBoundaryId();
            
            const auto& nodeIndices = elements[i].getNodeIndices();
            for (size_t nodeIndex : nodeIndices) {
                file << " " << nodeIndex + 1; // 转换为1-based索引
            }
            file << "\n";
        }
        file << "\n";
    }
    
    /**
     * @brief 读取Gmsh节点数据
     */
    static void readGmshNodes(std::ifstream& file, std::shared_ptr<Mesh> mesh) {
        std::string line;
        std::getline(file, line);
        
        size_t numberOfNodes;
        std::istringstream iss(line);
        iss >> numberOfNodes;
        
        for (size_t i = 0; i < numberOfNodes; ++i) {
            std::getline(file, line);
            std::istringstream nodeStream(line);
            
            size_t nodeId;
            double x, y, z;
            nodeStream >> nodeId >> x >> y >> z;
            
            mesh->getNodes().addNode(x, y, z);
        }
    }
    
    /**
     * @brief 读取Gmsh单元数据
     */
    static void readGmshElements(std::ifstream& file, std::shared_ptr<Mesh> mesh) {
        std::string line;
        std::getline(file, line);
        
        size_t numberOfElements;
        std::istringstream iss(line);
        iss >> numberOfElements;
        
        for (size_t i = 0; i < numberOfElements; ++i) {
            std::getline(file, line);
            std::istringstream elemStream(line);
            
            size_t elementId, elementType, numTags;
            elemStream >> elementId >> elementType >> numTags;
            
            // 读取标签
            std::vector<size_t> tags(numTags);
            for (size_t j = 0; j < numTags; ++j) {
                elemStream >> tags[j];
            }
            
            // 读取节点
            std::vector<size_t> nodeIndices;
            size_t nodeIndex;
            while (elemStream >> nodeIndex) {
                nodeIndices.push_back(nodeIndex - 1); // Gmsh使用1-based索引
            }
            
            // 根据单元类型创建相应单元
            ElementType elmerType = convertGmshType(elementType);
            if (elmerType != ElementType::LINEAR) { // 跳过未知类型
                Element element(elmerType);
                element.setNodeIndices(nodeIndices);
                
                if (tags.size() >= 2) {
                    if (isBoundaryElement(elementType)) {
                        element.setBoundaryId(static_cast<int>(tags[1]));
                        mesh->addBoundaryElement(element);
                    } else {
                        element.setBodyId(static_cast<int>(tags[1]));
                        mesh->addBulkElement(element);
                    }
                }
            }
        }
    }
    
    /**
     * @brief 将Gmsh单元类型转换为Elmer单元类型
     */
    static ElementType convertGmshType(size_t gmshType) {
        switch (gmshType) {
            case 1:  return ElementType::LINEAR;      // 2节点线
            case 2:  return ElementType::LINEAR;      // 3节点三角形
            case 3:  return ElementType::QUADRATIC;   // 4节点四边形
            case 4:  return ElementType::TETRAHEDRON; // 4节点四面体
            case 5:  return ElementType::HEXAHEDRON;  // 8节点六面体
            case 6:  return ElementType::PRISM;       // 6节点棱柱
            case 7:  return ElementType::PYRAMID;     // 5节点金字塔
            default: return ElementType::LINEAR;      // 默认类型
        }
    }
    
    /**
     * @brief 判断是否为边界单元
     */
    static bool isBoundaryElement(size_t gmshType) {
        return gmshType == 1 || gmshType == 2 || gmshType == 3; // 线、三角形、四边形
    }
    
    /**
     * @brief 写入VTK单元数据
     */
    static void writeVTKCells(std::ofstream& file, const std::shared_ptr<Mesh>& mesh) {
        // 计算总单元数和连接点总数
        size_t totalCells = mesh->numberOfBulkElements();
        size_t totalPoints = 0;
        
        const auto& elements = mesh->getBulkElements();
        for (const auto& element : elements) {
            totalPoints += element.numberOfNodes() + 1; // +1 for cell type
        }
        
        file << "CELLS " << totalCells << " " << totalPoints << "\n";
        
        // 写入单元连接性
        for (const auto& element : elements) {
            const auto& nodeIndices = element.getNodeIndices();
            file << nodeIndices.size();
            
            for (size_t nodeIndex : nodeIndices) {
                file << " " << nodeIndex;
            }
            file << "\n";
        }
        file << "\n";
        
        // 写入单元类型
        file << "CELL_TYPES " << totalCells << "\n";
        for (const auto& element : elements) {
            int vtkType = elmerToVTKElementType(element.getType());
            file << vtkType << "\n";
        }
        file << "\n";
    }
    
    /**
     * @brief 判断是否为体单元
     */
    static bool isBulkElement(ElementType type);
    
    /**
     * @brief 将Gmsh单元类型转换为Elmer单元类型
     */
    static ElementType gmshToElmerElementType(int gmshType);
    
    /**
     * @brief 将Elmer单元类型转换为VTK单元类型
     */
    static int elmerToVTKElementType(ElementType elmerType);
    
    /**
     * @brief 读取节点数据
     */
    static void readNodes(std::ifstream& file, std::shared_ptr<Mesh> mesh, size_t numberOfNodes);
    
    /**
     * @brief 读取体单元数据
     */
    static void readBulkElements(std::ifstream& file, std::shared_ptr<Mesh> mesh, size_t numberOfBulkElements);
    
    /**
     * @brief 读取边界单元数据
     */
    static void readBoundaryElements(std::ifstream& file, std::shared_ptr<Mesh> mesh, size_t numberOfBoundaryElements);
    

    /**
     * @brief 写入VTK节点数据
     */
    static void writeVTKNodes(std::ofstream& file, const std::shared_ptr<Mesh>& mesh);
    
    /**
     * @brief 写入VTK单元数据
     */
    static void writeVTKElements(std::ofstream& file, const std::shared_ptr<Mesh>& mesh);
};

} // namespace elmer