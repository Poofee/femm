#include "MeshIO.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <stdexcept>

namespace elmer {




// 私有方法实现
void MeshIO::readNodes(std::ifstream& file, std::shared_ptr<Mesh> mesh, size_t numberOfNodes) {
    std::string line;
    
    // 跳过空行和注释
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '!' && line[0] != '#') {
            break;
        }
    }
    
    // 读取节点数据
    for (size_t i = 0; i < numberOfNodes; ++i) {
        std::istringstream iss(line);
        size_t nodeId;
        double x, y, z;
        iss >> nodeId >> x >> y >> z;
        
        mesh->getNodes().addNode(x, y, z);
        
        if (i < numberOfNodes - 1) {
            std::getline(file, line);
        }
    }
}

void MeshIO::readBulkElements(std::ifstream& file, std::shared_ptr<Mesh> mesh, size_t numberOfBulkElements) {
    std::string line;
    
    // 跳过空行和注释
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '!' && line[0] != '#') {
            break;
        }
    }
    
    // 读取体单元数据
    for (size_t i = 0; i < numberOfBulkElements; ++i) {
        std::istringstream iss(line);
        size_t elementId, elementType, numNodes;
        iss >> elementId >> elementType >> numNodes;
        
        Element element(static_cast<ElementType>(elementType));
        
        // 读取节点索引
        for (size_t j = 0; j < numNodes; ++j) {
            size_t nodeIndex;
            iss >> nodeIndex;
            element.addNodeIndex(nodeIndex - 1); // Elmer网格节点索引从1开始
        }
        
        mesh->addBulkElement(element);
        
        if (i < numberOfBulkElements - 1) {
            std::getline(file, line);
        }
    }
}

void MeshIO::readBoundaryElements(std::ifstream& file, std::shared_ptr<Mesh> mesh, size_t numberOfBoundaryElements) {
    std::string line;
    
    // 跳过空行和注释
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '!' && line[0] != '#') {
            break;
        }
    }
    
    // 读取边界单元数据
    for (size_t i = 0; i < numberOfBoundaryElements; ++i) {
        std::istringstream iss(line);
        size_t elementId, boundaryId, numNodes;
        iss >> elementId >> boundaryId >> numNodes;
        
        Element element(ElementType::LINEAR);
        element.setBoundaryId(boundaryId);
        
        // 读取节点索引
        for (size_t j = 0; j < numNodes; ++j) {
            size_t nodeIndex;
            iss >> nodeIndex;
            element.addNodeIndex(nodeIndex - 1); // Elmer网格节点索引从1开始
        }
        
        mesh->addBoundaryElement(element);
        
        if (i < numberOfBoundaryElements - 1) {
            std::getline(file, line);
        }
    }
}





void MeshIO::writeVTKElements(std::ofstream& file, const std::shared_ptr<Mesh>& mesh) {
    // 计算总单元数和连接性数据大小
    size_t totalElements = mesh->totalElements();
    size_t connectivitySize = 0;
    
    for (const auto& element : mesh->getBulkElements()) {
        connectivitySize += element.numberOfNodes() + 1;
    }
    for (const auto& element : mesh->getBoundaryElements()) {
        connectivitySize += element.numberOfNodes() + 1;
    }
    
    file << "CELLS " << totalElements << " " << connectivitySize << "\n";
    
    // 写入体单元
    for (const auto& element : mesh->getBulkElements()) {
        const auto& indices = element.getNodeIndices();
        file << indices.size();
        for (size_t index : indices) {
            file << " " << index;
        }
        file << "\n";
    }
    
    // 写入边界单元
    for (const auto& element : mesh->getBoundaryElements()) {
        const auto& indices = element.getNodeIndices();
        file << indices.size();
        for (size_t index : indices) {
            file << " " << index;
        }
        file << "\n";
    }
    file << "\n";
    
    // 写入单元类型
    file << "CELL_TYPES " << totalElements << "\n";
    
    // 写入体单元类型
    for (const auto& element : mesh->getBulkElements()) {
        file << elmerToVTKElementType(element.getType()) << "\n";
    }
    
    // 写入边界单元类型
    for (const auto& element : mesh->getBoundaryElements()) {
        file << elmerToVTKElementType(element.getType()) << "\n";
    }
    file << "\n";
}

ElementType MeshIO::gmshToElmerElementType(int gmshType) {
    switch (gmshType) {
        case 1: return ElementType::LINEAR;        // 2节点线
        case 2: return ElementType::TRIANGLE;      // 3节点三角形
        case 3: return ElementType::QUADRILATERAL; // 4节点四边形
        case 4: return ElementType::TETRAHEDRON;   // 4节点四面体
        case 5: return ElementType::HEXAHEDRON;    // 8节点六面体
        case 6: return ElementType::PRISM;         // 6节点棱柱
        case 7: return ElementType::PYRAMID;       // 5节点金字塔
        default: return ElementType::UNKNOWN;
    }
}

bool MeshIO::isBulkElement(ElementType type) {
    switch (type) {
        case ElementType::TETRAHEDRON:
        case ElementType::HEXAHEDRON:
        case ElementType::PRISM:
        case ElementType::PYRAMID:
            return true;
        default:
            return false;
    }
}


} // namespace elmer