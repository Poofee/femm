#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace elmer {

/**
 * @brief 电磁材料基类
 */
class ElectromagneticMaterial {
public:
    std::string name;
    double relativePermittivity;    ///< 相对介电常数
    double relativePermeability;    ///< 相对磁导率
    double conductivity;            ///< 电导率 [S/m]
    double frequency;               ///< 频率 [Hz]
    
    ElectromagneticMaterial(double er = 1.0, double mur = 1.0, double sigma = 0.0, 
                           const std::string& matName = "Unknown")
        : name(matName), relativePermittivity(er), relativePermeability(mur), 
          conductivity(sigma), frequency(0.0) {}
    
    /**
     * @brief 计算介电常数
     */
    double permittivity() const {
        return relativePermittivity * 8.854187817e-12; // ε₀ = 8.854187817e-12 F/m
    }
    
    /**
     * @brief 计算磁导率
     */
    double permeability() const {
        return relativePermeability * 1.25663706212e-6; // μ₀ = 1.25663706212e-6 H/m
    }
    
    /**
     * @brief 计算角频率
     */
    double angularFrequency() const {
        return 2.0 * 3.141592653589793 * frequency;
    }
    
    /**
     * @brief 计算趋肤深度
     */
    double skinDepth() const {
        if (conductivity == 0.0 || frequency == 0.0) {
            return 1.0e9; // 非常大的值，表示无趋肤效应
        }
        double omega = angularFrequency();
        double mu = permeability();
        return std::sqrt(2.0 / (omega * mu * conductivity));
    }
    
    /**
     * @brief 计算波数
     */
    double waveNumber() const {
        double omega = angularFrequency();
        double epsilon = permittivity();
        double mu = permeability();
        return omega * std::sqrt(epsilon * mu);
    }
};

/**
 * @brief 材料数据库
 */
class MaterialDatabase {
private:
    std::unordered_map<std::string, ElectromagneticMaterial> materials;
    
public:
    MaterialDatabase() = default;
    
    /**
     * @brief 创建预定义材料
     */
    void createPredefinedMaterials() {
        // 空气
        materials["Air"] = ElectromagneticMaterial(1.0, 1.0, 0.0, "Air");
        
        // 铜
        materials["Copper"] = ElectromagneticMaterial(1.0, 1.0, 5.96e7, "Copper");
        
        // 铁（相对磁导率假设为1000）
        materials["Iron"] = ElectromagneticMaterial(1.0, 1000.0, 1.0e7, "Iron");
        
        // 铝
        materials["Aluminum"] = ElectromagneticMaterial(1.0, 1.0, 3.5e7, "Aluminum");
        
        // 钢
        materials["Steel"] = ElectromagneticMaterial(1.0, 500.0, 4.0e6, "Steel");
        
        // 真空
        materials["Vacuum"] = ElectromagneticMaterial(1.0, 1.0, 0.0, "Vacuum");
    }
    
    /**
     * @brief 获取材料
     */
    ElectromagneticMaterial getMaterial(const std::string& name) const {
        auto it = materials.find(name);
        if (it != materials.end()) {
            return it->second;
        }
        throw std::runtime_error("Material not found: " + name);
    }
    
    /**
     * @brief 添加材料
     */
    void addMaterial(const std::string& name, const ElectromagneticMaterial& material) {
        materials[name] = material;
    }
    
    /**
     * @brief 获取所有材料名称
     */
    std::vector<std::string> getMaterialNames() const {
        std::vector<std::string> names;
        for (const auto& pair : materials) {
            names.push_back(pair.first);
        }
        return names;
    }
    
    /**
     * @brief 检查材料是否存在
     */
    bool hasMaterial(const std::string& name) const {
        return materials.find(name) != materials.end();
    }
};

} // namespace elmer