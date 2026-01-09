#pragma once

#include <memory>
#include <vector>
#include <array>
#include <map>
#include <string>

namespace elmer {

/**
 * @brief Material properties for electromagnetic simulations
 */
class ElectromagneticMaterial {
public:
    // Basic electromagnetic properties
    double relativePermittivity;      ///< Relative permittivity (ε_r)
    double relativePermeability;      ///< Relative permeability (μ_r)
    double conductivity;              ///< Electrical conductivity (σ)
    
    // Vacuum constants
    static constexpr double VACUUM_PERMITTIVITY = 8.854187817e-12;  ///< ε₀ [F/m]
    static constexpr double VACUUM_PERMEABILITY = 4.0e-7 * M_PI;    ///< μ₀ [H/m]
    
    // Derived properties
    double permittivity() const { return relativePermittivity * VACUUM_PERMITTIVITY; }
    double permeability() const { return relativePermeability * VACUUM_PERMEABILITY; }
    
    // Nonlinear material properties
    bool isNonlinear;                 ///< Flag for nonlinear materials
    std::map<double, double> BHCurve; ///< B-H curve for nonlinear magnetic materials
    
    // Frequency-dependent properties (for harmonic analysis)
    double frequency;                 ///< Frequency [Hz]
    double angularFrequency() const { return 2.0 * M_PI * frequency; }
    
    // Complex properties for harmonic analysis
    double conductivityImag;          ///< Imaginary part of conductivity
    double relativePermittivityImag;  ///< Imaginary part of permittivity
    
    // Material name
    std::string name;
    
    // Default constructor
    ElectromagneticMaterial()
        : relativePermittivity(1.0),
          relativePermeability(1.0),
          conductivity(0.0),
          isNonlinear(false),
          frequency(0.0),
          conductivityImag(0.0),
          relativePermittivityImag(0.0),
          name("DefaultMaterial") {}
    
    // Constructor with basic properties
    ElectromagneticMaterial(double eps_r, double mu_r, double sigma, 
                           const std::string& matName = "Material")
        : relativePermittivity(eps_r),
          relativePermeability(mu_r),
          conductivity(sigma),
          isNonlinear(false),
          frequency(0.0),
          conductivityImag(0.0),
          relativePermittivityImag(0.0),
          name(matName) {}
    
    // Set B-H curve for nonlinear magnetic materials
    void setBHCurve(const std::vector<double>& H_values, 
                    const std::vector<double>& B_values) {
        if (H_values.size() != B_values.size()) {
            throw std::invalid_argument("H and B arrays must have the same size");
        }
        
        BHCurve.clear();
        for (size_t i = 0; i < H_values.size(); ++i) {
            BHCurve[H_values[i]] = B_values[i];
        }
        isNonlinear = true;
    }
    
    // Get permeability from B-H curve (for nonlinear materials)
    double getNonlinearPermeability(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            return permeability();
        }
        
        // Simple linear interpolation for now
        // In practice, this should use more sophisticated interpolation
        auto it = BHCurve.lower_bound(H);
        if (it == BHCurve.begin()) {
            return it->second / it->first;
        }
        if (it == BHCurve.end()) {
            --it;
            return it->second / it->first;
        }
        
        auto prev = it;
        --prev;
        double H1 = prev->first, B1 = prev->second;
        double H2 = it->first, B2 = it->second;
        
        double B = B1 + (B2 - B1) * (H - H1) / (H2 - H1);
        return B / H;
    }
    
    // Get differential permeability (differential permeability for nonlinear materials)
    double getDifferentialPermeability(double H) const {
        if (!isNonlinear || BHCurve.size() < 2) {
            return permeability();
        }
        
        // Calculate derivative of B with respect to H
        auto it = BHCurve.lower_bound(H);
        if (it == BHCurve.begin()) {
            auto next = it;
            ++next;
            return (next->second - it->second) / (next->first - it->first);
        }
        if (it == BHCurve.end()) {
            --it;
            auto prev = it;
            --prev;
            return (it->second - prev->second) / (it->first - prev->first);
        }
        
        auto prev = it;
        --prev;
        auto next = it;
        ++next;
        
        if (next == BHCurve.end()) {
            return (it->second - prev->second) / (it->first - prev->first);
        }
        
        // Central difference for better accuracy
        return (next->second - prev->second) / (next->first - prev->first);
    }
    
    // Get B from H using B-H curve
    double getBfromH(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            return H * permeability();
        }
        
        auto it = BHCurve.lower_bound(H);
        if (it == BHCurve.begin()) {
            return it->second;
        }
        if (it == BHCurve.end()) {
            --it;
            return it->second;
        }
        
        auto prev = it;
        --prev;
        double H1 = prev->first, B1 = prev->second;
        double H2 = it->first, B2 = it->second;
        
        return B1 + (B2 - B1) * (H - H1) / (H2 - H1);
    }
    
    // Get H from B using B-H curve (inverse lookup)
    double getHfromB(double B) const {
        if (!isNonlinear || BHCurve.empty()) {
            return B / permeability();
        }
        
        // Find the interval containing B
        auto it = BHCurve.begin();
        auto next = it;
        ++next;
        
        while (next != BHCurve.end()) {
            if (B >= it->second && B <= next->second) {
                double B1 = it->second, H1 = it->first;
                double B2 = next->second, H2 = next->first;
                return H1 + (H2 - H1) * (B - B1) / (B2 - B1);
            }
            ++it;
            ++next;
        }
        
        // If B is outside the range, use linear extrapolation
        if (B < BHCurve.begin()->second) {
            auto first = BHCurve.begin();
            auto second = first;
            ++second;
            return first->first + (second->first - first->first) * 
                   (B - first->second) / (second->second - first->second);
        } else {
            auto last = BHCurve.rbegin();
            auto prev = last;
            ++prev;
            return last->first + (last->first - prev->first) * 
                   (B - last->second) / (last->second - prev->second);
        }
    }
    
    // Get complex conductivity for harmonic analysis
    std::complex<double> getComplexConductivity() const {
        return std::complex<double>(conductivity, conductivityImag);
    }
    
    // Get complex permittivity for harmonic analysis
    std::complex<double> getComplexPermittivity() const {
        return std::complex<double>(permittivity(), relativePermittivityImag * VACUUM_PERMITTIVITY);
    }
    
    // Wave number calculation
    double waveNumber() const {
        if (frequency <= 0.0) return 0.0;
        double omega = angularFrequency();
        return omega * std::sqrt(permittivity() * permeability());
    }
    
    // Skin depth calculation
    double skinDepth() const {
        if (conductivity <= 0.0 || frequency <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        return 1.0 / std::sqrt(M_PI * frequency * permeability() * conductivity);
    }
    
    // Characteristic impedance
    double characteristicImpedance() const {
        return std::sqrt(permeability() / permittivity());
    }
};

/**
 * @brief Material database for electromagnetic simulations
 */
class MaterialDatabase {
private:
    std::map<std::string, ElectromagneticMaterial> materials;
    
public:
    // Add material to database
    void addMaterial(const std::string& name, const ElectromagneticMaterial& material) {
        materials[name] = material;
    }
    
    // Get material by name
    ElectromagneticMaterial& getMaterial(const std::string& name) {
        auto it = materials.find(name);
        if (it == materials.end()) {
            throw std::runtime_error("Material not found: " + name);
        }
        return it->second;
    }
    
    // Check if material exists
    bool hasMaterial(const std::string& name) const {
        return materials.find(name) != materials.end();
    }
    
    // Get all material names
    std::vector<std::string> getMaterialNames() const {
        std::vector<std::string> names;
        for (const auto& pair : materials) {
            names.push_back(pair.first);
        }
        return names;
    }
    
    // Create predefined materials
    void createPredefinedMaterials() {
        // Vacuum
        ElectromagneticMaterial vacuum(1.0, 1.0, 0.0, "Vacuum");
        addMaterial("Vacuum", vacuum);
        
        // Air
        ElectromagneticMaterial air(1.0006, 1.00000037, 0.0, "Air");
        addMaterial("Air", air);
        
        // Copper
        ElectromagneticMaterial copper(1.0, 1.0, 5.96e7, "Copper");
        addMaterial("Copper", copper);
        
        // Iron (linear approximation)
        ElectromagneticMaterial iron(1.0, 5000.0, 1.0e7, "Iron");
        addMaterial("Iron", iron);
        
        // Silicon steel (typical transformer core)
        ElectromagneticMaterial siliconSteel(1.0, 2000.0, 2.0e6, "SiliconSteel");
        addMaterial("SiliconSteel", siliconSteel);
        
        // Water (distilled)
        ElectromagneticMaterial water(80.1, 1.0, 5.5e-6, "Water");
        addMaterial("Water", water);
    }
};

} // namespace elmer