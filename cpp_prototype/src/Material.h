#pragma once

#include <memory>
#include <vector>
#include <array>
#include <map>
#include <string>

namespace elmer {

/**
 * @brief Comprehensive material properties for multiphysics simulations
 */
class Material {
public:
    // Material identification
    std::string name;
    
    // Electromagnetic properties
    double relativePermittivity;      ///< Relative permittivity (ε_r)
    double relativePermeability;      ///< Relative permeability (μ_r)
    double conductivity;              ///< Electrical conductivity (σ)
    
    // Vacuum constants
    static constexpr double VACUUM_PERMITTIVITY = 8.854187817e-12;  ///< ε₀ [F/m]
    static constexpr double VACUUM_PERMEABILITY = 4.0e-7 * M_PI;    ///< μ₀ [H/m]
    
    // Derived electromagnetic properties
    double permittivity() const { return relativePermittivity * VACUUM_PERMITTIVITY; }
    double permeability() const { return relativePermeability * VACUUM_PERMEABILITY; }
    
    // Thermal properties
    double thermalConductivity;       ///< Thermal conductivity (k) [W/(m·K)]
    double density;                   ///< Density (ρ) [kg/m³]
    double specificHeat;             ///< Specific heat capacity (c_p) [J/(kg·K)]
    double heatSource;                ///< Internal heat source (q) [W/m³]
    
    // Mechanical properties
    double youngsModulus;             ///< Young's modulus (E) [Pa]
    double poissonsRatio;             ///< Poisson's ratio (ν)
    double thermalExpansion;          ///< Thermal expansion coefficient (α) [1/K]
    
    // Nonlinear material properties
    bool isNonlinear;                 ///< Flag for nonlinear materials
    std::map<double, double> BHCurve; ///< B-H curve for nonlinear magnetic materials
    
    // Frequency-dependent properties (for harmonic analysis)
    double frequency;                 ///< Frequency [Hz]
    double angularFrequency() const { return 2.0 * M_PI * frequency; }
    
    // Complex properties for harmonic analysis
    double conductivityImag;          ///< Imaginary part of conductivity
    double relativePermittivityImag;  ///< Imaginary part of permittivity
    
    // Default constructor
    Material()
        : relativePermittivity(1.0),
          relativePermeability(1.0),
          conductivity(0.0),
          thermalConductivity(0.0),
          density(0.0),
          specificHeat(0.0),
          heatSource(0.0),
          youngsModulus(0.0),
          poissonsRatio(0.0),
          thermalExpansion(0.0),
          isNonlinear(false),
          frequency(0.0),
          conductivityImag(0.0),
          relativePermittivityImag(0.0),
          name("DefaultMaterial") {}
    
    // Constructor with basic electromagnetic properties
    Material(double eps_r, double mu_r, double sigma, 
             const std::string& matName = "Material")
        : relativePermittivity(eps_r),
          relativePermeability(mu_r),
          conductivity(sigma),
          thermalConductivity(0.0),
          density(0.0),
          specificHeat(0.0),
          heatSource(0.0),
          youngsModulus(0.0),
          poissonsRatio(0.0),
          thermalExpansion(0.0),
          isNonlinear(false),
          frequency(0.0),
          conductivityImag(0.0),
          relativePermittivityImag(0.0),
          name(matName) {}
    
    // Constructor with thermal properties
    Material(double k, double rho, double cp, double q,
             const std::string& matName = "Material")
        : relativePermittivity(1.0),
          relativePermeability(1.0),
          conductivity(0.0),
          thermalConductivity(k),
          density(rho),
          specificHeat(cp),
          heatSource(q),
          youngsModulus(0.0),
          poissonsRatio(0.0),
          thermalExpansion(0.0),
          isNonlinear(false),
          frequency(0.0),
          conductivityImag(0.0),
          relativePermittivityImag(0.0),
          name(matName) {}
    
    // Constructor with mechanical properties
    Material(double E, double nu, double alpha,
             const std::string& matName = "Material")
        : relativePermittivity(1.0),
          relativePermeability(1.0),
          conductivity(0.0),
          thermalConductivity(0.0),
          density(0.0),
          specificHeat(0.0),
          heatSource(0.0),
          youngsModulus(E),
          poissonsRatio(nu),
          thermalExpansion(alpha),
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
};

/**
 * @brief Material database for multiphysics simulations
 */
class MaterialDatabase {
private:
    std::map<std::string, Material> materials;
    
public:
    // Add material to database
    void addMaterial(const std::string& name, const Material& material) {
        materials[name] = material;
    }
    
    // Get material by name
    Material& getMaterial(const std::string& name) {
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
        Material vacuum(1.0, 1.0, 0.0, "Vacuum");
        vacuum.thermalConductivity = 0.0;
        vacuum.density = 0.0;
        vacuum.specificHeat = 0.0;
        addMaterial("Vacuum", vacuum);
        
        // Air
        Material air(1.0006, 1.00000037, 0.0, "Air");
        air.thermalConductivity = 0.026;     // W/(m·K)
        air.density = 1.2;                   // kg/m³
        air.specificHeat = 1005.0;           // J/(kg·K)
        addMaterial("Air", air);
        
        // Copper
        Material copper(1.0, 1.0, 5.96e7, "Copper");
        copper.thermalConductivity = 401.0;   // W/(m·K)
        copper.density = 8960.0;             // kg/m³
        copper.specificHeat = 385.0;         // J/(kg·K)
        copper.youngsModulus = 110e9;        // Pa
        copper.poissonsRatio = 0.34;
        copper.thermalExpansion = 17e-6;     // 1/K
        addMaterial("Copper", copper);
        
        // Iron (linear approximation)
        Material iron(1.0, 5000.0, 1.0e7, "Iron");
        iron.thermalConductivity = 80.0;     // W/(m·K)
        iron.density = 7870.0;               // kg/m³
        iron.specificHeat = 449.0;           // J/(kg·K)
        iron.youngsModulus = 200e9;          // Pa
        iron.poissonsRatio = 0.29;
        iron.thermalExpansion = 12e-6;        // 1/K
        addMaterial("Iron", iron);
        
        // Silicon steel (typical transformer core)
        Material siliconSteel(1.0, 2000.0, 2.0e6, "SiliconSteel");
        siliconSteel.thermalConductivity = 45.0;  // W/(m·K)
        siliconSteel.density = 7650.0;           // kg/m³
        siliconSteel.specificHeat = 480.0;      // J/(kg·K)
        siliconSteel.youngsModulus = 210e9;      // Pa
        siliconSteel.poissonsRatio = 0.30;
        siliconSteel.thermalExpansion = 11e-6;   // 1/K
        addMaterial("SiliconSteel", siliconSteel);
        
        // Water (distilled)
        Material water(80.1, 1.0, 5.5e-6, "Water");
        water.thermalConductivity = 0.6;        // W/(m·K)
        water.density = 1000.0;                 // kg/m³
        water.specificHeat = 4182.0;           // J/(kg·K)
        addMaterial("Water", water);
        
        // Aluminum
        Material aluminum(1.0, 1.0, 3.5e7, "Aluminum");
        aluminum.thermalConductivity = 237.0;   // W/(m·K)
        aluminum.density = 2700.0;             // kg/m³
        aluminum.specificHeat = 897.0;          // J/(kg·K)
        aluminum.youngsModulus = 69e9;          // Pa
        aluminum.poissonsRatio = 0.33;
        aluminum.thermalExpansion = 23e-6;      // 1/K
        addMaterial("Aluminum", aluminum);
        
        // Steel (structural)
        Material steel(1.0, 1.0, 1.0e6, "Steel");
        steel.thermalConductivity = 50.0;       // W/(m·K)
        steel.density = 7850.0;                // kg/m³
        steel.specificHeat = 466.0;            // J/(kg·K)
        steel.youngsModulus = 200e9;           // Pa
        steel.poissonsRatio = 0.30;
        steel.thermalExpansion = 12e-6;         // 1/K
        addMaterial("Steel", steel);
        
        // Concrete
        Material concrete(1.0, 1.0, 1.0e-3, "Concrete");
        concrete.thermalConductivity = 1.7;    // W/(m·K)
        concrete.density = 2400.0;              // kg/m³
        concrete.specificHeat = 880.0;          // J/(kg·K)
        concrete.youngsModulus = 30e9;          // Pa
        concrete.poissonsRatio = 0.20;
        concrete.thermalExpansion = 10e-6;      // 1/K
        addMaterial("Concrete", concrete);
    }
};

} // namespace elmer