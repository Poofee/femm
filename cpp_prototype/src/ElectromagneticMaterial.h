#pragma once

#include <memory>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <complex>

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
    static constexpr double VACUUM_PERMEABILITY = 4.0e-7 * 3.141592653589793;    ///< μ₀ [H/m]
    
    // Derived properties
    double permittivity() const { return relativePermittivity * VACUUM_PERMITTIVITY; }
    double permeability() const { return relativePermeability * VACUUM_PERMEABILITY; }
    
    // Nonlinear material properties
    bool isNonlinear;                 ///< Flag for nonlinear materials
    std::map<double, double> BHCurve; ///< B-H curve for nonlinear magnetic materials
    
    // Frequency-dependent properties (for harmonic analysis)
    double frequency;                 ///< Frequency [Hz]
    double angularFrequency() const { return 2.0 * 3.141592653589793 * frequency; }
    
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
        
        // Linear interpolation for B
        double B = B1 + (H - H1) * (B2 - B1) / (H2 - H1);
        return B / H;
    }
    
    // Get permeability and its derivative (for Newton-Raphson iteration)
    std::pair<double, double> getNonlinearPermeabilityWithDerivative(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            return {permeability(), 0.0}; // Linear material has zero derivative
        }
        
        auto it = BHCurve.lower_bound(H);
        if (it == BHCurve.begin()) {
            double mu = it->second / it->first;
            return {mu, 0.0}; // Assume constant at beginning
        }
        if (it == BHCurve.end()) {
            --it;
            double mu = it->second / it->first;
            return {mu, 0.0}; // Assume constant at end
        }
        
        auto prev = it;
        --prev;
        double H1 = prev->first, B1 = prev->second;
        double H2 = it->first, B2 = it->second;
        
        // Linear interpolation for B
        double B = B1 + (H - H1) * (B2 - B1) / (H2 - H1);
        double mu = B / H;
        
        // Derivative of mu with respect to H
        double dBdH = (B2 - B1) / (H2 - H1);
        double dMudH = (dBdH * H - B) / (H * H);
        
        return {mu, dMudH};
    }
    
    // Get reluctivity (inverse of permeability)
    double getReluctivity(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            return 1.0 / permeability();
        }
        
        double mu = getNonlinearPermeability(H);
        return 1.0 / mu;
    }
    
    // Get reluctivity and its derivative
    std::pair<double, double> getReluctivityWithDerivative(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            double nu = 1.0 / permeability();
            return {nu, 0.0}; // Linear material has zero derivative
        }
        
        auto [mu, dMudH] = getNonlinearPermeabilityWithDerivative(H);
        double nu = 1.0 / mu;
        double dNudH = -dMudH / (mu * mu);
        
        return {nu, dNudH};
    }
    
    // Get B from H using B-H curve
    double getBfromH(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            return permeability() * H;
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
    
    // Get differential permeability (dμ/dH) for nonlinear materials
    double getDifferentialPermeability(double H) const {
        if (!isNonlinear || BHCurve.empty()) {
            return 0.0; // Linear materials have constant permeability
        }
        
        // Find the interval containing H
        auto it = BHCurve.upper_bound(H);
        if (it == BHCurve.begin()) {
            // H is before the first point
            auto first = BHCurve.begin();
            auto second = first;
            ++second;
            return (second->second - first->second) / (second->first - first->first);
        } else if (it == BHCurve.end()) {
            // H is after the last point
            auto last = BHCurve.rbegin();
            auto prev = last;
            ++prev;
            return (last->second - prev->second) / (last->first - prev->first);
        } else {
            // H is within the curve
            auto prev = it;
            --prev;
            double H1 = prev->first, B1 = prev->second;
            double H2 = it->first, B2 = it->second;
            return (B2 - B1) / (H2 - H1);
        }
    }
    
    // Complex permeability for harmonic analysis
    std::complex<double> getComplexPermeability(double frequency = 0.0) const {
        double mu_real = permeability();
        double mu_imag = 0.0; // For now, assume no imaginary part
        
        // Frequency-dependent effects could be added here
        if (frequency > 0.0) {
            // Simple model: imaginary part proportional to conductivity
            mu_imag = conductivity / (2.0 * 3.141592653589793 * frequency);
        }
        
        return std::complex<double>(mu_real, mu_imag);
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
    double waveNumber(double frequency = 0.0) const {
        if (frequency <= 0.0) return 0.0;
        double omega = 2.0 * 3.141592653589793 * frequency;
        return omega * sqrt(permittivity() * permeability());
    }
    
    // Skin depth calculation
    double skinDepth(double frequency = 0.0) const {
        if (conductivity <= 0.0 || frequency <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        return 1.0 / sqrt(3.141592653589793 * frequency * permeability() * conductivity);
    }
    
    // Characteristic impedance
    double characteristicImpedance() const {
        return sqrt(permeability() / permittivity());
    }
};

/**
 * @brief Nonlinear material database for electromagnetic simulations
 */
class NonlinearMaterialDatabase {
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
    
    // Clear all materials
    void clear() {
        materials.clear();
    }
    
    // Get number of materials
    size_t size() const {
        return materials.size();
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