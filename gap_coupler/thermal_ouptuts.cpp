#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <iostream>

using namespace std;

struct ThermalStats {
    double T_min = 0, T_mean = 0, T_max = 0;
    double D_min = 0, D_max = 0;
    bool valid = false;
};

std::vector<double> extract_scalar_values(const std::string& filename, const std::string& scalar_name) {
    std::ifstream vtkfile(filename);
    std::vector<double> values;
    std::string line;

    while (std::getline(vtkfile, line)) {
        if (line.find("SCALARS " + scalar_name) != std::string::npos) {
            while (std::getline(vtkfile, line)) {
                if (line.find("LOOKUP_TABLE") != std::string::npos)
                    break;
            }
            while (std::getline(vtkfile, line)) {
                if (!line.empty() && std::isalpha(line[0])) break;
                std::istringstream iss(line);
                double val;
                while (iss >> val) {
                    values.push_back(val);
                }
            }
            break;
        }
    }
    return values;
}

ThermalStats parse_thermal_vtk(const std::string& filename) {
    ThermalStats stats;
    auto T_vals = extract_scalar_values(filename, "T_[C]");
    auto D_vals = extract_scalar_values(filename, "Total_Deformation_[um]");

    if (!T_vals.empty() && !D_vals.empty()) {
        stats.T_min = *std::min_element(T_vals.begin(), T_vals.end());
        stats.T_max = *std::max_element(T_vals.begin(), T_vals.end());
        stats.T_mean = std::accumulate(T_vals.begin(), T_vals.end(), 0.0) / T_vals.size();
        stats.D_min = *std::min_element(D_vals.begin(), D_vals.end());
        stats.D_max = *std::max_element(D_vals.begin(), D_vals.end());
        stats.valid = true;
    }
    return stats;
}

void log_thermal_summary(const std::string& label, int rev_num, const std::string& vtk_path, std::ostream& out = std::cout) {
    ThermalStats stats = parse_thermal_vtk(vtk_path);
    if (!stats.valid) return;

    out << "[" << rev_num << "] " << label << " Revolution " << rev_num << " complete." << std::endl;
    out << "Starting Thermal Analysis:\n" << std::endl;
    out << "T_" << label << "_Min:   " << fixed << setprecision(1) << stats.T_min << " °C" << std::endl;
    out << "T_" << label << "_Mean:  " << stats.T_mean << " °C" << std::endl;
    out << "T_" << label << "_Max:   " << stats.T_max << " °C" << std::endl << std::endl;
    out << "Thermal Deformation after Rev " << rev_num << ":" << std::endl;
    out << "DT_" << label << ":min:  " << stats.D_min << " µm" << std::endl;
    out << "DT_" << label << ":max:  " << stats.D_max << " µm" << std::endl << std::endl;
}

// Sample integration into your main loop (just before simulation ends):
// if (i == rev - 1 && interface_enabled[0]) {
//     log_thermal_summary("Piston", i + 1, "output/piston/vtk/piston_th_rev_" + to_string(i + 1) + ".vtk", logfile);
// }
// if (i == rev - 1 && interface_enabled[1]) {
//     log_thermal_summary("Block", i + 1, "output/block/vtk/cylinderblock.vtk", logfile);
// }
// if (i == rev - 1 && interface_enabled[2]) {
//     log_thermal_summary("Slipper", i + 1, "output/slipper/vtk/slipper_thermal.360.vtk", logfile);
// }