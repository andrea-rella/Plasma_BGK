#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include "ConfigData.hpp"
#include "SpaceMeshFV.hpp"
#include "VelocityMesh.hpp"
#include "utilities.hpp"
#include "SolverFV.hpp"
#include "metrics_utils.hpp"

void print_vector(const std::vector<std::pair<double, double>> &vec);
void print_vector(const std::vector<double> &vec);

int main(int argc, char *argv[])
{
    //-------------------------- USAGE OF CONFIG DATA --------------------------
    // -------------------------------------------------------------------------

    std::string configPath = "data/"; // default path

    if (argc > 1)
        configPath += std::string(argv[1]);
    else
        throw std::runtime_error("No config file specified. Please provide a config file name as a command line argument.");

    Bgk::ConfigData<double> Data(configPath);

    //-------------------------- USAGE OF SPACE MESH ---------------------------
    // -------------------------------------------------------------------------

    // Bgk::SpaceMeshFV<double> space_mesh(Data);
    // space_mesh.initialize_mesh();

    // Define uniform spacing
    // double x_start = 0.0;
    // double x_end = 10.0;

    // auto uniform_spacing = [x_start, x_end, &mesh](int i) -> double
    //{
    //     return x_start + (x_end - x_start) * static_cast<double>(i) / static_cast<double>(mesh.get_N());
    // };

    //// Initialize with uniform spacing
    // space_mesh.initialize_with_custom_spacing(uniform_spacing);

    // space_mesh.write_mesh_vtk("prova");
    //  space_mesh.write_mesh_txt("prova");

    //------------------------ USAGE OF VELOCITY MESH --------------------------
    // -------------------------------------------------------------------------

    // Bgk::VelocityMesh<double> velocity_mesh(Data);
    // velocity_mesh.initialize_mesh();

    // velocity_mesh.write_mesh_txt("prova");

    //-------------------- USAGE OF VELOCITY QUICK COEFF -----------------------
    // -------------------------------------------------------------------------

    // auto quick_coefficients = Bgk::QUICK_coefficients_p(space_mesh);
    // std::cout << "QUICK coefficients for Space Mesh (p):" << std::endl;
    // print_vector(quick_coefficients);
    //
    // auto quick_coefficients_n = Bgk::QUICK_coefficients_n(space_mesh);
    // std::cout << "QUICK coefficients for Space Mesh (n):" << std::endl;
    // print_vector(quick_coefficients_n);

    //--------------------------- USAGE OF SOLVER ------------------------------
    // -------------------------------------------------------------------------
    Bgk::SolverFV<double> solver(configPath);
    std::cout << "Solver initialized with config: " << configPath << "\n"
              << std::endl;

    solver.initialize();
    solver.write_initial_state_txt("first_test");
    solver.solve(Bgk::metrics::VectorNormType::L2, Bgk::metrics::RowAggregateType::Max);
    solver.write_all("first_test");

    return 0;
}

void print_vector(const std::vector<std::pair<double, double>> &vec)
{
    std::cout << "a_1: " << std::endl;
    for (const auto &val : vec)
    {
        std::cout << val.first << " ";
    }
    std::cout << std::endl;

    std::cout << "a_2: " << std::endl;
    for (const auto &val : vec)
    {
        std::cout << val.second << " ";
    }
    std::cout << std::endl;
};

void print_vector(const std::vector<double> &vec)
{
    for (const auto &val : vec)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
};