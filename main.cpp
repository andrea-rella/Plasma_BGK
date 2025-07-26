#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <muParser.h>

#include "ConfigData.hpp"
#include "SpaceMeshFV.hpp"

int main(int argc, char *argv[])
{
    std::string s;
    // This is a simple C++ program that prints "Hello, World!" to the console.
    std::cout << "Eigen, json and muparser work!!!!!!" << std::endl;

    //---- Example usage of ConfigData ----

    std::string configPath = "data/"; // default path

    if (argc > 1)
        configPath += std::string(argv[1]);
    else
        throw std::runtime_error("No config file specified. Please provide a config file name as a command line argument.");

    Bgk::ConfigData<double> Data(configPath);

    //---- Example usage of SpaceMesh ----

    Bgk::SpaceMeshFV<double> mesh(Data);
    mesh.initialize_mesh();

    // Define uniform spacing
    // double x_start = 0.0;
    // double x_end = 10.0;

    // auto uniform_spacing = [x_start, x_end, &mesh](int i) -> double
    //{
    //     return x_start + (x_end - x_start) * static_cast<double>(i) / static_cast<double>(mesh.get_N());
    // };

    //// Initialize with uniform spacing
    // mesh.initialize_with_custom_spacing(uniform_spacing);

    mesh.write_mesh_vtk("prova");
    mesh.write_mesh_txt("prova");

    return 0;
}
