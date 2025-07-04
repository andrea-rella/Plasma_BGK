#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <muParser.h>

#include "ConfigData.hpp"
#include "SpaceMesh.hpp"

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

    Bgk::SpaceMesh<double> mesh(Data);
    mesh.initialize_mesh();

    mesh.write_mesh_vtk("prova");

    return 0;
}
