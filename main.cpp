#include <iostream>
#include "ConfigData.hpp"
#include "SolverFV.hpp"

int main(int argc, char *argv[])
{
    //-------------------------- USAGE OF CONFIG DATA --------------------------
    // -------------------------------------------------------------------------

    std::string configPath = ""; // default path

    if (argc > 1)
        configPath += std::string(argv[1]);
    else
        throw std::runtime_error("No config file specified. Please provide a config file name as a command line argument.");
    Bgk::ConfigData<double> Data(configPath);

    //--------------------------- USAGE OF SOLVER ------------------------------
    // -------------------------------------------------------------------------

    Bgk::SolverFV<double> solver(Data);
    std::cout << "Solver initialized with config: " << configPath << "\n"
              << std::endl;

    solver.initialize();

    // solver.solve_parallel<Bgk::PlotStrategy::ONLYEND>(Bgk::metrics::VectorNormType::L2, Bgk::metrics::RowAggregateType::Max);
    solver.solve<Bgk::PlotStrategy::EACHSTEP>(Bgk::metrics::VectorNormType::L2,
                                              Bgk::metrics::RowAggregateType::Max);

    solver.write_all(Data.get_saving_folder_name());

    return 0;
}
