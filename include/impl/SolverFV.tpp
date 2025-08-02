// ██▓███   ██▓    ▄▄▄        ██████  ███▄ ▄███▓ ▄▄▄          ▄▄▄▄     ▄████  ██ ▄█▀
// ▓██░  ██▒▓██▒   ▒████▄    ▒██    ▒ ▓██▒▀█▀ ██▒▒████▄       ▓█████▄  ██▒ ▀█▒ ██▄█▒
// ▓██░ ██▓▒▒██░   ▒██  ▀█▄  ░ ▓██▄   ▓██    ▓██░▒██  ▀█▄     ▒██▒ ▄██▒██░▄▄▄░▓███▄░
// ▒██▄█▓▒ ▒▒██░   ░██▄▄▄▄██   ▒   ██▒▒██    ▒██ ░██▄▄▄▄██    ▒██░█▀  ░▓█  ██▓▓██ █▄
// ▒██▒ ░  ░░██████▒▓█   ▓██▒▒██████▒▒▒██▒   ░██▒ ▓█   ▓██▒   ░▓█  ▀█▓░▒▓███▀▒▒██▒ █▄
// ▒▓▒░ ░  ░░ ▒░▓  ░▒▒   ▓▒█░▒ ▒▓▒ ▒ ░░ ▒░   ░  ░ ▒▒   ▓▒█░   ░▒▓███▀▒ ░▒   ▒ ▒ ▒▒ ▓▒
// ░▒ ░     ░ ░ ▒  ░ ▒   ▒▒ ░░ ░▒  ░ ░░  ░      ░  ▒   ▒▒ ░   ▒░▒   ░   ░   ░ ░ ░▒ ▒░
// ░░         ░ ░    ░   ▒   ░  ░  ░  ░      ░     ░   ▒       ░    ░ ░ ░   ░ ░ ░░ ░
//              ░  ░     ░  ░      ░         ░         ░  ░    ░            ░ ░  ░
//
// Andrea Rella
// Politecnico di Milano
// https://github.com/andrea-rella/Plasma_BGK

#ifndef SOLVERFV_IMPL_CDF025BF_A5BC_40F3_8ED2_AED426E4A6ED
#define SOLVERFV_IMPL_CDF025BF_A5BC_40F3_8ED2_AED426E4A6ED

#include "../SolverFV.hpp"

namespace Bgk
{
    template <typename T>
    SolverFV<T>::SolverFV(const std::string &config_file_path)
        : Data(config_file_path), Space_mesh(Data), Velocity_mesh(Data)
    {
        A = Eigen::SparseMatrix<T>(Space_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        B = Eigen::SparseMatrix<T>(Space_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        R = Eigen::SparseMatrix<T>(Space_mesh.get_N() + 1, Space_mesh.get_N() + 1);

        g = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(Space_mesh.get_N() + 1, Velocity_mesh.get_N() + 1);
        h = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(Space_mesh.get_N() + 1, Velocity_mesh.get_N() + 1);

        SolVector = Eigen::Matrix<T, -1, 1>::Zero(Space_mesh.get_N() + 1);

        g0 = [](T)
        { return T(0); };
        g_infty = [](T)
        { return T(0); };
        h0 = [](T)
        { return T(0); };
        h_infty = [](T)
        { return T(0); };

        g_init = [](T)
        { return T(0); };
        h_init = [](T)
        { return T(0); };
    }
}

#endif /* SOLVERFV_IMPL_CDF025BF_A5BC_40F3_8ED2_AED426E4A6ED */
