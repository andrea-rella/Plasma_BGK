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

#ifndef SOLVERFV_A4D77DB0_B281_477E_A2CA_231E73195044
#define SOLVERFV_A4D77DB0_B281_477E_A2CA_231E73195044

#include "../SolverFV.hpp"
#include <numbers>
#include <cmath>

namespace Bgk
{

    // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    SolverFV<T>::SolverFV(const ConfigData<T> &InputData) : Data(InputData), Space_mesh(Data), Velocity_mesh(Data)
    {

        // numerical matrices
        A = Eigen::SparseMatrix<T>(Space_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        B = Eigen::SparseMatrix<T>(Space_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        R = Eigen::SparseMatrix<T>(Space_mesh.get_N() + 1, Space_mesh.get_N() + 1);

        // solution matrices
        g = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2 * Velocity_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        h = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2 * Velocity_mesh.get_N() + 1, Space_mesh.get_N() + 1);

        // vector used for the solution of the numerical systems
        SolVector = Eigen::Matrix<T, -1, 1>::Zero(Space_mesh.get_N());

        // boundary conditions for g and h
        T T_infty_w = Data.get_T_infty_w();
        T p_infty_w = Data.get_p_infty_w();
        T M_infty = Data.get_M_infty();

        g0 = [](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * std::exp(-z * z);
        };

        g_infty = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-1.5}) *
                   std::exp(-std::pow((z + std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        h0 = [](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * std::exp(-z * z);
        };

        h_infty = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-0.5}) *
                   std::exp(-std::pow((z + std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        // initial conditions for g and h
        g_init = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-1.5}) *
                   std::exp(-std::pow((z + std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        h_init = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-0.5}) *
                   std::exp(-std::pow((z + std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };
    }

    template <typename T>
    SolverFV<T>::SolverFV(const std::string &config_file_path) : SolverFV(ConfigData<T>(config_file_path)){};

    // ------ BUILD MATRICES / SETUP -----------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    void SolverFV<T>::initializeMeshes()
    {
        Space_mesh.initialize_mesh();
        Velocity_mesh.initialize_mesh();
    }

    template <typename T>
    void SolverFV<T>::setInitialState()
    {
        if (!Space_mesh.validate_mesh() || !Velocity_mesh.validate_mesh())
        {
            throw std::runtime_error(error_message("Meshes must be valid before initializing the initial state of the solver"));
        }

        size_t Space_N = Space_mesh.get_N();
        size_t Velocity_N = Velocity_mesh.get_N();

        // --- Negative velocities (j = 0 .. Velocity_N-1)
        for (size_t j = 0; j < Velocity_N; ++j)
        {
            T vj = Velocity_mesh[j];
            T g_val = g_init(vj);
            T h_val = h_init(vj);

            for (size_t i = 0; i < Space_N; ++i)
            {
                g(j, i) = g_val;
                h(j, i) = h_val;
            }

            g(j, Space_N) = g_infty(vj);
            h(j, Space_N) = h_infty(vj);
        }

        // --- Zero velocity (j = Velocity_N)
        {
            T v0 = Velocity_mesh[Velocity_N];
            T g_val = g_init(v0);
            T h_val = h_init(v0);

            for (size_t i = 0; i <= Space_N; ++i)
            {
                g(Velocity_N, i) = g_val;
                h(Velocity_N, i) = h_val;
            }
        }

        // --- Positive velocities (j = Velocity_N+1 .. 2*Velocity_N)
        for (size_t j = Velocity_N + 1; j <= 2 * Velocity_N; ++j)
        {
            T vj = Velocity_mesh[j];
            T g_val = g_init(vj);
            T h_val = h_init(vj);

            g(j, 0) = g0(vj);
            h(j, 0) = h0(vj);

            for (size_t i = 1; i <= Space_N; ++i)
            {
                g(j, i) = g_val;
                h(j, i) = h_val;
            }
        }
    }

    // ------ OUTPUT ---------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    void SolverFV<T>::write_sol_txt(const std::string &folder_name) const
    {
        std::filesystem::create_directories("output/" + folder_name);

        // Write g solution
        std::string filename_g = "output/" + folder_name + "/solution_g.txt";
        std::ofstream txt_file_g(filename_g);
        if (!txt_file_g.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_g << std::endl;
            return;
        }

        txt_file_g << "g solution in the phase space:\n";
        txt_file_g << "Rows: Velocity (j)  |  Columns Space (i)  |  g(j, i)\n";
        txt_file_g << "Dimensions: " << g.rows() << " x " << g.cols() << "\n";
        txt_file_g << "--------------------------------------------------------\n";

        for (Eigen::Index j = 0; j < g.rows(); ++j)
        {
            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                txt_file_g << g(j, i) << "  ";
            }
            txt_file_g << "\n";
        }

        // Write h solution
        std::string filename_h = "output/" + folder_name + "/solution_h.txt";
        std::ofstream txt_file_h(filename_h);
        if (!txt_file_h.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_h << std::endl;
            return;
        }

        txt_file_h << "h solution in the phase space:\n";
        txt_file_h << "Rows: Velocity (j)  |  Columns Space (i)  |  h(j, i)\n";
        txt_file_h << "Dimensions: " << h.rows() << " x " << h.cols() << "\n";
        txt_file_h << "--------------------------------------------------------\n";

        for (Eigen::Index j = 0; j < h.rows(); ++j)
        {
            for (Eigen::Index i = 0; i < h.cols(); ++i)
            {
                txt_file_h << h(j, i) << "  ";
            }
            txt_file_h << "\n";
        }
    }
}

#endif /* SOLVERFV_A4D77DB0_B281_477E_A2CA_231E73195044 */
