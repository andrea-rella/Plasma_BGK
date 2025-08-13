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

#ifndef SOLVERFV_F5CE4BC0_A11B_431A_A0BD_E176318E8B09
#define SOLVERFV_F5CE4BC0_A11B_431A_A0BD_E176318E8B09

#include "../SolverFV.hpp"
#include "phys_utils.hpp"
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
        A = Eigen::SparseMatrix<T>(Space_mesh.get_N(), Space_mesh.get_N());
        B = Eigen::SparseMatrix<T>(Space_mesh.get_N(), Space_mesh.get_N());
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

    // ------ BUILD NUMERICAL MATRICES / SETUP -------------------------------------------------------
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

    template <typename T>
    void SolverFV<T>::assemble_A()
    {
        const Eigen::Index N = A.rows();
        const std::vector<std::pair<T, T>> QUICK_a = QUICK_coefficients_p(Space_mesh);
        const std::vector<T> vol_sizes = Space_mesh.get_volume_sizes();

        // Reserve space efficiently
        const Eigen::Index nnz = N + 2 * std::max(Eigen::Index{0}, N - 1) + std::max(Eigen::Index{0}, N - 2);
        std::vector<Eigen::Triplet<T>> triplets;
        triplets.reserve(nnz);

        // Helper lambdas
        auto main_coeff = [&](Eigen::Index i) -> T
        {
            return (T{1} - QUICK_a[i + 2].first + QUICK_a[i + 2].second - QUICK_a[i + 1].first) / vol_sizes[i + 1];
        };

        auto sub1_coeff = [&](Eigen::Index i) -> T
        {
            return (QUICK_a[i + 2].second - T{1} + QUICK_a[i + 1].first - QUICK_a[i + 1].second) / vol_sizes[i + 1];
        };

        auto super1_coeff = [&](Eigen::Index i) -> T
        {
            return QUICK_a[i + 2].first / vol_sizes[i + 1];
        };

        auto sub2_coeff = [&](Eigen::Index i) -> T
        {
            return (-QUICK_a[i + 1].second) / vol_sizes[i + 1];
        };

        // First row (i = 0)
        triplets.emplace_back(0, 0, main_coeff(0));
        if (N > 1)
        {
            triplets.emplace_back(0, 1, super1_coeff(0));
        }

        // Second row (i = 1)
        if (N > 1)
        {
            triplets.emplace_back(1, 0, sub1_coeff(1));
            triplets.emplace_back(1, 1, main_coeff(1));
            if (N > 2)
            {
                triplets.emplace_back(1, 2, super1_coeff(1));
            }
        }

        // Interior rows (no conditionals in loop - maximum performance)
        // Loop from 2 to N-3 (inclusive)
        const Eigen::Index loop_end = std::max(Eigen::Index{2}, N - 2);
        for (Eigen::Index i = 2; i < loop_end; ++i)
        {
            triplets.emplace_back(i, i - 2, sub2_coeff(i));
            triplets.emplace_back(i, i - 1, sub1_coeff(i));
            triplets.emplace_back(i, i, main_coeff(i));
            triplets.emplace_back(i, i + 1, super1_coeff(i));
        }

        // Second-to-last row (i = N-2)
        if (N > 2)
        {
            const Eigen::Index i = N - 2;
            triplets.emplace_back(i, i - 2, sub2_coeff(i));
            triplets.emplace_back(i, i - 1, sub1_coeff(i));
            triplets.emplace_back(i, i, main_coeff(i));
            if (N > 3)
            {
                triplets.emplace_back(i, i + 1, super1_coeff(i));
            }
        }

        // Last row (i = N-1) (Custum handling)
        if (N > 1)
        {
            const Eigen::Index i = N - 1;
            if (N > 2)
            {
                triplets.emplace_back(i, i - 2, -QUICK_a[i + 1].second / vol_sizes[i + 1]);
            }
            triplets.emplace_back(i, i - 1, (QUICK_a[i + 1].first - T{1} - QUICK_a[i + 1].second) / vol_sizes[i + 1]);
            triplets.emplace_back(i, i, (T{1} - QUICK_a[i + 1].first) / vol_sizes[i + 1]);
        }

        A.setFromTriplets(triplets.begin(), triplets.end());

        std::cout << "Matrix A:\n"
                  << Eigen::MatrixXd(A) << "\n";

        std::cout << "Determinant: " << Eigen::MatrixXd(A).determinant() << "\n";
    }

    template <typename T>
    void SolverFV<T>::assemble_B()
    {
        const Eigen::Index N = B.rows();
        const std::vector<std::pair<T, T>> QUICK_b = QUICK_coefficients_n(Space_mesh);
        const std::vector<T> vol_sizes = Space_mesh.get_volume_sizes();

        const Eigen::Index nnz = N + 2 * std::max(Eigen::Index{0}, N - 1) + std::max(Eigen::Index{0}, N - 2);
        std::vector<Eigen::Triplet<T>> triplets;
        triplets.reserve(nnz);

        auto main_coeff = [&](Eigen::Index i) -> T
        {
            return (QUICK_b[i + 1].first - T{1} + QUICK_b[i].first - QUICK_b[i].second) / vol_sizes[i];
        };

        auto sub1_coeff = [&](Eigen::Index i) -> T
        {
            return -QUICK_b[i].first / vol_sizes[i];
        };

        auto super1_coeff = [&](Eigen::Index i) -> T
        {
            return (T{1} - QUICK_b[i + 1].first + QUICK_b[i + 1].second + QUICK_b[i].second) / vol_sizes[i];
        };

        auto super2_coeff = [&](Eigen::Index i) -> T
        {
            return -QUICK_b[i].second / vol_sizes[i];
        };

        // First row (i = 0)
        const Eigen::Index i = 0;
        triplets.emplace_back(0, 0, (QUICK_b[i + 1].first - T{1}) / vol_sizes[i]);
        if (N > 1)
        {
            triplets.emplace_back(0, 1, (T{1} - QUICK_b[i + 1].first + QUICK_b[i + 1].second) / vol_sizes[i]);
            if (N > 2)
            {
                triplets.emplace_back(0, 2, -QUICK_b[i + 1].second / vol_sizes[i]);
            }
        }

        // Second row (i = 1)
        if (N > 1)
        {
            triplets.emplace_back(1, 0, sub1_coeff(1));
            triplets.emplace_back(1, 1, main_coeff(1));
            if (N > 2)
            {
                triplets.emplace_back(1, 2, super1_coeff(1));
                if (N > 3)
                {
                    triplets.emplace_back(1, 3, super2_coeff(1));
                }
            }
        }

        // Interior rows. Loop from 2 to N-3 (inclusive)
        const Eigen::Index loop_end = std::max(Eigen::Index{2}, N - 2);
        for (Eigen::Index i = 2; i < loop_end; ++i)
        {
            triplets.emplace_back(i, i - 1, sub1_coeff(i));
            triplets.emplace_back(i, i, main_coeff(i));
            triplets.emplace_back(i, i + 1, super1_coeff(i));
            triplets.emplace_back(i, i + 2, super2_coeff(i));
        }

        // Second-to-last row (i = N-2)
        if (N > 2)
        {
            const Eigen::Index i = N - 2;
            triplets.emplace_back(i, i - 1, sub1_coeff(i));
            triplets.emplace_back(i, i, main_coeff(i));
            if (N > 3)
            {
                triplets.emplace_back(i, i + 1, super1_coeff(i));
            }
        }

        // Last row (i = N-1) (Custom handling)
        if (N > 1)
        {
            const Eigen::Index i = N - 1;
            triplets.emplace_back(i, i - 1, sub1_coeff(i));
            triplets.emplace_back(i, i, main_coeff(i));
        }

        B.setFromTriplets(triplets.begin(), triplets.end());

        std::cout << "Matrix B:\n"
                  << Eigen::MatrixXd(B) << "\n";

        std::cout << "Determinant: " << Eigen::MatrixXd(B).determinant() << "\n";
    }

    template <typename T>
    void SolverFV<T>::assemble_R()
    {
        const Eigen::Vector<T, Eigen::Dynamic> rho = phys::compute_density(g, Velocity_mesh);
        const std::vector<T> vol_sizes = Space_mesh.get_volume_sizes();
        std::vector<Eigen::Triplet<T>> triplets;

        for (Eigen::Index i = 0; i < rho.size(); ++i)
        {
            triplets.emplace_back(i, i, T{2} / std::numbers::pi_v<T> * rho[i] / vol_sizes[i]);
        }

        R.setFromTriplets(triplets.begin(), triplets.end());

        std::cout << "Matrix R:\n"
                  << Eigen::MatrixXd(R) << "\n";

        std::cout << "Determinant: " << Eigen::MatrixXd(R).determinant() << "\n";
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

#endif /* SOLVERFV_F5CE4BC0_A11B_431A_A0BD_E176318E8B09 */
