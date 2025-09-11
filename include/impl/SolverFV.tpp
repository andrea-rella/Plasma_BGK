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

#ifndef SOLVERFV_D5765103_3BE7_45C8_A8B9_1958FA21E0F4
#define SOLVERFV_D5765103_3BE7_45C8_A8B9_1958FA21E0F4

#include "../SolverFV.hpp"
#include "phys_utils.hpp"
#include "numerics_utils.hpp"
#include <numbers>
#include <cmath>
#include <Eigen/SparseLU>

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
        R = Eigen::Vector<T, Eigen::Dynamic>::Zero(Space_mesh.get_N() + 1);

        // solution matrices
        g = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2 * Velocity_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        h = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2 * Velocity_mesh.get_N() + 1, Space_mesh.get_N() + 1);

        // density, mean velocity and temperature vectors
        density = Eigen::Vector<T, Eigen::Dynamic>::Zero(Space_mesh.get_N() + 1);
        mean_velocity = Eigen::Vector<T, Eigen::Dynamic>::Zero(Space_mesh.get_N() + 1);
        temperature = Eigen::Vector<T, Eigen::Dynamic>::Zero(Space_mesh.get_N() + 1);

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

        // RHS
        G = [](T z, T rho, T v, T Temp) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * rho * std::pow(Temp, T{-0.5}) * std::exp(-((z - v) * (z - v)) / Temp) * T{0};
        };

        H = [](T z, T rho, T v, T Temp) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * rho * std::sqrt(Temp) * std::exp(-((z - v) * (z - v)) / Temp) * T{0};
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

        return;
    }

    template <typename T>
    void SolverFV<T>::set_physical_quantities()
    {
        density = phys::compute_density(g, Velocity_mesh);
        mean_velocity = phys::compute_meanGasVelocity(g, Velocity_mesh, density);
        temperature = phys::compute_temperature(g, h, Velocity_mesh, density, mean_velocity);

        return;
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

        set_physical_quantities();

        return;
    }

    template <typename T>
    void SolverFV<T>::assemble_A()
    {
        const Eigen::Index N = A.rows();
        const std::vector<std::pair<T, T>> QUICK_a = numerics::QUICKcoefficients_p<T>(Space_mesh);
        const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();

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
            return (-QUICK_a[i + 2].second - T{1} + QUICK_a[i + 1].first - QUICK_a[i + 1].second) / vol_sizes[i + 1];
        };

        auto super1_coeff = [&](Eigen::Index i) -> T
        {
            return QUICK_a[i + 2].first / vol_sizes[i + 1];
        };

        auto sub2_coeff = [&](Eigen::Index i) -> T
        {
            return (QUICK_a[i + 1].second) / vol_sizes[i + 1];
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
                triplets.emplace_back(i, i - 2, QUICK_a[i + 1].second / vol_sizes[i + 1]);
            }
            triplets.emplace_back(i, i - 1, (QUICK_a[i + 1].first - T{1} - QUICK_a[i + 1].second) / vol_sizes[i + 1]);
            triplets.emplace_back(i, i, (T{1} - QUICK_a[i + 1].first) / vol_sizes[i + 1]);
        }

        A.setFromTriplets(triplets.begin(), triplets.end());

        return;
    }

    template <typename T>
    void SolverFV<T>::assemble_B()
    {
        const Eigen::Index N = B.rows();
        const std::vector<std::pair<T, T>> QUICK_b = numerics::QUICKcoefficients_n<T>(Space_mesh);
        const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();

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
            return -QUICK_b[i + 1].second / vol_sizes[i];
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

        return;
    }

    template <typename T>
    void SolverFV<T>::assemble_R()
    {
        const Eigen::Index NN = R.size();
        for (Eigen::Index i = 0; i < NN; ++i)
        {
            R[i] = (T{2} / std::sqrt(std::numbers::pi_v<T>)) * density[i];
        }
        return;
    }

    template <typename T>
    void SolverFV<T>::initialize()
    {
        std::cout << "Initializing SolverFV..." << std::endl;
        initializeMeshes();
        std::cout << "Meshes initialized." << std::endl;
        setInitialState();
        std::cout << "Initial state set." << std::endl;
        assemble_A();
        assemble_B();
        assemble_R();
        std::cout << "Numerical matrices assembled." << std::endl;
        std::cout << "Assembly complete." << std::endl;

        // Eigen::Vector<T, Eigen::Dynamic> Ones = Eigen::Vector<T, Eigen::Dynamic>::Ones(Space_mesh.get_N());
        //
        // std::cout << "B*1: \n"
        //          << (B * Ones).transpose() << std::endl;
        //
        // const size_t Space_N = Space_mesh.get_N();
        // const std::pair<T, T> bN = numerics::QUICKcoefficients_n_at<T>(Space_mesh, Space_N);
        // const std::pair<T, T> bNm1 = numerics::QUICKcoefficients_n_at<T>(Space_mesh, Space_N - 1);
        // const T sigmaNm1 = T{1} - bN.first + bN.second + bNm1.second;
        // const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();
        //
        // std::cout << "Correction N-1: " << (T{1} / vol_sizes[Space_N - 1]) * (bN.second - sigmaNm1) << std::endl;
        // std::cout << "Correction N-2: " << (T{1} / vol_sizes[Space_N - 2]) * bNm1.second << std::endl;

        is_initialized = true;
    }

    // ------ SOLVE HELPERS  -------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    Eigen::Vector<T, Eigen::Dynamic> SolverFV<T>::assemble_U_pos(size_t j, T a1_2, T a2_2, T omega1) const
    {
        Eigen::Vector<T, Eigen::Dynamic> U_p(Space_mesh.get_N());
        const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();

        for (Eigen::Index i = 0; i < U_p.size(); ++i)
        {
            U_p[i] = T{2} / std::sqrt(std::numbers::pi_v<T>) * density[i + 1] *
                     G(Velocity_mesh[j], density[i + 1], mean_velocity[i + 1], temperature[i + 1]);
        }

        U_p[0] -= (Velocity_mesh[j] / vol_sizes[1]) * (a1_2 + omega1) * g0(Velocity_mesh[j]);
        U_p[1] -= (Velocity_mesh[j] / vol_sizes[2]) * a2_2 * g0(Velocity_mesh[j]);

        return U_p;
    }

    template <typename T>
    Eigen::Vector<T, Eigen::Dynamic> SolverFV<T>::assemble_W_pos(size_t j, T a1_2, T a2_2, T omega1) const
    {
        Eigen::Vector<T, Eigen::Dynamic> W_p(Space_mesh.get_N());
        const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();

        for (Eigen::Index i = 0; i < W_p.size(); ++i)
        {
            W_p[i] = T{2} / std::sqrt(std::numbers::pi_v<T>) * density[i + 1] *
                     H(Velocity_mesh[j], density[i + 1], mean_velocity[i + 1], temperature[i + 1]);
        }

        // boundary correction for the first two entries
        W_p[0] -= (Velocity_mesh[j] / vol_sizes[1]) * (a1_2 + omega1) * h0(Velocity_mesh[j]);
        W_p[1] -= (Velocity_mesh[j] / vol_sizes[2]) * a2_2 * h0(Velocity_mesh[j]);

        return W_p;
    }

    template <typename T>
    Eigen::Vector<T, Eigen::Dynamic> SolverFV<T>::assemble_U_zero() const
    {
        Eigen::Vector<T, Eigen::Dynamic> U_0(Space_mesh.get_N());

        for (Eigen::Index i = 0; i < U_0.size(); ++i)
        {
            U_0[i] = T{2} / std::sqrt(std::numbers::pi_v<T>) * density[i] *
                     G(T{0}, density[i], mean_velocity[i], temperature[i]);
        }

        return U_0;
    }

    template <typename T>
    Eigen::Vector<T, Eigen::Dynamic> SolverFV<T>::assemble_W_zero() const
    {
        Eigen::Vector<T, Eigen::Dynamic> W_0(Space_mesh.get_N());

        for (Eigen::Index i = 0; i < W_0.size(); ++i)
        {
            W_0[i] = T{2} / std::sqrt(std::numbers::pi_v<T>) * density[i] *
                     H(T{0}, density[i], mean_velocity[i], temperature[i]);
        }

        return W_0;
    }

    template <typename T>
    Eigen::Vector<T, Eigen::Dynamic> SolverFV<T>::assemble_U_neg(size_t j, T bN_2, T bNm1_2, T sigmaNm1) const
    {
        const size_t Space_N = Space_mesh.get_N();
        Eigen::Vector<T, Eigen::Dynamic> U_m(Space_N);
        const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();

        for (Eigen::Index i = 0; i < U_m.size(); ++i)
        {
            U_m[i] = T{2} / std::sqrt(std::numbers::pi_v<T>) * density[i] *
                     G(Velocity_mesh[j], density[i], mean_velocity[i], temperature[i]);
        }

        // boundary correction for the last two entries
        U_m[Space_N - 1] -= (Velocity_mesh[j] / vol_sizes[Space_N - 1]) * (-bN_2 + sigmaNm1) * g_infty(Velocity_mesh[j]);
        U_m[Space_N - 2] -= (Velocity_mesh[j] / vol_sizes[Space_N - 2]) * (-bNm1_2) * g_infty(Velocity_mesh[j]);

        return U_m;
    }

    template <typename T>
    Eigen::Vector<T, Eigen::Dynamic> SolverFV<T>::assemble_W_neg(size_t j, T bN_2, T bNm1_2, T sigmaNm1) const
    {
        const size_t Space_N = Space_mesh.get_N();
        Eigen::Vector<T, Eigen::Dynamic> W_m(Space_N);
        const std::vector<T> &vol_sizes = Space_mesh.get_volume_sizes();

        for (Eigen::Index i = 0; i < W_m.size(); ++i)
        {
            W_m[i] = T{2} / std::sqrt(std::numbers::pi_v<T>) * density[i] *
                     H(Velocity_mesh[j], density[i], mean_velocity[i], temperature[i]);
        }

        // boundary correction for the last two entries
        W_m[Space_N - 1] -= (Velocity_mesh[j] / vol_sizes[Space_N - 1]) * (-bN_2 + sigmaNm1) * h_infty(Velocity_mesh[j]);
        W_m[Space_N - 2] -= (Velocity_mesh[j] / vol_sizes[Space_N - 2]) * (-bNm1_2) * h_infty(Velocity_mesh[j]);

        return W_m;
    }

    template <typename T>
    void SolverFV<T>::solve_timestep_pos()
    {
        const size_t Velocity_N = Velocity_mesh.get_N();
        const size_t Space_N = Space_mesh.get_N();
        const T dt = Data.get_dt();

        // Positive velocities are indices Velocity_N+1 .. 2*Velocity_N
        const size_t j_begin = Velocity_N + 1;
        const size_t j_end = 2 * Velocity_N; // inclusive

        const std::pair<T, T> a1 = numerics::QUICKcoefficients_p_at<T>(Space_mesh, 1);
        const std::pair<T, T> a2 = numerics::QUICKcoefficients_p_at<T>(Space_mesh, 2);
        const T omega1 = -a2.second - T{1} + a1.first - a1.second;

        // Build diagonal (R part) once
        Eigen::Vector<T, Eigen::Dynamic> R_loc = R.tail(Space_N);

        // Ensure every diagonal entry exists (so adding identity is fast & pattern stable)
        // for (Eigen::Index i = 0; i < C.rows(); ++i)
        //  C.coeffRef(i, i) += T{0};
        // C.makeCompressed();

        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> solver;
        solver.analyzePattern(A);

        // Reusable buffers
        Eigen::Vector<T, Eigen::Dynamic> U(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> W(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> rhs_g(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> rhs_h(Space_N);
        Eigen::SparseMatrix<T> M; // will hold I + alpha * C each iteration

        for (size_t j = j_begin; j <= j_end; ++j)
        {
            U = assemble_U_pos(j, a1.second, a2.second, omega1);
            W = assemble_W_pos(j, a1.second, a2.second, omega1);

            Eigen::Vector<T, Eigen::Dynamic> g_j = g.row(j).tail(Space_N);
            Eigen::Vector<T, Eigen::Dynamic> h_j = h.row(j).tail(Space_N);

            // Rebuild numeric values of M: M = v_j*dt*A + dt*diag(R_loc) + I
            M = A;                        // copy A pattern & values
            M *= (Velocity_mesh[j] * dt); // scale A
            // Add dt*R and identity on diagonal
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(Space_N); ++i)
                // M.coeffRef(i, i) += dt * R_loc[i] + T{1}; //  You can avoid a pattern merge and one temporary by doing in‑place operations.
                M.coeffRef(i, i) += T{1};
            // (C already compressed; coeffRef keeps pattern stable)
            /**
             * @note Doing this for loop (modifying in place) is more efficient than summing directly the
             *       identity iff every diagonal entry already exists (which it does by contruction). Infact:
             *
             *       - Expression I + alpha*C triggers a sparse addition: it allocates a fresh matrix, merges
             *         two patterns (even if Identity is trivial), and copies all nnz.
             *
             *       - In‑place: copy C -> M (already compressed), scale (single pass over nnz), then increment
             *         N diagonal scalars. No pattern merge.
             *
             *       Caveat: Using coeffRef(i,i) does a (cheap) binary search per column; if the diagonal entry
             *       does NOT exist it INSERTS (costly, breaks pattern reuse). So ensure diagonal entries exist beforehand
             *       (e.g. once: C.coeffRef(i,i) += 0; C.makeCompressed()).
             */

            solver.factorize(M);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("LU factorization failed in solve_timestep_pos");

            rhs_g = g_j + dt * U;
            rhs_h = h_j + dt * W;

            auto x = solver.solve(rhs_g);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (g) failed in solve_timestep_pos");
            auto y = solver.solve(rhs_h);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (h) failed in solve_timestep_pos");

            g.row(j).tail(Space_N) = x.transpose();
            h.row(j).tail(Space_N) = y.transpose();
        }

        return;
    }

    template <typename T>
    void SolverFV<T>::solve_timestep_neg()
    {
        // Implementation for negative velocity timestep
        const size_t Velocity_N = Velocity_mesh.get_N();
        const size_t Space_N = Space_mesh.get_N();
        const T dt = Data.get_dt();

        const size_t j_begin = 0;
        const size_t j_end = Velocity_N - 1; // inclusive

        const std::pair<T, T> bN = numerics::QUICKcoefficients_n_at<T>(Space_mesh, Space_N);
        const std::pair<T, T> bNm1 = numerics::QUICKcoefficients_n_at<T>(Space_mesh, Space_N - 1);
        const T sigmaNm1 = T{1} - bN.first + bN.second + bNm1.second;

        // Precompute constant part C = A + R_mat
        Eigen::Vector<T, Eigen::Dynamic> R_loc = R.head(Space_N);

        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> solver;
        solver.analyzePattern(B);

        // Reusable buffers
        Eigen::Vector<T, Eigen::Dynamic> U(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> W(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> rhs_g(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> rhs_h(Space_N);
        Eigen::SparseMatrix<T> M; // will hold I + alpha * C each iteration

        for (size_t j = j_begin; j <= j_end; ++j)
        {
            U = assemble_U_neg(j, bN.second, bNm1.second, sigmaNm1);
            W = assemble_W_neg(j, bN.second, bNm1.second, sigmaNm1);

            Eigen::Vector<T, Eigen::Dynamic> g_j = g.row(j).head(Space_N);
            Eigen::Vector<T, Eigen::Dynamic> h_j = h.row(j).head(Space_N);

            M = B;                        // copy B pattern & values
            M *= (Velocity_mesh[j] * dt); // scale B
            // Add dt*R and identity on diagonal
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(Space_N); ++i)
                // M.coeffRef(i, i) += dt * R_loc[i] + T{1};
                M.coeffRef(i, i) += T{1};

            solver.factorize(M);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("LU factorization failed in solve_timestep_pos");

            rhs_g = g_j + dt * U;
            rhs_h = h_j + dt * W;

            auto x = solver.solve(rhs_g);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (g) failed in solve_timestep_pos");
            auto y = solver.solve(rhs_h);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (h) failed in solve_timestep_pos");

            g.row(j).head(Space_N) = x.transpose();
            h.row(j).head(Space_N) = y.transpose();
        }
    }

    template <typename T>
    void SolverFV<T>::solve_timestep_zero()
    {
        const size_t Velocity_N = Velocity_mesh.get_N();
        const size_t Space_N = Space_mesh.get_N();
        const T dt = Data.get_dt();

        // Assemble sources
        Eigen::Vector<T, Eigen::Dynamic> U = assemble_U_zero(); // size = Space_N
        Eigen::Vector<T, Eigen::Dynamic> W = assemble_W_zero(); // size = Space_N

        Eigen::Vector<T, Eigen::Dynamic> g_j = g.row(Velocity_N).head(Space_N); // copy (row vector -> column vector)
        Eigen::Vector<T, Eigen::Dynamic> h_j = h.row(Velocity_N).head(Space_N);

        Eigen::Vector<T, Eigen::Dynamic> rhs_g = g_j + dt * U;
        Eigen::Vector<T, Eigen::Dynamic> rhs_h = h_j + dt * W;

        // Diagonal solve: (I + dt * R) x = rhs
        // R has size Space_N + 1 (matches zero-velocity row length)
        Eigen::Vector<T, Eigen::Dynamic> denom = Eigen::Vector<T, Eigen::Dynamic>::Ones(Space_N);
        // denom.noalias() += dt * R.head(Space_N); // element-wise: 1 + dt*R_i

        Eigen::Vector<T, Eigen::Dynamic> x = rhs_g.cwiseQuotient(denom);
        Eigen::Vector<T, Eigen::Dynamic> y = rhs_h.cwiseQuotient(denom);

        g.row(Velocity_N).head(Space_N) = x.transpose();
        h.row(Velocity_N).head(Space_N) = y.transpose();
    }

    // ------ SOLVE  ---------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    void SolverFV<T>::solve(const metrics::VectorNormType vec_norm_type,
                            const metrics::RowAggregateType agg_type)
    {

        if (!is_initialized)
            initialize();

        size_t max_iter = Data.get_max_iter();
        size_t k = 0;
        T tol = Data.get_tol();
        T rel_err = std::numeric_limits<T>::max();
        auto mat_norm = metrics::MatrixNormFactory<T>::create(vec_norm_type, agg_type);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> g_old, h_old;

        std::cout << "Starting solver..." << std::endl;
        while (k < max_iter && rel_err > tol)
        {
            g_old = g;
            h_old = h;

            solve_timestep_neg();
            solve_timestep_zero();
            solve_timestep_pos();

            set_physical_quantities();

            assemble_R();

            rel_err = std::sqrt(std::pow(mat_norm->compute(g, g_old, Space_mesh.get_volume_sizes()), T{2}) +
                                std::pow(mat_norm->compute(h, h_old, Space_mesh.get_volume_sizes()), T{2}));
            ++k;
        }

        std::cout << "Solver finished after " << k << " iterations with relative error " << rel_err << std::endl;
        return;
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

        return;
    }

    template <typename T>
    void SolverFV<T>::write_phys_txt(const std::string &folder_name) const
    {
        std::filesystem::create_directories("output/" + folder_name);

        // Write physical quantities (four columns: index, density, mean_velocity, temperature)
        std::string filename_phys = "output/" + folder_name + "/physical_quantities.txt";
        std::ofstream txt_file_phys(filename_phys);
        if (!txt_file_phys.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_phys << std::endl;
            return;
        }

        const Eigen::Index n = density.size();
        txt_file_phys << "Physical quantities at spatial grid points:\n";
        txt_file_phys << "# index\tDensity\tMean Velocity\tTemperature\n";
        txt_file_phys << "--------------------------------------------------------\n";
        for (Eigen::Index i = 0; i < n; ++i)
        {
            txt_file_phys << i << '\t' << density[i] << '\t' << mean_velocity[i] << '\t' << temperature[i] << '\n';
        }

        return;
    }

    template <typename T>
    void SolverFV<T>::write_meshes_txt(const std::string &folder_name) const
    {
        Space_mesh.write_mesh_txt(folder_name);
        Velocity_mesh.write_mesh_txt(folder_name);
        return;
    }

    template <typename T>
    void SolverFV<T>::write_space_mesh_vtk(const std::string &folder_name) const
    {
        Space_mesh.write_mesh_vtk(folder_name);
        return;
    }

    template <typename T>
    void SolverFV<T>::write_initial_state_txt(const std::string &folder_name) const
    {
        std::filesystem::create_directories("output/" + folder_name);

        const size_t Space_N = Space_mesh.get_N();
        const size_t Velocity_N = Velocity_mesh.get_N();

        // Write initial g boundary values
        std::string filename_g0 = "output/" + folder_name + "/initial_state_g.txt";
        std::ofstream txt_file_g0(filename_g0);
        if (!txt_file_g0.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_g0 << std::endl;
            return;
        }

        txt_file_g0 << "Initial g boundary values:\n";
        txt_file_g0 << "Negative velocities (j = 0 .. " << Velocity_N - 1 << "): g(Space_N, j)\n";
        txt_file_g0 << "Positive velocities (j = " << Velocity_N + 1 << " .. " << 2 * Velocity_N << "): g(0, j)\n";
        txt_file_g0 << "--------------------------------------------------------\n";

        // Negative velocities
        for (size_t j = 0; j < Velocity_N; ++j)
        {
            txt_file_g0 << "g(" << Space_N << ", " << j << ") = " << g(j, Space_N) << "\n";
        }
        // zero velocity
        txt_file_g0 << "g(" << Space_N << ", " << Velocity_N << ") = " << g(Velocity_N, Space_N) << "\n";
        // Positive velocities
        for (size_t j = Velocity_N + 1; j <= 2 * Velocity_N; ++j)
        {
            txt_file_g0 << "g(0, " << j << ") = " << g(j, 0) << "\n";
        }

        // Write initial h boundary values
        std::string filename_h0 = "output/" + folder_name + "/initial_state_h.txt";
        std::ofstream txt_file_h0(filename_h0);
        if (!txt_file_h0.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_h0 << std::endl;
            return;
        }

        txt_file_h0 << "Initial h boundary values:\n";
        txt_file_h0 << "Negative velocities (j = 0 .. " << Velocity_N - 1 << "): h(Space_N, j)\n";
        txt_file_h0 << "Positive velocities (j = " << Velocity_N + 1 << " .. " << 2 * Velocity_N << "): h(0, j)\n";
        txt_file_h0 << "--------------------------------------------------------\n";

        // Negative velocities
        for (size_t j = 0; j < Velocity_N; ++j)
        {
            txt_file_h0 << "h(" << Space_N << ", " << j << ") = " << h(j, Space_N) << "\n";
        }
        // zero velocity
        txt_file_h0 << "h(" << Space_N << ", " << Velocity_N << ") = " << h(Velocity_N, Space_N) << "\n";
        // Positive velocities
        for (size_t j = Velocity_N + 1; j <= 2 * Velocity_N; ++j)
        {
            txt_file_h0 << "h(0, " << j << ") = " << h(j, 0) << "\n";
        }

        return;
    }

    template <typename T>
    void SolverFV<T>::write_all(const std::string &folder_name) const
    {
        write_sol_txt(folder_name);
        write_phys_txt(folder_name);
        write_meshes_txt(folder_name);
        write_space_mesh_vtk(folder_name);
        return;
    }
}

#endif /* SOLVERFV_D5765103_3BE7_45C8_A8B9_1958FA21E0F4 */
