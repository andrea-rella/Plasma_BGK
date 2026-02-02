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

#ifndef SOLVERFV_E0E66268_B084_4775_9396_92BB528BD61E
#define SOLVERFV_E0E66268_B084_4775_9396_92BB528BD61E

#include "../SolverFV.hpp"
#include "phys_utils.hpp"
#include "numerics_utils.hpp"
#include <numbers>
#include <cmath>
#include <Eigen/SparseLU>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP is not enabled. The parallel solver (if called) will run sequentially."
#endif

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
        g = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(2 * Velocity_mesh.get_N() + 1, Space_mesh.get_N() + 1);
        h = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(2 * Velocity_mesh.get_N() + 1, Space_mesh.get_N() + 1);

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
                   std::exp(-std::pow((z - std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        h0 = [](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * std::exp(-z * z);
        };

        h_infty = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-0.5}) *
                   std::exp(-std::pow((z - std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        // initial conditions for g and h
        g_init = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-1.5}) *
                   std::exp(-std::pow((z - std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        h_init = [T_infty_w, p_infty_w, M_infty](T z) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * p_infty_w * std::pow(T_infty_w, T{-0.5}) *
                   std::exp(-std::pow((z - std::sqrt(T{5} / T{6} * T_infty_w) * M_infty), T{2}) * (T{1} / T_infty_w));
        };

        // RHS
        G = [](T z, T rho, T v, T Temp) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * rho * std::pow(Temp, T{-0.5}) * std::exp(-((z - v) * (z - v)) / Temp);
        };

        H = [](T z, T rho, T v, T Temp) -> T
        {
            return T{1} / std::sqrt(std::numbers::pi_v<T>) * rho * std::sqrt(Temp) * std::exp(-((z - v) * (z - v)) / Temp);
        };
        return;
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
        auto get_jacobian = [&](size_t k) -> T
        {
            // Cast to signed integer to handle negative relative indices correctly
            long long j_signed = static_cast<long long>(k) - static_cast<long long>(Velocity_mesh.get_N());
            return Data.get_a1() + 3.0 * Data.get_a2() * (static_cast<T>(j_signed) * static_cast<T>(j_signed));
        };

        density = phys::compute_density(g, Velocity_mesh, g0(0), get_jacobian);
        mean_velocity = phys::compute_meanGasVelocity(g, Velocity_mesh, density, g0(0), get_jacobian);
        temperature = phys::compute_temperature(g, h, Velocity_mesh, density, mean_velocity, g0(0), h0(0), get_jacobian);

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

        A.makeCompressed();
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

        B.makeCompressed();

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
        std::cout << "Assembly complete. \n"
                  << std::endl;

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

        // Positive velocities indices
        const size_t j_begin = Velocity_N + 1;
        const size_t j_end = 2 * Velocity_N;

        // Pre-calculate coefficients once
        const std::pair<T, T> a1 = numerics::QUICKcoefficients_p_at<T>(Space_mesh, 1);
        const std::pair<T, T> a2 = numerics::QUICKcoefficients_p_at<T>(Space_mesh, 2);
        const T omega1 = -a2.second - T{1} + a1.first - a1.second;

        // --- Pre-compute diagonal offsets ---
        // This maps where the diagonal elements live in the raw value array.
        // Complexity: One-time O(NNZ), saves O(N log k) per timestep.
        std::vector<ptrdiff_t> diag_offsets(Space_N);

        // Ensure A is compressed to guarantee valuePtr safety
        // (Assuming A is already built; if not, call A.makeCompressed())
        if (!A.isCompressed())
            A.makeCompressed();

        for (int k = 0; k < A.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
            {
                if (it.row() == it.col())
                {
                    // Store the offset from the beginning of the value array
                    diag_offsets[it.row()] = &it.value() - A.valuePtr();
                }
            }
        }

        // --- Structure Reuse ---
        // Initialize M with A's pattern ONCE. We never change the pattern, only values.
        Eigen::SparseMatrix<T> M = A;

        // Analyze pattern once (Symbolic Factorization)
        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> solver;
        solver.analyzePattern(A);

        // Buffers
        Eigen::Vector<T, Eigen::Dynamic> U(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> W(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> rhs(Space_N); // Reused for both g and h

        // Pre-fetch constant vector parts
        const Eigen::Vector<T, Eigen::Dynamic> R_loc = R.tail(Space_N);
        const T *A_vals = A.valuePtr(); // Read-only source
        T *M_vals = M.valuePtr();       // Write target
        const Eigen::Index nnz = A.nonZeros();

        for (size_t j = j_begin; j <= j_end; ++j)
        {
            // 1. Matrix Assembly: M = (v_j * dt) * A + (I + dt * R)
            const T v_dt = Velocity_mesh[j] * dt;

            // Step A: Reset M values by scaling A (Vectorized copy)
            // This is much faster than M = A (which copies indices too) followed by M *= scalar
            for (Eigen::Index k = 0; k < nnz; ++k)
            {
                M_vals[k] = A_vals[k] * v_dt;
            }

            // Step B: Update diagonal (Direct access via offsets)
            for (size_t i = 0; i < Space_N; ++i)
            {
                M_vals[diag_offsets[i]] += (dt * R_loc[i] + T{1});
            }

            // 2. Numerical Factorization
            solver.factorize(M);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("LU factorization failed in solve_timestep_pos");

            // 3. Assemble RHS vectors
            U = assemble_U_pos(j, a1.second, a2.second, omega1);
            W = assemble_W_pos(j, a1.second, a2.second, omega1);

            // 4. Solve for g
            rhs = g.row(j).tail(Space_N).transpose();
            rhs += dt * U;

            auto x = solver.solve(rhs);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (g) failed");
            g.row(j).tail(Space_N) = x.transpose();

            // 5. Solve for h
            rhs = h.row(j).tail(Space_N).transpose();
            rhs += dt * W;

            auto y = solver.solve(rhs);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (h) failed");
            h.row(j).tail(Space_N) = y.transpose();
        }
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

        // --- Optimization 1: Pre-compute diagonal offsets for B ---
        std::vector<ptrdiff_t> diag_offsets(Space_N);

        // Ensure B is compressed to guarantee valuePtr safety
        if (!B.isCompressed())
            B.makeCompressed();

        for (int k = 0; k < B.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(B, k); it; ++it)
            {
                if (it.row() == it.col())
                {
                    // Correct syntax: address of value reference
                    diag_offsets[it.row()] = &it.value() - B.valuePtr();
                }
            }
        }

        // --- Optimization 2: Structure Reuse ---
        // Initialize M with B's pattern ONCE.
        Eigen::SparseMatrix<T> M = B;

        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> solver;
        solver.analyzePattern(B);

        // Reusable buffers
        Eigen::Vector<T, Eigen::Dynamic> U(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> W(Space_N);
        Eigen::Vector<T, Eigen::Dynamic> rhs(Space_N); // Reused for both g and h

        // Pre-fetch constant vector parts and pointers
        const Eigen::Vector<T, Eigen::Dynamic> R_loc = R.head(Space_N);
        const T *B_vals = B.valuePtr(); // Read-only source
        T *M_vals = M.valuePtr();       // Write target
        const Eigen::Index nnz = B.nonZeros();

        for (size_t j = j_begin; j <= j_end; ++j)
        {
            // 1. Matrix Assembly: M = (v_j * dt) * B + (I + dt * R)
            const T v_dt = Velocity_mesh[j] * dt;

            // Step A: Reset M values by scaling B (Vectorized copy)
            for (Eigen::Index k = 0; k < nnz; ++k)
            {
                M_vals[k] = B_vals[k] * v_dt;
            }

            // Step B: Update diagonal (Direct access via offsets)
            for (size_t i = 0; i < Space_N; ++i)
            {
                M_vals[diag_offsets[i]] += (dt * R_loc[i] + T{1});
            }

            // 2. Numerical Factorization
            solver.factorize(M);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("LU factorization failed in solve_timestep_neg");

            // 3. Assemble RHS vectors
            U = assemble_U_neg(j, bN.second, bNm1.second, sigmaNm1);
            W = assemble_W_neg(j, bN.second, bNm1.second, sigmaNm1);

            // 4. Solve for g
            // rhs = g_j + dt * U (Avoiding extra allocations)
            rhs = g.row(j).head(Space_N).transpose();
            rhs += dt * U;

            auto x = solver.solve(rhs);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (g) failed in solve_timestep_neg");
            g.row(j).head(Space_N) = x.transpose();

            // 5. Solve for h
            rhs = h.row(j).head(Space_N).transpose();
            rhs += dt * W;

            auto y = solver.solve(rhs);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Solve (h) failed in solve_timestep_neg");
            h.row(j).head(Space_N) = y.transpose();
        }
    }

    template <typename T>
    void SolverFV<T>::solve_timestep_zero()
    {
        const size_t Velocity_N = Velocity_mesh.get_N();
        const size_t Space_N = Space_mesh.get_N();
        const T dt = Data.get_dt();

        // 1. Assemble sources
        Eigen::Vector<T, Eigen::Dynamic> U = assemble_U_zero();
        Eigen::Vector<T, Eigen::Dynamic> W = assemble_W_zero();

        // 2. Pre-calculate Inverse Denominator
        Eigen::Vector<T, Eigen::Dynamic> inv_denom =
            (Eigen::Vector<T, Eigen::Dynamic>::Ones(Space_N) + dt * R.head(Space_N)).cwiseInverse();

        // 3. Update 'g' in-place (No intermediate allocations)
        g.row(Velocity_N).head(Space_N) = (g.row(Velocity_N).head(Space_N).transpose() + dt * U)
                                              .cwiseProduct(inv_denom)
                                              .transpose();

        // 4. Update 'h' in-place
        h.row(Velocity_N).head(Space_N) = (h.row(Velocity_N).head(Space_N).transpose() + dt * W)
                                              .cwiseProduct(inv_denom)
                                              .transpose();
    }

    // ---- PARALLEL -----
    template <typename T>
    void SolverFV<T>::solve_timestep_pos_parallel()
    {
        const size_t Velocity_N = Velocity_mesh.get_N();
        const size_t Space_N = Space_mesh.get_N();
        const T dt = Data.get_dt();

        // Positive velocities indices
        const size_t j_begin = Velocity_N + 1;
        const size_t j_end = 2 * Velocity_N;

        // Pre-calculate coefficients once (Shared)
        const std::pair<T, T> a1 = numerics::QUICKcoefficients_p_at<T>(Space_mesh, 1);
        const std::pair<T, T> a2 = numerics::QUICKcoefficients_p_at<T>(Space_mesh, 2);
        const T omega1 = -a2.second - T{1} + a1.first - a1.second;

        // Ensure A is compressed for raw pointer access
        if (!A.isCompressed())
            A.makeCompressed();

        // --- Pre-compute diagonal offsets (Shared) ---
        // This map is valid for any matrix with A's pattern (including thread-local M)
        std::vector<ptrdiff_t> diag_offsets(Space_N);
        const T *A_global_ptr = A.valuePtr(); // Pointer to shared A values

        for (int k = 0; k < A.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
            {
                if (it.row() == it.col())
                {
                    diag_offsets[it.row()] = &it.value() - A_global_ptr;
                }
            }
        }

        // Shared pointers/refs for loop access
        const Eigen::Vector<T, Eigen::Dynamic> &R_ref = R;
        const Eigen::Index nnz = A.nonZeros();

        // Exception capture for OpenMP
        std::atomic<bool> error_occurred{false};

// --- Parallel Region ---
// Note: We create thread-local resources HERE, outside the loop
#pragma omp parallel
        {
            // 1. Thread-Local Allocations
            Eigen::SparseMatrix<T> M_loc = A;

            // Solver instance for this thread
            Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> solver_loc;
            solver_loc.analyzePattern(A); // Analyze pattern ONCE per thread

            // Thread-local Buffers
            Eigen::Vector<T, Eigen::Dynamic> U_loc(Space_N);
            Eigen::Vector<T, Eigen::Dynamic> W_loc(Space_N);
            Eigen::Vector<T, Eigen::Dynamic> rhs_loc(Space_N);

            // Pointers for fast access
            T *M_vals = M_loc.valuePtr();
            const Eigen::Vector<T, Eigen::Dynamic> R_loc = R_ref.tail(Space_N);

// 2. Parallel Loop
// schedule(dynamic) is usually better if solve times vary,
// otherwise schedule(static) has lower overhead.
#pragma omp for schedule(static)
            for (size_t j = j_begin; j <= j_end; ++j)
            {
                if (error_occurred)
                    continue; // Skip work if another thread failed

                try
                {
                    // --- Matrix Assembly ---
                    const T v_dt = Velocity_mesh[j] * dt;

                    // Step A: Vectorized reset using raw pointers
                    // We use A_global_ptr (Shared read-only) to reset M_vals (Thread-local)
                    for (Eigen::Index k = 0; k < nnz; ++k)
                    {
                        M_vals[k] = A_global_ptr[k] * v_dt;
                    }

                    // Step B: Update diagonal using pre-calc offsets
                    for (size_t i = 0; i < Space_N; ++i)
                    {
                        M_vals[diag_offsets[i]] += (dt * R_loc[i] + T{1});
                    }

                    // --- Factorization ---
                    solver_loc.factorize(M_loc);
                    if (solver_loc.info() != Eigen::Success)
                        throw std::runtime_error("LU factorization failed");

                    // --- RHS Assembly ---
                    U_loc = assemble_U_pos(j, a1.second, a2.second, omega1);
                    W_loc = assemble_W_pos(j, a1.second, a2.second, omega1);

                    // --- Solve for g ---
                    // @note: .row(j) write is thread-safe as j is unique per thread
                    rhs_loc = g.row(j).tail(Space_N).transpose();
                    rhs_loc += dt * U_loc;

                    auto x = solver_loc.solve(rhs_loc);
                    if (solver_loc.info() != Eigen::Success)
                        throw std::runtime_error("Solve g failed");
                    g.row(j).tail(Space_N) = x.transpose();

                    // --- Solve for h ---
                    rhs_loc = h.row(j).tail(Space_N).transpose();
                    rhs_loc += dt * W_loc;

                    auto y = solver_loc.solve(rhs_loc);
                    if (solver_loc.info() != Eigen::Success)
                        throw std::runtime_error("Solve h failed");
                    h.row(j).tail(Space_N) = y.transpose();
                }
                catch (...)
                {
                    error_occurred = true;
                }
            }
        } // End of parallel region

        if (error_occurred)
            throw std::runtime_error(error_message("An error occurred during parallel execution of solve_timestep_pos."));

        return;
    }

    template <typename T>
    void SolverFV<T>::solve_timestep_neg_parallel()
    {
        // Implementation for negative velocity timestep
        const size_t Velocity_N = Velocity_mesh.get_N();
        const size_t Space_N = Space_mesh.get_N();
        const T dt = Data.get_dt();

        const size_t j_begin = 0;
        const size_t j_end = Velocity_N - 1; // inclusive

        // Pre-calculate coefficients once (Shared)
        // This map is valid for any matrix with A's pattern (including thread-local M)
        const std::pair<T, T> bN = numerics::QUICKcoefficients_n_at<T>(Space_mesh, Space_N);
        const std::pair<T, T> bNm1 = numerics::QUICKcoefficients_n_at<T>(Space_mesh, Space_N - 1);
        const T sigmaNm1 = T{1} - bN.first + bN.second + bNm1.second;

        // Ensure B is compressed to guarantee valuePtr safety
        if (!B.isCompressed())
            B.makeCompressed();

        // --- Pre-compute diagonal offsets for B (Shared) ---
        std::vector<ptrdiff_t> diag_offsets(Space_N);
        const T *B_global_ptr = B.valuePtr(); // Pointer to shared B values

        for (int k = 0; k < B.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(B, k); it; ++it)
            {
                if (it.row() == it.col())
                {
                    diag_offsets[it.row()] = &it.value() - B_global_ptr;
                }
            }
        }

        // Shared references for loop access
        const Eigen::Vector<T, Eigen::Dynamic> &R_ref = R;
        const Eigen::Index nnz = B.nonZeros();

        // Exception capture for OpenMP
        std::atomic<bool> error_occurred{false};

// --- Parallel Region ---
// Allocating heavy resources (Matrix copy, Solver, Buffers) ONCE per thread
#pragma omp parallel
        {
            // 1. Thread-Local Allocations
            Eigen::SparseMatrix<T> M_loc = B;

            // Solver instance for this thread
            Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> solver_loc;
            solver_loc.analyzePattern(B); // Analyze pattern ONCE per thread

            // Thread-local Buffers
            Eigen::Vector<T, Eigen::Dynamic> U_loc(Space_N);
            Eigen::Vector<T, Eigen::Dynamic> W_loc(Space_N);
            Eigen::Vector<T, Eigen::Dynamic> rhs_loc(Space_N);

            // Pointers for fast access
            T *M_vals = M_loc.valuePtr();                                       // Thread-local write target
            const Eigen::Vector<T, Eigen::Dynamic> R_loc = R_ref.head(Space_N); // Pre-fetch R head

// 2. Parallel Loop
#pragma omp for schedule(static)
            for (size_t j = j_begin; j <= j_end; ++j)
            {
                if (error_occurred)
                    continue;

                try
                {
                    // --- Matrix Assembly: M = (v_j * dt) * B + (I + dt * R) ---
                    const T v_dt = Velocity_mesh[j] * dt;

                    // Step A: Reset M values by scaling B (Vectorized copy using shared source)
                    for (Eigen::Index k = 0; k < nnz; ++k)
                    {
                        M_vals[k] = B_global_ptr[k] * v_dt;
                    }

                    // Step B: Update diagonal (Direct access via shared offsets)
                    for (size_t i = 0; i < Space_N; ++i)
                    {
                        M_vals[diag_offsets[i]] += (dt * R_loc[i] + T{1});
                    }

                    // --- Numerical Factorization ---
                    solver_loc.factorize(M_loc);
                    if (solver_loc.info() != Eigen::Success)
                        throw std::runtime_error("LU factorization failed in solve_timestep_neg");

                    // --- Assemble RHS vectors ---
                    U_loc = assemble_U_neg(j, bN.second, bNm1.second, sigmaNm1);
                    W_loc = assemble_W_neg(j, bN.second, bNm1.second, sigmaNm1);

                    // --- Solve for g ---
                    // rhs = g_j + dt * U
                    rhs_loc = g.row(j).head(Space_N).transpose();
                    rhs_loc += dt * U_loc;

                    auto x = solver_loc.solve(rhs_loc);
                    if (solver_loc.info() != Eigen::Success)
                        throw std::runtime_error("Solve (g) failed in solve_timestep_neg");
                    g.row(j).head(Space_N) = x.transpose();

                    // --- Solve for h ---
                    rhs_loc = h.row(j).head(Space_N).transpose();
                    rhs_loc += dt * W_loc;

                    auto y = solver_loc.solve(rhs_loc);
                    if (solver_loc.info() != Eigen::Success)
                        throw std::runtime_error("Solve (h) failed in solve_timestep_neg");
                    h.row(j).head(Space_N) = y.transpose();
                }
                catch (...)
                {
                    error_occurred = true;
                }
            }
        } // End of parallel region

        if (error_occurred)
            throw std::runtime_error(error_message("An error occurred during parallel execution of solve_timestep_neg."));

        return;
    }

    // ------ SOLVE  ---------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    template <PlotStrategy Strategy>
    void SolverFV<T>::solve(const metrics::VectorNormType vec_norm_type,
                            const metrics::RowAggregateType agg_type)
    {
        if (!is_initialized)
            initialize();

        size_t max_iter = Data.get_max_iter();
        size_t k = 0;
        T tol = Data.get_tol();
        T rel_err = std::numeric_limits<T>::max();

        size_t plot_every_k_steps = 0;
        if constexpr (Strategy == PlotStrategy::EACHSTEP)
        {
            plot_every_k_steps = Data.get_plot_every_k_steps();
            std::cout << "Plotting every " << plot_every_k_steps << " steps." << std::endl;
            write_phys_instant(Data.get_saving_folder_name(), 0);
        }

        auto mat_norm = metrics::MatrixNormFactory<T>::create(vec_norm_type, agg_type);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> g_old, h_old;

        while (k < max_iter && rel_err > tol)
        {
            g_old = g;
            h_old = h;

            solve_timestep_neg();
            solve_timestep_zero();
            solve_timestep_pos();

            set_physical_quantities();
            assemble_R();

            // Error calculation logic
            rel_err = std::sqrt(std::pow(mat_norm->compute(g, g_old, Space_mesh.get_volume_sizes()), T{2}) +
                                std::pow(mat_norm->compute(h, h_old, Space_mesh.get_volume_sizes()), T{2}));
            ++k;

            // Compile-time branching: Zero overhead if Strategy is OnlyEnd
            if constexpr (Strategy == PlotStrategy::EACHSTEP)
            {
                if (k % plot_every_k_steps == 0 || k == 1 || k == 20 || k == 40 || k == 60 || k == 80 || k == 100)
                {
                    write_phys_instant(Data.get_saving_folder_name(), k);
                }
            }
        }

        // Always plot at the end regardless of strategy
        write_phys_txt(Data.get_saving_folder_name());

        std::cout << "Solver finished after " << k << " iterations with relative error " << rel_err << std::endl;

        return;
    }

    template <typename T>
    template <PlotStrategy Strategy>
    void SolverFV<T>::solve_parallel(const metrics::VectorNormType vec_norm_type,
                                     const metrics::RowAggregateType agg_type)
    {
        if (!is_initialized)
            initialize();

        size_t max_iter = Data.get_max_iter();
        size_t k = 0;
        T tol = Data.get_tol();
        T rel_err = std::numeric_limits<T>::max();

        size_t plot_every_k_steps = 0;
        if constexpr (Strategy == PlotStrategy::EACHSTEP)
        {
            plot_every_k_steps = Data.get_plot_every_k_steps();
            std::cout << "Plotting every " << plot_every_k_steps << " steps." << std::endl;
            write_phys_instant(Data.get_saving_folder_name(), 0);
        }

        auto mat_norm = metrics::MatrixNormFactory<T>::create(vec_norm_type, agg_type);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> g_old, h_old;

        auto start_time = std::chrono::high_resolution_clock::now();

        while (k < max_iter && rel_err > tol)
        {
            g_old = g;
            h_old = h;

            solve_timestep_neg_parallel();
            solve_timestep_zero();
            solve_timestep_pos_parallel();

            set_physical_quantities();
            assemble_R();

            // Error calculation logic
            rel_err = std::sqrt(std::pow(mat_norm->compute(g, g_old, Space_mesh.get_volume_sizes()), T{2}) +
                                std::pow(mat_norm->compute(h, h_old, Space_mesh.get_volume_sizes()), T{2}));
            ++k;

            // Compile-time branching: Zero overhead if Strategy is OnlyEnd
            if constexpr (Strategy == PlotStrategy::EACHSTEP)
            {
                if (k % plot_every_k_steps == 0 || k == 1 || k == 20 || k == 40 || k == 60 || k == 80 || k == 100)
                {
                    write_phys_instant(Data.get_saving_folder_name(), k);
                }
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();

        // Always plot at the end regardless of strategy
        write_phys_txt(Data.get_saving_folder_name());

        std::cout << "Solver finished after " << k << " iterations with relative error " << rel_err << std::endl;

        std::chrono::duration<double> elapsed = end_time - start_time;
        std::cout << "Total elapsed time: " << elapsed.count() << " seconds." << std::endl;

        return;
    }
    // ------ OUTPUT ---------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    void SolverFV<T>::write_sol_txt(const std::string &folder_path) const
    {
        std::filesystem::create_directories(folder_path);

        // Write g solution
        std::string filename_g = folder_path + "/solution_g.txt";
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
        std::string filename_h = folder_path + "/solution_h.txt";
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
    void SolverFV<T>::write_sol_instant_txt(const std::string &folder_path, size_t iter) const
    {
        std::filesystem::create_directories(folder_path);

        // Write g solution snapshot
        std::string filename_g = folder_path + "/g_iter_" + std::to_string(iter) + ".txt";
        std::ofstream txt_file_g(filename_g);
        if (!txt_file_g.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_g << std::endl;
            return;
        }

        txt_file_g << "g solution in the phase space (iteration " << iter << "):\n";
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

        // Write h solution snapshot
        std::string filename_h = folder_path + "/h_iter_" + std::to_string(iter) + ".txt";
        std::ofstream txt_file_h(filename_h);
        if (!txt_file_h.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_h << std::endl;
            return;
        }

        txt_file_h << "h solution in the phase space (iteration " << iter << "):\n";
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
    void SolverFV<T>::write_phys_txt(const std::string &folder_path) const
    {
        std::filesystem::create_directories(folder_path);

        // Write physical quantities (four columns: index, density, mean_velocity, temperature)
        std::string filename_phys = folder_path + "/physical_quantities.txt";
        std::ofstream txt_file_phys(filename_phys);
        if (!txt_file_phys.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_phys << std::endl;
            return;
        }

        const Eigen::Index n = density.size();
        txt_file_phys << "Physical quantities at spatial grid points:\n";
        txt_file_phys << "# index -- Density -- Mean Velocity -- Temperature\n";
        txt_file_phys << "--------------------------------------------------------\n";
        for (Eigen::Index i = 0; i < n; ++i)
        {
            txt_file_phys << i << '\t' << density[i] << '\t' << mean_velocity[i] << '\t' << temperature[i] << '\n';
        }

        return;
    }

    template <typename T>
    void SolverFV<T>::write_meshes_txt(const std::string &folder_path) const
    {
        Space_mesh.write_mesh_txt(folder_path);
        Velocity_mesh.write_mesh_txt(folder_path);
        return;
    }

    template <typename T>
    void SolverFV<T>::write_space_mesh_vtk(const std::string &folder_path) const
    {
        Space_mesh.write_mesh_vtk(folder_path);
        return;
    }

    template <typename T>
    void SolverFV<T>::write_initial_state_txt(const std::string &folder_path) const
    {
        std::filesystem::create_directories(folder_path);

        const size_t Space_N = Space_mesh.get_N();
        const size_t Velocity_N = Velocity_mesh.get_N();

        // Write initial g boundary values
        std::string filename_g0 = folder_path + "/initial_state_g.txt";
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
        std::string filename_h0 = folder_path + "/initial_state_h.txt";
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
    void SolverFV<T>::write_phys_instant(const std::string &folder_path, size_t iter) const
    {
        std::filesystem::create_directories(folder_path);

        // Write temperature at current instant
        std::string filename_phys = folder_path + "/phys_iter_" + std::to_string(iter) + ".txt";
        std::ofstream txt_file_phys(filename_phys);
        if (!txt_file_phys.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename_phys << std::endl;
            return;
        }

        const Eigen::Index n = density.size();
        txt_file_phys << "Physical quantities at spatial grid points:\n";
        txt_file_phys << "# index -- Density -- Mean Velocity -- Temperature\n";
        txt_file_phys << "--------------------------------------------------------\n";
        for (Eigen::Index i = 0; i < n; ++i)
        {
            txt_file_phys << i << '\t' << density[i] << '\t' << mean_velocity[i] << '\t' << temperature[i] << '\n';
        }

        return;
    }

    template <typename T>
    void SolverFV<T>::write_all(const std::string &folder_path) const
    {
        write_sol_txt(folder_path);
        write_phys_txt(folder_path);
        write_meshes_txt(folder_path);
        write_space_mesh_vtk(folder_path);
        return;
    }
}

#endif /* SOLVERFV_E0E66268_B084_4775_9396_92BB528BD61E */
