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

#ifndef PHYS_UTILS_FB91EE41_406D_4600_A84B_EB367B78C2CB
#define PHYS_UTILS_FB91EE41_406D_4600_A84B_EB367B78C2CB

#include "../phys_utils.hpp"
#include "../utilities.hpp"

namespace Bgk
{
    namespace phys
    {

        // ------ GAS DENSITY----------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                         const VelocityMesh<T> &velocity_mesh,
                                                         T wall_g_val,
                                                         JacobianFunc get_jacobian)
        {
            Eigen::Vector<T, Eigen::Dynamic> densities(g.cols());
            if (g.cols() == 0)
            {
                return densities;
            }

            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            // --- Case i = 0 (Wall) ---
            {
                constexpr Eigen::Index i = 0;
                T integral_neg = T{0};
                T integral_pos = T{0};

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    integral_neg += weight * g(k, i) * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing) with discontinuity fix
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);
                    T f_val = (k == zero_idx) ? wall_g_val : g(k, i);
                    integral_pos += weight * f_val * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);

                densities.coeffRef(i) = integral_neg + integral_pos;
            }

            // --- Cases i > 0 (Away from Wall) ---
            for (Eigen::Index i = 1; i < g.cols(); ++i)
            {
                T integral_neg = T{0};
                T integral_pos = T{0};

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    integral_neg += weight * g(k, i) * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing)
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);
                    integral_pos += weight * g(k, i) * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);

                densities.coeffRef(i) = integral_neg + integral_pos;
            }

            return densities;
        }

        template <typename T, typename JacobianFunc>
        T compute_density_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                             const VelocityMesh<T> &velocity_mesh,
                             const Eigen::Index i,
                             T wall_g_val,
                             JacobianFunc get_jacobian)
        {
            if (i < 0 || i >= g.cols())
            {
                // This check from your old code is still good practice
                throw std::out_of_range("Index out of bounds when computing the density");
            }

            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            T integral_neg = T{0};
            T integral_pos = T{0};

            // 1. Negative Velocity Integral (Incoming)
            // This part is the same for i=0 and i>0
            for (size_t k = 0; k <= zero_idx; ++k)
            {
                T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                integral_neg += weight * g(k, i) * get_jacobian(k);
            }
            integral_neg *= (1.0 / 3.0);

            // 2. Positive Velocity Integral (Outgoing)
            // This part handles the discontinuity at the wall
            for (size_t k = zero_idx; k <= max_idx; ++k)
            {
                T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                // Apply the fix: use wall_maxwellian_val ONLY at i=0 and k=zero_idx
                T f_val = (i == 0 && k == zero_idx) ? wall_g_val : g(k, i);

                integral_pos += weight * f_val * get_jacobian(k);
            }
            integral_pos *= (1.0 / 3.0);

            return integral_neg + integral_pos;
        }

        // ------ MEAN VELOCITY   -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                                 T wall_g_val, // Value of g at wall for zeta=0+
                                                                 JacobianFunc get_jacobian)
        {
            Eigen::Vector<T, Eigen::Dynamic> velocities(g.cols());
            if (g.cols() == 0)
            {
                return velocities;
            }
            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            // --- Case i = 0 (Wall) ---
            {
                constexpr Eigen::Index i = 0;
                T integral_neg = T{0};
                T integral_pos = T{0};

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    integral_neg += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing) with discontinuity fix
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                    // Apply the fix: use wall_g_val at k=zero_idx
                    T g_val = (k == zero_idx) ? wall_g_val : g(k, i);

                    integral_pos += weight * (g_val * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                velocities.coeffRef(i) = (integral_neg + integral_pos) / densities(i);
            }

            // --- Cases i > 0 (Away from Wall) ---
            for (Eigen::Index i = 1; i < g.cols(); ++i)
            {
                T integral_neg = T{0};
                T integral_pos = T{0};

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    integral_neg += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing)
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);
                    integral_pos += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                velocities.coeffRef(i) = (integral_neg + integral_pos) / densities(i);
            }

            return velocities;
        }

        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 T wall_g_val, // For velocity
                                                                 JacobianFunc get_jacobian)
        {
            Eigen::Vector<T, Eigen::Dynamic> velocities(g.cols());
            if (g.cols() == 0)
            {
                return velocities;
            }

            // --- 1. COMPUTE DENSITIES (using vector function) ---
            Eigen::Vector<T, Eigen::Dynamic> densities =
                compute_density(g, velocity_mesh, wall_g_val, get_jacobian);

            // --- 2. COMPUTE MEAN VELOCITIES (full loop, using new logic) ---
            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            // --- Case i = 0 (Wall) ---
            {
                constexpr Eigen::Index i = 0;
                T integral_neg = T{0};
                T integral_pos = T{0};

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    integral_neg += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing) with discontinuity fix
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                    // Apply the fix: use wall_g_val at k=zero_idx
                    T g_val = (k == zero_idx) ? wall_g_val : g(k, i);

                    integral_pos += weight * (g_val * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                velocities.coeffRef(i) = (integral_neg + integral_pos) / densities(i);
            }

            // --- Cases i > 0 (Away from Wall) ---
            for (Eigen::Index i = 1; i < g.cols(); ++i)
            {
                T integral_neg = T{0};
                T integral_pos = T{0};

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    integral_neg += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing)
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);
                    integral_pos += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                velocities.coeffRef(i) = (integral_neg + integral_pos) / densities(i);
            }

            return velocities;
        }

        template <typename T, typename JacobianFunc>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     T wall_g_val, // For velocity calc
                                     JacobianFunc get_jacobian)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the mean gas velocity"));
            }

            // Call the "new" density function, which requires wall_maxwellian_val
            T density_i = compute_density_at(g, velocity_mesh, i, wall_g_val, get_jacobian);

            if (density_i == T{0})
            {
                // Avoid division by zero, return 0 or throw
                return T{0};
            }

            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            T integral_neg = T{0};
            T integral_pos = T{0};

            // 1. Negative Velocity Integral (Incoming)
            for (size_t k = 0; k <= zero_idx; ++k)
            {
                T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                integral_neg += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
            }
            integral_neg *= (1.0 / 3.0);

            // 2. Positive Velocity Integral (Outgoing)
            for (size_t k = zero_idx; k <= max_idx; ++k)
            {
                T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                // Apply the fix: use wall_g_val ONLY at i=0 and k=zero_idx
                T g_val = (i == 0 && k == zero_idx) ? wall_g_val : g(k, i);

                integral_pos += weight * (g_val * velocity_mesh[k]) * get_jacobian(k);
            }
            integral_pos *= (1.0 / 3.0);

            return (integral_neg + integral_pos) / density_i;
        }

        template <typename T, typename JacobianFunc>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     const T density_i, // Pre-computed
                                     T wall_g_val,      // For velocity calc
                                     JacobianFunc get_jacobian)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the mean gas velocity"));
            }

            if (density_i == T{0})
            {
                // Avoid division by zero, return 0 or throw
                return T{0};
            }

            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            T integral_neg = T{0};
            T integral_pos = T{0};

            // 1. Negative Velocity Integral (Incoming)
            for (size_t k = 0; k <= zero_idx; ++k)
            {
                T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                integral_neg += weight * (g(k, i) * velocity_mesh[k]) * get_jacobian(k);
            }
            integral_neg *= (1.0 / 3.0);

            // 2. Positive Velocity Integral (Outgoing)
            for (size_t k = zero_idx; k <= max_idx; ++k)
            {
                T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                // Apply the fix: use wall_g_val ONLY at i=0 and k=zero_idx
                T g_val = (i == 0 && k == zero_idx) ? wall_g_val : g(k, i);

                integral_pos += weight * (g_val * velocity_mesh[k]) * get_jacobian(k);
            }
            integral_pos *= (1.0 / 3.0);

            return (integral_neg + integral_pos) / density_i;
        }

        // ------ GAS TEMPERATURE -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &mean_velocities,
                                                             T wall_g_val, // Value of g at wall for zeta=0+
                                                             T wall_h_val, // Value of h at wall for zeta=0+
                                                             JacobianFunc get_jacobian)
        {
            Eigen::Vector<T, Eigen::Dynamic> temperatures(g.cols());
            if (g.cols() == 0)
            {
                return temperatures;
            }
            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            // --- Case i = 0 (Wall) ---
            {
                constexpr Eigen::Index i = 0;
                T integral_neg = T{0};
                T integral_pos = T{0};
                T v_bar = mean_velocities(i);

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                    integral_neg += weight * integrand * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing) with discontinuity fix
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                    // Apply the fix: use wall values at k=zero_idx
                    T g_val = (k == zero_idx) ? wall_g_val : g(k, i);
                    T h_val = (k == zero_idx) ? wall_h_val : h(k, i);

                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g_val) + h_val;
                    integral_pos += weight * integrand * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                temperatures.coeffRef(i) = (T{2} / T{3}) * (integral_neg + integral_pos) / densities(i);
            }

            // --- Cases i > 0 (Away from Wall) ---
            for (Eigen::Index i = 1; i < g.cols(); ++i)
            {
                T integral_neg = T{0};
                T integral_pos = T{0};
                T v_bar = mean_velocities(i);

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                    integral_neg += weight * integrand * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing)
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);
                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                    integral_pos += weight * integrand * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                temperatures.coeffRef(i) = (T{2} / T{3}) * (integral_neg + integral_pos) / densities(i);
            }

            return temperatures;
        }

        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             T wall_g_val, // For velocity & temp
                                                             T wall_h_val, // For temp
                                                             JacobianFunc get_jacobian)
        {
            Eigen::Vector<T, Eigen::Dynamic> temperatures(g.cols());
            if (g.cols() == 0)
            {
                return temperatures;
            }

            // --- 1. COMPUTE DENSITIES (using vector function) ---
            Eigen::Vector<T, Eigen::Dynamic> densities =
                compute_density(g, velocity_mesh, wall_g_val, get_jacobian);

            // --- 2. COMPUTE MEAN VELOCITIES (using vector function) ---
            Eigen::Vector<T, Eigen::Dynamic> mean_velocities =
                compute_meanGasVelocity(g, velocity_mesh, densities, wall_g_val, get_jacobian);

            // --- 3. COMPUTE TEMPERATURES (full loop, using new logic) ---
            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            // --- Case i = 0 (Wall) ---
            {
                constexpr Eigen::Index i = 0;
                T integral_neg = T{0};
                T integral_pos = T{0};
                T v_bar = mean_velocities(i);

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                    integral_neg += weight * integrand * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing) with discontinuity fix
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                    // Apply the fix: use wall values at k=zero_idx
                    T g_val = (k == zero_idx) ? wall_g_val : g(k, i);
                    T h_val = (k == zero_idx) ? wall_h_val : h(k, i);

                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g_val) + h_val;
                    integral_pos += weight * integrand * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                temperatures.coeffRef(i) = (T{2} / T{3}) * (integral_neg + integral_pos) / densities(i);
            }

            // --- Cases i > 0 (Away from Wall) ---
            for (Eigen::Index i = 1; i < g.cols(); ++i)
            {
                T integral_neg = T{0};
                T integral_pos = T{0};
                T v_bar = mean_velocities(i);

                // 1. Negative Velocity Integral (Incoming)
                for (size_t k = 0; k <= zero_idx; ++k)
                {
                    T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                    integral_neg += weight * integrand * get_jacobian(k);
                }
                integral_neg *= (1.0 / 3.0);

                // 2. Positive Velocity Integral (Outgoing)
                for (size_t k = zero_idx; k <= max_idx; ++k)
                {
                    T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);
                    T v_diff = velocity_mesh[k] - v_bar;
                    T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                    integral_pos += weight * integrand * get_jacobian(k);
                }
                integral_pos *= (1.0 / 3.0);
                temperatures.coeffRef(i) = (T{2} / T{3}) * (integral_neg + integral_pos) / densities(i);
            }

            return temperatures;
        }

        //----------------------------------------------------------------------------

        template <typename T, typename JacobianFunc>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 T wall_g_val, // For internal velocity & temp
                                 T wall_h_val, // For temp
                                 JacobianFunc get_jacobian)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the temperature"));
            }

            // Call "new" helper functions
            T density_i = compute_density_at(g, velocity_mesh, i, wall_g_val, get_jacobian);
            T mean_velocity_i = compute_meanGasVelocity_at(g, velocity_mesh, i, density_i, wall_g_val, get_jacobian);

            if (density_i == T{0})
            {
                return T{0}; // Avoid division by zero
            }

            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            T integral_neg = T{0};
            T integral_pos = T{0};
            T v_bar = mean_velocity_i;

            // 1. Negative Velocity Integral (Incoming)
            for (size_t k = 0; k <= zero_idx; ++k)
            {
                T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                T v_diff = velocity_mesh[k] - v_bar;
                T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                integral_neg += weight * integrand * get_jacobian(k);
            }
            integral_neg *= (1.0 / 3.0);

            // 2. Positive Velocity Integral (Outgoing)
            for (size_t k = zero_idx; k <= max_idx; ++k)
            {
                T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                // Apply the fix: use wall values ONLY at i=0 and k=zero_idx
                T g_val = (i == 0 && k == zero_idx) ? wall_g_val : g(k, i);
                T h_val = (i == 0 && k == zero_idx) ? wall_h_val : h(k, i);

                T v_diff = velocity_mesh[k] - v_bar;
                T integrand = (v_diff * v_diff * g_val) + h_val;
                integral_pos += weight * integrand * get_jacobian(k);
            }
            integral_pos *= (1.0 / 3.0);

            return (T{2} / T{3}) * (integral_neg + integral_pos) / density_i;
        }

        template <typename T, typename JacobianFunc>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 const T density_i,       // Pre-computed
                                 const T mean_velocity_i, // Pre-computed
                                 T wall_g_val,            // For temp
                                 T wall_h_val,            // For temp
                                 JacobianFunc get_jacobian)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the temperature"));
            }

            if (density_i == T{0})
            {
                return T{0}; // Avoid division by zero
            }

            const size_t zero_idx = velocity_mesh.get_N();
            const size_t max_idx = velocity_mesh.size() - 1;

            T integral_neg = T{0};
            T integral_pos = T{0};
            T v_bar = mean_velocity_i; // Use the provided value

            // 1. Negative Velocity Integral (Incoming)
            for (size_t k = 0; k <= zero_idx; ++k)
            {
                T weight = (k == 0 || k == zero_idx) ? 1.0 : ((zero_idx - k) % 2 != 0 ? 4.0 : 2.0);
                T v_diff = velocity_mesh[k] - v_bar;
                T integrand = (v_diff * v_diff * g(k, i)) + h(k, i);
                integral_neg += weight * integrand * get_jacobian(k);
            }
            integral_neg *= (1.0 / 3.0);

            // 2. Positive Velocity Integral (Outgoing)
            for (size_t k = zero_idx; k <= max_idx; ++k)
            {
                T weight = (k == zero_idx || k == max_idx) ? 1.0 : ((k - zero_idx) % 2 != 0 ? 4.0 : 2.0);

                // Apply the fix: use wall values ONLY at i=0 and k=zero_idx
                T g_val = (i == 0 && k == zero_idx) ? wall_g_val : g(k, i);
                T h_val = (i == 0 && k == zero_idx) ? wall_h_val : h(k, i);

                T v_diff = velocity_mesh[k] - v_bar;
                T integrand = (v_diff * v_diff * g_val) + h_val;
                integral_pos += weight * integrand * get_jacobian(k);
            }
            integral_pos *= (1.0 / 3.0);

            return (T{2} / T{3}) * (integral_neg + integral_pos) / density_i;
        }
    }

};

#endif /* PHYS_UTILS_FB91EE41_406D_4600_A84B_EB367B78C2CB */
