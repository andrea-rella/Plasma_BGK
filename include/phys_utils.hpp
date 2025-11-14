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

#ifndef PHYS_UTILS_EDCE3BAA_66FA_48CD_B277_A15BDAD8AFB0
#define PHYS_UTILS_EDCE3BAA_66FA_48CD_B277_A15BDAD8AFB0

#include <Eigen/Sparse>
#include <VelocityMesh.hpp>

namespace Bgk
{
    namespace phys
    {
        // ------ GAS DENSITY----------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /**
         * @brief Computes the normalised density for a specific computational point in space
         *
         * Computes the normalised density @f$ \overline{\rho} = \frac{\rho}{\rho_w} = \int_{-\infty}^\infty g \,d\zeta @f$
         * for a specific computational point in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param i Index of the space computational point
         * @return The value of the normalised density at the specified point
         *
         * @throw std::out_of_range if i exceeds the limits
         *
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                         const VelocityMesh<T> &velocity_mesh,
                                                         T wall_g_val,
                                                         JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic gas density at a single spatial point.
         *
         * This function implements the density calculation (Eq. 11) by performing
         * two separate integrations using Simpson's 1/3 rule: one for the
         * incoming molecules (zeta < 0) and one for the outgoing molecules (zeta > 0).
         *
         * It correctly handles the discontinuity at the wall (i=0) for the
         * outgoing part of the distribution.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param i The spatial index at which to compute the density.
         * @param wall_maxwellian_val The value of the Maxwellian at the wall (for i=0, k=zero_idx).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return T The computed density at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_density_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                             const VelocityMesh<T> &velocity_mesh,
                             const Eigen::Index i,
                             T wall_g_val,
                             JacobianFunc get_jacobian);
        // ------ MEAN GAS VELOCITY ---------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /**
         * @brief Computes the macroscopic mean gas velocity at each spatial point.
         *
         * Implements Eq. 12 by separately integrating for zeta < 0 and zeta > 0
         * using Simpson's 1/3 rule. Correctly handles the wall discontinuity.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param densities A pre-computed vector of densities at each spatial point.
         * @param wall_g_val Value of g at wall for zeta=0+ (for discontinuity fix).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return Eigen::Vector<T, Eigen::Dynamic> A vector of mean gas velocities.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                                 T wall_g_val, // Value of g at wall for zeta=0+
                                                                 JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic mean gas velocity at each spatial point.
         *
         * This version first calls the vector-based `compute_density` function,
         * then uses its results to compute the mean velocities. All calculations
         * use the "new" integration philosophy.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param wall_maxwellian_val Value for internal density comp. at wall.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity comp.).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return Eigen::Vector<T, Eigen::Dynamic> A vector of mean gas velocities.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 T wall_g_val, // For velocity
                                                                 JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic mean gas velocity at a single spatial point.
         *
         * This version computes the density internally using the "new" philosophy.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param i The spatial index at which to compute the velocity.
         * @param wall_maxwellian_val Value of Maxwellian for density comp. at wall.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity comp.).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return T The computed mean gas velocity at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     T wall_g_val, // For velocity calc
                                     JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic mean gas velocity at a single spatial point.
         *
         * This overload accepts a pre-computed density for efficiency.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param i The spatial index at which to compute the velocity.
         * @param density_i The pre-computed density at point i.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity comp.).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return T The computed mean gas velocity at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     const T density_i, // Pre-computed
                                     T wall_g_val,      // For velocity calc
                                     JacobianFunc get_jacobian);

        // ------ GAS TEMPERATURE -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------
        /**
         * @brief Computes the macroscopic gas temperature (performant overload).
         *
         * Implements Eq. 13 by separately integrating for zeta < 0 and zeta > 0
         * using Simpson's 1/3 rule. Correctly handles the wall discontinuity.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param h The second reduced velocity distribution function, h(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param densities Pre-computed densities (must be correct).
         * @param mean_velocities Pre-computed velocities (must be correct).
         * @param wall_g_val Value of g at wall for zeta=0+ (for discontinuity fix).
         * @param wall_h_val Value of h at wall for zeta=0+ (for discontinuity fix).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return Eigen::Vector<T, Eigen::Dynamic> A vector of temperatures.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &mean_velocities,
                                                             T wall_g_val, // Value of g at wall for zeta=0+
                                                             T wall_h_val, // Value of h at wall for zeta=0+
                                                             JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic gas temperature at each spatial point.
         *
         * This version first calls the vector-based `compute_density` and
         * `compute_meanGasVelocity` functions, then uses their results to
         * compute the temperatures. All calculations use the "new"
         * integration philosophy.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param h The second reduced velocity distribution function, h(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param wall_maxwellian_val Value for internal density comp. at wall.
         * @param wall_g_val Value of g at wall for zeta=0+ (for internal velocity/temp comp.).
         * @param wall_h_val Value of h at wall for zeta=0+ (for internal temp comp.).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return Eigen::Vector<T, Eigen::Dynamic> A vector of temperatures.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             T wall_g_val, // For velocity & temp
                                                             T wall_h_val, // For temp
                                                             JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic gas temperature at a single spatial point.
         *
         * This version computes the density and mean velocity internally
         * using the "new" philosophy.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param h The second reduced velocity distribution function, h(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param i The spatial index at which to compute the temperature.
         * @param wall_maxwellian_val Value for density comp. at wall.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity/temp comp.).
         * @param wall_h_val Value of h at wall for zeta=0+ (for temp comp.).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return T The computed temperature at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 T wall_g_val, // For internal velocity & temp
                                 T wall_h_val, // For temp
                                 JacobianFunc get_jacobian);

        /**
         * @brief Computes the macroscopic gas temperature at a single spatial point.
         *
         * This overload accepts pre-computed density and mean velocity for efficiency.
         *
         * @param g The reduced velocity distribution function, g(j, i).
         * @param h The second reduced velocity distribution function, h(j, i).
         * @param velocity_mesh The non-uniform velocity mesh object.
         * @param i The spatial index at which to compute the temperature.
         * @param density_i The pre-computed density at point i.
         * @param mean_velocity_i The pre-computed mean velocity at point i.
         * @param wall_g_val Value of g at wall for zeta=0+ (for temp comp.).
         * @param wall_h_val Value of h at wall for zeta=0+ (for temp comp.).
         * @param get_jacobian A function to get the Jacobian for the integration.
         * @return T The computed temperature at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 const T density_i,       // Pre-computed
                                 const T mean_velocity_i, // Pre-computed
                                 T wall_g_val,            // For temp
                                 T wall_h_val,            // For temp
                                 JacobianFunc get_jacobian);
    }
}

#include "impl/phys_utils.tpp"

#endif /* PHYS_UTILS_EDCE3BAA_66FA_48CD_B277_A15BDAD8AFB0 */
