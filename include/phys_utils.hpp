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
         * @brief Computes the normalised density for all the computational points in space
         *
         * Computes the normalised density @f$ \overline{\rho} = \frac{\rho}{\rho_w} = \int_{-\infty}^\infty g \,d\zeta @f$
         * for all the computational points in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @return Eigen::Vector<T, Eigen::Dynamic> containing the normalized densities in each point in space
         */
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                         const VelocityMesh<T> &velocity_mesh);

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
        template <typename T>
        T compute_density_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                             const VelocityMesh<T> &velocity_mesh,
                             const Eigen::Index i);

        // ------ MEAN GAS VELOCITY ---------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /**
         * @brief Computes the mean gas velocity for all the computational points in space
         *
         * Computes the normalized mean gas velocity @f$ \overline{u} = \frac{1}{\overline{\rho}} \int_{-\infty}^\infty \zeta g \,d\zeta @f$
         * for all the computational points in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @return Eigen::Vector<T, Eigen::Dynamic> containing the mean gas velocities in each point in space
         */
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh);

        /**
         * @overload compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &,
         *                                   const VelocityMesh<T> &);
         *
         * Computes the normalized mean gas velocity @f$ \overline{u} = \frac{1}{\overline{\rho}} \int_{-\infty}^\infty \zeta g \,d\zeta @f$
         * for all the computational points in space using the trapezoidal rule. It accepts the precomputed density values as an input.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param densities Eigen vector containing the densities for each spatial point
         * @return Eigen::Vector<T, Eigen::Dynamic>
         */
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 const Eigen::Vector<T, Eigen::Dynamic> &densities);

        /**
         * @brief Computes the mean gas velocity at a specific spatial point
         *
         * Computes the mean gas velocity @f$ \overline{u} = \frac{1}{\overline{\rho}} \int_{-\infty}^\infty \zeta g \,d\zeta @f$
         * at a specific spatial point using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param i Index of the spatial point
         * @return The mean gas velocity at the specified spatial point
         *
         * @throw std::out_of_range if i exceeds the limits
         *
         */
        template <typename T>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i);

        /**
         * @overload compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &,
         *                                     const VelocityMesh<T> &,
         *                                     const Eigen::Index);
         *
         * Computes the mean gas velocity @f$ \overline{u} = \frac{1}{\overline{\rho}} \int_{-\infty}^\infty \zeta g \,d\zeta @f$
         * at a specific spatial point using the trapezoidal rule. Accepts the precomputed density value as an input.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param i Index of the spatial point
         * @param density_i Density at the spatial point
         * @return The mean gas velocity at the specified spatial point
         *
         * @throw std::out_of_range if i exceeds the limits
         *
         */
        template <typename T>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     const T density_i);

        // ------ GAS TEMPERATURE -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /**
         * @brief Computes the normalized gas temperature for all points in space
         *
         * Computes the normalized gas temperature @f$ \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left(
         * \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)@f$ for
         * all the computational points in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @return Eigen::Vector<T, Eigen::Dynamic> containing the normalized temperatures in each point in space
         */
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh);

        /**
         * @overload compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &,
         *                               const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &,
         *                               const VelocityMesh<T> &)
         *
         * Computes the normalized gas temperature @f$ \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left(
         * \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)@f$ for
         * all the computational points in space using the trapezoidal rule. Accepts the precomputed normalized density
         * and normalized mean velocity values as input.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param densities Eigen vector containing the precomputed normalized densities
         * @param mean_velocities Eigen vector containing the precomputed normalized mean velocities
         * @return Eigen::Vector<T, Eigen::Dynamic> containing the normalized temperatures in each point in space
         *
         */
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &mean_velocities);

        /**
         * @brief Computes the normalized gas temperature at a specific point in space.
         *
         * Computes the normalized gas temperature at a specific point @f$ \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left(
         * \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)@f$ for
         * all the computational points in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param i Index of the spatial point
         * @return T containing the normalized temperature at the specified point in space
         *
         * @throw std::out_of_range if i exceeds the limits
         */
        template <typename T>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i);

        /**
         * @overload compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &,
         *                       const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &,
         *                       const VelocityMesh<T> &,
         *                       const Eigen::Index);
         *
         * Computes the normalized gas temperature at a specific point @f$ \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left(
         * \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)@f$ for
         * all the computational points in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points)
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols -> spatial points)
         * @param velocity_mesh Bgk::VelocityMesh
         * @param i Index of the spatial point
         * @param density_i Density at the spatial point
         * @param mean_velocity_i Mean velocity at the spatial point
         * @return T containing the normalized temperature at the specified point in space
         *
         * @throw std::out_of_range if i exceeds the limits
         */
        template <typename T>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 const T density_i,
                                 const T mean_velocity_i);
    }
}

#include "impl/phys_utils.tpp"

#endif /* PHYS_UTILS_EDCE3BAA_66FA_48CD_B277_A15BDAD8AFB0 */
