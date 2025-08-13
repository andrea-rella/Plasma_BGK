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
        /**
         * @brief Computes the normalised density for all the computational points in space
         *
         * Computes the normalised density @f$ \overline{\rho} = \frac{\rho}{\rho_w} = \int_\infty^\infty g \,d\zeta @f$
         * for all the computational points in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocities, cols -> spatial points)
         * @param velocities Bgk::VelocityMesh
         * @return Eigen::Vector<T, Eigen::Dynamic>
         */
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                         const VelocityMesh<T> &velocities);

        //----------------------------------------------------------------------------

        /**
         * @brief Computes the normalised density for a specific computational point in space
         *
         * Computes the normalised density @f$ \overline{\rho} = \frac{\rho}{\rho_w} = \int_\infty^\infty g \,d\zeta @f$
         * for a specific computational point in space using the trapezoidal rule.
         *
         * @tparam T floating type precision
         * @param g Eigen matrix containing the g quantity (rows -> velocities, cols -> spatial points)
         * @param velocities Bgk::VelocityMesh
         * @param i Index of the space computational point
         * @return T
         *
         * @throw std::out_of_range if i exceeds the limits
         *
         * @overload
         */
        template <typename T>
        T compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g, const VelocityMesh<T> &velocities,
                          const Eigen::Index i);
    }
}

#include "impl/phys_utils.tpp"

#endif /* PHYS_UTILS_EDCE3BAA_66FA_48CD_B277_A15BDAD8AFB0 */
