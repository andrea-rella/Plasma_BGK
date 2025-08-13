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

#ifndef PHYS_UTILS_BD38F0A0_5655_4640_BAA9_4F7E62E0ADA7
#define PHYS_UTILS_BD38F0A0_5655_4640_BAA9_4F7E62E0ADA7

#include "../phys_utils.hpp"
#include "../utilities.hpp"

namespace Bgk
{
    namespace phys
    {
        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                         const VelocityMesh<T> &velocities)
        {
            Eigen::Vector<T, Eigen::Dynamic> densities(g.cols());

            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                T integral = T{0};
                for (size_t j = 0; j < velocities.size() - 1; ++j)
                {
                    integral += (g(j, i) + g(j + 1, i)) * (velocities[j + 1] - velocities[j]) * T{0.5};
                }
                densities.coeffRef(i) = integral;
            }

            return densities;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        T compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g, const VelocityMesh<T> &velocities,
                          const Eigen::Index i)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the density"));
            }

            T density = T{0};

            for (size_t j = 0; j < velocities.size() - 1; ++j)
            {
                density += (g(j, i) + g(j + 1, i)) * (velocities[j + 1] - velocities[j]) * T{0.5};
            }
            return density;
        }
    }
}

#endif /* PHYS_UTILS_BD38F0A0_5655_4640_BAA9_4F7E62E0ADA7 */
