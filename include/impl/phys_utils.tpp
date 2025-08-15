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

        // ------ GAS DENSITY----------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                         const VelocityMesh<T> &velocity_mesh)
        {
            Eigen::Vector<T, Eigen::Dynamic> densities(g.cols());

            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                T integral = T{0};
                for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
                {
                    integral += (g(j, i) + g(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
                }
                densities.coeffRef(i) = integral;
            }

            return densities;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        T compute_density_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                             const VelocityMesh<T> &velocity_mesh,
                             const Eigen::Index i)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the density"));
            }

            T density = T{0};

            for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
            {
                density += (g(j, i) + g(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
            }
            return density;
        }

        // ------ MEAN GAS VELOCITY ---------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh)
        {
            Eigen::Vector<T, Eigen::Dynamic> velocities(g.cols());
            Eigen::Vector<T, Eigen::Dynamic> densities = compute_density_at(g, velocity_mesh);

            auto dv = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return g(j, i) * velocity_mesh[j];
            };

            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                T integral = T{0};
                for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
                {
                    integral += (dv(j, i) + dv(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
                }
                velocities.coeffRef(i) = integral / densities(i);
            }
            return velocities;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 const Eigen::Vector<T, Eigen::Dynamic> &densities)
        {
            Eigen::Vector<T, Eigen::Dynamic> velocities(g.cols());

            auto dv = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return g(j, i) * velocity_mesh[j];
            };

            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                T integral = T{0};
                for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
                {
                    integral += (dv(j, i) + dv(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
                }
                velocities.coeffRef(i) = integral / densities(i);
            }
            return velocities;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the mean gas velocity"));
            }

            T density_i = compute_density_at(g, velocity_mesh, i);

            auto dv = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return g(j, i) * velocity_mesh[j];
            };

            T integral = T{0};
            for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
            {
                integral += (dv(j, i) + dv(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
            }
            return integral / density_i;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     const T density_i)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the mean gas velocity"));
            }

            auto dv = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return g(j, i) * velocity_mesh[j];
            };

            T integral = T{0};
            for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
            {
                integral += (dv(j, i) + dv(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
            }
            return integral / density_i;
        }

        // ------ GAS TEMPERATURE -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh)
        {
            Eigen::Vector<T, Eigen::Dynamic> temperatures(g.cols());
            Eigen::Vector<T, Eigen::Dynamic> densities = compute_density_at(g, velocity_mesh);
            Eigen::Vector<T, Eigen::Dynamic> mean_velocities = compute_meanGasVelocity_at(g, velocity_mesh, densities);

            auto dT = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return (velocity_mesh[j] - mean_velocities(i)) * (velocity_mesh[j] - mean_velocities(i)) * g(j, i) + h(j, i);
            };

            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                T integral = T{0};
                for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
                {
                    integral += (dT(j, i) + dT(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
                }
                temperatures.coeffRef(i) = T{2} / T{3} * integral / densities(i);
            }
            return temperatures;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &mean_velocities)
        {
            Eigen::Vector<T, Eigen::Dynamic> temperatures(g.cols());

            auto dT = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return (velocity_mesh[j] - mean_velocities(i)) * (velocity_mesh[j] - mean_velocities(i)) * g(j, i) + h(j, i);
            };

            for (Eigen::Index i = 0; i < g.cols(); ++i)
            {
                T integral = T{0};
                for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
                {
                    integral += (dT(j, i) + dT(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
                }
                temperatures.coeffRef(i) = T{2} / T{3} * integral / densities(i);
            }
            return temperatures;
        }

        //----------------------------------------------------------------------------

        template <typename T>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i)
        {

            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the mean gas velocity"));
            }

            T density_i = compute_density_at(g, velocity_mesh, i);
            T mean_velocity_i = compute_meanGasVelocity_at(g, velocity_mesh, i, density_i);

            auto dT = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return (velocity_mesh[j] - mean_velocity_i) * (velocity_mesh[j] - mean_velocity_i) * g(j, i) + h(j, i);
            };

            T integral = T{0};

            for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
            {
                integral += (dT(j, i) + dT(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
            }
            return (T{2} / T{3}) * (integral / density_i);
        }

        //----------------------------------------------------------------------------

        template <typename T>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 const T density_i,
                                 const T mean_velocity_i)
        {
            if (i < 0 || i >= g.cols())
            {
                throw std::out_of_range(error_message("Index out of bounds when computing the mean gas velocity"));
            }

            auto dT = [&](Eigen::Index j, Eigen::Index i) -> T
            {
                return (velocity_mesh[j] - mean_velocity_i) * (velocity_mesh[j] - mean_velocity_i) * g(j, i) + h(j, i);
            };

            T integral = T{0};

            for (size_t j = 0; j < velocity_mesh.size() - 1; ++j)
            {
                integral += (dT(j, i) + dT(j + 1, i)) * (velocity_mesh[j + 1] - velocity_mesh[j]) * T{0.5};
            }
            return (T{2} / T{3}) * (integral / density_i);
        }

    }
}

#endif /* PHYS_UTILS_BD38F0A0_5655_4640_BAA9_4F7E62E0ADA7 */
