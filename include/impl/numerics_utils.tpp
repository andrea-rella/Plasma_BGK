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

#ifndef NUMERICS_UTILS_AB72451F_1E13_4923_A66C_EA9C225414E5
#define NUMERICS_UTILS_AB72451F_1E13_4923_A66C_EA9C225414E5

#include "../numerics_utils.hpp"

namespace Bgk
{
    namespace numerics
    {
        template <typename T>
        std::vector<std::pair<T, T>> QUICKcoefficients_p(const SpaceMeshFV<T> &mesh)
        {
            if (!mesh.validate_mesh())
                throw std::invalid_argument(error_message(
                    "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh."));

            std::vector<std::pair<T, T>> coefficients;
            const std::vector<T> &vol_boundaries = mesh.get_volume_boundaries();
            const std::vector<T> &x_comp = mesh.get_XComp();

            size_t V = vol_boundaries.size();
            T a1, a2;

            coefficients.reserve(V);

            // Case i = 0
            coefficients.emplace_back(T{0}, T{0});

            // Case i = 1: Account for the symmetric ghost node x^{(-1)} (x_comp[-1])
            a1 = ((vol_boundaries[1] - x_comp[0]) * (vol_boundaries[1] + x_comp[1])) /
                 ((x_comp[1] - x_comp[0]) * (x_comp[1] + x_comp[1]));

            a2 = ((vol_boundaries[1] - x_comp[0]) * (x_comp[1] - vol_boundaries[1])) /
                 ((x_comp[0] + x_comp[1]) * (x_comp[1] + x_comp[1]));

            coefficients.emplace_back(a1, a2);

            for (size_t i = 2; i < V - 1; ++i)
            {
                a1 = ((vol_boundaries[i] - x_comp[i - 1]) * (vol_boundaries[i] - x_comp[i - 2])) /
                     ((x_comp[i] - x_comp[i - 1]) * (x_comp[i] - x_comp[i - 2]));

                a2 = ((vol_boundaries[i] - x_comp[i - 1]) * (x_comp[i] - vol_boundaries[i])) /
                     ((x_comp[i - 1] - x_comp[i - 2]) * (x_comp[i] - x_comp[i - 2]));

                coefficients.emplace_back(a1, a2);
            }

            // Case i = V - 1
            coefficients.emplace_back(T{0}, T{0});

            return coefficients;
        }

        //--------------------------------------------------------------------------
        template <typename T>
        std::pair<T, T> QUICKcoefficients_p_at(const SpaceMeshFV<T> &mesh, const size_t i)
        {
            if (!mesh.validate_mesh())
                throw std::invalid_argument(error_message(
                    "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh."));

            if (i > mesh.get_N() + 2)
                throw std::out_of_range(error_message(
                    "Invalid index: trying to compute QUICK coefficients for an out-of-bounds volume boundary."));

            const std::vector<T> &vol_boundaries = mesh.get_volume_boundaries();
            const std::vector<T> &x_comp = mesh.get_XComp();

            if (i == 0 || i == vol_boundaries.size() - 1)
                return {T{0}, T{0}};

            if (i == 1)
            {
                T a1 = ((vol_boundaries[1] - x_comp[0]) * (vol_boundaries[1] + x_comp[1])) /
                       ((x_comp[1] - x_comp[0]) * (x_comp[1] + x_comp[1]));

                T a2 = ((vol_boundaries[1] - x_comp[0]) * (x_comp[1] - vol_boundaries[1])) /
                       ((x_comp[0] + x_comp[1]) * (x_comp[1] + x_comp[1]));
                return {a1, a2};
            }

            T a1 = ((vol_boundaries[i] - x_comp[i - 1]) * (vol_boundaries[i] - x_comp[i - 2])) /
                   ((x_comp[i] - x_comp[i - 1]) * (x_comp[i] - x_comp[i - 2]));

            T a2 = ((vol_boundaries[i] - x_comp[i - 1]) * (x_comp[i] - vol_boundaries[i])) /
                   ((x_comp[i - 1] - x_comp[i - 2]) * (x_comp[i] - x_comp[i - 2]));

            return {a1, a2};
        }

        //--------------------------------------------------------------------------

        template <typename T>
        std::vector<std::pair<T, T>> QUICKcoefficients_n(const SpaceMeshFV<T> &mesh)
        {

            if (!mesh.validate_mesh())
                throw std::invalid_argument(
                    "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh.");

            std::vector<std::pair<T, T>> coefficients;
            const std::vector<T> &vol_boundaries = mesh.get_volume_boundaries();
            const std::vector<T> &x_comp = mesh.get_XComp();

            size_t V = vol_boundaries.size();
            T b1, b2;

            coefficients.reserve(V);

            // Case i = 0:
            coefficients.emplace_back(T{0}, T{0});

            // General case:
            for (size_t i = 1; i < V - 2; ++i)
            {
                b1 = ((vol_boundaries[i] - x_comp[i]) * (vol_boundaries[i] - x_comp[i + 1])) /
                     ((x_comp[i - 1] - x_comp[i]) * (x_comp[i - 1] - x_comp[i + 1]));

                b2 = ((vol_boundaries[i] - x_comp[i]) * (x_comp[i - 1] - vol_boundaries[i])) /
                     ((x_comp[i] - x_comp[i + 1]) * (x_comp[i - 1] - x_comp[i + 1]));

                coefficients.emplace_back(b1, b2);
            }

            // Case i = V - 2: Account for the symmetric ghost node x^{(N+1)} (x_comp[V - 1])
            b1 = ((vol_boundaries[V - 2] - x_comp[V - 2]) * (vol_boundaries[V - 2] - (2. * x_comp[V - 2] - x_comp[V - 3]))) /
                 ((x_comp[V - 3] - x_comp[V - 2]) * (x_comp[V - 3] - (2. * x_comp[V - 2] - x_comp[V - 3])));

            b2 = ((vol_boundaries[V - 2] - x_comp[V - 2]) * (x_comp[V - 3] - vol_boundaries[V - 2])) /
                 ((x_comp[V - 2] - (2. * x_comp[V - 2] - x_comp[V - 3])) * (x_comp[V - 3] - (2. * x_comp[V - 2] - x_comp[V - 3])));

            coefficients.emplace_back(b1, b2);

            // Case i = V - 1:
            coefficients.emplace_back(T{0}, T{0});

            return coefficients;
        }

        //--------------------------------------------------------------------------

        template <typename T>
        std::pair<T, T> QUICKcoefficients_n_at(const SpaceMeshFV<T> &mesh, const size_t i)
        {
            if (!mesh.validate_mesh())
                throw std::invalid_argument(
                    "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh.");

            if (i > mesh.get_N() + 2)
                throw std::out_of_range(error_message(
                    "Invalid index: trying to compute QUICK coefficients for an out-of-bounds volume boundary."));

            const std::vector<T> &vol_boundaries = mesh.get_volume_boundaries();
            const size_t V = vol_boundaries.size();
            const std::vector<T> &x_comp = mesh.get_XComp();

            if (i == 0 || i == vol_boundaries.size() - 1)
                return {T{0}, T{0}};

            if (i == vol_boundaries.size() - 2)
            {
                T b1 = ((vol_boundaries[V - 2] - x_comp[V - 2]) * (vol_boundaries[V - 2] - (2. * x_comp[V - 2] - x_comp[V - 3]))) /
                       ((x_comp[V - 3] - x_comp[V - 2]) * (x_comp[V - 3] - (2. * x_comp[V - 2] - x_comp[V - 3])));

                T b2 = ((vol_boundaries[V - 2] - x_comp[V - 2]) * (x_comp[V - 3] - vol_boundaries[V - 2])) /
                       ((x_comp[V - 2] - (2. * x_comp[V - 2] - x_comp[V - 3])) * (x_comp[V - 3] - (2. * x_comp[V - 2] - x_comp[V - 3])));

                return {b1, b2};
            }

            T b1 = ((vol_boundaries[i] - x_comp[i]) * (vol_boundaries[i] - x_comp[i + 1])) /
                   ((x_comp[i - 1] - x_comp[i]) * (x_comp[i - 1] - x_comp[i + 1]));

            T b2 = ((vol_boundaries[i] - x_comp[i]) * (x_comp[i - 1] - vol_boundaries[i])) /
                   ((x_comp[i] - x_comp[i + 1]) * (x_comp[i - 1] - x_comp[i + 1]));

            return {b1, b2};
        }

        template <typename T>
        T CDScoefficients_at(const SpaceMeshFV<T> &mesh, const size_t i)
        {
            if (!mesh.validate_mesh())
                throw std::invalid_argument(
                    "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh.");

            if (i == 0 || i > mesh.get_N())
                throw std::out_of_range(error_message(
                    "Invalid index: trying to compute CDS coefficients for either a boundary CV volume face or \
                    an out-of-bounds volume face."));

            const std::vector<T> &vol_boundaries = mesh.get_volume_boundaries();
            const std::vector<T> &x_comp = mesh.get_XComp();

            return (vol_boundaries[i] - x_comp[i - 1]) / (x_comp[i] - x_comp[i - 1]);
        }

    }
}

#endif /* NUMERICS_UTILS_AB72451F_1E13_4923_A66C_EA9C225414E5 */
