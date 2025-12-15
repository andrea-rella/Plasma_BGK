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

#ifndef NUMERICS_UTILS_ADAF7230_7005_496D_B00B_61A2C1DAB1D7
#define NUMERICS_UTILS_ADAF7230_7005_496D_B00B_61A2C1DAB1D7

#include "SpaceMeshFV.hpp"

namespace Bgk
{

    // Forward declaration of SpaceMeshFV to avoid circular dependencies.
    template <typename T>
    class SpaceMeshFV;

    namespace numerics
    {
        /** @brief Computes QUICK coefficients for a given SpaceMeshFV and a POSITIVE velocity.
         *
         * Computes the QUICK (quadratic upwind interpolation) coefficients @f$ a_1, a_2 @f$ for a given SpaceMeshFV
         * mesh such that @f$ \phi_e = \phi_U + a_1 (\phi_D-\phi_U) + a_2 (\phi_U-\phi_{UU}) @f$, where:
         *
         * - @f$ a_1 = \frac{ (x_e-x_U) (x_e-x_{UU}) }{ (x_D-x_U) (x_D-x_{UU}) } @f$
         *
         * - @f$ a_2 = \frac{ (x_e-x_U) (x_D-x_e) }{ (x_U-x_{UU}) (x_D-x_{UU}) } @f$
         *
         * with U = Upstream, D = Downstream, UU = Upstream-Upstream, e = Edge.
         *
         * @tparam T precision type of the mesh components.
         * @param mesh Space mesh to compute the coefficients for.
         * @return std::vector<std::pair<T, T>> containing the coefficients for all volume boundaries.
         *
         * @throws std::runtime_error if the mesh is not initialized or constructed.
         *
         * @note The first and last elements of the vector are set to (0., 0.) as they are not used in this specific
         *       finite volume formulation of the problem. The reason is that for a positive velocity the first boundary
         *       coincide with the emission boundary condition and the last boundary is not used since the it's identified
         *       with the last computational cell.
         *
         * @relatesalso SpaceMeshFV.
         *
         * @cite ferziger2019computational
         */
        template <typename T>
        std::vector<std::pair<T, T>> QUICKcoefficients_p(const SpaceMeshFV<T> &mesh);

        /** @brief Computes QUICK coefficient for a given volume boundary of a SpaceMeshFV and a POSITIVE velocity.
         *
         * Computes the QUICK (quadratic upwind interpolation) coefficient @f$ a_1, a_2 @f$ for a given SpaceMeshFV mesh
         * volume boundary such that @f$ \phi_e = \phi_U + a_1 (\phi_D-\phi_U) + a_2 (\phi_U-\phi_{UU}) @f$, where:
         *
         * - @f$ a_1 = \frac{ (x_e-x_U) (x_e-x_{UU}) }{ (x_D-x_U) (x_D-x_{UU}) } @f$
         *
         * - @f$ a_2 = \frac{ (x_e-x_U) (x_D-x_e) }{ (x_U-x_{UU}) (x_D-x_{UU}) } @f$
         *
         * with U = Upstream, D = Downstream, UU = Upstream-Upstream, e = Edge.
         *
         * @tparam T precision type of the mesh components.
         * @param mesh Space mesh to compute the coefficients for.
         * @param i Index of the volume boundary for which to compute the coefficients.
         * @return std::pair<T, T> containing the coefficients for the specified volume boundary.
         *
         * @throws std::runtime_error if the mesh is not initialized or constructed.
         * @throws std::out_of_range if the index @p i is out of bounds of the mesh volume boundaries.
         *
         * @note The first and last elements of the vector are set to (0., 0.) as they are not used in this specific
         *       finite volume formulation of the problem. The reason is that for a positive velocity the first boundary
         *       coincide with the emission boundary condition and the last boundary is not used since the it's identified
         *       with the last computational cell.
         *
         * @relatesalso SpaceMeshFV.
         *
         * @cite ferziger2019computational
         */
        template <typename T>
        std::pair<T, T> QUICKcoefficients_p_at(const SpaceMeshFV<T> &mesh, const size_t i);

        /** @brief Computes QUICK coefficients for a given SpaceMeshFV and a NEGATIVE velocity.
         *
         * Computes the QUICK (quadratic upwind interpolation) coefficients @f$ a_1, a_2 @f$ for a given SpaceMeshFV
         * mesh such that @f$ \phi_e = \phi_U + a_1 (\phi_D-\phi_U) + a_2 (\phi_U-\phi_{UU}) @f$, where:
         *
         * - @f$ a_1 = \frac{ (x_e-x_U) (x_e-x_{UU}) }{ (x_D-x_U) (x_D-x_{UU}) } @f$
         *
         * - @f$ a_2 = \frac{ (x_e-x_U) (x_D-x_e) }{ (x_U-x_{UU}) (x_D-x_{UU}) } @f$
         *
         * with U = Upstream, D = Downstream, UU = Upstream-Upstream, e = Edge.
         *
         * @tparam T precision type of the mesh components.
         * @param mesh Space mesh to compute the coefficients for.
         * @return std::vector<std::pair<T, T>> containing the coefficients for all volume boundaries.
         *
         * @throws std::runtime_error if the mesh is not initialized or constructed.
         *
         * @note The first and last elements of the vector are set to (0., 0.) as they are not used in this specific
         *       finite volume formulation of the problem. Since if the velocity is negative the first boundary is identified
         *       with the first computational cell and the last boundary coincide with the inflow boundary condition.
         *
         * @relatesalso SpaceMeshFV
         *
         * @cite ferziger2019computational
         */
        template <typename T>
        std::vector<std::pair<T, T>> QUICKcoefficients_n(const SpaceMeshFV<T> &mesh);

        /** @brief Computes QUICK coefficients for a given SpaceMeshFV and a NEGATIVE velocity.
         *
         * Computes the QUICK (quadratic upwind interpolation) coefficients @f$ a_1, a_2 @f$ for a given SpaceMeshFV mesh
         * volume boundary such that @f$ \phi_e = \phi_U + a_1 (\phi_D-\phi_U) + a_2 (\phi_U-\phi_{UU}) @f$, where:
         *
         * - @f$ a_1 = \frac{ (x_e-x_U) (x_e-x_{UU}) }{ (x_D-x_U) (x_D-x_{UU}) } @f$
         *
         * - @f$ a_2 = \frac{ (x_e-x_U) (x_D-x_e) }{ (x_U-x_{UU}) (x_D-x_{UU}) } @f$
         *
         * with U = Upstream, D = Downstream, UU = Upstream-Upstream, e = Edge.
         *
         * @tparam T precision type of the mesh components.
         * @param mesh Space mesh to compute the coefficients for.
         * @param i Index of the volume boundary for which to compute the coefficients.
         * @return std::pair<T, T> containing the coefficients for the specified volume boundary.
         *
         * @throws std::runtime_error if the mesh is not initialized or constructed.
         * @throws std::out_of_range if the index @p i is out of bounds of the mesh volume boundaries.
         *
         * @note The first and last elements of the vector are set to (0., 0.) as they are not used in this specific
         *       finite volume formulation of the problem. Since if the velocity is negative the first boundary is identified
         *       with the first computational cell and the last boundary coincide with the inflow boundary condition.
         *
         * @relatesalso SpaceMeshFV
         *
         * @cite ferziger2019computational
         */
        template <typename T>
        std::pair<T, T> QUICKcoefficients_n_at(const SpaceMeshFV<T> &mesh, const size_t i);

    }
}

#include "impl/numerics_utils.tpp"

#endif /* NUMERICS_UTILS_ADAF7230_7005_496D_B00B_61A2C1DAB1D7 */
