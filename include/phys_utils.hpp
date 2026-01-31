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

#ifndef PHYS_UTILS_B3D31D31_524A_4009_AA77_ABC07533CF00
#define PHYS_UTILS_B3D31D31_524A_4009_AA77_ABC07533CF00

#include <Eigen/Sparse>
#include <VelocityMesh.hpp>

namespace Bgk
{
    namespace phys
    {
        // ------ GAS DENSITY----------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /** @brief Computes the macroscopic normalised gas density for all computational points in space.
         *
         * Computes the normalised density:
         * \f[
         * \overline{\rho} = \frac{\rho}{\rho_w} = \int_{-\infty}^\infty g \,d\zeta
         * \f]
         * for all computational points in space. The integral is computed numerically by splitting it into
         * negative and positive velocity domains:
         * \f[
         * \overline{\rho}_i = \int_{-\infty}^{\zeta=0} g(\zeta, x_i) \,d\zeta + \int_{\zeta=0}^{\infty} g(\zeta, x_i) \,d\zeta
         * \f]
         *
         * The numerical integration is performed over the discrete computational index $j$ using a change of variables:
         * \f[
         * \int g(\zeta) \,d\zeta \approx \sum w_j g(\zeta_j) \cdot \left| \frac{d\zeta}{dj} \right|_j \cdot \Delta j
         * \f]
         *
         * The function is implemented using **Simpson's \f$\frac{1}{3}\f$ Rule** where the computational spacing $\Delta j$ is
         * constant ($\Delta j=1$) and the sum of weights includes the \f$\frac{1}{3}\f$ factor. The term returned by
         * `get_jacobian(k)` corresponds to the **Jacobian of the transformation** \f$ \left| \frac{d\zeta}{dj} \right| \cdot \Delta j \f$.
         * This represents the derivative of the physical velocity $\zeta$ with respect to the computational index $j$ (which is equivalent
         * to the velocity index $k$ in the loop), scaled by the computational step size.
         *
         * A special treatment is applied to the **Case i = 0 (Wall Boundary)** where a potential
         * **discontinuity** in the distribution function is handled at the \f$ \zeta=0 \f$ velocity point.
         * Specifically, for the positive velocity integral (outgoing), the value \f$ g(k, i) \f$ at \f$ k = \text{zero\_idx} \f$
         * is replaced by a prescribed **wall value** (`wall_g_val`) to correctly capture the boundary
         * condition and numerical flux.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param wall_g_val Value of g at wall for zeta=0+ (for density calc)
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return Eigen::Vector<T, Eigen::Dynamic> The vector of normalised densities for all spatial points.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_density(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                                         const VelocityMesh<T> &velocity_mesh,
                                                         T wall_g_val,
                                                         JacobianFunc get_jacobian);

        /** @brief Computes the macroscopicnormalized gas density at a single spatial point.
         *
         * Computes the normalised density:
         * \f[
         * \overline{\rho}_i = \frac{\rho_i}{\rho_w} = \int_{-\infty}^\infty g(\zeta, x_i) \,d\zeta
         * \f]
         * at a specific spatial index `i`. The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * For further details on the implementation and mathematical formulation, please refer to the vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param i The spatial index at which to compute the density.
         * @param wall_g_val Value of g at wall for zeta=0+ (for density calc).
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return T The computed density at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_density_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                             const VelocityMesh<T> &velocity_mesh,
                             const Eigen::Index i,
                             T wall_g_val,
                             JacobianFunc get_jacobian);

        // ------ MEAN GAS VELOCITY ---------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /** @brief Computes the macroscopic normalized mean gas velocity at each spatial point.
         *
         * Computes the normalised mean gas velocity:
         * \f[
         * \bar{v} = \bar{\rho}^{-1} \int_{-\infty}^{\infty} \zeta g \, d\zeta
         * \f]
         * for all computational points in space. The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param densities Pre-computed densities.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity calc)
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return Eigen::Vector<T, Eigen::Dynamic> The vector of normalised mean velocities for all spatial points.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                                 T wall_g_val, // Value of g at wall for zeta=0+
                                                                 JacobianFunc get_jacobian);

        /** @brief Computes the macroscopic normalized mean gas velocity at each spatial point.
         *
         * Computes the normalised mean gas velocity:
         * \f[
         * \bar{v} = \bar{\rho}^{-1} \int_{-\infty}^{\infty} \zeta g \, d\zeta
         * \f]
         * for all computational points in space. The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * This version first calls the vector-based `compute_density` function to obtain the densities,
         * then uses these results to compute the mean velocities.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity calc)
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return Eigen::Vector<T, Eigen::Dynamic> The vector of normalised mean velocities for all spatial points.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_meanGasVelocity(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                                                 const VelocityMesh<T> &velocity_mesh,
                                                                 T wall_g_val, // For velocity
                                                                 JacobianFunc get_jacobian);

        /** @brief Computes the macroscopic normalized mean gas velocity at a single spatial point.
         *
         * Computes the normalised mean gas velocity at a specific spatial index `i`:
         * \f[
         * \bar{v}_i = \bar{\rho}_i^{-1} \int_{-\infty}^{\infty} \zeta g(\zeta, x_i) \, d\zeta
         * \f]
         *
         * The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * Accepts a pre-computed density for efficiency.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols ->
         * spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param i The spatial index at which to compute the mean gas velocity.
         * @param density_i The pre-computed density at point i.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity calc).
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return T The computed mean gas velocity at the i-th spatial point.
         *
         */
        template <typename T, typename JacobianFunc>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     const T density_i, // Pre-computed
                                     T wall_g_val,      // For velocity calc
                                     JacobianFunc get_jacobian);

        /** @brief Computes the macroscopic normalizedmean gas velocity at a single spatial point.
         *
         * Computes the normalised mean gas velocity at a specific spatial index `i`:
         * \f[
         * \bar{v}_i = \bar{\rho}_i^{-1} \int_{-\infty}^{\infty} \zeta g(\zeta, x_i) \, d\zeta
         * \f]
         *
         * The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * This version computes the density internally calling the `compute_density_at` function.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols ->
         * spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param i The spatial index at which to compute the mean gas velocity.
         * @param wall_g_val Value of g at wall for zeta=0+ (for velocity calc).
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return T The computed mean gas velocity at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_meanGasVelocity_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                     const VelocityMesh<T> &velocity_mesh,
                                     const Eigen::Index i,
                                     T wall_g_val, // For velocity calc
                                     JacobianFunc get_jacobian);

        // ------ GAS TEMPERATURE -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /** @brief Computes the macroscopic normalized mean gas temperature at each spatial point.
         *
         * Computes the normalised mean gas temperature:
         * \f[
         * \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left( \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)
         * \f]
         * for all computational points in space. The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points).
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols -> spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param densities Pre-computed densities.
         * @param mean_velocities Pre-computed mean velocities.
         * @param wall_g_val Value of g at wall for zeta=0+ (for temperature calc).
         * @param wall_h_val Value of h at wall for zeta=0+ (for temperature calc)..
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return Eigen::Vector<T, Eigen::Dynamic> The vector of normalised temperatures for all spatial points.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &densities,
                                                             const Eigen::Vector<T, Eigen::Dynamic> &mean_velocities,
                                                             T wall_g_val, // Value of g at wall for zeta=0+
                                                             T wall_h_val, // Value of h at wall for zeta=0+
                                                             JacobianFunc get_jacobian);

        /** @brief Computes the macroscopic normalized mean gas temperature at each spatial point.
         *
         * Computes the normalised mean gas temperature:
         * \f[
         * \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left( \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)
         * \f]
         * for all computational points in space. The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * This version first calls the vector-based `compute_density` and `compute_meanGasVelocity` functions to obtain the densities,
         * then uses these results to compute the mean velocities.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols -> spatial points).
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols -> spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param wall_g_val Value of g at wall for zeta=0+ (for temperature calc).
         * @param wall_h_val Value of h at wall for zeta=0+ (for temperature calc).
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return Eigen::Vector<T, Eigen::Dynamic> The vector of normalised temperatures for all spatial points.
         */
        template <typename T, typename JacobianFunc>
        Eigen::Vector<T, Eigen::Dynamic> compute_temperature(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                                             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &h,
                                                             const VelocityMesh<T> &velocity_mesh,
                                                             T wall_g_val, // For velocity & temp
                                                             T wall_h_val, // For temp
                                                             JacobianFunc get_jacobian);

        /** @brief Computes the macroscopic mean gas temperature at a single spatial point.
         *
         * Computes the normalised mean gas temperature at a specific spatial index `i`:
         * \f[
         * \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left( \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)
         * \f]
         *
         * The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * Accepts pre-computed density and mean velocity for efficiency.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols ->
         * spatial points).
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols ->
         * spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param i The spatial index at which to compute the mean gas temperature.
         * @param density_i The pre-computed density at point i.
         * @param mean_velocity_i The pre-computed mean velocity at point i.
         * @param wall_g_val Value of g at wall for zeta=0+ (for temperature calc).
         * @param wall_h_val Value of h at wall for zeta=0+ (for temperature calc).
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return T The computed mean gas temperature at the i-th spatial point.
         *
         */
        template <typename T, typename JacobianFunc>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 const T density_i,       // Pre-computed
                                 const T mean_velocity_i, // Pre-computed
                                 T wall_g_val,            // For temp
                                 T wall_h_val,            // For temp
                                 JacobianFunc get_jacobian);

        /** @brief Computes the macroscopic mean gas temperature at a single spatial point.
         *
         * Computes the normalised mean gas temperature at a specific spatial index `i`:
         * \f[
         * \bar{T} = \frac{2}{3} \bar{\rho}^{-1} \left( \int_{-\infty}^{\infty} (\zeta - \bar{v})^2 g \, d\zeta + \int_{-\infty}^{\infty} h \, d\zeta \right)
         * \f]
         *
         * The integral is computed numerically by splitting it into
         * negative and positive velocity domains. The numerical integration is performed using Simpson's 1/3 Rule and a
         * change of variables. A special treatment is applied for the wall boundary condition at `i = 0` to handle
         * the discontinuity in the distribution function at \f$ \zeta=0 \f$.
         *
         * This version computes the density internally calling the `compute_density_at` and `compute_mean_velocity_at` functions.
         *
         * For further details on the implementation and mathematical formulation, please refer to the density
         * vectorized version.
         *
         * @tparam T floating type precision (e.g., float, double)
         * @tparam JacobianFunc The type of the callable that computes the Jacobian term for the integral.
         * @param g Eigen matrix containing the g quantity (rows -> velocity points, cols ->
         * spatial points).
         * @param h Eigen matrix containing the h quantity (rows -> velocity points, cols ->
         * spatial points).
         * @param velocity_mesh A VelocityMesh object providing the discrete velocity points and grid information.
         * @param i The spatial index at which to compute the mean gas temperature.
         * @param wall_g_val Value of g at wall for zeta=0+ (for temperature calc).
         * @param wall_h_val Value of h at wall for zeta=0+ (for temperature calc).
         * @param get_jacobian A callable (function/lambda) that takes the velocity index $k$ and returns the
         * Jacobian term \f$ J_k = \left| \frac{d\zeta}{dj} \right|_k \cdot \Delta j \f$.
         * @return T The computed mean gas temperature at the i-th spatial point.
         */
        template <typename T, typename JacobianFunc>
        T compute_temperature_at(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &g,
                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &h,
                                 const VelocityMesh<T> &velocity_mesh,
                                 const Eigen::Index i,
                                 T wall_g_val, // For internal velocity & temp
                                 T wall_h_val, // For temp
                                 JacobianFunc get_jacobian);
    }
}

#include "impl/phys_utils.tpp"

#endif /* PHYS_UTILS_B3D31D31_524A_4009_AA77_ABC07533CF00 */
