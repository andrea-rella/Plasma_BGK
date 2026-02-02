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

#ifndef VELOCITYMESH_FA26FD6E_B376_43B5_B21B_CDE75DDA7C7F
#define VELOCITYMESH_FA26FD6E_B376_43B5_B21B_CDE75DDA7C7F

#include "ConfigData.hpp"
#include "BaseMesh1D.hpp"
#include "utilities.hpp"

namespace Bgk
{
    /** @brief Class representing a one-dimensional velocity mesh (extends BaseMesh1D<T, std::vector<T>>, MeshNature::VELOCITY).
     *
     * This class implements a 1D velocity mesh for the BGK model. It defines the computational points
     * based on the configuration data provided in a symmetric way around zero. It offers methods for
     * initializing the mesh with custom or default spacing.
     *
     * @tparam T The data type for the mesh points (e.g., float, double).
     */
    template <typename T>
    class VelocityMesh : public BaseMesh1D<T, std::vector<T>, MeshNature::VELOCITY>
    {
    private:
        // Spacing parameter
        T a1;
        // Spacing parameter
        T a2;

        std::function<T(size_t)> jacobian_func;

    public:
        // ----- CONSTRUCTORS AND DESTRUCTORS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Default constructor.
        VelocityMesh() = default;

        /** @brief Constructor based on configuration data.
         *
         * @param config Configuration data object containing mesh parameters.
         */
        VelocityMesh(const ConfigData<T> &config);

        /// @brief Virtual destructor.
        ~VelocityMesh() = default;

        // ------ SETTERS --------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Sets the jacobian function associated with the mesh.
         *
         * @param jacobian A callable object that defines the jacobian function (see SpacingFunction concept).
         *
         * @note This method relies on universal references to allow for flexibility in the type
         *       of jacobian function provided.
         */
        template <SpacingFunction<T> Jacobian>
        void set_jacobian_function(Jacobian &&jacobian);

        // ------ GETTERS --------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Returns the jacobian function associated with the mesh.
        /// @return The jacobian function. (std::function<T(size_t)>)
        std::function<T(size_t)> get_jacobian_function() const { return jacobian_func; }

        // ------ OPERATORS ------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Access operator to get the computational point at a specific index. Read only.
         *
         * The method is overridden w.r.t. BaseMesh1D to account for the specific size of the velocity mesh
         * which in this case, given the symmetry around zero, contains 2 * N + 1 points.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @pre index must be within the dimension of the mesh (<= 2 * N + 1).
         */
        T operator[](size_t index) const override;

        /** @brief Access operator to get the computational point at a specific index.
         *
         * The method is overridden w.r.t. BaseMesh1D to account for the specific size of the velocity mesh
         * which in this case, given the symmetry around zero, contains 2 * N + 1 points.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @pre index must be within the dimension of the mesh (<= 2 * N + 1).
         *
         * @overload
         */
        T &operator[](size_t index) override;

        /** @brief Bounds-checked access to the computational point at a specific index. Read only.
         *
         * The method is overridden w.r.t. BaseMesh1D to account for the specific size of the velocity mesh
         * which in this case, given the symmetry around zero, contains 2 * N + 1 points.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @throws std::out_of_range if index is greater than 2 * N + 1.
         *
         * @note This method provides safe access with bounds checking, unlike operator[].
         */
        T at(size_t index) const override;

        /** @brief Bounds-checked access to the computational point at a specific index.
         *
         * The method is overridden w.r.t. BaseMesh1D to account for the specific size of the velocity mesh
         * which in this case, given the symmetry around zero, contains 2 * N + 1 points.
         *
         * @param index The index of the computational point to access.
         * @return Reference to the computational point at the specified index.
         *
         * @throws std::out_of_range if index is greater than 2 * N + 1.
         *
         * @note This method provides safe access with bounds checking, unlike operator[].
         *
         * @overload
         */
        T &at(size_t index) override;

        // ------ INITIALIZATION -------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Initializes the mesh with custom spacing.
         *
         * This method allows the user to define a custom spacing function (i.e. a function that given i return x_i)
         * for the mesh then the computational points are computed based on the configuration data. The mesh is
         * symmetric around zero, with points at negative indices mirroring those at positive indices; the spacing
         * is supposed to be defined for the positive indices only.
         *
         * @param spacing_func A callable object that defines the custom spacing (see SpacingFunction concept).
         *
         * @note This method relies on universal references to allow for flexibility in the type
         *       of spacing function provided.
         * @throws std::runtime_error if the VelocityMeshFV object is not constructed.
         */
        template <SpacingFunction<T> Spacing>
        void initialize_with_custom_spacing(Spacing &&spacing_func);

        /** @brief Initializes the mesh the mesh by filling x_comp with the default spacing (see note).
         *
         * This method computes the computational points based on the configuration data calling
         * initialize_with_custom_spacing. It uses the polynomial spacing following the schedule
         * defined in the note
         *
         * @throws std::runtime_error if the VelocityMeshFV object is not constructed.
         *
         * @note The mesh is initialized with the following formula @cite aoki1990numerical:
         *       @f{align*}
         *       zeta_j = a_1 \cdot j + a_2 \cdot j^3, \quad j = -2\overline{N}\dots,2\overline{N} \\
         *       @f}
         *
         */
        void initialize_mesh() override;
    };
}

#include "impl/VelocityMesh.tpp"

#endif /* VELOCITYMESH_FA26FD6E_B376_43B5_B21B_CDE75DDA7C7F */
