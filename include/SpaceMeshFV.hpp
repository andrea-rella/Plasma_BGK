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

#ifndef SPACEMESHFV_ADAC74E0_9DD7_4843_93D3_6435A47FB6E7
#define SPACEMESHFV_ADAC74E0_9DD7_4843_93D3_6435A47FB6E7

#include <vector>
#include "ConfigData.hpp"
#include "BaseMesh1D.hpp"
#include "utilities.hpp"

namespace Bgk
{
    /** * @brief 1D Finite Volume Spache Mesh Class (extends BaseMesh1D<T, std::vector<T>>, MeshNature::SPACE).
     *
     * This class implements a 1D finite volume space mesh for the BGK model. It defines the computational
     * points and volume boundaries based on the configuration data provided.
     *
     * @tparam T Precision type for the mesh (e.g., float, double).
     *
     * @extends BaseMesh1D<T, std::vector<T>, MeshNature::SPACE>
     * @see BaseMesh1D
     *
     */
    template <typename T>
    class SpaceMeshFV : public BaseMesh1D<T, std::vector<T>, MeshNature::SPACE>
    {
    private:
        /// Number of points in the polynomially stretched region.
        size_t N0;
        /// Spacing parameters.
        T d1;
        T d2;

        /// Volume boundaries
        std::vector<T> x_vol;
        /// Volume sizes
        std::vector<T> vol_sizes;

    public:
        // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Default constructor.
        SpaceMeshFV() = default;

        /** * @brief Constructor based on configuration data.
         *
         * @param config Configuration data object containing mesh parameters.
         */
        SpaceMeshFV(const ConfigData<T> &);

        /// @brief Virtual destructor.
        ~SpaceMeshFV() = default;

        // ---------- MESH INITIALIZATION METHODS --------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /**
         * @brief Initializes the mesh with custom spacing.
         *
         * This method allows the user to define a custom spacing function for the mesh then the computational
         * points, volume boundaries and volume sizesare computed based on the configuration data.
         *
         * @tparam SpacingFunc A callable type that takes an integer index and returns a T value
         *         representing the spacing at that index.
         * @param spacing_func A callable object that defines the custom spacing.
         *
         * @note This method relies on universal references to allow for flexibility in the type
         *       of spacing function provided.
         * @throws std::runtime_error if the SpaceMeshFV object is not constructed.
         */
        template <SpacingFunction<T> Spacing>
        void initialize_with_custom_spacing(Spacing &&spacing_func);

        /**
         * @brief Initializes the mesh the mesh by filling x_comp and x_vol vectors with the default spacing (see note).
         *
         * This method computes the computational points and volume boundaries based on the configuration data
         * calling initialize_with_custom_spacing. It uses the polynomial spacing for the first N0 points and
         * uniform spacing for the remaining points.
         *
         * @throws std::runtime_error if the SpaceMeshFV object is not constructed.
         *
         * @note The mesh is initialized with the following formula @cite aoki1990numerical
         *       @f{align*}
         *       x_i = d_1 i + 0.25 (d_2-d_1)(i^4/N_0^3), \quad i = 0,\dots,N_0 \\
         *       x_i = x_{N_0} + (i-N_0) d_2, \quad i = N_0+1,\dots,N \\
         *       x_{i+1/2} = 0.5 (x_i + x_{i+1}), \quad i = 0,\dots,N-1
         *       @f}
         *
         */
        void initialize_mesh() override;

        // ---------- GETTERS ----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Returns the computational points vector.
        const std::vector<T> &get_volume_boundaries() const { return x_vol; }

        /// @brief Returns the volume sizes vector.
        const std::vector<T> &get_volume_sizes() const { return vol_sizes; }

        // ---------- OUTPUT METHODS ---------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /**
         * @brief Writes the mesh to a text file.
         *
         * It creates a .txt file containing the computational points and volume boundaries of the mesh.
         * The file will be saved as output/<folder_name>/space_mesh.txt
         *
         * @param folder_name The name of the directory.
         *
         * @throws std::runtime_error if the mesh is not initialized.
         *
         * @note If not already present, the output directory will be created.
         */
        void write_mesh_txt(const std::string &folder_name) const override;

        /**
         * @brief Writes the mesh to a VTK file.
         *
         * It creates a .vtk file for visualization in VTK-compatible software such as ParaView or VisIt.
         * The file will contain the computational points and volume boundaries of the mesh and it will be
         * saved as output/<folder_name>/space_mesh.vtk
         *
         * @param folder_name The name of the directory.
         *
         * @throws std::runtime_error if the mesh is not initialized.
         *
         * @note The VTK file format is structured to allow visualization of the mesh in 3D space.
         * @note If not already present, the output directory will be created.
         */
        void write_mesh_vtk(const std::string &folder_name) const;
    };
}

#include "impl/SpaceMeshFV.tpp"

#endif /* SPACEMESHFV_ADAC74E0_9DD7_4843_93D3_6435A47FB6E7 */
