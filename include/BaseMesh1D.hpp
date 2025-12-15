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

#ifndef BASEMESH1D_AA966284_DA7A_4358_8C05_DB9CDA497B65
#define BASEMESH1D_AA966284_DA7A_4358_8C05_DB9CDA497B65

#include <concepts>
#include <ranges>
#include <string>
#include "ConfigData.hpp"
#include "utilities.hpp"
#include <cassert>

namespace Bgk
{

    /** @brief Pure virtual base class for 1D meshes.
     *
     * This class serves as a base for different types of 1D meshes, providing common functionality such as
     * initialization, validation, and access to mesh parameters and components. As base members it stores the
     * number of points N and a container (defined with a concept that requires basic iterable functionality)
     * for the computational points x_comp.
     *
     * @tparam T A floating-point type for mesh coordinates.
     * @tparam Container A container type for storing mesh points.
     * @tparam Nature The nature of the mesh, either SPACE or VELOCITY.
     */
    template <FloatingPoint T, MeshContainer1D<T> Container = std::vector<T>, MeshNature Nature = MeshNature::SPACE>
    class BaseMesh1D
    {
    protected:
        /// @brief Number of points in the mesh.
        size_t N;
        /// @brief Container for computational points of the mesh.
        Container x_comp;

        /// @brief Flag indicating whether the mesh has been constructed.
        bool is_constructed = false;
        /// @brief Flags indicating whether the mesh container(s) have been filled.
        bool is_initialized = false;

    public:
        // ----- CONSTRUCTORS AND DESTRUCTORS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Default constructor.
        BaseMesh1D() = default;

        /// @brief Constructor based on configuration data.
        BaseMesh1D(ConfigData<T> config) : is_initialized(false), is_constructed(true)
        {
            if constexpr (Nature == MeshNature::SPACE)
                N = config.get_N();
            else if constexpr (Nature == MeshNature::VELOCITY)
                N = config.get_BarN();
            else
                static_assert(Nature == MeshNature::SPACE || Nature == MeshNature::VELOCITY, "Unsupported MeshNature type");
        };

        /// @brief Virtual destructor.
        virtual ~BaseMesh1D() = default;

        // ------ INITIALIZATION -------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Pure virtual method to initialize the mesh.
         *
         * This method must be implemented by derived classes to set up the mesh points based on specific
         * criteria or algorithms.
         */
        virtual void initialize_mesh() = 0;

        /** @brief Validates the mesh initialization.
         *
         * This method checks whether the mesh has been properly initialized by verifying that the
         * initialization flag is set and that the computational points container is not empty.
         *
         * @return true if the mesh is initialized and valid; false otherwise.
         */
        virtual bool validate_mesh() const
        {
            return is_initialized && !x_comp.empty();
        }

        /** @brief Resets the mesh components and initialization status.
         *
         * This method clears the computational points and sets the initialization flag to false.
         * It can be used to reinitialize the mesh with new parameters or configurations.
         *
         * @note - This method does not reset the mesh parameters and neither the constructed status.
         * @note - This method does not deallocate memory
         */
        virtual void reset_mesh()
        {
            x_comp.clear();
            is_initialized = false;
        }

        /// @brief  Returns boolean indicating whether the mesh is initialized.
        bool isInitialized() const { return is_initialized; }

        /// @brief  Returns boolean indicating whether the mesh is constructed.
        bool isConstructed() const { return is_constructed; }

        // ------ GETTERS FOR MESH PARAMETERS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Returns the total number of computational points in the mesh.
        /// @return The number of points, as a size_t value.
        size_t get_N() const { return N; }

        // ------ SETTERS FOR MESH PARAMETERS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Sets the total number of computational points in the mesh.
         *
         * @param n The number of points to set.
         */
        void set_N(size_t n) { N = n; }

        // ------ GETTERS FOR MESH COMPONENTS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Returns a constant reference to the container holding the computational points of the mesh.
         *
         * @return const Container&
         */
        const Container &get_XComp() const { return x_comp; }

        /** @brief Returns an iterator to the beginning of the computational points container.
         *
         * @return auto
         */
        auto begin() const { return std::ranges::begin(x_comp); }

        /** @brief Returns an iterator to the end of the computational points container.
         *
         * @return auto
         */
        auto end() const { return std::ranges::end(x_comp); }

        /** @brief Returns the size of the computational points container.
         *
         * @return auto
         */
        auto size() const { return std::ranges::size(x_comp); }

        // ------ OPERATORS ------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Access operator to get the computational point at a specific index.
         *
         * This version returns a non-const reference, allowing the point's value to be modified.
         *
         * @param index The index of the computational point to access.
         * @return A writable reference to the computational point at the specified index.
         *
         * @note index must be within the dimension of the mesh (<= N).
         */
        virtual T &operator[](size_t index)
        {
            assert(index <= N && "Index out of range in BaseMesh1D::operator[]");
            return x_comp[index];
        }

        /** @brief Access operator to get the computational point at a specific index (Read-Only).
         *
         * This is the const-qualified overload. It returns a copy of the value, preventing modification
         * of the mesh data.
         *
         * @param index The index of the computational point to access.
         * @return A copy of the computational point value at the specified index.
         *
         * @note index must be within the dimension of the mesh (<= N).
         *
         * @overload
         */
        virtual T operator[](size_t index) const
        {
            assert(index <= N && "Index out of range in BaseMesh1D::operator[]");
            return x_comp[index];
        }

        /** @brief Bounds-checked access to the computational point at a specific index.
         *
         * This version returns a non-const reference, allowing the point's value to be modified.
         *
         * @param index The index of the computational point to access.
         * @return Reference to the computational point at the specified index.
         *
         * @throws std::out_of_range if index is greater than N.
         *
         * @note This method provides safe access with bounds checking, unlike operator[].
         *
         */
        virtual T &at(size_t index)
        {
            if (index > N)
                throw std::out_of_range("Index out of range in BaseMesh1D::at()");

            return x_comp[index];
        }

        /** @brief Bounds-checked access to the computational point at a specific index. Read only.
         *
         * This is the const-qualified overload. It returns a copy of the value, preventing modification.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @throws std::out_of_range if index is greater than N.
         *
         * @note This method provides safe access with bounds checking, unlike operator[].
         *
         * @overload
         */
        virtual T at(size_t index) const
        {
            if (index > N)
                throw std::out_of_range("Index out of range in BaseMesh1D::at()");

            return x_comp[index];
        }

        // ------ OUTPUT MESH ----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Writes the mesh to a text file.
         *
         * It creates a .txt file containing the computational points of the mesh. The file will be saved as
         * <folder_path>/<Nature>_mesh.txt and it will contain the index and the corresponding
         * computational point.
         *
         * @param folder_path The path of the directory where the mesh file will be saved.
         *
         * @throws std::runtime_error if the mesh is not initialized.
         *
         * @note If not already present, the output directory will be created.
         *
         * @example write_mesh_txt("./test1")
         */
        virtual void write_mesh_txt(const std::string &folder_path) const
        {
            if (!this->is_initialized)
                throw std::runtime_error(
                    error_message("Trying to write to .txt a mesh wich is not initialized. Call initialize_mesh() first."));

            std::filesystem::create_directories(folder_path);

            std::string filename;
            if constexpr (Nature == MeshNature::SPACE)
                filename = folder_path + "/space_mesh.txt";
            else if constexpr (Nature == MeshNature::VELOCITY)
                filename = folder_path + "/velocity_mesh.txt";
            else
                // Fallback for any other nature, if applicable
                static_assert(Nature == MeshNature::SPACE || Nature == MeshNature::VELOCITY, "Unsupported MeshNature type");

            std::ofstream txt_file(filename);
            if (!txt_file.is_open())
            {
                std::cerr << "Failed to open file for writing: " << filename << std::endl;
                return;
            }

            txt_file << "Computational Points (index - x_comp):\n";
            for (std::size_t i = 0; i < x_comp.size(); ++i)
                txt_file << i << "    " << this->x_comp[i] << "\n";
        }
    };
}

#endif /* BASEMESH1D_AA966284_DA7A_4358_8C05_DB9CDA497B65 */
