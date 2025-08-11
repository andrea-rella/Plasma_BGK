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

#ifndef BASEMESH1D_D265AE3B_FC27_492A_A8A8_6B79C689F620
#define BASEMESH1D_D265AE3B_FC27_492A_A8A8_6B79C689F620

#include <concepts>
#include <ranges>
#include <string>
#include "ConfigData.hpp"
#include "utilities.hpp"
#include <cassert>

namespace Bgk
{

    /**
     * @brief Pure virtual class representing a base 1D mesh.
     *
     * @tparam Container for the mesh components, needs to be a range type. Defaulted to std::vector<T>.
     * @tparam T Type of the mesh components / values precision. Needs to be a floating point type.
     * @tparam Nature The nature of the mesh, either SPACE or TIME. Defaulted to SPACE.
     */
    template <FloatingPoint T, MeshContainer<T> Container = std::vector<T>, MeshNature Nature = MeshNature::SPACE>
    class BaseMesh1D
    {
    protected:
        size_t N;
        Container x_comp;

        bool is_initialized = false;
        bool is_constructed = false;

    public:
        // ----- CONSTRUCTORS AND DESTRUCTORS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        BaseMesh1D() = default;
        BaseMesh1D(ConfigData<T> config) : is_initialized(false), is_constructed(true)
        {
            if constexpr (Nature == MeshNature::SPACE)
                N = config.get_N();
            else if constexpr (Nature == MeshNature::VELOCITY)
                N = config.get_BarN();
            else
                static_assert(Nature == MeshNature::SPACE || Nature == MeshNature::VELOCITY, "Unsupported MeshNature type");
        };
        virtual ~BaseMesh1D() = default;

        // ------ INITIALIZATION -------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        virtual void initialize_mesh() = 0;

        virtual bool validate_mesh() const
        {
            return is_initialized && !x_comp.empty();
        }

        /** * @brief Resets the mesh components and initialization status.
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

        size_t get_N() const { return N; }

        // ------ GETTERS FOR MESH COMPONENTS ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------
        const Container &get_XComp() const { return x_comp; }

        auto begin() const { return std::ranges::begin(x_comp); }
        auto end() const { return std::ranges::end(x_comp); }
        auto size() const { return std::ranges::size(x_comp); }

        // ------ OPERATORS ------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /**
         * @brief Access operator to get the computational point at a specific index. Read only.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @pre index must be within the dimension of the mesh (<= N).
         */
        virtual T operator[](size_t index) const
        {
            assert(index <= N && "Index out of range in BaseMesh1D::operator[]");
            return x_comp[index];
        }

        /**
         * @brief Access operator to get the computational point at a specific index.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @pre index must be within the dimension of the mesh (<= N).
         *
         * @overload
         */
        virtual T &operator[](size_t index)
        {
            assert(index <= N && "Index out of range in BaseMesh1D::operator[]");
            return x_comp[index];
        }

        /**
         * @brief Bounds-checked access to the computational point at a specific index. Read only.
         *
         * @param index The index of the computational point to access.
         * @return The computational point at the specified index.
         *
         * @throws std::out_of_range if index is greater than N.
         *
         * @note This method provides safe access with bounds checking, unlike operator[].
         */
        virtual T at(size_t index) const
        {
            if (index > N)
                throw std::out_of_range("Index out of range in BaseMesh1D::at()");

            return x_comp[index];
        }

        /**
         * @brief Bounds-checked access to the computational point at a specific index.
         *
         * @param index The index of the computational point to access.
         * @return Reference to the computational point at the specified index.
         *
         * @throws std::out_of_range if index is greater than N.
         *
         * @note This method provides safe access with bounds checking, unlike operator[].
         *
         * @overload
         */
        virtual T &at(size_t index)
        {
            if (index > N)
                throw std::out_of_range("Index out of range in BaseMesh1D::at()");

            return x_comp[index];
        }

        // ------ OUTPUT MESH ---------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        /**
         * @brief Writes the mesh to a text file.
         *
         * It creates a .txt file containing the computational points of the mesh. The file will be saved as
         * output/<folder_name>/<Nature>_mesh.txt and it will contain the index and the corresponding
         * computational point.
         *
         * @param folder_name The name of the directory.
         *
         *
         * @throws std::runtime_error if the mesh is not initialized.
         *
         * @note If not already present, the output directory will be created.
         */
        virtual void write_mesh_txt(const std::string &folder_name) const
        {
            if (!this->is_initialized)
                throw std::runtime_error("Trying to write to .txt a mesh wich is not initialized. Call initialize_mesh() first. [" +
                                         std::string(__FILE__) + ":" + std::to_string(__LINE__) + "]");

            std::filesystem::create_directories("output/" + folder_name);

            std::string filename;
            if constexpr (Nature == MeshNature::SPACE)
                filename = "output/" + folder_name + "/space_mesh.txt";
            else if constexpr (Nature == MeshNature::VELOCITY)
                filename = "output/" + folder_name + "/velocity_mesh.txt";
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

#endif /* BASEMESH1D_D265AE3B_FC27_492A_A8A8_6B79C689F620 */
