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

    template <typename T>
    class VelocityMesh : public BaseMesh1D<T, std::vector<T>, MeshNature::VELOCITY>
    {
    private:
        T a1;
        T a2;

    public:
        // ----- Constructors and destructors -----------------------------------------------
        // ----------------------------------------------------------------------------------

        VelocityMesh() = default;
        VelocityMesh(const ConfigData<T> &config);
        ~VelocityMesh() = default;

        // ------ Initialization and validation methods -------------------------------------
        // ----------------------------------------------------------------------------------

        /**
         * @brief Initializes the mesh with custom spacing.
         *
         * This method allows the user to define a custom spacing function for the mesh then the computational
         * points are computed based on the configuration data. The mesh is symmetric around zero, with points
         * at negative indices mirroring those at positive indices; the spacing is supposed to be defined for
         * the positive indices only.
         *
         * @tparam SpacingFunc A callable type that takes an integer index and returns a T value
         *         representing the spacing at that index.
         * @param spacing_func A callable object that defines the custom spacing.
         *
         * @note This method relies on universal references to allow for flexibility in the type
         *       of spacing function provided.
         * @throws std::runtime_error if the VelocityMeshFV object is not constructed.
         */
        template <SpacingFunction<T> Spacing>
        void initialize_with_custom_spacing(Spacing &&spacing_func);

        /**
         * @brief Initializes the mesh the mesh by filling x_comp with the default spacing (see note).
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

#include "impl/VelocityMesh_impl.hpp"

#endif /* VELOCITYMESH_FA26FD6E_B376_43B5_B21B_CDE75DDA7C7F */
