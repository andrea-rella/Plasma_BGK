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

#ifndef BASEMESH_E6CE4ADC_7736_4D8B_95AF_44F071A5FB5A
#define BASEMESH_E6CE4ADC_7736_4D8B_95AF_44F071A5FB5A

#include <concepts>
#include <ranges>
#include <string>
#include "ConfigData.hpp"
#include "utilities.hpp"

namespace Bgk
{

    /**
     * @brief Pure virtual class representing a base 1D mesh.
     *
     * @tparam Container for the mesh components, needs to be a range type.
     * @tparam T Type of the mesh components / values precision. Needs to be a floating point type.
     */
    template <FloatingPoint T, MeshContainer<T> Container = std::vector<T>>
    class BaseMesh
    {
    protected:
        int N;
        Container x_comp;

        bool is_initialized = false;
        bool is_constructed = false;

    public:
        // ----- Constructors and destructors -----------------------------------------------
        // ----------------------------------------------------------------------------------

        BaseMesh() = default;
        BaseMesh(ConfigData<T> config) : N(config.get_N()), is_initialized(false), is_constructed(true) {};
        virtual ~BaseMesh() = default;

        // ------ Initialization and validation methods -------------------------------------
        // ----------------------------------------------------------------------------------

        virtual void initialize_mesh() = 0;

        virtual bool validate_mesh() const
        {
            return is_initialized && !x_comp.empty();
        }

        virtual void reset_mesh()
        {
            x_comp.clear();
            is_initialized = false;
        }

        bool isInitialized() const { return is_initialized; }
        bool isConstructed() const { return is_constructed; }

        // ------ Getters for mesh parameters -----------------------------------------------
        // ----------------------------------------------------------------------------------

        int get_N() const { return N; }

        // ------ Getters for mesh components -----------------------------------------------
        // ----------------------------------------------------------------------------------
        const Container &get_XComp() const { return x_comp; }

        auto begin() const { return std::ranges::begin(x_comp); }
        auto end() const { return std::ranges::end(x_comp); }
        auto size() const { return std::ranges::size(x_comp); }
    };
}

#endif /* BASEMESH_E6CE4ADC_7736_4D8B_95AF_44F071A5FB5A */
