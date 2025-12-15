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

#ifndef VELOCITYMESH_IMPL_D96EB33E_0366_4848_A105_C24157A1D6BB
#define VELOCITYMESH_IMPL_D96EB33E_0366_4848_A105_C24157A1D6BB

#include "../VelocityMesh.hpp"
#include <ranges>
#include <numeric>
#include <iomanip>
#include "../utilities.hpp"

namespace Bgk
{
    template <typename T>
    VelocityMesh<T>::VelocityMesh(const ConfigData<T> &config) : BaseMesh1D<T, std::vector<T>, MeshNature::VELOCITY>(config),
                                                                 a1(config.get_a1()), a2(config.get_a2()){};

    // ------ OPERATORS ------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    T VelocityMesh<T>::operator[](size_t index) const
    {
        assert(index <= 2 * this->N && "Index out of range in VelocityMesh::operator[]");
        return this->x_comp[index];
    }

    template <typename T>
    T &VelocityMesh<T>::operator[](size_t index)
    {
        assert(index <= 2 * this->N && "Index out of range in VelocityMesh::operator[]");
        return this->x_comp[index];
    }

    template <typename T>
    T VelocityMesh<T>::at(size_t index) const
    {
        if (index > 2 * this->N)
            throw std::out_of_range("Index out of range in VelocityMesh::at");
        return this->x_comp[index];
    }

    template <typename T>
    T &VelocityMesh<T>::at(size_t index)
    {
        if (index > 2 * this->N)
            throw std::out_of_range("Index out of range in VelocityMesh::at");
        return this->x_comp[index];
    }

    // ------ INITIALIZATION -------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------

    template <typename T>
    template <SpacingFunction<T> Spacing>
    void VelocityMesh<T>::initialize_with_custom_spacing(Spacing &&spacing_func)
    {
        if (!this->is_constructed)
            throw std::runtime_error(error_message(
                "Velocity Mesh object not constructed. Call the constructor with ConfigData first."));

        if (this->is_initialized)
        {
            std::cerr << "Velocity Mesh already initialized. Skipping initialization." << std::endl;
            return;
        }

        // Computational points x^{(i)}
        this->x_comp.resize(2 * this->N + 1);

        this->x_comp[this->N] = T{0}; // Center point at zero

        T value;

        // initialize positive side and mirror to negative side
        for (size_t i = this->N + 1; i <= 2 * this->N; ++i)
        {
            value = spacing_func(i - this->N);
            this->x_comp[i] = value;
            this->x_comp[2 * this->N - i] = -value; // Symmetric point
        }

        this->is_initialized = true;
    }

    // ------------------------------------------------------------------------------

    template <typename T>
    void VelocityMesh<T>::initialize_mesh()
    {
        auto default_spacing = [this](size_t i) -> T
        {
            return a1 * i + a2 * std::pow<T>(i, 3);
        };

        initialize_with_custom_spacing(default_spacing);
    }
}

#endif /* VELOCITYMESH_IMPL_D96EB33E_0366_4848_A105_C24157A1D6BB */
