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

#ifndef SPACEMESH_IMPL_B6DAD095_AFD5_4569_BA1B_F73464BEC0AA
#define SPACEMESH_IMPL_B6DAD095_AFD5_4569_BA1B_F73464BEC0AA

#include "../SpaceMesh.hpp"
#include <ranges>
#include <numeric>
#include <algorithm>
#include <iomanip>

namespace Bgk
{
    template <typename T>
    SpaceMesh<T>::SpaceMesh(const ConfigData<T> &config) : D(config.get_D()), N(config.get_N()), N0(config.get_N0()),
                                                           d1(config.get_d1()), d2(config.get_d2()), is_constructed{true} {};

    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------

    template <typename T>
    void SpaceMesh<T>::initialize_mesh()
    {
        if (!is_constructed)
            throw std::runtime_error("SpaceMesh not constructed. Call the constructor with ConfigData first.");

        if (is_initialized)
            std::cerr << "Mesh already initialized. Skipping initialization." << std::endl;

        // Computational points x^{(i)}

        x_comp.resize(N + 1);

        if (N <= N0)
        {
            for (int i = 0; i <= N; ++i)
                x_comp[i] = i * d1 + T{0.25} * (d2 - d1) * std::pow<T>(i, 4) / std::pow<T>(N0, 3);
        }
        else
        {
            // First part: polynomial spacing
            for (int i = 0; i <= N0; ++i)
                x_comp[i] = i * d1 + T{0.25} * (d2 - d1) * std::pow<T>(i, 4) / std::pow<T>(N0, 3);

            // Second part: uniform spacing
            for (int i = N0 + 1; i <= N; ++i)
                x_comp[i] = x_comp[N0] + (i - N0) * d2;
        }

        // Volumes boundaries x^{(i+1/2)}
        x_vol.resize(N);
        for (int i = 0; i < N; ++i)
            x_vol[i] = T{0.5} * (x_comp[i] + x_comp[i + 1]);

        is_initialized = true;
        std::cout << "Mesh initialized with " << N + 1 << " computational points and " << N << " volume boundaries." << std::endl;
    }

    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------

    template <typename T>
    void SpaceMesh<T>::write_mesh_vtk(const std::string &name) const
    {
        if (!is_initialized)
            throw std::runtime_error("Mesh not initialized. Call initialize_mesh() first.");

        std::filesystem::create_directories("output/" + name);

        const std::string filename = "output/" + name + "/space_mesh.vtk";
        std::ofstream vtk_file(filename);
        if (!vtk_file.is_open())
        {
            std::cerr << "Failed to open file for writing: " << filename << std::endl;
            return;
        }

        int num_cells = N + 1;        // Number of computational cells (N+1 total)
        int num_points = 4 * (N + 1); // 4 corners per quad
        T height = 0.3;               // Height of each block

        vtk_file << "# vtk DataFile Version 3.0\n";
        vtk_file << "1D Computational Mesh\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET UNSTRUCTURED_GRID\n";

        // === Write all points ===
        vtk_file << "POINTS " << num_points << " float\n";
        for (int i = 0; i <= N; ++i) // Loop from 0 to N (inclusive) for N+1 cells
        {
            T x_left, x_right;
            if (i == 0)
            {
                x_left = x_comp[0];
                x_right = x_vol[0];
            }
            else if (i == N)
            {
                x_left = x_vol[N - 1];
                x_right = x_comp[N];
            }
            else
            {
                x_left = x_vol[i - 1];
                x_right = x_vol[i];
            }

            T y_bottom = -height / 2;
            T y_top = height / 2;

            // Each quad has 4 points: bottom-left, bottom-right, top-right, top-left
            vtk_file << std::fixed << std::setprecision(6)
                     << x_left << " " << y_bottom << " 0.0\n";
            vtk_file << std::fixed << std::setprecision(6)
                     << x_right << " " << y_bottom << " 0.0\n";
            vtk_file << std::fixed << std::setprecision(6)
                     << x_right << " " << y_top << " 0.0\n";
            vtk_file << std::fixed << std::setprecision(6)
                     << x_left << " " << y_top << " 0.0\n";
        }

        // === Write each quad (cell) ===
        vtk_file << "\nCELLS " << num_cells << " " << (5 * num_cells) << "\n";
        for (int i = 0; i <= N; ++i)
        {
            int idx = i * 4;
            vtk_file << "4 " << idx << " " << idx + 1 << " " << idx + 2 << " " << idx + 3 << "\n";
        }

        // === Cell types ===
        vtk_file << "\nCELL_TYPES " << num_cells << "\n";
        for (int i = 0; i < num_cells; ++i)
            vtk_file << "9\n"; // VTK_QUAD

        // === Add cell data to visualize mesh structure ===

        // Add cell width data for better visualization
        vtk_file << "\nCELL_DATA " << num_cells << "\n";
        vtk_file << "SCALARS cell_width float 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (int i = 0; i <= N; ++i)
        {
            T width;
            if (i == 0)
                width = x_vol[0] - x_comp[0];
            else if (i == N)
                width = x_comp[N] - x_vol[N - 1];
            else
                width = x_vol[i] - x_vol[i - 1];
            vtk_file << std::fixed << std::setprecision(6) << width << "\n";
        }

        vtk_file << "SCALARS cell_id int 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (int i = 0; i < num_cells; ++i)
            vtk_file << i << "\n";

        vtk_file.close();
        std::cout << "2D block mesh written to: " << filename << std::endl;
    }

}

#endif /* SPACEMESH_IMPL_B6DAD095_AFD5_4569_BA1B_F73464BEC0AA */
