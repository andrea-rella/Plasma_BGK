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

#include "../SpaceMeshFV.hpp"
#include <ranges>
#include <numeric>
#include <algorithm>
#include <iomanip>

namespace Bgk
{
	template <typename T>
	SpaceMeshFV<T>::SpaceMeshFV(const ConfigData<T> &config) : BaseMesh<T, std::vector<T>>(config), N0(config.get_N0()),
															   d1(config.get_d1()), d2(config.get_d2()) {}

	// ------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------

	template <typename T>
	template <SpacingFunction<T> Spacing>
	void SpaceMeshFV<T>::initialize_with_custom_spacing(Spacing &&spacing_func)
	{
		if (!this->is_constructed)
			throw std::runtime_error("SpaceMesh object not constructed. Call the constructor with ConfigData first.");

		if (this->is_initialized)
			std::cerr << "Mesh already initialized. Skipping initialization." << std::endl;

		// Computational points x^{(i)}
		this->x_comp.resize(this->N + 1);

		for (int i = 0; i <= this->N; ++i)
			this->x_comp[i] = spacing_func(i);

		// ...existing volume boundary code...
		x_vol.resize(this->N + 2);
		x_vol[0] = this->x_comp[0];
		for (int i = 0; i < this->N; ++i)
			x_vol[i + 1] = T{0.5} * (this->x_comp[i] + this->x_comp[i + 1]);
		x_vol[this->N + 1] = this->x_comp[this->N];

		this->is_initialized = true;
		std::cout << "Mesh initialized with " << this->x_comp.size() << " computational points and " << this->x_vol.size() << " volume boundaries." << std::endl;
	}

	// ------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------

	template <typename T>
	void SpaceMeshFV<T>::initialize_mesh()
	{
		auto default_spacing = [this](int i) -> T
		{
			if (this->N <= N0)
			{
				return i * d1 + T{0.25} * (d2 - d1) * std::pow<T>(i, 4) / std::pow<T>(N0, 3);
			}
			else
			{
				if (i <= N0)
					return i * d1 + T{0.25} * (d2 - d1) * std::pow<T>(i, 4) / std::pow<T>(N0, 3);
				else
					return this->x_comp[N0] + (i - N0) * d2;
			}
		};

		initialize_with_custom_spacing(default_spacing);
	}

	// ------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------

	template <typename T>
	void SpaceMeshFV<T>::write_mesh_vtk(const std::string &folder) const
	{
		if (!this->is_initialized)
			throw std::runtime_error("Mesh not initialized. Call initialize_mesh() first.");

		std::filesystem::create_directories("output/" + folder);

		const std::string filename = "output/" + folder + "/space_mesh.vtk";
		std::ofstream vtk_file(filename);
		if (!vtk_file.is_open())
		{
			std::cerr << "Failed to open file for writing: " << filename << std::endl;
			return;
		}

		int num_cells = this->N + 1;		// Number of computational cells (N+1 total)
		int num_points = 4 * (this->N + 1); // 4 corners per quad
		T height = 0.3;						// Height of each block

		vtk_file << "# vtk DataFile Version 3.0\n";
		vtk_file << "1D Computational Mesh\n";
		vtk_file << "ASCII\n";
		vtk_file << "DATASET UNSTRUCTURED_GRID\n";

		// === Write all points ===
		vtk_file << "POINTS " << num_points << " float\n";
		for (int i = 0; i <= this->N; ++i) // Loop from 0 to N (inclusive) for N+1 cells
		{
			T x_left, x_right;
			x_left = x_vol[i];
			x_right = x_vol[i + 1];

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
		for (int i = 0; i <= this->N; ++i)
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
		for (int i = 0; i <= this->N; ++i)
		{
			T width;
			width = x_vol[i + 1] - x_vol[i];
			vtk_file << std::fixed << std::setprecision(6) << width << "\n";
		}

		vtk_file << "SCALARS cell_id int 1\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int i = 0; i < num_cells; ++i)
			vtk_file << i << "\n";

		vtk_file.close();
		std::cout << "2D block mesh written to: " << filename << std::endl;
	}

	// ------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------

	template <typename T>
	void SpaceMeshFV<T>::write_mesh_txt(const std::string &folder) const
	{
		if (!this->is_initialized)
			throw std::runtime_error("Mesh not initialized. Call initialize_mesh() first.");

		std::filesystem::create_directories("output/" + folder);

		const std::string filename = "output/" + folder + "/space_mesh.txt";
		std::ofstream txt_file(filename);
		if (!txt_file.is_open())
		{
			std::cerr << "Failed to open file for writing: " << filename << std::endl;
			return;
		}

		txt_file << "Computational Points (x_comp):\n";
		for (const auto &x : this->x_comp)
			txt_file << x << "\n";

		txt_file << "\nVolume Boundaries (x_vol):\n";
		for (const auto &x : x_vol)
			txt_file << x << "\n";

		txt_file.close();
		std::cout << "Mesh written to: " << filename << std::endl;
	}
}

#endif /* SPACEMESH_IMPL_B6DAD095_AFD5_4569_BA1B_F73464BEC0AA */
