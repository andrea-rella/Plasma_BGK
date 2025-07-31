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
	SpaceMeshFV<T>::SpaceMeshFV(const ConfigData<T> &config) : BaseMesh1D<T, std::vector<T>>(config), N0(config.get_N0()),
															   d1(config.get_d1()), d2(config.get_d2()) {}

	// ---------- MESH INITIALIZATION METHODS --------------------------------------------------------
	// -----------------------------------------------------------------------------------------------

	template <typename T>
	template <SpacingFunction<T> Spacing>
	void SpaceMeshFV<T>::initialize_with_custom_spacing(Spacing &&spacing_func)
	{
		if (!this->is_constructed)
			throw std::runtime_error("SpaceMesh object not constructed. Call the constructor with ConfigData first.");

		if (this->is_initialized)
			std::cerr << "SpaceMesh already initialized. Skipping initialization." << std::endl;

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
	}

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

	// ---------- GETTERS ----------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------

	template <typename T>
	std::vector<T> SpaceMeshFV<T>::get_volume_sizes() const
	{
		if (!this->is_initialized)
			throw std::runtime_error("Mesh not initialized. Call initialize_mesh() first.");

		std::vector<T> volume_sizes(this->N + 1);
		for (int i = 0; i <= this->N; ++i)
			volume_sizes[i] = x_vol[i + 1] - x_vol[i];

		return volume_sizes;
	}

	// ---------- OUTPUT METHODS ---------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------

	template <typename T>
	void SpaceMeshFV<T>::write_mesh_txt(const std::string &folder_name) const
	{
		BaseMesh1D<T, std::vector<T>, MeshNature::SPACE>::write_mesh_txt(folder_name);

		const std::string filename = "output/" + folder_name + "/volume_boundaries_space.txt";
		std::ofstream txt_file(filename);
		if (!txt_file.is_open())
		{
			std::cerr << "Failed to open file for writing: " << filename << std::endl;
			return;
		};
		txt_file << "Volume Boundaries (index - x_vol):\n";
		for (std::size_t i = 0; i < x_vol.size(); ++i)
			txt_file << i << "    " << x_vol[i] << "\n";
	}

	// ------------------------------------------------------------------------------

	template <typename T>
	void SpaceMeshFV<T>::write_mesh_vtk(const std::string &folder) const
	{
		if (!this->is_initialized)
			throw std::runtime_error("Mesh not initialized. Call initialize_mesh() first.");

		constexpr T BLOCK_HEIGHT = 0.3;
		constexpr int VTK_QUAD_TYPE = 9;

		std::filesystem::create_directories("output/" + folder);
		const std::string filename = "output/" + folder + "/space_mesh.vtk";
		std::ofstream vtk_file(filename);

		if (!vtk_file)
		{
			std::cerr << "Failed to open file for writing: " << filename << std::endl;
			return;
		}

		vtk_file << std::fixed << std::setprecision(6);
		int num_cells = this->N + 1;
		int num_points = 4 * num_cells;
		vtk_file << "# vtk DataFile Version 3.0\n"
				 << "1D Computational Mesh\n"
				 << "ASCII\n"
				 << "DATASET UNSTRUCTURED_GRID\n"
				 << "POINTS " << num_points << " float\n";

		for (int i = 0; i < num_cells; ++i)
		{
			T x_left = x_vol[i];
			T x_right = x_vol[i + 1];
			T y_bottom = -BLOCK_HEIGHT / 2;
			T y_top = BLOCK_HEIGHT / 2;

			vtk_file << x_left << " " << y_bottom << " 0.0\n"
					 << x_right << " " << y_bottom << " 0.0\n"
					 << x_right << " " << y_top << " 0.0\n"
					 << x_left << " " << y_top << " 0.0\n";
		}

		vtk_file << "\nCELLS " << num_cells << " " << (5 * num_cells) << "\n";

		for (int i = 0; i < num_cells; ++i)
		{
			int idx = i * 4;
			vtk_file << "4 " << idx << " " << idx + 1 << " " << idx + 2 << " " << idx + 3 << "\n";
		}

		vtk_file << "\nCELL_TYPES " << num_cells << "\n";

		for (int i = 0; i < num_cells; ++i)
			vtk_file << VTK_QUAD_TYPE << "\n";

		vtk_file << "\nCELL_DATA " << num_cells << "\n"
				 << "SCALARS cell_width float 1\n"
				 << "LOOKUP_TABLE default\n";

		for (int i = 0; i < num_cells; ++i)
			vtk_file << x_vol[i + 1] - x_vol[i] << "\n";

		vtk_file << "SCALARS cell_id int 1\n"
				 << "LOOKUP_TABLE default\n";

		for (int i = 0; i < num_cells; ++i)
			vtk_file << i << "\n";
	}

}

#endif /* SPACEMESH_IMPL_B6DAD095_AFD5_4569_BA1B_F73464BEC0AA */
