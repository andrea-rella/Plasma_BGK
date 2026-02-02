# `Plasma_BGK`: A BGK Kinetic Model for the study of Gas - Condensed-Phase Interaction

`Plasma_BGK` is a 1D Finite Volume solver for the Bhatnagar-Gross-Krook (BGK) kinetic equation. Developed at **Politecnico di Milano**, this project focuses on simulating rarefied gas dynamics, specifically targeting problems involving evaporation, condensation, and boundary layer discontinuities. The problem setting is the one described as:

<p align="center">
  <img src="docs/images/problem_setting.png" alt="Problem Setting" width="700"/>
</p>

$$
\frac{\partial f}{\partial t} + \xi_1 \frac{\partial f}{\partial X_1} = A_c \rho (f^{(0)} - f)
$$

$$
f^{(0)} = \frac{\rho}{(2\pi R T)^{3/2}} \exp \left( -\frac{(\xi_1 - v_1)^2 + \xi_2^2 + \xi_3^2}{2RT} \right)
$$

Intially adressed with a finite difference method in [[1](#ref-1)] and [[2](#ref-2)]; here we propose a finite volume (Q.U.I.C.K.) approach integrated with the fast `Eigen` library. For a full theoretical description and implementation details please refer to the project report in the [`docs/`](docs/) folder.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Building the Project](#building-the-project)
- [Main Features](#main-features)
- [Usage & Examples](#usage--examples)
- [Project Structure](#project-structure)
- [Authors](#authors)

## Requirements

To build and run the C++ solver, you will need:

- **C++ Compiler**: A compiler supporting **C++20**.
- **Make**: GNU Make build system.
- [**Eigen3**](https://libeigen.gitlab.io): A template library for linear algebra.
- [**nlohmann/json**](https://github.com/nlohmann/json): A library for parsing JSON configuration files.
- [**OpenMP**](https://www.openmp.org): For parallel solver execution. (_Optional_)

To run the post-processing scripts, you will need **Python 3.x** and the following packages:

- `numpy`
- `matplotlib`
- `scipy` (optional, depending on specific analysis scripts)

If needed a virtual enviroment called `Bgk-venv` with all the necessary packages is provided in the projects's [Google Drive Folder](https://drive.google.com/drive/folders/14fjwMq2cQmRHNADEw7G7EoBn30CHxR80?usp=sharing)

## Installation

You can download the code by cloning the repository from GitHub:

```bash
git clone git@github.com:andrea-rella/Plasma_BGK.git
cd Plasma_BGK
```

## Building the Project

The project uses a standard `Makefile`. Before building, ensure the generic include paths in the Makefile match your system's library locations (specifically for Eigen and nlohmann-json).

To compile the executables:

```bash
make
```

To clean build artifacts:

```bash
make clean
```

To clean everything:

```bash
make distclean
```

## Main Features

- **Mathematical Model**: Solves the 1D BGK kinetic equation.
- **Numerical Method**: Finite Volume Method (FVM) for spatial discretization with specific support for velocity grid meshes (`SpaceMeshFV`, `VelocityMesh`).
- **Configuration**: Full simulation control via JSON input files (time steps, grid size, physical parameters, output settings).
- **Boundary Conditions**: Supports specific boundary conditions for evaporation/condensation problems and discontinuities.
- **Post-Processing**: Includes Python scripts (`postprocessing/`) for plotting density profiles, temperature, velocity, and comparing results with theoretical models (e.g., Ytrehus model).
- **Documentation**: Integrated Doxygen documentation support.

## Usage & Examples

### Running the Solver

The main executable requires a path to a valid JSON configuration file.

**Basic Syntax:**

```bash
./main.exe <path_to_config_file>
```

**Example:**
To run a simulation using the data in `data/speed.json`:

```bash
./main.exe data/speed.json
```

### Configuration (JSON)

Input files (located in `data/`) control the simulation. An example structure includes:

```json
{
  "physical": {
    "p_infty_w": 3,
    "M_infty": 4.0
  },
  "simulation": {
    "time_step": 0.01,
    "max_iterations": 1000,
    "plot_every_k_steps": 10
  },
  "general": {
    "saving_folder_name": "output/results"
  }
}
```

### Post-Processing

After the simulation finishes, results are stored in the folder specified in the JSON (e.g., `output/`). You can visualize them using the Python scripts:

```bash
# Example: Plotting results
python3 postprocessing/plot.py
```

## Project Structure

- 📁 [`src/`](src/): Source files (`.cpp`).
- 📁 [`include/`](include/): Header files (`.hpp`) and template implementations (`impl/`).
- 📁 [`data/`](data/): JSON configuration files for different test cases.
- 📁 [`docs/`](docs/): Doxygen documentation configuration and output.
- 📁 [`postprocessing/`](postprocessing/): Python scripts for data analysis and plotting.
- 📁 [`output/`](output/): Directory where simulation results are generated.

## Authors

- **Andrea Rella** - _Politecnico di Milano_

## References

<a id="ref-1"></a>[1] Aoki, Kazuo and Sone, Yoshio and Yamada, Tatsuo "Numerical analysis of gas flows condensing on its plane condensed phase on the basis of kinetic theory" Physics of Fluids A: Fluid Dynamics, 1990, American Institute of Physics.

<a id="ref-2"></a>[2] Sone, Yoshio and Sugimoto, Hiroshi "Strong evaporation from a plane condensed phase" Shinku, 1988, The Vacuum Society of Japan.
