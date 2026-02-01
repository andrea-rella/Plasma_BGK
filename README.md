# `Plasma_BGK`: A BGK kinetic model for the study of gas - condensed-phase interaction

**Plasma BGK** is a 1D Finite Volume solver for the Bhatnagar-Gross-Krook (BGK) kinetic equation. Developed at **Politecnico di Milano**, this project focuses on simulating rarefied gas dynamics and plasma flows, specifically targeting problems involving evaporation, condensation, and boundary layer discontinuities.

The solver utilizes a velocity grid approach and is configurable via JSON files, with a dedicated Python suite for post-processing and visualization.

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
- **Eigen3**: A template library for linear algebra.
- **nlohmann/json**: A library for parsing JSON configuration files.

To run the post-processing scripts, you will need **Python 3.x** and the following packages:

- `numpy`
- `matplotlib`
- `scipy` (optional, depending on specific analysis scripts)

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
