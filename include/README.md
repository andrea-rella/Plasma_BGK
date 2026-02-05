## `include/` Folder Overview

The `include/`folder contains both the template headers (`.hpp`) and implementation files (`.tpp`) of the core functionalities of the code. It is organized in the following way:

- ` utils`: A collection of utilities used through the whole code. Contains concept definition, enumerations and error diagnostics.

- ` ConfigData`: Class that stores the essential simulation data. Reads from a `.json`file.

- ` BaseMesh1D`: Pure virtual class defining the key structures of a general 1D mesh.

- ` SpaceMeshFV`: Specialization of ` BaseMesh1D` to implement Finite Volume logic for space.

- ` VelocityMesh`: Specialization of ` BaseMesh1D` for symmetric velocity space.

- ` numeric_utils`: Collection of functions for the computation of Finite Volume QUICK coefficients.

- ` phys_utils`: Collection of function for the computation of physical quantities $\overline{\rho}$, $\overline{v}$ and $\overline{T}$.

- ` metric_utils`: Object factory for the error computation during solve procedure.

- ` SolverFV`: Solver Class. Incorporates setup, iteration (both serial and parallel) and results output.
