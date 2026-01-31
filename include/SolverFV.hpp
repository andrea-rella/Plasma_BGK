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

#ifndef SOLVERFV_CF3739E4_0179_477E_9C04_0D0ED2512066
#define SOLVERFV_CF3739E4_0179_477E_9C04_0D0ED2512066

#include "utilities.hpp"
#include "phys_utils.hpp"
#include "ConfigData.hpp"
#include "SpaceMeshFV.hpp"
#include "VelocityMesh.hpp"
#include "metrics_utils.hpp"

#include "Eigen/Sparse"
#include "math.h"
#include "functional"

namespace Bgk
{
    /** @brief Finite Volume Solver For the BGK Boltzmann Equation
     *
     * This class implements a Finite Volume solver for the BGK Boltzmann equation, it manages the meshes,
     * physical quantities, numerical matrices assembly, and the solution process over time. It also provides
     * methods for outputting results in various formats.
     *
     * @tparam T precision of the computationss
     */
    template <typename T>
    class SolverFV
    {
    private:
        /// Flag indicating whether the solver has been initialized (i.e., meshes and matrices are set up).
        bool is_initialized = false;

        /// Configuration data object containing simulation parameters.
        ConfigData<T> Data;
        /// Space mesh object for spatial discretization.
        SpaceMeshFV<T> Space_mesh;
        /// Velocity mesh object for velocity discretization.
        VelocityMesh<T> Velocity_mesh;

        // ------ PHYSICAL QUANTITIES -------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// Grid of g distribution function values
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> g;
        /// Grid of h distribution function values
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h;

        /// Macroscopic normalized density vector
        Eigen::Vector<T, Eigen::Dynamic> density;
        /// Macroscopic normalized mean velocity vector
        Eigen::Vector<T, Eigen::Dynamic> mean_velocity;
        /// Macroscopic normalized temperature vector
        Eigen::Vector<T, Eigen::Dynamic> temperature;

        // ------ NUMERICAL MATRICES --------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// Positive Velocity Advection Matrix
        Eigen::SparseMatrix<T> A;
        /// Negative Velocity Advection Matrix
        Eigen::SparseMatrix<T> B;
        /// Right-hand side Matrix
        Eigen::Vector<T, Eigen::Dynamic> R;

        // ------ BOUNDARY AND INITIAL CONDITIONS --------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// Boundary condition functions for g at the left and right boundaries
        std::function<T(T z)> g0, g_infty;
        /// Boundary condition functions for h at the left and right boundaries
        std::function<T(T z)> h0, h_infty;
        /// Initial condition functions for g and h
        std::function<T(T z)> g_init, h_init;

        /// ------ COLLISION TERMS -----------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// Collision term functions G and H
        std::function<T(T z, T rho, T v, T Temp)> G, H;

        // ------ SOLVE HELPERS --------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Assembles the reaction vector for the g equation for positive velocities
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @param j velocity index
         * @param a1_2 QUICK coefficient @f$ \alpha^{1}_2 @f$
         * @param a2_2 QUICK coefficient @f$ \alpha^{2}_2 @f$
         * @param omega1 QUICK coefficient @f$ \omega^{1} @f$
         * @return Eigen::Vector<T, Eigen::Dynamic> reaction vector
         */
        Eigen::Vector<T, Eigen::Dynamic> assemble_U_pos(size_t j, T a1_2, T a2_2, T omega1) const;
        /** @brief Assembles the reaction vector for the h equation for positive velocities
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @param j velocity index
         * @param a1_2 QUICK coefficient @f$ \alpha^{1}_2 @f$
         * @param a2_2 QUICK coefficient @f$ \alpha^{2}_2 @f$
         * @param omega1 QUICK coefficient @f$ \omega^{1} @f$
         * @return Eigen::Vector<T, Eigen::Dynamic> reaction vector
         */
        Eigen::Vector<T, Eigen::Dynamic> assemble_W_pos(size_t j, T a1_2, T a2_2, T omega1) const;
        /** @brief Assembles the reaction vector for the g equation for zero velocity
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @return Eigen::Vector<T, Eigen::Dynamic> reaction vector
         */
        Eigen::Vector<T, Eigen::Dynamic> assemble_U_zero() const;
        /** @brief Assembles the reaction vector for the h equation for zero velocity
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @return Eigen::Vector<T, Eigen::Dynamic> reaction vector
         */
        Eigen::Vector<T, Eigen::Dynamic> assemble_W_zero() const;
        /** @brief Assembles the reaction vector for the g equation for negative velocities
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @param j velocity index
         * @param bN1_2 QUICK coefficient @f$ \beta^{N-1}_2 @f$
         * @param bN2_2 QUICK coefficient @f$ \beta^{N-2}_2 @f$
         * @param sigmaN1 QUICK coefficient @f$ \sigma^{N-1} @f$
         * @return Eigen::Vector<T, Eigen::Dynamic> reaction vector
         */
        Eigen::Vector<T, Eigen::Dynamic> assemble_U_neg(size_t j, T bN1_2, T bN2_2, T sigmaN1) const;
        /** @brief Assembles the reaction vector for the h equation for negative velocities
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @param j velocity index
         * @param bN1_2 QUICK coefficient @f$ \beta^{N-1}_2 @f$
         * @param bN2_2 QUICK coefficient @f$ \beta^{N-1}_2 @f$
         * @param sigmaN1 QUICK coefficient @f$ \sigma^{N-1} @f$
         * @return Eigen::Vector<T, Eigen::Dynamic> reaction vector
         */
        Eigen::Vector<T, Eigen::Dynamic> assemble_W_neg(size_t j, T bN1_2, T bN2_2, T sigmaN1) const;

        // ------ SOLVE  ---------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Solve a single time step for positive velocities
         *
         * This function solves the discretized BGK Boltzmann equations of g and h for positive velocities
         * for a single time step. It uses the pre-assembled QUICK advection matrix A, reaction matrix R and
         * the reaction vectors following a semi-implicit time integration scheme.
         *
         * For more details on the numerical scheme and implementation, refer to the library report.
         *
         */
        void solve_timestep_pos();
        /** @brief Solve a single time step for negative velocities
         *
         * This function solves the discretized BGK Boltzmann equations of g and h for negative velocities
         * for a single time step. It uses the pre-assembled QUICK advection matrix B, reaction matrix R
         * and the reaction vectors following a semi-implicit time integration scheme.
         *
         * For more details on the numerical scheme and implementation, refer to the library report.
         *
         */
        void solve_timestep_neg();
        /** @brief Solve a single time step for zero velocity
         *
         * This function solves the discretized BGK Boltzmann equations of g and h for zero velocity
         * for a single time step. It uses the reaction matrix R and the reaction vectors
         * following a semi-implicit time integration scheme.
         *
         * For more details on the numerical scheme and implementation, refer to the library report.
         *
         */
        void solve_timestep_zero();

    public:
        // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Default constructor
        SolverFV() = default;
        /** @brief Construct a new Solver FV object through a ConfigData object
         *
         * This constructor initializes the SolverFV instance using a provided ConfigData object. It sets up
         * the simulation parameters, meshes and boundary conditions based on configuration data.
         *
         * @note It DOES NOT initalize the meshes neither the numerical matrices. The user must call
         * the global initialize() method or the individual methods initializeMeshes(), set_physical_quantities(),
         * setInitialState() and the various matrix assemble methods to complete the setup before solving.
         *
         *
         * @param InputData configuration data object
         */
        SolverFV(const ConfigData<T> &InputData);
        /** @brief Construct a new Solver FV object from a configuration file. .json configuration file
         *
         * This constructor initializes the SolverFV instance by reading simulation parameters,
         * meshes and boundary conditions from a specified JSON configuration file.
         *
         * @note It DOES NOT initalize the meshes neither the numerical matrices. The user must call
         * the global initialize() method or the individual methods initializeMeshes(), set_physical_quantities(),
         * setInitialState() and the various matrix assemble methods to complete the setup before solving.
         *
         * @param config_file_path path to the configuration file
         */
        SolverFV(const std::string &config_file_path);
        /// @brief Default destructor
        ~SolverFV() = default;

        // ------ INITIALIZE METHODS --------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Initializes the space and velocity meshes based on configuration data
         *
         * This method sets up the spatial and velocity meshes using the parameters specified
         * in the ConfigData object. It configures the mesh sizes, spacing, and distribution
         * according to the simulation requirements.
         *
         * @note It calls the SpaceMeshFV and VelocityMesh initialization methods.
         *
         */
        void initializeMeshes();
        /** @brief Sets physical quantities and boundary conditions based on configuration data
         *
         * This method initializes the physical quantities such as density, mean velocity,
         * and temperature.
         *
         * @note It calls the phys_utils functions
         *
         */
        void set_physical_quantities();
        /** @brief Sets the initial state of the distribution functions g and h
         *
         * This method initializes the distribution function grids g and h using the
         * provided initial and boundary condition functions. It populates the g and h matrices
         * based on the spatial mesh and the initial conditions.
         *
         */
        void setInitialState();

        // ------ BUILD NUMERICAL MATRICES / SETUP -------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Assemble the numerical advection matrix @f$ A / \zeta^{(j)} @f$ for the case velocity > 0
         *
         * This function assembles the numerical advection matrix A for the case where the velocity is positive.
         * It computes the QUICK coefficients from Space_mesh and populates the matrix according
         * to the finite volume discretization. (Refer to the library report for further details).
         *
         * Such matrix will be used to solve for the points from 1 to N given that the solution in 0 is given by
         * the boundary condition.
         *
         * @see Bgk::QUICK_coefficients_p
         *
         */
        void assemble_A();
        /** @brief Assemble the numerical advection matrix @f$ B / \zeta^{(j)} @f$ for the case velocity < 0
         *
         * This function assembles the numerical advection matrix A for the case where the velocity is positive.
         * It computes the QUICK coefficients from Space_mesh and populates the matrix according
         * to the finite volume discretization. (Refer to the library report for further details).
         *
         * Such matrix will be used to solve for the points from 0 to N-1 given that the solution in N is given by
         * the boundary condition.
         *
         * @see Bgk::QUICK_coefficients_p
         *
         */
        void assemble_B();
        /** @brief Assembles the reaction diagonal matrix R
         *
         * This function assembles the reaction diagonal matrix R used in the semi-implicit time integration
         * scheme for solving the BGK Boltzmann equations.
         *
         * For more details on the assembly process, refer to the library report.
         *
         * @note since the matrix is diagonal, it is stored as a vector for efficiency.
         *
         */
        void assemble_R();

        // ------ GLOBAL INITIALIZE AND RESET ------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Global initialization method that sets up meshes, physical quantities, initial state and numerical matrices.
        void initialize();
        /// @brief Resets the solver to its initial uninitialized state.
        void reset();

        // ------ SOLVE  ---------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Solves the BGK Boltzmann equations over time until convergence or maximum iterations
         *
         * This method iteratively solves the discretized BGK Boltzmann equations for the distribution
         * functions g and h over time. It continues the simulation until either convergence is achieved
         * based on the specified vector norm and aggregation type, or the maximum number of iterations
         * is reached.
         *
         * @param vec_norm_type Type of vector norm to use for convergence checking (e.g., L2, L∞).
         * @param agg_type Type of aggregation method to use for convergence checking (e.g., max, average).
         */
        void solve(const metrics::VectorNormType vec_norm_type,
                   const metrics::RowAggregateType agg_type);

        // ------ OUTPUT ---------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Writes the solution distribution functions g and h to text files
         *
         * This method outputs the current state of the distribution functions g and h
         * to text files in the specified folder path.
         *
         * @param folder_path Path to the folder where the text files will be saved. (e.g., "./results")
         *
         *
         */
        void write_sol_txt(const std::string &folder_path) const;
        /** @brief Writes the instantaneous distribution functions g and h at a specific iteration to text files
         *
         * Outputs the full g and h matrices snapshot at a given iteration
         * into two files in the specified folder path.
         *
         * @param folder_path Path to the folder where the text files will be saved.
         * @param iter Current iteration number
         */
        void write_sol_instant_txt(const std::string &folder_path, size_t iter) const;
        /** @brief Writes the physical quantities to text files
         *
         * This method outputs the current state of the physical quantities
         * to text files in the specified folder path.
         *
         * @param folder_path Path to the folder where the text files will be saved. (e.g., "./results")
         *
         */
        void write_phys_txt(const std::string &folder_path) const;
        /** @brief Writes the space mesh to a text file
         *
         * This method outputs the spatial mesh data
         * to a text file in the specified folder path.
         *
         * @param folder_path Path to the folder where the text file will be saved. (e.g., "./results")
         *
         * @note calls SpaceMeshFV::write_mesh_txt() and VelocityMesh::write_mesh_txt()
         */
        void write_meshes_txt(const std::string &folder_path) const;
        /** @brief Writes the space mesh to a VTK file
         *
         * This method outputs the spatial mesh data
         * to a VTK file in the specified folder path.
         *
         * @param folder_path Path to the folder where the VTK file will be saved. (e.g., "./results")
         *
         * @note calls SpaceMeshFV::write_mesh_vtk()
         */
        void write_space_mesh_vtk(const std::string &folder_path) const;
        /** @brief Writes the initial state of the simulation to text files
         *
         * This method outputs the initial state of the distribution functions g and h,
         * as well as the physical quantities, to text files in the specified folder path.
         *
         * @param folder_path Path to the folder where the text files will be saved. (e.g., "./results")
         *
         */
        void write_initial_state_txt(const std::string &folder_path) const;
        /** @brief Writes the instantaneous physical quantities at a specific iteration to text files
         *
         * This method outputs the physical quantities at a given iteration
         * to text files in the specified folder path.
         *
         * @param folder_path Path to the folder where the text files will be saved. (e.g., "./results")
         * @param iter Current iteration number
         *
         */
        void write_phys_instant(const std::string &folder_path, size_t iter) const;
        /** @brief Writes all relevant simulation data to files
         *
         * This method outputs the distribution functions, physical quantities,
         * and mesh data to files in the specified folder path.
         *
         * @param folder_path Path to the folder where the files will be saved. (e.g., "./results")
         *
         */
        void write_all(const std::string &folder_path) const;
    };
}

#include "impl/SolverFV.tpp"

#endif /* SOLVERFV_CF3739E4_0179_477E_9C04_0D0ED2512066 */
