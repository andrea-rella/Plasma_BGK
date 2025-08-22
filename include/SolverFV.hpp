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

#ifndef SOLVERFV_A020EBB4_6BD0_412F_8B28_FB5C8A3B4D66
#define SOLVERFV_A020EBB4_6BD0_412F_8B28_FB5C8A3B4D66

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
    /**
     * @brief Finite Volume Solver For the BGK Boltzmann Equation
     *
     * @tparam T precision of the computations
     * @param Data configuration data @see ConfigData
     * @param Space_mesh spatial finite volume mesh @see SpaceMeshFV
     * @param Velocity_mesh velocity mesh @see VelocityMesh
     * @param g, h solution matrices in the phase space @f$ (x,z) @f$ (space, velocity)
     * @param A, B, R numerical Eigen matrices used to solve the system
     * @param SolVector vector used to store the solution of the numerical systems
     * @param g0, g_infty, h0, h_infty boundary conditions
     * @param g_init, h_init initial conditions
     */
    template <typename T>
    class SolverFV
    {
        bool is_initialized = false;

        ConfigData<T> Data;
        SpaceMeshFV<T> Space_mesh;
        VelocityMesh<T> Velocity_mesh;

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> g;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h;
        Eigen::Vector<T, Eigen::Dynamic> density;
        Eigen::Vector<T, Eigen::Dynamic> mean_velocity;
        Eigen::Vector<T, Eigen::Dynamic> temperature;

        Eigen::SparseMatrix<T> A;
        Eigen::SparseMatrix<T> B;
        Eigen::Vector<T, Eigen::Dynamic> R;

        Eigen::Matrix<T, -1, 1> SolVector;

        std::function<T(T z)> g0, g_infty;
        std::function<T(T z)> h0, h_infty;

        std::function<T(T z)> g_init, h_init;

        std::function<T(T z, T rho, T v, T Temp)> G, H;

        // ------ SOLVE HELPERS --------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        Eigen::Vector<T, Eigen::Dynamic> assemble_U_pos(size_t j, T a1_2, T a2_2, T omega1) const;
        Eigen::Vector<T, Eigen::Dynamic> assemble_W_pos(size_t j, T a1_2, T a2_2, T omega1) const;
        Eigen::Vector<T, Eigen::Dynamic> assemble_U_zero() const;
        Eigen::Vector<T, Eigen::Dynamic> assemble_W_zero() const;
        Eigen::Vector<T, Eigen::Dynamic> assemble_U_neg(size_t j, T bN1_2, T bN2_2, T sigmaN1) const;
        Eigen::Vector<T, Eigen::Dynamic> assemble_W_neg(size_t j, T bN1_2, T bN2_2, T sigmaN1) const;

        // ------ SOLVE  ---------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        void solve_timestep_pos();
        void solve_timestep_neg();
        void solve_timestep_zero();

    public:
        // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        SolverFV() = default;
        SolverFV(const ConfigData<T> &InputData);
        SolverFV(const std::string &config_file_path);
        ~SolverFV() = default;

        // ------ BUILD NUMERICAL MATRICES / SETUP -------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        void initializeMeshes();
        void set_physical_quantities();
        void setInitialState();

        /**
         * @brief Assemble the numerical advection matrix @f$ A / \zeta^{(j)} @f$ for the case velocity > 0
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

        /**
         * @brief Assemble the numerical advection matrix @f$ B / \zeta^{(j)} @f$ for the case velocity < 0
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
        void assemble_R();
        void initialize();
        void reset();

        // ------ SOLVE  ---------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        void solve(const metrics::VectorNormType vec_norm_type,
                   const metrics::RowAggregateType agg_type);
        // solve
        // getters in general
        // setters in general
        // compute Tinf and all the other needed quantities

        // ------ OUTPUT ---------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        void write_sol_txt(const std::string &folder_name) const;
        void write_phys_txt(const std::string &folder_name) const;
        void write_meshes_txt(const std::string &folder_name) const;
        void write_space_mesh_vtk(const std::string &folder_name) const;
        void write_all(const std::string &folder_name) const;
    };
}

#include "impl/SolverFV.tpp"

#endif /* SOLVERFV_A020EBB4_6BD0_412F_8B28_FB5C8A3B4D66 */
