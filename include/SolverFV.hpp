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

#ifndef SOLVERFV_E5B713FB_A80D_4161_86AB_0F68175D69D1
#define SOLVERFV_E5B713FB_A80D_4161_86AB_0F68175D69D1

#include "utilities.hpp"
#include "ConfigData.hpp"
#include "SpaceMeshFV.hpp"
#include "VelocityMesh.hpp"

#include "Eigen/Sparse"
#include "math.h"
#include "functional"

namespace Bgk
{
    template <typename T>
    class SolverFV
    {
        bool is_initialized = false;

        ConfigData<T> Data;

        SpaceMeshFV<T> Space_mesh;

        VelocityMesh<T> Velocity_mesh;

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> g;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h;

        Eigen::SparseMatrix<T> A;
        Eigen::SparseMatrix<T> B;
        Eigen::SparseMatrix<T> R;

        Eigen::Matrix<T, -1, 1> SolVector;

        std::function<T(T)> g0, g_infty;
        std::function<T(T)> h0, h_infty;

        std::function<T(T)> g_init, h_init;

    public:
        // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        SolverFV() = default;
        SolverFV(const std::string &config_file_path);
        ~SolverFV() = default;

        // ------ BUILD MATRICES / SETUP -----------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        void build_A();
        void build_B();
        void build_R();
        void setup();
        void reset();

        // solve
        // getters in general
        // setters in general
        // compute Tinf and all the other needed quantities
    };
}

#include "impl/SolverFV.tpp"

#endif /* SOLVERFV_E5B713FB_A80D_4161_86AB_0F68175D69D1 */
