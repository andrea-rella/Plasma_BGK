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

#ifndef CONFIGDATA_F996EE8C_E537_4E04_90A6_8289EAEBDA9C
#define CONFIGDATA_F996EE8C_E537_4E04_90A6_8289EAEBDA9C

#include <string>

namespace Bgk
{
    /**
     * @brief Configuration data for the simulation. It holds all the user defined parameters
     *       for the simulation, such as mesh parameters and physical parameters.
     *
     * @tparam T Data type for configuration parameters / computational precision.
     *
     * @param N Number of spatial points in the simulation.
     * @param barN Number of velocity points in the simulation.
     * @param N0 Space mesh discritization parameter
     * @param d1, d2, a1, a2 Meshes dicretization parameters. (d1, d2 refer to space, a1, a2 refer to velocity).
     * @param v_infty Asymptotic velocity for the simulation.
     *
     */
    template <typename T>
    class ConfigData
    {
    private:
        // Mesh parameters
        size_t N;
        size_t N0;
        size_t barN;
        T d1;
        T d2;
        T a1;
        T a2;

        // Physical parameters
        T v_infty;
        T T_w_infty;
        T p_w_infty;
        T M_infty;

        // Simulation parameters
        T dt;
        T tol;
        size_t max_iter;

    public:
        // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
        // -----------------------------------------------------------------------------------------------
        ConfigData() = default;
        ConfigData(std::string filename);
        ~ConfigData() = default;

        // ------ GETTERS --------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------
        size_t get_N() const { return N; }
        size_t get_N0() const { return N0; }
        size_t get_BarN() const { return barN; }
        T get_d1() const { return d1; }
        T get_d2() const { return d2; }
        T get_a1() const { return a1; }
        T get_a2() const { return a2; }
        T get_v_infty() const { return v_infty; }
        T get_T_w_infty() const { return T_w_infty; }
        T get_p_w_infty() const { return p_w_infty; }
        T get_M_infty() const { return M_infty; }
        T get_dt() const { return dt; }
        T get_tol() const { return tol; }
        size_t get_max_iter() const { return max_iter; }
    };

}

#include "impl/ConfigData.tpp"

#endif /* CONFIGDATA_F996EE8C_E537_4E04_90A6_8289EAEBDA9C */
