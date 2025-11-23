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
        T T_infty_w;
        T p_infty_w;
        T M_infty;

        // Simulation parameters
        T dt;
        T tol;
        size_t max_iter;
        size_t plot_every_k_steps;

        // General
        std::string saving_folder_name;

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
        T get_T_infty_w() const { return T_infty_w; }
        T get_p_infty_w() const { return p_infty_w; }
        T get_M_infty() const { return M_infty; }
        T get_dt() const { return dt; }
        T get_tol() const { return tol; }
        size_t get_max_iter() const { return max_iter; }
        size_t get_plot_every_k_steps() const { return plot_every_k_steps; }
        std::string get_saving_folder_name() const { return saving_folder_name; }
    };

}

#include "impl/ConfigData.tpp"

#endif /* CONFIGDATA_F996EE8C_E537_4E04_90A6_8289EAEBDA9C */
