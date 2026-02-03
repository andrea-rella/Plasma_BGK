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
    /** @brief Configuration data for the simulation.
     *
     * It holds all the user defined parameters for the simulation, such as mesh parameters and
     * physical parameters. It offers getters to access these parameters and the possibility to
     * load them from a configuration file (.json)
     *
     * @tparam T Data type for configuration parameters / computational precision.
     *
     */
    template <typename T>
    class ConfigData
    {
    private:
        // ------ MESH PARAMETERS -----------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /// Space mesh total number of points
        size_t N;
        /// Space mesh number of polynomially spaced points
        size_t N0;
        /// Velocity mesh total number of points
        size_t barN;
        /// Space mesh spacing parameter
        T d1;
        /// Space mesh spacing parameter
        T d2;
        /// Velocity mesh spacing parameter
        T a1;
        /// Velocity mesh spacing parameter
        T a2;

        // ------ PHYSICAL PARAMETERS -------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /// Ratio of infinity to wall temperature
        T T_infty_w;
        /// Ratio of infinity to wall pressure
        T p_infty_w;
        /// Mach number at infinity
        T M_infty;

        // ------ SIMULATION PARAMETERS -----------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /// Time stepping size
        T dt;
        /// Convergence relative tolerance
        T tol;
        /// Maximum number of iterations
        size_t max_iter;
        /// Plot every k steps
        size_t save_every_k_steps;

        // ------ SAVING PARAMETERS ---------------------------------------------------------------------
        // ----------------------------------------------------------------------------------------------

        /// Folder name to save results
        std::string saving_folder_name;

    public:
        // ------ CONSTRUCTORS AND DESTRUCTORS -----------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /// @brief Default constructor
        ConfigData() = default;
        /** @brief Construct a new Config Data object from a configuration file. .json configuration file
         *
         * @param filename path to the configuration file
         */
        ConfigData(std::string filename);
        /// @brief Default destructor
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
        size_t get_save_every_k_steps() const { return save_every_k_steps; }
        std::string get_saving_folder_name() const { return saving_folder_name; }
    };

}

#include "impl/ConfigData.tpp"

#endif /* CONFIGDATA_F996EE8C_E537_4E04_90A6_8289EAEBDA9C */
