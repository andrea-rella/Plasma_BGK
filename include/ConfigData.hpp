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
     * @param Z Length of the simulation space.
     * @param D Velocity interval for the simulation.
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
        T Z;
        T D;
        int N;
        int N0;
        int barN;
        T d1;
        T d2;
        T a1;
        T a2;

        // Physical parameters
        T v_infty;

    public:
        // Constructors and destructors
        ConfigData() = default;
        ConfigData(std::string filename);
        ~ConfigData() = default;

        // Getters
        T get_Z() const { return Z; }
        T get_D() const { return D; }
        int get_N() const { return N; }
        int get_N0() const { return N0; }
        int get_BarN() const { return barN; }
        T get_d1() const { return d1; }
        T get_d2() const { return d2; }
        T get_a1() const { return a1; }
        T get_a2() const { return a2; }
        T get_v_infty() const { return v_infty; }
    };

}

#include "impl/ConfigData_impl.hpp"

#endif /* CONFIGDATA_F996EE8C_E537_4E04_90A6_8289EAEBDA9C */
