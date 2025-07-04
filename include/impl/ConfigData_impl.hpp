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

#ifndef CONFIGDATA_IMPL_EB256154_4013_44BB_AFC8_449639C989E9
#define CONFIGDATA_IMPL_EB256154_4013_44BB_AFC8_449639C989E9

#include "../ConfigData.hpp"
#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

namespace Bgk
{
    template <typename T>
    ConfigData<T>::ConfigData(std::string filename)
    {
        // Parameters will be read from a .json file

        std::ifstream f(filename);
        json data = json::parse(f);

        // ---- Mesh parameters ----

        D = data["mesh"].value("velocity_interval", static_cast<T>(5.0));
        Z = data["mesh"].value("space_length", static_cast<T>(10.0));
        N = data["mesh"].value("space_points", 100);
        N0 = data["mesh"].value("N0", 90);
        barN = data["mesh"].value("velocity_points", 10);
        d1 = data["mesh"].value("d1", static_cast<T>(0.1));
        d2 = data["mesh"].value("d2", static_cast<T>(0.1));
        a1 = data["mesh"].value("a1", static_cast<T>(0.1));
        a2 = data["mesh"].value("a2", static_cast<T>(0.1));

        // ---- Physical parameters ----

        v_infty = data["physical_parameters"].value("v_infty", static_cast<T>(1.0));
    };

}

#endif /* CONFIGDATA_IMPL_EB256154_4013_44BB_AFC8_449639C989E9 */
