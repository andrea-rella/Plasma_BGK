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
        Z = data["mesh"]["space_length"].get<T>();
        D = data["mesh"]["velocity_interval"].get<T>();

        std::cout << "Z: " << Z << ", D: " << D << std::endl;
    };

}

#endif /* CONFIGDATA_IMPL_EB256154_4013_44BB_AFC8_449639C989E9 */
