#ifndef CONFIGDATA_B50FF598_5E29_4C64_B14A_21660E232A71
#define CONFIGDATA_B50FF598_5E29_4C64_B14A_21660E232A71

#include <string>

namespace Bgk
{

    template <typename T>
    class ConfigData
    {
    private:
        T Z;
        T D;

    public:
        ConfigData() = default;
        ConfigData(std::string filename);
        ~ConfigData() = default;
    };

}

#include "impl/ConfigData_impl.hpp"

#endif /* CONFIGDATA_B50FF598_5E29_4C64_B14A_21660E232A71 */
