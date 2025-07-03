//__________.__                                   __________  ________ ____  __.
//\______   \  | _____    ______ _____ _____      \______   \/  _____/|    |/ _|
// |     ___/  | \__  \  /  ___//     \\__  \      |    |  _/   \  ___|      <
// |    |   |  |__/ __ \_\___ \|  Y Y  \/ __ \_    |    |   \    \_\  \    |  \
// |____|   |____(____  /____  >__|_|  (____  /____|______  /\______  /____|__\
//                    \/     \/      \/     \/_____/      \/        \/        \/
//
// Andrea Rella
// Politecnico di Milano
// https://github.com/andrea-rella/Plasma_BGK

#ifndef SPACEMESH_A76B3E2F_9D34_41BA_B9B7_4464981210C3
#define SPACEMESH_A76B3E2F_9D34_41BA_B9B7_4464981210C3

#include <vector>
#include "ConfigData.hpp"

namespace Bgk
{
    template <typename T>
    class SpaceMesh
    {
    private:
        bool is_initialized = false;
        bool is_constructed = false;

        int N;
        int N0;
        T d1;
        T d2;
        std::vector<T> x_comp;
        std::vector<T> x_vol;

    public:
        // Constructors and destructors
        SpaceMesh() = default;
        SpaceMesh(const ConfigData<T> &);
        ~SpaceMesh() = default;
    };
}

#include "impl/SpaceMesh_impl.hpp"

#endif /* SPACEMESH_A76B3E2F_9D34_41BA_B9B7_4464981210C3 */
