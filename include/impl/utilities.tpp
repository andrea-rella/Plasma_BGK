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

#ifndef UTILITIES_IMPL_BF745DF9_95F2_42D6_BF01_CA4D3185119B
#define UTILITIES_IMPL_BF745DF9_95F2_42D6_BF01_CA4D3185119B

#include "../utilities.hpp"
#include "../SpaceMeshFV.hpp"

namespace Bgk
{

     std::string error_message(const std::string &message, const std::source_location &loc)
     {
          std::string error_msg;
          error_msg += "Error at " + std::string(loc.file_name()) + ":" + std::to_string(loc.line()) +
                       " in " + loc.function_name() + ": " + message + "\n";

          return error_msg;
     };

     template <typename T>
     std::vector<std::pair<T, T>> QUICK_coefficients_p(const SpaceMeshFV<T> &mesh)
     {
          if (!mesh.validate_mesh())
               throw std::invalid_argument(
                   "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh.");

          std::vector<std::pair<T, T>> coefficients;
          std::vector<T> vol_boundaries = mesh.get_volume_boundaries();
          std::vector<T> x_comp = mesh.get_XComp();

          size_t V = vol_boundaries.size();
          T a1, a2;

          coefficients.reserve(V);

          // Case i = 0
          coefficients.emplace_back(T{0}, T{0});

          // Case i = 1: Account for the symmetric ghost node x^{(-1)} (x_comp[-1])
          a1 = ((vol_boundaries[1] - x_comp[0]) * (vol_boundaries[1] + x_comp[1])) /
               ((x_comp[1] - x_comp[0]) * (x_comp[1] + x_comp[1]));

          a2 = ((vol_boundaries[1] - x_comp[0]) * (x_comp[1] - vol_boundaries[1])) /
               ((x_comp[0] + x_comp[1]) * (x_comp[1] + x_comp[1]));

          coefficients.emplace_back(a1, a2);

          for (size_t i = 2; i < V - 1; ++i)
          {
               a1 = ((vol_boundaries[i] - x_comp[i - 1]) * (vol_boundaries[i] - x_comp[i - 2])) /
                    ((x_comp[i] - x_comp[i - 1]) * (x_comp[i] - x_comp[i - 2]));

               a2 = ((vol_boundaries[i] - x_comp[i - 1]) * (x_comp[i] - vol_boundaries[i])) /
                    ((x_comp[i - 1] - x_comp[i - 2]) * (x_comp[i] - x_comp[i - 2]));

               coefficients.emplace_back(a1, a2);
          }

          // Case i = V - 1
          coefficients.emplace_back(T{0}, T{0});

          return coefficients;
     }

     //--------------------------------------------------------------------------

     template <typename T>
     std::vector<std::pair<T, T>> QUICK_coefficients_n(const SpaceMeshFV<T> &mesh)
     {

          if (!mesh.validate_mesh())
               throw std::invalid_argument(
                   "Invalid mesh: trying to compute QUICK coefficients for an uninitialized or unconstructed mesh.");

          std::vector<std::pair<T, T>> coefficients;
          std::vector<T> vol_boundaries = mesh.get_volume_boundaries();
          std::vector<T> x_comp = mesh.get_XComp();

          size_t V = vol_boundaries.size();
          T b1, b2;

          coefficients.reserve(V);

          // Case i = 0:
          coefficients.emplace_back(T{0}, T{0});

          // General case:
          for (size_t i = 1; i < V - 2; ++i)
          {
               b1 = ((vol_boundaries[i] - x_comp[i]) * (vol_boundaries[i] - x_comp[i + 1])) /
                    ((x_comp[i - 1] - x_comp[i]) * (x_comp[i - 1] - x_comp[i + 1]));

               b2 = ((vol_boundaries[i] - x_comp[i]) * (x_comp[i - 1] - vol_boundaries[i])) /
                    ((x_comp[i] - x_comp[i + 1]) * (x_comp[i - 1] - x_comp[i + 1]));

               coefficients.emplace_back(b1, b2);
          }

          // Case i = V - 2: Account for the symmetric ghost node x^{(N+1)} (x_comp[V - 1])
          b1 = ((vol_boundaries[V - 2] - x_comp[V - 2]) * (vol_boundaries[V - 2] - (2. * x_comp[V - 2] - x_comp[V - 3]))) /
               ((x_comp[V - 3] - x_comp[V - 2]) * (x_comp[V - 3] - (2. * x_comp[V - 2] - x_comp[V - 3])));

          b2 = ((vol_boundaries[V - 2] - x_comp[V - 2]) * (x_comp[V - 3] - vol_boundaries[V - 2])) /
               ((x_comp[V - 2] - (2. * x_comp[V - 2] - x_comp[V - 3])) * (x_comp[V - 3] - (2. * x_comp[V - 2] - x_comp[V - 3])));

          coefficients.emplace_back(b1, b2);

          // Case i = V - 1:
          coefficients.emplace_back(T{0}, T{0});

          return coefficients;
     }

}

#endif /* UTILITIES_IMPL_BF745DF9_95F2_42D6_BF01_CA4D3185119B */
