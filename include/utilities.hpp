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

/**
 * @brief Storage of utility functions and concepts for the Plasma_BGK project.
 *
 * In this file, we define various utility functions and concepts that are used throughout the Plasma_BGK project.
 * These utilities include
 *
 * - @see MeshContainer concept to check if a type is a range and holds mesh components.
 * - @see FloatingPoint concept to check if a type is a floating point type.
 * - @see SpacingFunction concept to define a callable type for mesh spacing.
 * - @see MeshNature enum class to define the nature of a mesh (space or velocity).
 * - @see QUICK_coefficients_p Function to compute QUICK coefficients for finite volume meshes.
 * - @see QUICK_coefficients_n Function to compute QUICK coefficients for finite volume meshes.
 *
 * These utilities are essential for ensuring type safety and providing common functionality across different mesh implementations.
 */

#ifndef UTILITIES_F9FED3EC_E850_404A_8753_23CF1EDC3674
#define UTILITIES_F9FED3EC_E850_404A_8753_23CF1EDC3674

#include <concepts>
#include <ranges>
#include <string>
#include <source_location>
#include <Eigen/Sparse>

/**
 * @brief Throws a std::runtime_error in debug builds if the condition is false.
 *        In release builds, logs a warning message instead.
 *
 * This macro checks the expression `expr`.
 * - If `expr` evaluates to false in debug mode (`NDEBUG` **not** defined), it throws an exception with the provided message `msg`.
 * - If in release mode (`NDEBUG` defined), it prints a warning to `std::cerr` but does not throw.
 *
 * @param expr The condition to evaluate.
 * @param msg  The error message used for the exception or warning.
 *
 * @note Use this macro to enforce debug-only sanity checks without crashing release builds.
 *
 * @example
 * DEBUG_THROW_MSG(x >= 0, "x must be non-negative");
 */
#ifndef NDEBUG
#define DEBUG_THROW_MSG(expr, msg)         \
    do                                     \
    {                                      \
        if (!(expr))                       \
            throw std::runtime_error(msg); \
    } while (0)
#else
#define DEBUG_THROW_MSG(expr, msg)                                                    \
    do                                                                                \
    {                                                                                 \
        if (!(expr))                                                                  \
        {                                                                             \
            std::cerr << "Warning: skipped debug check failed: " << msg << std::endl; \
            /* optionally: log error, notify, etc. */                                 \
        }                                                                             \
    } while (0)
#endif

/**
 * @brief Wrapper for more informative error messages
 *
 * @param message error message
 * @param loc source code location of the error. Defaulted to std::source_location::current()
 * @return std::string with the complete error message
 */
std::string error_message(const std::string &message, const std::source_location &loc = std::source_location::current());

namespace Bgk
{

    // Forward declaration of SpaceMeshFV to avoid circular dependencies.
    template <typename T>
    class SpaceMeshFV;

    /**
     * @brief Concept for a mesh container.
     *
     * This concept checks if a type is a range and if its value type matches the specified type T.
     *
     * @tparam Container container type that should hold mesh components.
     * @tparam T precision type of the mesh components.
     */
    template <typename Container, typename T>
    concept MeshContainer = std::ranges::range<Container> &&
                            std::same_as<std::ranges::range_value_t<Container>, T>;

    /**
     * @brief Concept for floating point types.
     *
     * @tparam T type to check.
     */
    template <typename T>
    concept FloatingPoint = std::floating_point<T>;

    /**
     * @brief Concept for a spacing function.
     *
     * This concept checks if a callable type can be invoked with an integer index and returns a value
     * convertible to T.
     *
     * @tparam SpacingFunc callable type that defines the spacing.
     * @tparam T precision type of the mesh components.
     */
    template <typename SpacingFunc, typename T>
    concept SpacingFunction = requires(SpacingFunc func, size_t index) {
        { func(index) } -> std::convertible_to<T>;
    };

    /**
     * @brief Possible nature of a mesh.
     *
     * This enum class defines the nature of a mesh, which can either be in space or in velocity.
     * It's used to define the context in which the mesh is being used to then apply the appropriate
     * versions of the methods.
     *
     */
    enum class MeshNature
    {
        /// Space mesh identifier.
        SPACE,
        /// Velocity mesh identifier.
        VELOCITY
    };

    /**
     * @brief Concept for a function that takes a single argument and returns a value of type T.
     *
     * This concept checks if a callable type can be invoked with an argument of type T and returns a value
     * convertible to T.
     *
     * @tparam Func callable type to check.
     * @tparam T type of the argument and return value.
     */
    template <typename Func, typename T>
    concept SingleArgFunction = requires(Func func, T arg) {
        { func(arg) } -> std::convertible_to<T>;
    };

}

#include "impl/utilities.tpp"

#endif /* UTILITIES_F9FED3EC_E850_404A_8753_23CF1EDC3674 */
