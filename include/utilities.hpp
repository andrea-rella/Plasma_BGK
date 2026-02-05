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

#ifndef UTILITIES_D1960633_B6A3_4A73_81CA_EEFCDE02F6EE
#define UTILITIES_D1960633_B6A3_4A73_81CA_EEFCDE02F6EE

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
     * To define a mesh container this concept checks:
     *
     * 1. The container must be a standard range (provides begin(), end(), and supports range-for loops)
     *
     * 2. The container must specify the exact type of the elements (the value type)
     *
     * 3. The container must have a size() member function
     *
     * 4. The container must support operator[] for element access
     *
     * 5. The container should provide random access iterators: it guarantees operator[] is O(1) and that
     *   the container is a contiguous block of data or similar
     *
     * @tparam Container container type that should hold mesh components.
     * @tparam T precision type of the mesh components.
     */
    template <typename Container, typename T>
    concept MeshContainer1D =
        // 1.
        std::ranges::range<Container> &&

        // 2.
        std::same_as<std::ranges::range_value_t<Container>, T> &&

        // 3.
        requires(const Container &c) {
            // Calls the size() member function. We don't strictly require
            // the return type to be size_t, but it should be an integral type.
            { c.size() } -> std::integral;

            // Ensure size() is non-zero when the container is not empty,
            // although std::ranges::range_value_t<Container> is sufficient
            // to prove that the container can hold objects of type T.
        } &&

        // 4.
        requires(const Container &c, std::size_t n) {
            // Calls operator[] with an index.
            // The return type should be convertible or the same as T (const reference for const container).
            { c[n] } -> std::convertible_to<const T &>;
        } &&

        // 5.
        std::random_access_iterator<std::ranges::iterator_t<Container>>;

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

    /**
     * @brief Strategy for plotting during the solve process.
     *
     * This enum class defines the strategy for plotting the solution during the solve process.
     * It can either plot at every step or only at the end of the simulation.
     */
    enum class PlotStrategy
    {
        /// Plot at every step.
        EACHSTEP,
        /// Plot at the end only.
        ONLYEND
    };

}

#include "impl/utilities.tpp"

#endif /* UTILITIES_D1960633_B6A3_4A73_81CA_EEFCDE02F6EE */
