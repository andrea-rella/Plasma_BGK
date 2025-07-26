#ifndef UTILITIES_B121BCD3_6455_43E7_9630_E0343B18BD8C
#define UTILITIES_B121BCD3_6455_43E7_9630_E0343B18BD8C

#include <concepts>
#include <ranges>
#include <string>

namespace Bgk
{

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
     * This concept checks if a callable type can be invoked with an integer index and returns a value convertible to T.
     *
     * @tparam SpacingFunc callable type that defines the spacing.
     * @tparam T precision type of the mesh components.
     */
    template <typename SpacingFunc, typename T>
    concept SpacingFunction = requires(SpacingFunc func, int index) {
        { func(index) } -> std::convertible_to<T>;
    };

}

#endif /* UTILITIES_B121BCD3_6455_43E7_9630_E0343B18BD8C */
