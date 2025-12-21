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

#ifndef METRICS_UTILS_CDB44BDA_5E7D_48F7_84A7_B2F8EBB4B9E8
#define METRICS_UTILS_CDB44BDA_5E7D_48F7_84A7_B2F8EBB4B9E8

#include <Eigen/Dense>
#include <memory>
#include <cmath>
#include "utilities.hpp"
#include "SpaceMeshFV.hpp"

namespace Bgk
{
    namespace metrics
    {

        enum class VectorNormType
        {
            L1,
            L2,
            Linf
        };

        // ----- VECTOR NORMS ----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Abstract base class for calculating vector norms.
         *
         * Provides an interface for computing norms on Eigen vectors. Supports both standard
         * and weighted norm calculations.
         *
         * @tparam T The scalar type (e.g., double, float).
         */
        template <typename T>
        class VectorNorm
        {
        public:
            /// @brief Virtual destructor.
            virtual ~VectorNorm() = default;

            /** @brief Computes the norm of a vector.
             * @param v The Eigen vector to process.
             * @return T The computed norm value.
             */
            virtual T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const = 0;

            /** @brief Computes the weighted norm of a vector.
             *
             * This is useful for simulations on non-uniform grids where elements
             * must be weighted by cell volume or length.
             *
             * @param v The Eigen vector to process.
             * @param weights A vector of weights corresponding to each element in v.
             * @return T The computed weighted norm value.
             */
            virtual T compute(const Eigen::Vector<T, Eigen::Dynamic> &v,
                              const std::vector<T> &) const { return compute(v); }
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Implementation of the L2 (Euclidean) norm.
         *
         * Calculates the square root of the sum of squares. For weighted norms:
         * \f[ \sqrt{\sum |v_i|^2 \cdot w_i} \f]
         *
         * @example
         * @code
         * auto l2 = std::make_unique<L2VectorNorm<double>>();
         * Eigen::VectorXd v = {3.0, 4.0};
         * double res = l2->compute(v); // Returns 5.0
         * @endcode
         */
        template <typename T>
        class L2VectorNorm : public VectorNorm<T>
        {
        public:
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const override
            {
                return v.norm();
            }
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v, const std::vector<T> &weights) const
            {
                if (static_cast<size_t>(v.size()) != weights.size())
                    throw std::invalid_argument(error_message("Incompatible mesh and weights sizes"));
                T sum = T{0};
                for (size_t i = 0; i < weights.size(); ++i)
                {
                    sum += v[i] * v[i] * weights[i];
                }
                return std::sqrt(sum);
            }
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Implementation of the L1 (Manhattan) norm.
         *
         * Calculates the sum of absolute values. For weighted norms:
         * \f[ \sum |v_i| \cdot w_i \f]
         *
         * @example
         * @code
         * auto l1 = std::make_unique<L1VectorNorm<double>>();
         * Eigen::VectorXd v = {3.0, -4.0};
         * double res = l1->compute(v); // Returns 7.0
         * @endcode
         */
        template <typename T>
        class L1VectorNorm : public VectorNorm<T>
        {
        public:
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const override
            {
                // Use Eigen's built-in lpNorm for L1. ".template" is required because lpNorm is a template
                // member function being called on a template-dependent type
                return v.template lpNorm<1>();
            }
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v, const std::vector<T> &weights) const
            {
                if (static_cast<size_t>(v.size()) != weights.size())
                    throw std::invalid_argument(error_message("Incompatible vector and weights sizes"));
                T sum = T{0};

                for (size_t i = 0; i < weights.size(); ++i)
                {
                    sum += std::abs(v[i]) * weights[i];
                }
                return sum;
            }
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Implementation of the L∞ (Maximum) norm.
         *
         * Calculates the maximum absolute value in the vector.
         *
         * @example
         * @code
         * auto linf = std::make_unique<LinfVectorNorm<double>>();
         * Eigen::VectorXd v = {3.0, -4.0, 2.0};
         * double res = linf->compute(v); // Returns 4.0
         * @endcode
         */
        template <typename T>
        class LinfVectorNorm : public VectorNorm<T>
        {
        public:
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const override
            {
                // Use Eigen's built-in lpNorm for L∞. ".template" is required because lpNorm is a template
                // member function being called on a template-dependent type
                return v.template lpNorm<Eigen::Infinity>();
            }
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &, const std::vector<T> &weights) const override
            {
                if (static_cast<size_t>(v.size()) != weights.size())
                    throw std::invalid_argument(error_message("Incompatible vector and weights sizes"));
                T max_val = T{0};
                for (size_t i = 0; i < weights.size(); ++i)
                {
                    T val = std::abs(v[i]) * weights[i];
                    if (val > max_val)
                        max_val = val;
                }
                return max_val;
            }
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Factory for creating VectorNorm instances.
         *
         * Provides a static method to return a unique pointer to a VectorNorm instance based on the specified type.
         *
         * @tparam T The scalar type (e.g., double, float).
         */
        template <typename T>
        class VectorNormFactory
        {
        public:
            /** @brief Creates a VectorNorm instance based on the specified type.
             *
             * @param type The type of vector norm to create (L1, L2, Linf).
             * @return std::unique_ptr<VectorNorm<T>> A unique pointer to the created VectorNorm instance.
             * @throws std::invalid_argument If an unknown vector norm type is specified.
             */
            static std::unique_ptr<VectorNorm<T>> create(VectorNormType type)
            {
                switch (type)
                {
                case VectorNormType::L1:
                    return std::make_unique<L1VectorNorm<T>>();
                case VectorNormType::L2:
                    return std::make_unique<L2VectorNorm<T>>();
                case VectorNormType::Linf:
                    return std::make_unique<LinfVectorNorm<T>>();
                default:
                    throw std::invalid_argument("Unknown vector norm type");
                }
            }
        };

        // ----- ROW AGGREGATORS -------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        enum class RowAggregateType
        {
            Average,
            Max
        };

        /** @brief Abstract base class for aggregating row-wise metrics.
         *
         * Provides an interface for aggregating a collection of values derived from
         * row-wise computations.
         *
         * @tparam T The scalar type (e.g., double, float).
         */
        template <typename T>
        class RowAggregator
        {
        public:
            /// @brief Virtual destructor.
            virtual ~RowAggregator() = default;
            /** @brief Aggregates a collection of values.
             * @param values The collection of values to aggregate.
             * @return T The aggregated result.
             */
            virtual T aggregate(const std::vector<T> &values) const = 0;
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Implementation of the average row aggregator.
         *
         * Computes the arithmetic mean of a collection of values.
         *
         * @example
         * @code
         * auto avg_agg = std::make_unique<AverageRowAggregator<double>>();
         * std::vector<double> values = {2.0, 4.0, 6.0};
         * double res = avg_agg->aggregate(values); // Returns 4.0
         * @endcode
         */
        template <typename T>
        class AverageRowAggregator : public RowAggregator<T>
        {
        public:
            T aggregate(const std::vector<T> &values) const override
            {
                if (values.empty())
                    throw std::invalid_argument("AverageRowAggregator: no values to aggregate");
                const T sum = std::accumulate(values.begin(), values.end(), T{0});
                return sum / static_cast<T>(values.size());
            }
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Implementation of the maximum row aggregator.
         *
         * Identifies the maximum value from a collection of values.
         *
         * @example
         * @code
         * auto max_agg = std::make_unique<MaxRowAggregator<double>>();
         * std::vector<double> values = {2.0, 4.0, 6.0};
         * double res = max_agg->aggregate(values); // Returns 6.0
         * @endcode
         */
        template <typename T>
        class MaxRowAggregator : public RowAggregator<T>
        {
        public:
            T aggregate(const std::vector<T> &values) const override
            {
                if (values.empty())
                    throw std::invalid_argument("MaxRowAggregator: no values to aggregate");
                return *std::max_element(values.begin(), values.end());
            }
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Factory for creating RowAggregator instances.
         *
         * Provides a static method to return a unique pointer to a RowAggregator instance based on the specified type.
         *
         * @tparam T The scalar type (e.g., double, float).
         */
        template <typename T>
        class RowAggregatorFactory
        {
        public:
            /** @brief Creates a RowAggregator instance based on the specified type.
             *
             * @param type The type of row aggregator to create (Average, Max).
             * @return std::unique_ptr<RowAggregator<T>> A unique pointer to the created RowAggregator instance.
             * @throws std::invalid_argument If an unknown row aggregate type is specified.
             */
            static std::unique_ptr<RowAggregator<T>> create(RowAggregateType type)
            {
                switch (type)
                {
                case RowAggregateType::Average:
                    return std::make_unique<AverageRowAggregator<T>>();
                case RowAggregateType::Max:
                    return std::make_unique<MaxRowAggregator<T>>();
                default:
                    throw std::invalid_argument("Unknown row aggregate type");
                }
            }
        };

        // ----- MATRIX NORMS ----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Abstract base class for calculating matrix norms.
         *
         * Provides an interface for computing norms on Eigen matrices. Supports both standard
         * and weighted norm calculations.
         *
         * @tparam T The scalar type (e.g., double, float).
         */
        template <typename T>
        class MatrixNorm
        {
        public:
            /// @brief Virtual destructor.
            virtual ~MatrixNorm() = default;

            /** @brief Computes the norm of the difference between two matrices.
             * @param current The current Eigen matrix.
             * @param prev The previous Eigen matrix.
             * @return T The computed norm value.
             */
            virtual T compute(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &current,
                              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &prev) const = 0;

            /** @brief Computes the weighted norm of the difference between two matrices.
            @param current The current Eigen matrix.
            @param prev The previous Eigen matrix.
            @param weights A vector of weights corresponding to each row in the matrices.
            @return T The computed weighted norm value.
             */
            virtual T compute(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &current,
                              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &prev,
                              const std::vector<T> &) const
            {
                return compute(current, prev);
            }
            // needed because overload resolution happens on the static type, not the dynamic one. So if i want
            // the factory to be as generic as possible i have to declare it even for methods that do not need it
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Implementation of a row-wise matrix norm.
         *
         * Calculates the norm of the difference between two matrices on a row-wise basis,
         * normalizing each row's difference by the norm of the current row, and then
         * aggregating the results using a specified row aggregator.
         *
         * @example
         * @code
         * auto vec_norm = std::make_unique<L2VectorNorm<double>>();
         * auto aggregator = std::make_unique<AverageRowAggregator<double>>();
         * auto row_wise_norm = std::make_unique<RowWiseMatrixNorm<double>>(std::move(vec_norm), std::move(aggregator));
         * Eigen::MatrixXd current = ...;
         * Eigen::MatrixXd prev = ...;
         * double res = row_wise_norm->compute(current, prev); // Computes the row-wise norm
         * @endcode
         */
        template <typename T>
        class RowWiseMatrixNorm : public MatrixNorm<T>
        {
        public:
            /// @brief Constructor that initializes the vector norm and row aggregator.
            RowWiseMatrixNorm(std::unique_ptr<VectorNorm<T>> vec_norm,
                              std::unique_ptr<RowAggregator<T>> aggregator)
                : vec_norm_(std::move(vec_norm)), aggregator_(std::move(aggregator)) {}

            T compute(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &current,
                      const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &prev) const override
            {
                const auto rows = static_cast<int>(current.rows());
                if (rows == 0)
                    throw std::invalid_argument("RowWiseMatrixNorm: matrix has zero rows");

                std::vector<T> row_norms;
                row_norms.reserve(static_cast<size_t>(rows));
                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> diff = current - prev;

                for (int i = 0; i < rows; ++i)
                {
                    // Copy row into a column vector for the VectorNorm interface
                    Eigen::Vector<T, Eigen::Dynamic> diff_v = diff.row(i).transpose();
                    Eigen::Vector<T, Eigen::Dynamic> current_v = current.row(i).transpose();
                    row_norms.push_back(vec_norm_->compute(diff_v) / vec_norm_->compute(current_v));
                }
                return aggregator_->aggregate(row_norms);
            }

            T compute(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &current,
                      const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &prev,
                      const std::vector<T> &weights) const override
            {
                const auto rows = static_cast<int>(current.rows());
                if (rows == 0)
                    throw std::invalid_argument("RowWiseMatrixNorm: matrix has zero rows");

                std::vector<T> row_norms;
                row_norms.reserve(static_cast<size_t>(rows));
                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> diff = current - prev;

                for (int i = 0; i < rows; ++i)
                {
                    // Copy row into a column vector for the VectorNorm interface
                    Eigen::Vector<T, Eigen::Dynamic> diff_v = diff.row(i).transpose();
                    Eigen::Vector<T, Eigen::Dynamic> current_v = current.row(i).transpose();
                    row_norms.push_back(vec_norm_->compute(diff_v, weights) / vec_norm_->compute(current_v, weights));
                }
                return aggregator_->aggregate(row_norms);
            }

        private:
            std::unique_ptr<VectorNorm<T>> vec_norm_;
            std::unique_ptr<RowAggregator<T>> aggregator_;
        };

        // -----------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------

        /** @brief Factory for creating MatrixNorm instances.
         *
         * Provides a static method to return a unique pointer to a MatrixNorm instance
         * based on the specified vector norm type and row aggregate type.
         *
         * @tparam T The scalar type (e.g., double, float).
         */
        template <typename T>
        class MatrixNormFactory
        {
        public:
            /** @brief Creates a MatrixNorm instance based on the specified types.
             *
             * @param vec_norm_type The type of vector norm to use (L1, L2, Linf).
             * @param agg_type The type of row aggregator to use (Average, Max).
             * @return std::unique_ptr<MatrixNorm<T>> A unique pointer to the created MatrixNorm instance.
             */
            static std::unique_ptr<MatrixNorm<T>> create(VectorNormType vec_norm_type,
                                                         RowAggregateType agg_type)
            {
                auto vec_norm = VectorNormFactory<T>::create(vec_norm_type);
                auto aggregator = RowAggregatorFactory<T>::create(agg_type);
                return std::make_unique<RowWiseMatrixNorm<T>>(std::move(vec_norm), std::move(aggregator));
            }
        };

    }
}

#endif /* METRICS_UTILS_CDB44BDA_5E7D_48F7_84A7_B2F8EBB4B9E8 */
