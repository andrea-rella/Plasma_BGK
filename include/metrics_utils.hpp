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

        template <typename T>
        class VectorNorm
        {
        public:
            virtual ~VectorNorm() = default;
            virtual T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const = 0;
            virtual T compute(const Eigen::Vector<T, Eigen::Dynamic> &v,
                              const std::vector<T> &) const { return compute(v); }
        };

        template <typename T>
        class L2VectorNorm : public VectorNorm<T>
        {
        public:
            using VectorNorm<T>::compute;
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

        template <typename T>
        class L1VectorNorm : public VectorNorm<T>
        {
        public:
            using VectorNorm<T>::compute;
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const override
            {
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

        template <typename T>
        class LinfVectorNorm : public VectorNorm<T>
        {
        public:
            using VectorNorm<T>::compute;
            T compute(const Eigen::Vector<T, Eigen::Dynamic> &v) const override
            {
                return v.template lpNorm<Eigen::Infinity>();
            }
        };

        template <typename T>
        class VectorNormFactory
        {
        public:
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

        enum class RowAggregateType
        {
            Average,
            Max
        };

        template <typename T>
        class RowAggregator
        {
        public:
            virtual ~RowAggregator() = default;
            virtual T aggregate(const std::vector<T> &values) const = 0;
        };

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

        template <typename T>
        class RowAggregatorFactory
        {
        public:
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

        // ---- NEW: Matrix norm that applies a row-wise vector norm + aggregation

        template <typename T>
        class MatrixNorm
        {
        public:
            virtual ~MatrixNorm() = default;
            virtual T compute(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &current,
                              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &prev) const = 0;

            virtual T compute(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &current,
                              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &prev,
                              const std::vector<T> &) const
            {
                return compute(current, prev);
            }
            // needed because overload resolution happens on the static type, not the dynamic one. So if i want
            // the factory to be as generic as possible i have to declare it even for methods that do not need it
        };

        template <typename T>
        class RowWiseMatrixNorm : public MatrixNorm<T>
        {
        public:
            RowWiseMatrixNorm(std::unique_ptr<VectorNorm<T>> vec_norm,
                              std::unique_ptr<RowAggregator<T>> aggregator)
                : vec_norm_(std::move(vec_norm)), aggregator_(std::move(aggregator)) {}

            using MatrixNorm<T>::compute;

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

        template <typename T>
        class MatrixNormFactory
        {
        public:
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
