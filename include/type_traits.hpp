#pragma once

#include <Eigen/Dense>

namespace type_traits{
// For calculating the error estimator
template <typename T>
struct is_plain_double {
    static constexpr bool value = std::is_same<T, double>::value;
};

template <typename T>
struct is_real_double_eigen {
    static constexpr bool value = false;
};

template <typename T>
struct is_real_double_eigen<Eigen::Matrix<T, Eigen::Dynamic, 1>> {
    static constexpr bool value = std::is_same<T, double>::value;
};

template <typename T, int Size>
struct is_real_double_eigen<Eigen::Matrix<T, Size, 1>> {
    static constexpr bool value = std::is_same<T, double>::value;
};

template <typename T>
struct is_real_double_eigen<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> {
    static constexpr bool value = std::is_same<T, double>::value;
};

template <typename T, int Rows, int Cols>
struct is_real_double_eigen<Eigen::Matrix<T, Rows, Cols>> {
    static constexpr bool value = std::is_same<T, double>::value;
};

// For initializing the the ks
template <typename T>
struct is_eigen_matrix {
    static constexpr bool value = false;
};

template <typename T, int Rows, int Cols>
struct is_eigen_matrix<Eigen::Matrix<T, Rows, Cols>> {
    static constexpr bool value = true;
};

template <typename T>
using enable_if_eigen_matrix = std::enable_if_t<is_eigen_matrix<std::decay_t<T>>::value, int>;


}
