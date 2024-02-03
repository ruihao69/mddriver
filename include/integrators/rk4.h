// rk45.h
#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "type_traits.hpp"

using namespace type_traits;

namespace rhbi {

// Function to get a zero-initialized object of type Y
template <typename Y>
Y get_zero(size_t size = 1) {
    if constexpr (is_eigen_matrix<Y>::value) {
        if constexpr (std::decay_t<Y>::RowsAtCompileTime == Eigen::Dynamic) {
            if constexpr (std::decay_t<Y>::ColsAtCompileTime == Eigen::Dynamic) {
                // Dynamically sized Eigen matrix
                return Y::Zero(size, size);
            } else {
                // Dynamically sized Eigen vector (column vector)
                return Y::Zero(size, 1);
            }
        } else {
            // Fixed-sized Eigen matrix or vector
            return Y::Zero();
        }
    } else {
        // For other types, use the default constructor (assumes the type has a default constructor)
        return Y{};
    }
}

// Function to get a vector of ks
template <typename Y>
std::vector<Y> get_ks(size_t stage, size_t size = 1) {
    std::vector<Y> ks(stage);

    for (size_t i = 0; i < stage; ++i) {
        ks[i] = get_zero<Y>(size);
    }

    return ks;
}

template <typename X, typename Y>
struct RungeKutta4 {
  std::function<Y(X, Y)> f;
  double dt;
  size_t size; // size if the integrand is some Vector or Matrix
  size_t steps;
  size_t i;
  std::vector<Y> k;
  Y r;
  RungeKutta4(std::function<Y(X, Y)> func) : f(func), size(1), steps(0), i(0), r(get_zero<Y>(), k(get_ks<Y>(4, 1))){}
  RungeKutta4(std::function<Y(X, Y)> func, size_t size) : f(func), size(size), steps(0), i(0), r(get_zero<Y>(size), k(get_ks<Y>(4, size))){}
  void operator()(const X& xthis, const Y& ythis, X& xnext, Y& ynext, double h) {
    k[0] = h * f(xthis + h * 0.5, ythis);

    r = ythis + k[0] * 0.5;
    k[1] = h * f(xthis + h * 0.5, r);

    r = ythis + k[1] * 0.5;
    k[2] = h * f(xthis + h * 0.5, r);

    r = ythis + k[2];
    k[3] = h * f(xthis + h, r);
  }
};  // struct RungeKutta4

}  // namespace rhbi
