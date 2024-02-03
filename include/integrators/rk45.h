// rk45.h
#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "type_traits.hpp"

using namespace type_traits;

namespace rhbi {
// Vector A
static constexpr auto A =
    std::array<double, 6>({0, 2.0 / 9, 1.0 / 3, 3.0 / 4, 1.0, 5.0 / 6});

// Matrix B
static constexpr auto B = std::array<std::array<double, 5>, 6>(
    {std::array<double, 5>{0, 0, 0, 0, 0},
     std::array<double, 5>{2.0 / 9, 0, 0, 0, 0},
     std::array<double, 5>{1.0 / 12, 1.0 / 4, 0, 0, 0},
     std::array<double, 5>{69.0 / 128, -243.0 / 128, 135.0 / 64, 0, 0},
     std::array<double, 5>{-17.0 / 12, 27.0 / 4, -27.0 / 5, 16.0 / 15, 0},
     std::array<double, 5>{65.0 / 432, -5.0 / 16, 13.0 / 16, 4.0 / 27,
                           5.0 / 144}});

// Vector C
static constexpr auto C =
    std::array<double, 6>({1.0 / 9, 0, 9.0 / 20, 16.0 / 45, 1.0 / 12, 0});

// Vector CH
static constexpr auto CH = std::array<double, 6>(
    {47.0 / 450, 0, 12.0 / 25, 32.0 / 225, 1.0 / 30, 6.0 / 25});

// Vector CT
static constexpr auto CT = std::array<double, 6>(
    {1.0 / 150, 0, -3.0 / 100, 16.0 / 75, 1.0 / 20, -6.0 / 25});

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


double adapt_h(double h, double tol, double r) {
  double hnew = h;
  double h_scale = 0.84 * std::pow((tol / r), 0.25);
  if (h_scale <= 0.1) {
    hnew *= 0.1;
  } else if (h_scale >= 4.0) {
    hnew *= 4.0;
  } else {
    hnew *= h_scale;
  }
  return hnew;
}

template <typename Y>
double error_estimator(const std::vector<Y>& k, Y& r, double h){

  r *= 0;

  for (auto ii = 0; ii < 6; ++ii) {
    r += CT[ii] * k[ii];
  }

  double err;

  if constexpr (is_plain_double<Y>::value) {
        err = std::abs(r) / h;
    } else if constexpr (is_real_double_eigen<Y>::value) {
        // Handle Eigen vectors and matrices
        // You may need to implement appropriate norms or other methods for error estimation
        err = r.norm() / h;
    } else {
        // Throw a "not implemented" error for unsupported types
        throw std::logic_error("Error estimator not implemented for this type");
    }

  return err;
}

template <typename X, typename Y>
struct RungeKutta45 {
  std::function<Y(X, Y)> f;
  size_t size; // size if the integrand is some Vector or Matrix
  size_t steps;
  size_t i;
  Y r;
  RungeKutta45(std::function<Y(X, Y)> func) : f(func), size(1), steps(0), i(0), r(get_zero<Y>()){}
  RungeKutta45(std::function<Y(X, Y)> func, size_t size) : f(func), size(size), steps(0), i(0), r(get_zero<Y>(size)){}
  void operator()(const X& xthis, const Y& ythis, X& xnext, Y& ynext, double& h,
                  double tol = 1e-5, double hmin = 1e-4, double hmax = 1.0) {
    bool flag = true;
    while (flag) {
      steps++;
      // std::cout << "is k alright?" << std::endl;
      auto k = get_ks<Y>(6, size);
      k[0] = h * f(xthis + h * A[0], ythis);
      k[1] = h * f(xthis + h * A[1], ythis + B[1][0] * k[0]);
      k[2] = h * f(xthis + h * A[2], ythis + B[2][0] * k[0] + B[2][1] * k[1]);
      k[3] = h * f(xthis + h * A[3],  //
                   ythis + B[3][0] * k[0] + B[3][1] * k[1] + B[3][2] * k[2]);
      k[4] = h * f(xthis + h * A[4],  //
                   ythis + B[4][0] * k[0] + B[4][1] * k[1] + B[4][2] * k[2] +
                       B[4][3] * k[3]);
      k[5] = h * f(xthis + h * A[5],  //
                   ythis + B[5][0] * k[0] + B[5][1] * k[1] + B[5][2] * k[2] +
                       B[5][3] * k[3] + B[5][4] * k[4]);

      double err = error_estimator<Y>(k, r, h);
      // std::cout << "is k alright? err: "<< err << std::endl;

      if (err <= tol) {
        xnext = xthis + h;
        ynext = ythis + (CH[0] * k[0] + CH[1] * k[1] + CH[2] * k[2] +
                         CH[3] * k[3] + CH[4] * k[4] + CH[5] * k[5]);
        i++;
      }

      h = adapt_h(h, tol, err);

      if (h > hmax) {
        h = hmax;
      } else if (h < hmin) {
        h = hmin;
      }

      if (err <= tol) {
        flag = false;
      }
    }
  }
};  // struct RungeKutta45

}  // namespace rhbi
