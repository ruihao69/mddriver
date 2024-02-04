#include "pulses/gaussian_envolope.h"

#include <cmath>  // Add missing include
#include <iostream>

namespace rhbi {
namespace pulses {

GaussianEnvolope::GaussianEnvolope(double A, double t0, double tau, double omega, double phi, double tol)
    : A(A), t0(t0), tau(tau), omega(omega), phi(phi) {
  estimate_pulse_range(tol);
}

double GaussianEnvolope::operator()(double t) const {
  // approximate the pulse by zero when t is out of range
  if (t < tmin || t > tmax) {
    return 0;
  }
  // evaluate the gaussian envolope pulse
  return A * gaussian(t) * cosine(t);
}

double GaussianEnvolope::gaussian(double t) const {
  return exp(-pow(0.5 * (t - t0) / tau, 2));
}

double GaussianEnvolope::cosine(double t) const {
  return cos(omega * t + phi);
}

void GaussianEnvolope::estimate_pulse_range(double tol) {
  double _tmax = t0, _tmin = t0;
  double step = 0.1 * tau;  // Initial a step size

  while (gaussian(_tmin) > tol) {
    _tmin -= step;
  }

  // Perform gradient descent to estimate pulse range
  // Use adaptive step size to avoid overshooting
  // while (gaussian(_tmin) > tol && step > STEP_SIZE_THRESHOLD) {
  //   double gradient = -2 * (_tmin - t0) / pow(tau, 2) * gaussian(_tmin);
  //   double prevTmin = _tmin;

  //   _tmin -= gradient * step;

  //   // Apply Armijo rule to adjust step size
  //   while (gaussian(_tmin) > gaussian(prevTmin) - 0.5 * step * gradient * gradient && step > STEP_SIZE_THRESHOLD) {
  //     step *= 0.5;  // Reduce step size by half
  //   }

  //   // Update _tmin outside the inner while loop
  //   _tmin = prevTmin - gradient * step;

  //   std::cout << _tmin << std::endl;
  // }

  double Dt = t0 - _tmin;
  _tmax = t0 + Dt;

  tmin = _tmin;
  tmax = _tmax;
}

// void GaussianEnvolope::estimate_pulse_range(double tol) {
//   double _tmax, _tmin = t0;
//   constexpr double STEP_SIZE_THRESHOLD = 1e-6;
//   double step = 0.1 * tau;  // Initial a step size
//
//   // Perform gradient descent to estimate pulse range
//   // Use adaptive step size to avoid overshooting
//   while (gaussian(_tmin) > tol && step > STEP_SIZE_THRESHOLD) {  // Modify the condition to include a check for step size
//     double gradient = 2 * (_tmin - t0) / pow(tau, 2) * gaussian(_tmin);
//     double prevTmin = _tmin;
//     _tmin -= gradient * step;
//
//     // Apply Armijo rule to adjust step size
//     while (gaussian(_tmin) > gaussian(prevTmin) - 0.5 * step * gradient * gradient && step > STEP_SIZE_THRESHOLD) {  // Modify the condition to include a check for step size
//       step *= 0.5;                                                                                                   // Reduce step size by half
//       _tmin = prevTmin - gradient * step;
//     }
//   }
//
//   double Dt = t0 - _tmin;
//   _tmax = t0 + Dt;
//
//   tmin = _tmin;
//   tmax = _tmax;
// }

}  // namespace pulses
}  // namespace rhbi
