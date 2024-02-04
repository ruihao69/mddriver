#include <cmath>
#include <iomanip>
#include <iostream>

#include "pulses/gaussian_envolope.h"

constexpr size_t WIDTH = 15;
struct DoubleGaussianEnvolope {
  rhbi::pulses::GaussianEnvolope pump;
  rhbi::pulses::GaussianEnvolope prob;

  DoubleGaussianEnvolope(double A, double t0, double tau, double omega, double delay, double phi = 0.0, double tol = 1e-10)
      : pump(A, t0, tau, omega, phi, tol), prob(A, t0 + delay, tau, omega, phi, tol) {}
  double operator()(double t) const { return pump(t) + prob(t); }
};

int main(int argc, char** argv) {
  // GaussianEnvolope(double A, double t0, double tau, double omega, double phi = 0.0, double tol=1e-10);
  double A = 1;
  double t0 = 0;
  double tau = 1;
  double omega = 20;
  double delay = 15;

  // rhbi::pulses::GaussianEnvolope pulse(A, t0, tau, omega, 0.0, 1e-6);
  auto pulse = DoubleGaussianEnvolope(A, t0, tau, omega, delay, 0, 1e-8);

  double dt = 0.03 * 2.0 * M_PI / omega;

  double t = -15;
  std::cout << "# t" << std::setw(WIDTH) << "pulse" << std::setw(WIDTH) << std::endl;
  while (t < 30) {
    std::cout << t << std::setw(WIDTH) << std::setprecision(6) << std::fixed
              << pulse(t) << std::setw(WIDTH) << std::setprecision(6) << std::endl;
    // << pulse(t) << std::endl;
    t += dt;
  }
  return 0;
}