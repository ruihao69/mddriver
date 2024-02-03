#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>

#include "integrators/rk45.h"

using namespace rhbi;

struct Lorenz{
  double sigma;
  double rho;
  double beta;
  Eigen::Vector3d dY;

  Lorenz(double sigma, double rho, double beta) : sigma(sigma), rho(rho), beta(beta){}

  Eigen::Vector3d operator()(double X, Eigen::Vector3d Y){
    // Eigen::Vector3d dY = Eigen::Vector3d::Zero();
    dY(0) = sigma * (Y(1) - Y(0)); 
    dY(1) = Y(0) * (rho - Y(2)); 
    dY(2) = Y(0) * Y(1) - beta * Y(2); 

    return dY;
  }
};

int main() {
  std::cout << "Hello world!" << std::endl;

  auto lorenz = Lorenz(10.0, 28.0, 8.0/3);

  auto rk45 = RungeKutta45<double, Eigen::Vector3d>(lorenz);

  double x0 = 0.0;
  auto y0 = Eigen::Vector3d({2.0, 1.0, 1.0});

  double xf = 10.;

  double xthis = x0;
  double xnext = xthis;
  Eigen::Vector3d ythis = y0;
  Eigen::Vector3d ynext = y0;

  double hmax = 1.0e-1;
  double hmin = 1.0e-4;
  double tol = 1.0e-8;

  double h = hmax;

  while (xthis < xf) {
    // std::cout << "just checking" << std::endl;
    rk45(xthis, ythis, xnext, ynext, h, tol, hmin, hmax);
    // std::cout << "just checking" << std::endl;
    xthis = xnext;
    ythis = ynext;

    std::cout << xthis << std::setw(15) << ythis[0] << std::setw(15)
              << ythis[1] << std::setw(15) << ythis[2] << std::endl;
  }

  return 0;
}
