#include "pulse.h"

namespace rhbi {
namespace pulses {

class GaussianEnvolope : public Pulse<double> {
private:
  double A;          // amplitude
  double t0;         // center
  double tau;        // width
  double omega;      // frequency
  double phi;        // phase

  double tmin;       // the estimated time to kick off the pulse
  double tmax;       // the estimated time to end the pulse
public:
  GaussianEnvolope(double A, double t0, double tau, double omega, double phi = 0.0, double tol=1e-8);
  double operator()(double t) const override;
private:
  double gaussian(double t) const;
  double cosine(double t) const;
  void estimate_pulse_range(double tol);
  // enum estimate_pulse_range() const;
};

} // namespace pulses
} // namespace rhbi