#ifndef PULSE_H
#define PULSE_H

namespace rhbi {

template <typename T>
class Pulse {
 public:
  virtual T operator()(double t) const = 0;
};

}  // namespace rhbi

#endif  // PULSE_H
