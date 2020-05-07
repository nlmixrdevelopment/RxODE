#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/policies/error_handling.hpp>

extern "C" double gamma_p(double a, double z) {
  return boost::math::gamma_p<double, double>(a, z);
}

extern "C" double gamma_q(double a, double z) {
  return boost::math::gamma_q<double, double>(a, z);
}

extern "C" double tgamma_lower(double a, double z) {
  return boost::math::tgamma_lower<double, double>(a, z);
}

extern "C" double tgamma_upper(double a, double z) {
  return boost::math::tgamma<double, double>(a, z);
}

extern "C" double gamma_p_derivative(double a, double x) {
  return boost::math::gamma_p_derivative<double, double>(a, x);
}

extern "C" double gamma_q_inv(double a, double q) {
  return boost::math::gamma_q_inv<double, double>(a, q);
}

extern "C" double gamma_p_inv(double a, double p) {
  return boost::math::gamma_p_inv<double, double>(a, p);
}

extern "C" double gamma_q_inva(double a, double q) {
  return boost::math::gamma_q_inv<double, double>(a, q);
}

extern "C" double gamma_p_inva(double a, double p) {
  return boost::math::gamma_p_inv<double, double>(a, p);
}
