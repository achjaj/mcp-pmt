// wrapper for a function from boost c++ library

#include <boost/math/special_functions/gamma.hpp>

#include "boost_wrapper.hpp"

extern "C" {
    double gamma_p_inv(double a, double p) {
        return boost::math::gamma_p_inv(a, p);
    }
    
    double erf_inv(double x) {
        return boost::math::erf_inv(x);
    }
}
