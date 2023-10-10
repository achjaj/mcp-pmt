// wrapper for a function from boost c++ library

#ifndef BOOST_WRAPPER_H_INCLUDED
#define BOOST_WRAPPER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

    double gamma_p_inv(double a, double p);
    double erf_inv(double x);
    
#ifdef __cplusplus
}
#endif

#endif // BOOST_WRAPPER_H_INCLUDED
