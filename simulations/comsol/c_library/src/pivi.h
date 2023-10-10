/**
 * Implementaion of the Furman-Pivi model of secondary electron emission:
 *          Furman, M. and Pivi, M. Probabilistic model for the simulation of secondary electron emission.
 *          Physical review special topics-accelerators and beams, 5(12):124404, 2002.
 *
 * Based on the Python Implementaion done for PyECLOUD: https://github.com/PyCOMPLETE/PyECLOUD/blob/master/sec_emission_model_furman_pivi.py
 */

#ifndef PIVI_H_INCLUDED
#define PIVI_H_INCLUDED

#include <gsl/gsl_rng.h>

#define PI 3.14159265358979323846

struct PiviParams {
    int maxNoe;
    
    double e1,
          e2,
          r1,
          r2;

    double P1e_inf,
          P1e,
          Ee,
          W,
          p;

    double P1r_inf,
          Er,
          r;

    double delta,
          t1,
          t2,
          t3,
          t4;

    double sigma,
          q,
          E;

    double *epn,
            *pn;
            
    double s; // this is the parameter needed for the function D, it can be also measured; see pivi paper
};

struct PiviParams *init_pivi(const char* path); // parse parameters from *path into *params
int see(double E0, double theta, struct PiviParams *params, int M, double *out);
double *getRndVec();
double getElasticEnergy(double E0, double sigma);

double f1e(double E, double E0, double de, double sigma);
double f1r(double E, double E0, double dr, double q);
double fnts(double E, double p, double epsilon);

#endif // PIVI_H_INCLUDED
