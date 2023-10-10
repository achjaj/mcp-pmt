/**
 * Implementaion of the Furman-Pivi model of secondary electron emission:
 *          Furman, M. and Pivi, M. Probabilistic model for the simulation of secondary electron emission.
 *          Physical review special topics-accelerators and beams, 5(12):124404, 2002.
 *
 * Based on the Python Implementaion done for PyECLOUD: https://github.com/PyCOMPLETE/PyECLOUD/blob/master/sec_emission_model_furman_pivi.py
 */

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <string.h>

#include "pivi.h"
#include "boost_wrapper.hpp"

const gsl_rng_type *rng_t;
const gsl_rng *rng;

double *deltas = NULL;

void init_rng() {
    gsl_rng_env_setup();
    
    rng_t = gsl_rng_default;
    rng = gsl_rng_alloc(rng_t);
    
    FILE *fp = fopen("/dev/urandom", "r");
    
    size_t nbytes = sizeof(unsigned long int);
    unsigned char bytes[nbytes];
    
    fread(bytes, sizeof(*bytes), sizeof(bytes), fp);
    fclose(fp);
    
    unsigned long int seed = bytes[0];
    for (int i = 1; i < nbytes; i++) {
        seed += bytes[i] << i*8;
    }
    
    gsl_rng_set(rng, seed);
}

double *parse_arr(char *str, const char *delim, int *n_items) {
    int idx = 0;
    //fprintf(stderr, "parse_arr malloc...");
    double *arr = malloc(20*sizeof(double));
    //fprintf(stderr, "ok\n");
    char *token = strtok(str, delim);
    
    if (token == NULL)
        return NULL;
    
    arr[idx++] = atof(token);
    while ((token = strtok(NULL, delim)) != NULL) {
        arr[idx++] = atof(token);
    }
    
    *n_items = idx;
    
    return arr;
}

struct PiviParams *init_pivi(const char* path) {
    init_rng();
    
    //fprintf(stderr, "init_pivi free deltas\n");
    free(deltas);
    
    //fprintf(stderr, "init_pivi malloc deltas...");
    deltas = malloc(3*sizeof(double));
    //fprintf(stderr, "ok\n");
    
    static struct PiviParams params;
    
    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open pivi params file\n");
        return NULL;
    }
    size_t nbytes = 0;
    size_t line_len;
    char *line = NULL;
    
    // read pn; values are delimited by \t
    line_len = getline(&line, &nbytes, fp);
    if (line_len == -1) {
        fprintf(stderr, "Load pivi params: faile reading pn line\n");
        return NULL;
    }
    params.pn = parse_arr(line, "\t", &params.maxNoe);
    
    // read epn; values are delimited by \t
    line_len = getline(&line, &nbytes, fp);
    if (line_len == -1) {
        fprintf(stderr, "Load pivi params: faile reading pn line\n");
        return NULL;
    }
    params.epn = parse_arr(line, "\t", &params.maxNoe);
    
    params.maxNoe -= 1;
    
    // read the rest
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.e1 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.e2 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.r1 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.r2 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params. P1e_inf = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.P1e = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.Ee = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.W = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.p = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.P1r_inf = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.Er = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.r= atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.delta = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.t1 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.t2 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.t3 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.t4 = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.sigma = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.q = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.E = atof(line);
    
    line_len = getline(&line, &nbytes, fp);
    line[line_len - 1] = '\0';
    params.s= atof(line);
    
    // cleanup
    fclose(fp);
    //fprintf(stderr, "pivi_init: free line\n");
    free(line);
    
    return &params;
}

int min(int a, int b, int c) {
    if (a < b && a < c)
        return a;
    else if (b < a && b < c)
        return b;
    
    return c;
}

double f1e(double E, double E0, double de, double sigma) {
    return (2*exp(-(pow(E-E0, 2))/(2*pow(sigma, 2))))/(sqrt(2*PI)*sigma*erf(E0/(sqrt(2)*sigma)));
}

double f1r(double E, double E0, double dr, double q) {
    return (q + 1)*(pow(E, q))/pow(E0, q+1);
}

double fnts(double E, double p, double epsilon) {
    return pow(E, p - 1)*exp(-E/epsilon);
}

double D(double x, double s) {
    return s*x/(s - 1 + pow(x, s));
}

double Ets(double theta, double E, double t3, double t4) {
    return E*(1 + t3*(1 - pow(cos(theta), t4)));
}

double cdf_ts(double x, double eps, double p) {
    return gsl_sf_gamma_inc_P(p, x/eps);
}

double inv_cdf_ts(double x, double eps, double p, double norm) {
    return eps*gamma_p_inv(p, x/norm);
}

double inv_cdf_r(double E0, double x, double q) {
    return E0*pow(x, 1/(q+1));
}

double *getDeltas0(double E, double theta, struct PiviParams *params) {
    deltas[0] = params->P1e_inf + (params->P1e - params->P1e_inf)*exp(-(pow(fabs(E - params->Ee)/params->W, params->p))/params->p);  
    deltas[1] = params->P1r_inf*(1 - exp(-pow(E/params->Er, params->r)));
    deltas[2] = params->delta*(1 + params->t1*(1 - pow(cos(theta), params->t2)));
    
    return deltas;
}

double *getDeltas(double E, double theta, struct PiviParams *params) {
   double *ds0 = getDeltas0(E, theta, params); 

   ds0[0] *= (1 + params->e1*(1 - pow(cos(theta), params->e2)));
   ds0[1] *= (1 + params->r1*(1 - pow(cos(theta), params->r2)));
   ds0[2] *= D(E/Ets(theta, params->E, params->t3, params->t4), params->s);

   ds0[2] = ds0[2]/(1 - ds0[0] - ds0[1]);

   return ds0;
}

void getTSEnergies(double E0, int n, double *epses, double *ps, double *energies) {
   double eps = epses[n];
   double p = ps[n];
   //fprintf(stderr, "getTSEnergies: E0: %12.5e, n: %d, eps: %f, p: %f\n", E0, n, eps, p);
   double A = cdf_ts(E0/n, eps, p);
   double norm = 1/A;

   for (int i =0; i < n; i++) {
        energies[i] = inv_cdf_ts(gsl_rng_uniform(rng), eps, p, norm);
   }
}

double getRedifusedEnergy(double E0, double q) {
    return inv_cdf_r(E0, gsl_rng_uniform(rng), q);
}

double getElasticEnergy(double E0, double sigma) {
    double sqrt2 = sqrt(2);
    double u = gsl_rng_uniform(rng);
    double erfv = erf(E0/(sqrt2*sigma));
    double erfinvv = erf_inv((u - 1)*erfv);
    //fprintf(stderr, "EE: erfv = %12.5e u = %12.5e sqrt2 = %12.5e erfinvv = %12.5e\n", erfv, u, sqrt2, erfinvv);
    if (isnan(u))
        fprintf(stderr, "WARNING: gsl_rng_uniform returned nan!");
    
    return E0 + sqrt2*sigma*erfinvv;
}

int checkForNAN(double *arr, size_t len) {
    for (int i = 0; i < len; i++) {
        if (isnan(arr[i]))
            return i;
    }
    
    return -1;
}

// calculates number of secondary electrons and their energies; return number of new electrons and puts the energies into 'out' array
int see(double E0, double theta, struct PiviParams *params, int M, double *out) {
    //fprintf(stderr, "see: E0: %12.5e theta: %12.5e\n", E0, theta);
    
    double de, dr, dts_tild;
    double *deltas = getDeltas(E0, theta, params);
    int try = 0; // sometimes this function returns NAN, I want to try to prevent this by generating the energies again, but maximali 5 times
    
    de = deltas[0];
    dr = deltas[1];
    dts_tild = deltas[2];

    double u = gsl_rng_uniform(rng);

    if (u > de + dr) {
        unsigned int pv = gsl_ran_poisson(rng, dts_tild);
        int N = pv < M ? pv : M;
        
        int nanPos = -1;
        do {
            getTSEnergies(E0, N, params->epn, params->pn, out);
            nanPos = checkForNAN(out, N);
        } while(nanPos != -1 && try++ < 5);
        
        if (nanPos != -1) {
            fprintf(stderr, "WARNING: see: ts energies with nan: ");
            for (int i = 0; i < N; i++) {
                fprintf(stderr, "%12.5e ", out[i]);
            }
            fprintf(stderr, "\n");
        }
        
        /*printf("t,");
        for (int i = 0; i < N-1; i++) {
            printf("%2.5e,", out[i]);
        }
        printf("%2.5e\n", out[N-1]);*/
        
        return N;
    } else if (u < de) {
        double eE = -1;
        do {
            eE = getElasticEnergy(E0, params->sigma);
            out[0] = eE;
        } while(isnan(out[0]) && try++ < 5);
        
        if (isnan(out[0]))
            fprintf(stderr, "WARNING: see: elastic energy is nan! eE = %12.5e E0 = %12.5e sigma = %f\n", eE, E0, params->sigma);
        
        //printf("e,%2.5e\n", eE);
    } else {
        do {
            out[0] = getRedifusedEnergy(E0, params->q);
        } while(isnan(out[0]) && try++ < 5);
        
        if (isnan(out[0]))
            fprintf(stderr, "WARNING: see: redifused energy is nan!\n");
        
        //printf("r,%2.5e\n", out[0]);
    }
    
    return 1;
}

/*
   Based on "The New Method" from:
       Marsaglia, George. “Choosing a Point from the Surface of a Sphere.” The Annals of Mathematical Statistics 43, no. 2 (1972): 645–46. http://www.jstor.org/stable/2240001.
   However we need only one quadrant with z and y >= 0 so the folowing changes were made:
       1. The variables are picked from interval [0, 1]
       2. The resulting vector is rotated by -90° around y axis
*/
double *getRndVec() {
    double s, v1, v2;

    do {
        v1 = gsl_rng_uniform(rng);
        v2 = gsl_rng_uniform(rng);
        
        s = pow(v1, 2) + pow(v2, 2);
    } while (s >= 1);
    
    double *res = malloc(2*sizeof(double));
    res[0] = v1;
    res[1] = v2;
    
    return res;
}

