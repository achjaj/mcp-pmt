#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <stdarg.h>

#include "pivi.h"

#define max(a, b) a > b ? a : b
#define EVJFACTOR 1.602176634e-19

#define LOG_ON 0 // do not log: 0; do log: 1

pthread_mutex_t eval_lock = PTHREAD_MUTEX_INITIALIZER;

struct PiviParams *params = NULL;

double *newEnergies = NULL;

double *reinitializationEnergies = NULL;
int *reinitializationFlag = NULL;
double *reinitV1 = NULL;
double *reinitV2 = NULL;

int *removeParticle = NULL;

int nextPidx = -1;
int maxParticles = -1;

unsigned long N = -1;
unsigned long maxN = -1;

char initialized = 0;

void printlog(const char *fmt, ...) {
    if (LOG_ON == 1) {
        va_list arg_ptr;
        va_start(arg_ptr, fmt);
        vfprintf(stderr, fmt, arg_ptr);
        va_end(arg_ptr);
    }
}

void dumpArray(double *arr) {
    for (int i = 0; i < maxParticles; i++) {
        fprintf(stderr, "%12.5e ", arr[i]);
    }
    fprintf(stderr, "\n");
}

void dumpArrayI(int *arr) {
    for (int i = 0; i < maxParticles; i++) {
        fprintf(stderr, "%d ", arr[i]);
    }
    
    fprintf(stderr, "\n");
}

void dump() {
    if (LOG_ON == 0)
        return;
    
    fprintf(stderr, "===================================================================\n");
    
    fprintf(stderr, "newEnergies: ");
    dumpArray(newEnergies);
    
    fprintf(stderr, "reinitializationEnergies: ");
    dumpArray(reinitializationEnergies);
    
    fprintf(stderr, "reinitializationFlag: ");
    dumpArrayI(reinitializationFlag);
    
    fprintf(stderr, "reinitV1: ");
    dumpArray(reinitV1);
    
    fprintf(stderr, "reinitV2: ");
    dumpArray(reinitV2);
    
    fprintf(stderr, "removeParticle: ");
    dumpArrayI(removeParticle);
    
    fprintf(stderr, "nextPidx: %d\n", nextPidx);
    
    fprintf(stderr, "maxParticles: %d\n", maxParticles);
    
    fprintf(stderr, "N: %lu\n", N);
    
    fprintf(stderr, "maxN: %lu\n", maxN);
    
    fprintf(stderr, "===================================================================\n");
}

void lock_calls() {
    pthread_mutex_lock(&eval_lock);
}

void unlock_calls() {
    pthread_mutex_unlock(&eval_lock);
}

double speedFromEV(double ev, double m) {
    return sqrt(2*ev*EVJFACTOR/m);
}

double sum(double *arr, size_t len) {
    double res = 0;
    for (int i = 0; i < len; i++) {
        res += arr[i];
    }
    
    return res;
}

void cleanup() {
    printlog("Cleannig up\n");
    
    free(newEnergies);
    free(reinitializationEnergies);
    free(reinitializationFlag);
    free(reinitV1);
    free(reinitV2);
    free(removeParticle);
    
}

int init(const char *str) {
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
    
    lock_calls();
    if (initialized == 1) {
        printlog("Already initialized\n");
        unlock_calls();
        return 1;
    }
    initialized = 1;
    
    if (params != NULL) {
        free(params->pn);
        free(params->epn);
    }
    free(params);
    cleanup();
    nextPidx = 1;
    printlog("================================ %s =====================================\n", __TIME__);
    
    params = init_pivi(str);
    
    if (params == NULL) {
        //error = "init_pivi returned NULL";
        printlog("ERROR: init_pivi returned NULL\n");
        unlock_calls();
        return 0;
    }
    
    printlog("INIT DONE\n");
    unlock_calls();
    return 1;
}

const char *getLastError() {
    return "";
}

int checkArgs(int nArgs, int expected, const char *func) {
    if (nArgs != expected) {
            printlog("ERROR: %s: wrong number of arguments, expecting %d", func, expected);
            return 0;
        }
        
        return 1;
}

void genReinitVec(int pidx) {
    double *v1v2 = getRndVec();

    reinitV1[pidx] = v1v2[0];
    reinitV2[pidx] = v1v2[1];
    
    free(v1v2);
}

int initPiviArrays(int maxParticles) {
    cleanup();
    
    size_t das = sizeof(double[maxParticles]);
    size_t ias = sizeof(int[maxParticles]);
    
    newEnergies = malloc(das);
    reinitializationEnergies = malloc(das);
    removeParticle = malloc(ias);
    reinitializationFlag = malloc(ias);
    reinitV1 = malloc(das);
    reinitV2 = malloc(das);
        
    if (newEnergies == NULL || reinitializationEnergies == NULL || removeParticle == NULL || reinitializationFlag == NULL || reinitV1 == NULL || reinitV2 == NULL) {
        printlog("ERROR: initPiviArrays: could not initialize arrays\n");
        return 0;
    }
        
    for (int i = 0; i < maxParticles; i++) {
        newEnergies[i] = NAN;
        reinitializationEnergies[i] = NAN;
        removeParticle[i] = 0;
        reinitializationFlag[i] = 0;
        reinitV1[i] = NAN;
        reinitV2[i] = NAN;
    }
    
    nextPidx = 1;
    
    return 1;
}

int eval(const char *func, int nArgs, const double **inReal, const double **inImag, int blockSize, double *outReal, double *outImag) {
    lock_calls();
    printlog("%s\n", func);
    dump();
    
    if (strcmp("initPiviArrays", func) == 0) {
        
        if (!checkArgs(nArgs, 1, func)) {
            unlock_calls();
            return 0;
        }
        
        maxParticles = (int) inReal[0][0] + 10;
        
        outReal[0] = maxParticles;
        int res = initPiviArrays(maxParticles);
        unlock_calls();
        return res;
    }
    else if (strcmp("shouldBeRemoved", func) == 0) {
        if (!checkArgs(nArgs, 1, func)) {
            unlock_calls();
            return 0;
        }
        
        for (int i = 0; i < blockSize; i++) {
            int pidx = ((int) inReal[0][i]) - 1; // in this library particles are indexed from 0, but COMSOL indexes particles from 1
            
            if (pidx >= maxParticles) {
                printlog("pidx out of range\n");
                unlock_calls();
                return 0;
            }
            
            outReal[i] = removeParticle[pidx];
        }
    } 
    else if (strcmp("shouldBeReinitialized", func) == 0) {
        if (!checkArgs(nArgs, 1, func)) {
            unlock_calls();
            return 0;
        }
        
        for (int i = 0; i < blockSize; i++) {
            int pidx = ((int) inReal[0][i]) - 1; // in this library particles are indexed from 0, but COMSOL indexes particles from 1
            
            if (pidx >= maxParticles) {
                printlog("pidx out of range\n");
                unlock_calls();
                return 0;
            }
            
            outReal[i] = reinitializationFlag[pidx];
        }
    } 
    else if (strcmp("reinitializeVelocity", func) == 0) {
        if (!checkArgs(nArgs, 3, func)) {
            unlock_calls();
            return 0;
        }
        
        for (int i = 0; i < blockSize; i++) {
            int pidx = ((int) inReal[0][i]) - 1; // in this library particles are indexed from 0, but COMSOL indexes particles from 1
            double m = inReal[1][i];
            int vcomponent = (int) inReal[2][i];
            
            if (pidx >= maxParticles) {
                printlog("pidx out of range\n");
                unlock_calls();
                return 0;
            }
            
            reinitializationFlag[pidx] = max(0, reinitializationFlag[pidx] - 1);
            
            double v1 = reinitV1[pidx];
            double v2 = reinitV2[pidx];
            double s = pow(v1, 2) + pow(v2, 2);
            double speed = speedFromEV(reinitializationEnergies[pidx], m);
            
            switch(vcomponent) {
                case 1: // x component
                    outReal[i] = (2*s-1)*speed;
                    break;
                case 2: // y component
                    outReal[i] = 2*sqrt(1-s)*speed*v2;
                    break;
                case 3: // z component
                    outReal[i] = 2*sqrt(1-s)*speed*v1;
                    break;
            }
        }
    }
    else if (strcmp("pivi", func) == 0) {
        if (!checkArgs(nArgs, 3, func)){
            unlock_calls();
            return 0;
        }
        
        for (int i = 0; i < blockSize; i++) {
            if (nextPidx >= maxParticles) {
                printlog("nextPidx is at maximum, no new electrons should be produced!\n");
                outReal[i] = 0;
                continue;
            }
            
            double Ek = inReal[0][i]; // ASSUMES TO BE IN eV !!!!!
            double theta = inReal[1][i];
            int pidx = ((int) inReal[2][i]) - 1; // in this library particles are indexed from 0, but COMSOL indexes particles from 1
            
            if (pidx >= maxParticles) {
                printlog("ERROR: pidx out of range\n");
                return 0;
            }
            
            double resEng[params->maxNoe+5];
            int noe = see(Ek, theta, params, N < params->maxNoe ? N : params->maxNoe, resEng);
            printlog("NOE: %d\n", noe);
            printlog("New eneries: ");
            
            for (int j = 0; j < noe; j++) {
                if (isnan(resEng[j])) {
                    printlog("ERROR: One of new energies is NAN!\n");
                    unlock_calls();
                    return 0;
                }
                
                newEnergies[nextPidx] = resEng[j];
                if (++nextPidx >= maxParticles) {
                    printlog("nextPidx is at maximum, no new electrons should be produced!\n");
                    noe = j + 1;
                    break;
                }
                
                printlog("%12.5e ", resEng[j]);
            }
            printlog("\n");
            
            if (noe == 0) {
                removeParticle[pidx] = 1;
                outReal[i] = 0;
            } else if (noe == 1) {
                reinitializationEnergies[pidx] = resEng[0];
                reinitializationFlag[pidx] = 3;
                genReinitVec(pidx);
                outReal[i] = 0;
                nextPidx--;
            } else {
                double s = sum(resEng, noe);
                if (s > Ek) {
                    printlog("WARNING: resEng sum > Ek; taking value from resEng to prevent crashing!\n");
                    s = -resEng[0] - Ek;
                }
                reinitializationEnergies[pidx] = Ek - s;
                reinitializationFlag[pidx] = 3;
                genReinitVec(pidx);
                N = max(0, N - noe + 1);
                outReal[i] = noe;
            }
        }
    }
    else if (strcmp("initSpeed", func) == 0) {
        if (!checkArgs(nArgs, 2, func)) {
            unlock_calls();
            return 0;
        }
        
        for (int i  = 0; i < blockSize; i++) {
            int pidx = ((int) inReal[0][i]) - 1; // in this library particles are indexed from 0, but COMSOL indexes particles from 1
            
            
            
            if (pidx >= maxParticles) {
                printlog("pidx out of range\n");
                unlock_calls();
                return 0;
            }
            
            if (isnan(newEnergies[pidx])) {
                printlog("ERROR: New energy is nan, pidx: %d\n", pidx);
                return 0;
            }
            
            double mass = inReal[1][i];
            double v = speedFromEV(newEnergies[pidx], mass);
            
            printlog("New energy/speed for %d: %12.5e/%12.5e\n", pidx, newEnergies[pidx], v);
            
            outReal[i] = v;
        }
    }
    else if (strcmp("chargeWall", func) == 0) {
        if (!checkArgs(nArgs, 2, func)) {
            unlock_calls();
            return 0;
        }
        
        if (N < maxN) {
            double dNodt = inReal[0][0];
            double dt = inReal[1][0];
            
            unsigned long newN = round(N + dNodt*dt);
            
            N = newN < maxN ? newN : maxN;
        }
        
        outReal[0] = N;
    }
    else if (strcmp("initWallCharge", func) == 0) {
        if (!checkArgs(nArgs, 3, func)) {
            unlock_calls();
            return 0;
        }
        
        double C = inReal[0][0];
        double V = inReal[1][0];
        double q = inReal[2][0];
        
        N = (unsigned long) (C*V/q);
        maxN = N;
        
        outReal[0] = N;
    }
    else {
        printlog("ERROR: libcfun: Unknown function name\n");
        
        unlock_calls();
        return 0;
    }
    
    unlock_calls();
    return 1;
} 
