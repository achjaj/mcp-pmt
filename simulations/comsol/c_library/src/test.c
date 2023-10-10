// This program was used only for basic testing of the library

#include "pivi.h"

#include <stdio.h>
#include <dlfcn.h>
#include <math.h>

#include "erfinv.h"
#include "pivi.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }
    
    if (argc == 2) {
        int N = atoi(argv[1]);
        struct PiviParams* params = init_pivi("/home/jakub/pivi_params.txt");
        double out[100];
    
        for (int i = 0; i < N; i++) {
            see(50, 0, params, 10, out);
        }
    } else {
        char t = argv[1][0];
        double E = atof(argv[2]);
        
        if (t == 'e') {
            double E0 = atof(argv[3]);
            double de = atof(argv[4]);
            double sigma = atof(argv[5]);
            
            printf("%2.5e\n", f1e(E, E0, de, sigma));
        } else if (t == 'r') {
            double E0 = atof(argv[3]);
            double dr = atof(argv[4]);
            double q = atof(argv[5]);
            
            printf("%2.5e\n", f1r(E, E0, dr, q));
        } else {
            double p = atof(argv[3]);
            double epsilon = atof(argv[4]);
            
            printf("%2.5e\n", fnts(E, p, epsilon));
        }
    }
    
    return 0;
}

/*int main() {
    void *lib = dlopen("./libcfunc.so", RTLD_NOW);
    int (*init)(const char *str) = dlsym(lib, "init");
    
    init("/home/jakub/pivi_params.txt");

    double E = 1e3;
    double t = 1.36;
    
    int (*eval)(const char *func, int nArgs, const double **inReal, const double **inImag, int blockSize, double *outReal, double *outImag) = dlsym(lib, "eval");
    double **in = malloc(3*sizeof(*in));
    for (int i = 0; i < 3; i++) {
        in[i] = malloc(sizeof(*in[i]));
    }
    
    in[0][0] = 100;
    
    double out[1] = {0};
    eval("initPiviArrays", 1, in, NULL, 1, out, NULL);
    
    
    in[0][0] = 2.86e-13;
    in[1][0] = 2e3;
    in[2][0] = 1.6e-19;
    
    eval("initWallCharge", 3, in, NULL, 1, out, NULL);
    
    for (int i = 0; i < 1000; i ++) {
        in[0][0] = 1000;
        in[1][0] = 1.24210;
        in[2][0] = 1;
        int res = eval("pivi", 3, in, NULL, 1, out, NULL);
        printf("%f\n", out[0]);
        
        if (res == 0) {
            printf("ERROR!!!!!!!!!!!!!!!!!!!!!!!!\n");
            break;
        }
        
        /*in[0][0] = 1;
        in[1][0] = 1;
        eval("chargeWall", 2, in, NULL, 1, out, NULL);
        printf("%f\n", out[0]);
        
        in[0][0] = 1;
        eval("shouldBeRemoved", 1, in, NULL, 1, out, NULL);
        printf("%f\n", out[0]);
        
        in[0][0] = 2;
        in[1][0] = 9.81e-32;
        eval("initSpeed", 2, in, NULL, 1, out, NULL);
        printf("%12.5e\n", out[0]);*/
    /*}
    
    return 0;
}*/




/*static struct PiviParams *params = NULL;

int main() {
    for (int i = 0; i < 10; i++) {
    struct PiviParams *params = init_pivi("/home/jakub/pivi_params.txt");
    
    if (params == NULL) {
        return 1;
    }
    
    int m = 10;
    double out[m];
    
    /*for (int i = 0; i < 10; i++) {
        for (double theta = 0; theta <= PI/2; theta += PI/20) {
            int generated = see(205, PI/2, params, m, out);
            printf("[");
            for (int i = 0; i < generated; i++) {
                printf("%f ", out[i]);
            }
            printf("]\n");
       // }
        
        printf("==========================================\n");
    //}
    }
    
    return 0;
}*/
