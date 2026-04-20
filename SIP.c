#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Define grid dimensions (excluding boundaries)
#define NX 50
#define NY 50
#define MAX_ITER 1000
#define TOL 1e-6

void solve_sip(double b[NX][NY], double x[NX][NY], double alpha) {
    // Coefficients for a 5-point stencil: ap*x[i][j] + aw*x[i-1][j] + ae*x[i+1][j] + as*x[i][j-1] + an*x[i][j+1] = b[i][j]
    static double ap[NX][NY], aw[NX][NY], ae[NX][NY], as[NX][NY], an[NX][NY];
    static double Lw[NX][NY], Ls[NX][NY], Ld[NX][NY], Ue[NX][NY], Un[NX][NY];
    static double res[NX][NY], v[NX][NY];
    
    // 1. Incomplete LU Factorization (Stone's Method)
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            double lw = (i > 0) ? aw[i][j] / (1.0 + alpha * Un[i-1][j]) : 0;
            double ls = (j > 0) ? as[i][j] / (1.0 + alpha * Ue[i][j-1]) : 0;
            
            Lw[i][j] = lw;
            Ls[i][j] = ls;
            Ld[i][j] = ap[i][j] + alpha * (lw * Un[i-1][j] + ls * Ue[i][j-1]) - lw * Ue[i-1][j] - ls * Un[i][j-1];
            Ue[i][j] = (ae[i][j] - alpha * lw * Un[i-1][j]) / Ld[i][j];
            Un[i][j] = (an[i][j] - alpha * ls * Ue[i][j-1]) / Ld[i][j];
        }
    }

    // 2. Iterative Loop
    for (int iter = 0; iter < MAX_ITER; iter++) {
        double max_res = 0;

        // Calculate Residual and Forward Substitution: L*v = r
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                double r = b[i][j] - (ap[i][j]*x[i][j] + 
                           (i > 0 ? aw[i][j]*x[i-1][j] : 0) + (i < NX-1 ? ae[i][j]*x[i+1][j] : 0) +
                           (j > 0 ? as[i][j]*x[i][j-1] : 0) + (j < NY-1 ? an[i][j]*x[i][j+1] : 0));
                
                v[i][j] = (r - (i > 0 ? Lw[i][j]*v[i-1][j] : 0) - (j > 0 ? Ls[i][j]*v[i][j-1] : 0)) / Ld[i][j];
                max_res = fmax(max_res, fabs(r));
            }
        }

        if (max_res < TOL) break;

        // Backward Substitution: U*dx = v and update x
        for (int j = NY - 1; j >= 0; j--) {
            for (int i = NX - 1; i >= 0; i--) {
                double dx = v[i][j] - (i < NX-1 ? Ue[i][j]*v[i+1][j] : 0) - (j < NY-1 ? Un[i][j]*v[i][j+1] : 0);
                x[i][j] += dx;
                v[i][j] = dx; // Reuse v to store dx for next pass
            }
        }
    }
}
