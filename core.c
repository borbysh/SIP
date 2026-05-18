#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "heat.h" 


#define MAX_ITER 1000
#define TOL 1e-1
#define alpha 0.94

/* Update the temperature values using five-point stencil with the strongly implicit procedure (Stone's method)*/ 
void evolve(field *curr, double a, double dt) {
    int i, j;
    double dx2, dy2;
    dx2 = curr->dx * curr->dx;
    dy2 = curr->dy * curr->dy; 
    
    int NX = curr->nx + 1;
    int NY = curr->ny + 1;
    
    //double **ac, **aw, **ae, **as, **an;
    double **Lw, **Ls, **Lc, **Ue, **Un;
    double **dlt, **v, **temp;
    
/*
    ac = malloc_2d(NX, NY);
    aw = malloc_2d(NX, NY);
    ae = malloc_2d(NX, NY);
    as = malloc_2d(NX, NY);
    an = malloc_2d(NX, NY);
*/
    Lw = malloc_2d(NX, NY);
    Ls = malloc_2d(NX, NY);
    Lc = malloc_2d(NX, NY);
    Ue = malloc_2d(NX, NY);
    Un = malloc_2d(NX, NY);

    dlt = malloc_2d(NX, NY);
    v = malloc_2d(NX, NY);

// Initialize coefficient matrices (5-point stencil for Laplacian):
    // ac*x[i][j] + aw*x[i-1][j] + ae*x[i+1][j] + as*x[i][j-1] + an*x[i][j+1] = b[i][j]
    double ac = 2.0 * a * dt * ((1.0 / dx2) + (1.0 / dy2)) + 1.0;
    double aw = -a * dt / dx2;
    double ae = -a * dt / dx2;
    double as = -a * dt / dy2;
    double an = -a * dt / dy2;

    temp = curr->data;

    //#pragma omp for private(i, j)
    /*
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            ac[i][j] = 2.0 * a * dt * ((1.0 / dx2) + (1.0 / dy2)) + 1.0;
            aw[i][j] = (i > 0) ? -a*dt/dx2 : 0.0;
            ae[i][j] = (i < NX - 1) ? -a*dt/dx2 : 0.0;
            as[i][j] = (j > 0) ? - a * dt / dy2 : 0.0;
            an[i][j] = (j < NY - 1) ? - a * dt / dy2 : 0.0;
        }
    }
    */
    
    // 1. Incomplete LU Factorization (Stone's Method)
    //#pragma omp single
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            double lw = (i > 0) ? aw / (1.0 + alpha * Un[i - 1][j]) : 0.0;
            double ls = (j > 0) ? as / (1.0 + alpha * Ue[i][j - 1]) : 0.0;

            double Un_i_minus = (i > 0) ? Un[i - 1][j] : 0.0;
            double Ue_i_minus = (i > 0) ? Ue[i - 1][j] : 0.0;

            double Un_j_minus = (j > 0) ? Un[i][j - 1] : 0.0;
            double Ue_j_minus = (j > 0) ? Ue[i][j - 1] : 0.0; 

            Lw[i][j] = lw;
            Ls[i][j] = ls;
            Lc[i][j] = ac + alpha * (lw * Un_i_minus + ls * Ue_j_minus) - lw * Ue_i_minus - ls * Un_j_minus;

            if (fabs(Lc[i][j]) < 1e-14) {
                fprintf(stderr, "Warning: Near-zero diagonal element at (%d, %d)\n", i, j);
                Lc[i][j] = 1e-14;
            }

            Ue[i][j] = (ae - alpha * lw * Un_i_minus) / Lc[i][j];
            Un[i][j] = (an - alpha * ls * Ue_j_minus) / Lc[i][j];

        }
    }

    // 2. Initialize v and dx 
    // Calculate Residual and Forward Substitution: L*v = r
    for (int j = 1; j < NY; j++) {
        for (int i = 1; i < NX; i++) {
            double r = curr->data[i][j] - (ac * temp[i][j] +
                (i > 0 ? aw * temp[i - 1][j] : 0.0) + (i < NX ? ae * temp[i + 1][j] : 0.0) +
                (j > 0 ? as * temp[i][j - 1] : 0.0) + (j < NY ? an * temp[i][j + 1] : 0.0));

            v[i][j] = (r - (i > 1 ? Lw[i][j] * v[i - 1][j] : 0.0) - (j > 1 ? Ls[i][j] * v[i][j - 1] : 0.0)) / Lc[i][j];
            
        }
    }
        
    //#pragma omp single
    // Backward Substitution: U*dx = v
    for (int j = NY - 1; j > 0; j--) {
        for (int i = NX - 1; i > 0; i--) {
            dlt[i][j] = v[i][j] - (i < NX - 1 ? Ue[i][j] * dlt[i + 1][j] : 0.0) - (j < NY - 1 ? Un[i][j] * dlt[i][j + 1] : 0.0);
            //temp[i][j] += dlt[i][j];

        }
    }
    /* Determine the temperature field at next time step
     * As we have fixed boundary conditions, the outermost gridpoints
     * are not updated. */
    // 3. Iterative Loop
    for (int it = 1; it < MAX_ITER; it++) {
        double max_res_1 = 0.0;
        double max_res_2 = 0.0;

        //#pragma omp single
        // 3.1 Process black cells
        #pragma omp for private(j, i)
        for (j = 1; j < NY; j++) {
            //#pragma omp parallel for
            //#pragma omp for private(i)
            for (i = (j % 2) + 1; i < NX; i += 2) {
                // Calculate Residual and Forward Substitution: L*v = r
                double r = curr->data[i][j] - (ac * temp[i][j] +
                    (i > 0 ? aw * temp[i - 1][j] : 0.0) + (i < NX ? ae * temp[i + 1][j] : 0.0) +
                    (j > 0 ? as * temp[i][j - 1] : 0.0) + (j < NY ? an * temp[i][j + 1] : 0.0));

                v[i][j] = (r - (i > 1 ? Lw[i][j] * v[i - 1][j] : 0.0) - (j > 1 ? Ls[i][j] * v[i][j - 1] : 0.0)) / Lc[i][j];
                max_res_1 = fmax(max_res_1, fabs(r));
                // Backward Substitution: U*dx = v and update x
                dlt[i][j] = v[i][j] - (i < NX - 1 ? Ue[i][j] * dlt[i + 1][j] : 0.0) - (j < NY - 1 ? Un[i][j] * dlt[i][j + 1] : 0.0);
                temp[i][j] += dlt[i][j];

            }
        }
        // 3.2 Process white cells
        #pragma omp for private(j, i)
        for (j = 1; j < NY; j++) {
            //#pragma omp parallel for
            //#pragma omp for private(i)
            for (i = ((j + 1) % 2) + 1; i < NX; i += 2) {
                // Calculate Residual and Forward Substitution: L*v = r
                double r = curr->data[i][j] - (ac * temp[i][j] +
                    (i > 0 ? aw * temp[i - 1][j] : 0.0) + (i < NX ? ae * temp[i + 1][j] : 0.0) +
                    (j > 0 ? as * temp[i][j - 1] : 0.0) + (j < NY ? an * temp[i][j + 1] : 0.0));

                v[i][j] = (r - (i > 1 ? Lw[i][j] * v[i - 1][j] : 0.0) - (j > 1 ? Ls[i][j] * v[i][j - 1] : 0.0)) / Lc[i][j];
                max_res_2 = fmax(max_res_2, fabs(r));
                // Backward Substitution: U*dx = v and update x
                dlt[i][j] = v[i][j] - (i < NX - 1 ? Ue[i][j] * dlt[i + 1][j] : 0.0) - (j < NY - 1 ? Un[i][j] * dlt[i][j + 1] : 0.0);
                temp[i][j] += dlt[i][j];

            }
        }

        if (fmax(max_res_1, max_res_2) < TOL) break;
        
           
    }

    curr->data = temp;
    
    /*
    free_2d(ac);
    free_2d(ae);
    free_2d(aw);
    free_2d(an);
    free_2d(as);
    */
    free_2d(Lw);
    free_2d(Ls);
    free_2d(Lc);
    free_2d(Ue);
    free_2d(Un);

    free_2d(dlt);
    free_2d(v);
    //free_2d(temp);
}

/* Update the temperature values using five-point stencil */
void evolve_not(field* curr, double a, double dt)
{
    int i, j;
    double dx2, dy2;
    double max_res;
    double **temp;

    temp = curr->data;

    /* Determine the temperature field at next time step
     * As we have fixed boundary conditions, the outermost gridpoints
     * are not updated. */
    dx2 = curr->dx * curr->dx;
    dy2 = curr->dy * curr->dy;

    for (int it = 0; it < MAX_ITER; it++) {

        max_res = 0.0;
        #pragma omp for private(i, j)
        for (i = 1; i < curr->nx + 1; i++) {
            for (j = 1; j < curr->ny + 1; j++) {
                double r = curr->data[i][j] - temp[i][j] + a * dt *
                    ((temp[i + 1][j] -
                        2.0 * temp[i][j] +
                        temp[i - 1][j]) / dx2 +
                        (temp[i][j + 1] -
                            2.0 * temp[i][j] +
                            temp[i][j - 1]) / dy2);
                temp[i][j] += alpha * r;
                max_res = fmax(max_res, fabs(r));
            }
        }
        if (max_res < TOL) break;
    }
    curr->data = temp;
    //fprintf(stderr, "Numerical misfit: %f\n", max_res);
}
