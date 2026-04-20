#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "heat_equation.h"

#define N 100            // Grid size
#define T 1000           // Time steps
#define DELTA_T 0.01     // Time increment
#define DX 0.1           // Space increment
#define ALPHA 0.1        // Thermal diffusivity

void initialize(double u[N][N]);
void apply_boundary_conditions(double u[N][N]);
void laplace_solver(double u[N][N], double u_new[N][N]);

void initialize(double u[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            u[i][j] = 0.0;  // Initialize to zero
        }
    }
    // Set Dirac delta in the center
    u[N/2][N/2] = 1.0;  // Example value for Dirac delta spike
}

void apply_boundary_conditions(double u[N][N]) {
    // Apply Dirichlet boundary conditions (u=0)
    for (int i = 0; i < N; i++) {
        u[i][0] = 0.0;
        u[i][N-1] = 0.0;
        u[0][i] = 0.0;
        u[N-1][i] = 0.0;
    }
}

void laplace_solver(double u[N][N], double u_new[N][N]) {
    double factor = (DELTA_T / (DX * DX)) * ALPHA;
    
    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            u_new[i][j] = u[i][j] + factor * 
                          (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4.0 * u[i][j]);
        }
    }
}

void solve_heat_equation(void) {
    double (*u)[N] = (double (*)[N])malloc(sizeof(double) * N * N);
    double (*u_new)[N] = (double (*)[N])malloc(sizeof(double) * N * N);
    
    if (u == NULL || u_new == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return;
    }
    
    initialize(u);
    
    for (int t = 0; t < T; t++) {
        apply_boundary_conditions(u);
        laplace_solver(u, u_new);
        
        // Swap pointers
        double (*temp)[N] = u;
        u = u_new;
        u_new = temp;
    }
    
    free(u);
    free(u_new);
}
