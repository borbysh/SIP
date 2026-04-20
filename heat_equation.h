#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#define N 100
#define T 1000
#define DELTA_T 0.01
#define DX 0.1
#define ALPHA 0.1

// Function to apply boundary conditions
void apply_boundary_conditions(double u[N][N]);

// Function to solve the heat equation given initial conditions
void solve_heat_equation(void);

// Function to initialize the solution field
void initialize(double u[N][N]);

// Function to perform Laplacian discretization
void laplace_solver(double u[N][N], double u_new[N][N]);

#endif // HEAT_EQUATION_H
