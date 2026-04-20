#ifndef ANALYTICAL_H
#define ANALYTICAL_H

// Function to compute the fundamental solution (Green's function) for the 2D heat equation
double green_function(double x, double y, double t);

// Function to compute the analytical solution of the 2D heat equation using the Green's function
void analytical_solution(double *result, double x_min, double x_max, double y_min, double y_max, double t, int n);

#endif // ANALYTICAL_H