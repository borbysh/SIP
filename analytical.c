#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "analytical.h"

// Function to compute the fundamental solution (Green's function) for the 2D heat equation
double green_function(double x, double y, double t) {
    if (t <= 0.0) return 0.0;
    double coeff = 1.0 / (4.0 * M_PI * t);
    return coeff * exp(-(x * x + y * y) / (4.0 * t)); 
}

// Function to compute the analytical solution of the 2D heat equation
// with Dirac delta initial condition at the center
void analytical_solution(double *result, double x_min, double x_max, double y_min, double y_max, double t, int n) {
    if (result == NULL || n < 0) {
        fprintf(stderr, "Error: Invalid parameters in analytical_solution\n");
        return;
    }
    
    double dx = (x_max - x_min) / n;
    double dy = (y_max - y_min) / n;
    int index = 0;

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            result[index++] = green_function(x, y, t);
        }
    }
}
