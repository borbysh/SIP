#include <stdio.h>
#include <math.h>

// Function to compute the fundamental solution (Green's function) for the 2D heat equation
double green_function(double x, double y, double t) {
    double coeff = 1.0 / (4.0 * M_PI * t);
    return coeff * exp(-(x * x + y * y) / (4.0 * t)); 
}

// Function to compute the analytical solution of the 2D heat equation
// with Dirac delta initial condition at the center
void analytical_solution(double *result, double x_min, double x_max, double y_min, double y_max, double t, int n) {
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

int main() {
    int n = 100; // Number of grid points along each axis
    double x_min = -5.0;
    double x_max = 5.0;
    double y_min = -5.0;
    double y_max = 5.0;
    double t = 1.0; // Time
    double *result = (double *)malloc((n + 1) * (n + 1) * sizeof(double));

    analytical_solution(result, x_min, x_max, y_min, y_max, t, n);

    // Printing the results
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%.6f ", result[i * (n + 1) + j]);
        }
        printf("\n");
    }

    free(result);
    return 0;
}