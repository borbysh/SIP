#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SIP.h"
#include "heat_equation.h"
#include "analytical.h"
#include "comparison.h"
#include "vtk_output.h"

#define GRID_SIZE 100
#define NUM_POINTS (GRID_SIZE * GRID_SIZE)

// Function prototypes
void initialize_simulation(void);
void run_sip_solver(void);
void compute_analytical_solutions(void);
void generate_vtk_files(void);
void print_error_statistics(void);

int main(void) {
    printf("=== SIP Solver for Heat Equation ===\n\n");
    
    // Initialize the simulation
    initialize_simulation();

    // Run the SIP solver
    run_sip_solver();

    // Compute analytical solutions
    compute_analytical_solutions();

    // Generate VTK output for visualization
    generate_vtk_files();

    // Print error statistics
    print_error_statistics();

    printf("\n=== Simulation completed successfully ===\n");
    return 0;
}

void initialize_simulation(void) {
    printf("Initializing simulation parameters...\n");
    initialize_vtk_visualization();
    printf("Initialization completed.\n\n");
}

void run_sip_solver(void) {
    printf("Running SIP solver...\n");
    
    // Allocate memory for the system matrix and solution vectors
    double (*b)[NX] = (double (*)[NX])malloc(sizeof(double) * NX * NY);
    double (*x)[NX] = (double (*)[NX])malloc(sizeof(double) * NX * NY);
    
    if (b == NULL || x == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return;
    }
    
    // Initialize right-hand side and solution
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            b[i][j] = 1.0;  // Example RHS
            x[i][j] = 0.0;  // Initial guess
        }
    }
    
    // Solve the system using SIP
    solve_sip(b, x, 0.9);  // alpha = 0.9 for Stone's method
    
    printf("SIP solver completed.\n\n");
    
    free(b);
    free(x);
}

void compute_analytical_solutions(void) {
    printf("Computing analytical solutions...\n");
    
    double *result = (double *)malloc((101) * (101) * sizeof(double));
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return;
    }
    
    analytical_solution(result, -5.0, 5.0, -5.0, 5.0, 1.0, 100);
    
    printf("Analytical solutions computed.\n\n");
    
    free(result);
}

void generate_vtk_files(void) {
    printf("Generating VTK output...\n");
    
    int num_points = GRID_SIZE * GRID_SIZE;
    double *numerical = (double *)malloc(num_points * sizeof(double));
    double *analytical = (double *)malloc(num_points * sizeof(double));
    
    if (numerical == NULL || analytical == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return;
    }
    
    // Initialize with dummy data for demonstration
    for (int i = 0; i < num_points; i++) {
        numerical[i] = (double)i / num_points;
        analytical[i] = sin((double)i / num_points * 3.14159);
    }
    
    generate_vtk_output("solution.vtk", numerical, analytical, GRID_SIZE, GRID_SIZE);
    
    printf("VTK output generated (solution.vtk).\n\n");
    
    free(numerical);
    free(analytical);
}

void print_error_statistics(void) {
    printf("Computing error statistics...\n");
    
    int num_points = GRID_SIZE * GRID_SIZE;
    double *numerical = (double *)malloc(num_points * sizeof(double));
    double *analytical = (double *)malloc(num_points * sizeof(double));
    
    if (numerical == NULL || analytical == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return;
    }
    
    // Initialize with dummy data for demonstration
    for (int i = 0; i < num_points; i++) {
        numerical[i] = 1.0 + 0.01 * i;
        analytical[i] = 1.0 + 0.009 * i;
    }
    
    printf("\nError Metrics:\n");
    printf("  MAE:        %.6e\n", compute_mae(numerical, analytical, num_points));
    printf("  MSE:        %.6e\n", compute_mse(numerical, analytical, num_points));
    printf("  RMSE:       %.6e\n", compute_rmse(numerical, analytical, num_points));
    printf("  L2 norm:    %.6e\n", l2_norm_error(numerical, analytical, num_points));
    printf("  L-infinity: %.6e\n", l_inf_norm_error(numerical, analytical, num_points));
    printf("  Relative_error: %.6e\n\n", relative_error(numerical, analytical, num_points));
    
    free(numerical);
    free(analytical);
}
