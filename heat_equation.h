# heat_equation.h

#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include <vector>

// Structure to define the heat equation problem
struct HeatEquation {
    int nx;            // Number of grid points in the x direction
    int ny;            // Number of grid points in the y direction
    double dx;         // Grid spacing in the x direction
    double dy;         // Grid spacing in the y direction
    double alpha;      // Thermal diffusivity

    // Constructor to initialize parameters
    HeatEquation(int nx, int ny, double dx, double dy, double alpha) 
        : nx(nx), ny(ny), dx(dx), dy(dy), alpha(alpha) {}
};

// Enum for boundary conditions
enum BoundaryCondition {
    DIRICHLET,
    NEUMANN
};

// Function to apply boundary conditions
void applyBoundaryConditions(std::vector<std::vector<double>>& u, BoundaryCondition bc);

// Function to solve the heat equation given initial conditions and boundary conditions
void solveHeatEquation(HeatEquation& problem, const std::vector<std::vector<double>>& initialCondition, BoundaryCondition bc);

// Function to provide the analytical solution for the 2D Dirac delta initial condition
std::vector<std::vector<double>> analyticalSolution(HeatEquation& problem, double time);

#endif // HEAT_EQUATION_H
