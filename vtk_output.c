#include <stdio.h>
#include <stdlib.h>
#include "vtk_output.h"

void write_vtk_header(FILE *file, const char *title, int nx, int ny) {
    if (file == NULL || title == NULL) {
        fprintf(stderr, "Error: Invalid parameters in write_vtk_header\n");
        return;
    }
    
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "%s\n", title);
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d %d\n", nx, ny, 1);
    fprintf(file, "ORIGIN 0 0 0\n");
    fprintf(file, "SPACING 1 1 1\n");
}

void write_data(FILE *file, const double *data, int num_points) {
    if (file == NULL || data == NULL || num_points <= 0) {
        fprintf(stderr, "Error: Invalid parameters in write_data\n");
        return;
    }
    
    fprintf(file, "POINT_DATA %d\n", num_points);
    fprintf(file, "SCALARS solution float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    
    for (int i = 0; i < num_points; i++) {
        fprintf(file, "%f\n", data[i]);
    }
}

void generate_vtk_output(const char *filename, const double *numerical_solution, 
                        const double *analytical_solution, int nx, int ny) {
    if (filename == NULL || numerical_solution == NULL) {
        fprintf(stderr, "Error: Invalid parameters in generate_vtk_output\n");
        return;
    }
    
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    
    int num_points = nx * ny;
    write_vtk_header(file, "VTK output for heat equation solutions", nx, ny);
    write_data(file, numerical_solution, num_points);
    
    if (analytical_solution != NULL) {
        fprintf(file, "\nSCALARS analytical_solution float 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        write_data(file, analytical_solution, num_points);
    }
    
    if (fclose(file) != 0) {
        perror("Error closing file");
    }
}

void initialize_vtk_visualization(void) {
    printf("VTK visualization initialized.\n");
}

void finalize_vtk_visualization(void) {
    printf("VTK visualization finalized.\n");
}
