#include <stdio.h>
#include <stdlib.h>

void write_vtk_header(FILE *file, const char *title, int num_points) {
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "%s\n", title);
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d %d\n", num_points, 1, 1);
    fprintf(file, "ORIGIN 0 0 0\n");
    fprintf(file, "SPACING 1 1 1\n");
}

void write_data(FILE *file, double *data, int num_points) {
    fprintf(file, "POINT_DATA %d\n", num_points);
    fprintf(file, "SCALARS solution float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (int i = 0; i < num_points; i++) {
        fprintf(file, "%f\n", data[i]);
    }
}

void generate_vtk(const char *filename, double *numerical_solution, double *analytical_solution, int num_points) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }
    
    write_vtk_header(file, "VTK output for solutions", num_points);
    write_data(file, numerical_solution, num_points);
    
    fclose(file);
}

int main() {
    int num_points = 100; // Example number
    double numerical_solution[100];
    for (int i = 0; i < num_points; i++) {
        numerical_solution[i] = (double)i; // Sample data
    }
    
    generate_vtk("output.vtk", numerical_solution, NULL, num_points);
    
    return 0;
}