#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

// Function to write VTK header information
void write_vtk_header(FILE *file, const char *title, int nx, int ny);

// Function to write data to VTK file
void write_data(FILE *file, const double *data, int num_points);

// Function to generate VTK output files
void generate_vtk_output(const char *filename, const double *numerical_solution, 
                        const double *analytical_solution, int nx, int ny);

// Function to initialize VTK settings
void initialize_vtk_visualization(void);

// Function to finalize VTK visualization
void finalize_vtk_visualization(void);

#endif // VTK_OUTPUT_H
