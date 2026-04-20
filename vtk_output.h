#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

// Function to generate VTK output files
void generateVTKOutput(const char* filename, const double* temperatureField, int fieldSize);

// Function to initialize VTK settings
void initializeVTKVisualization();

// Function to finalize VTK visualization
void finalizeVTKVisualization();

#endif // VTK_OUTPUT_H
