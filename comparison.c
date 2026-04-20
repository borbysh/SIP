#include <stdio.h>
#include <math.h>

// Function to compute L2 norm error
double l2_norm_error(double *numerical, double *analytical, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(numerical[i] - analytical[i], 2);
    }
    return sqrt(sum / size);
}

// Function to compute L-infinity norm error
double l_inf_norm_error(double *numerical, double *analytical, int size) {
    double max_error = 0.0;
    for (int i = 0; i < size; i++) {
        double error = fabs(numerical[i] - analytical[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}

// Function to compute relative error
double relative_error(double numerical, double analytical) {
    return fabs((numerical - analytical) / analytical);
}

int main() {
    // Example usage (with dummy values)
    double numerical[] = {1.0, 2.0, 3.0};
    double analytical[] = {1.1, 1.9, 3.0};
    int size = sizeof(numerical) / sizeof(numerical[0]);
    
    printf("L2 norm error: %f\n", l2_norm_error(numerical, analytical, size));
    printf("L-infinity norm error: %f\n", l_inf_norm_error(numerical, analytical, size));
    printf("Relative error for first element: %f\n", relative_error(numerical[0], analytical[0]));

    return 0;
}