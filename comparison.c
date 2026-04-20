#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "comparison.h"

// Function to compute Mean Absolute Error (MAE)
double compute_mae(const double *numerical, const double *analytical, int size) {
    if (numerical == NULL || analytical == NULL || size <= 0) {
        fprintf(stderr, "Error: Invalid parameters in compute_mae\n");
        return -1.0;
    }
    
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(numerical[i] - analytical[i]);
    }
    return sum / size;
}

// Function to compute Mean Squared Error (MSE)
double compute_mse(const double *numerical, const double *analytical, int size) {
    if (numerical == NULL || analytical == NULL || size <= 0) {
        fprintf(stderr, "Error: Invalid parameters in compute_mse\n");
        return -1.0;
    }
    
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = numerical[i] - analytical[i];
        sum += diff * diff;
    }
    return sum / size;
}

// Function to compute Root Mean Squared Error (RMSE)
double compute_rmse(const double *numerical, const double *analytical, int size) {
    double mse = compute_mse(numerical, analytical, size);
    if (mse < 0.0) return -1.0;
    return sqrt(mse);
}

// Function to compute L2 norm error
double l2_norm_error(double *numerical, double *analytical, int size) {
    if (numerical == NULL || analytical == NULL || size <= 0) {
        fprintf(stderr, "Error: Invalid parameters in l2_norm_error\n");
        return -1.0;
    }
    
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(numerical[i] - analytical[i], 2);
    }
    return sqrt(sum / size);
}

// Function to compute L-infinity norm error
double l_inf_norm_error(double *numerical, double *analytical, int size) {
    if (numerical == NULL || analytical == NULL || size <= 0) {
        fprintf(stderr, "Error: Invalid parameters in l_inf_norm_error\n");
        return -1.0;
    }
    
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
    if (fabs(analytical) < 1e-14) {
        fprintf(stderr, "Warning: Analytical solution near zero, relative error undefined\n");
        return fabs(numerical);
    }
    return fabs((numerical - analytical) / analytical);
}
