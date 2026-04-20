#ifndef COMPARISON_H
#define COMPARISON_H

#include <cmath>

// Function to compute Mean Absolute Error (MAE)
double compute_mae(const double* numerical, const double* analytical, int size);

// Function to compute Mean Squared Error (MSE)
double compute_mse(const double* numerical, const double* analytical, int size);

// Function to compute Root Mean Squared Error (RMSE)
double compute_rmse(const double* numerical, const double* analytical, int size);

// Function to compute Relative Error
double compute_relative_error(const double* numerical, const double* analytical, int size);

#endif // COMPARISON_H
