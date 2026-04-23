#ifndef COMPARISON_H
#define COMPARISON_H

// Function to compute Mean Absolute Error (MAE)
double compute_mae(const double *numerical, const double *analytical, int size);

// Function to compute Mean Squared Error (MSE)
double compute_mse(const double *numerical, const double *analytical, int size);

// Function to compute Root Mean Squared Error (RMSE)
double compute_rmse(const double *numerical, const double *analytical, int size);

// Function to compute L2 norm error
double l2_norm_error(double *numerical, double *analytical, int size);

// Function to compute L-infinity norm error
double l_inf_norm_error(double *numerical, double *analytical, int size);

// Function to compute relative error
double relative_error(double *numerical, double *analytical, int size);

#endif // COMPARISON_H
