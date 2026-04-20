#ifndef SIP_H
#define SIP_H

#define NX 50
#define NY 50

// Function to solve the linear system using the Strongly Implicit Procedure
void solve_sip(double b[NX][NY], double x[NX][NY], double alpha);

// Function to initialize coefficient matrices for the SIP solver
void initialize_coefficients(double ap[NX][NY], double aw[NX][NY], double ae[NX][NY], 
                             double as[NX][NY], double an[NX][NY]);

#endif // SIP_H