#ifndef FINAL_JACOBI_H
#define FINAL_JACOBI_H

double **jacobi_alg(int n, double **matrix);
double calcOff(double **matrix, int n);
double *get_c_and_s_at_ij(int n,double **matrix, int i, int j);
double **get_rotation_mat(int n, int i, int j, double c, double s);
#endif
