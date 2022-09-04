#ifndef FINAL_USEFUL_H
#define FINAL_USEFUL_H

double euclidean_norm_powered(int d, double *p1, double *p2);
double euclidean_norm(int d, double *p1, double *p2);
int* allocate_memory_array_of_size(int k);
double** allocate_memory_array_of_points(int columns, int rows);
double* allocate_memory_array_of_doubles_of_size(int k);
int sign(double num);
void free_matrix(double **matrix);
void print_mat(int columns, int rows, double **mat);
void invalid_input();
void error();


#endif 
