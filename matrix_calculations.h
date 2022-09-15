/*
   Created by Noa Caspi on 02/08/2022.
*/

#ifndef FINAL_MATRIX_CALCULATIONS_H
#define FINAL_MATRIX_CALCULATIONS_H

double **createWeightedAdjacencyMatrix(int n, int d, double** data_points);
double **createDiagonalDegreeMatrix(int n, double ** weightedMatrix);
double ** getDiagonalMatrixPoweredByMinusHalf(int n, double ** matrix);
double **multiplyMatrix(int n, double **matrixA, double **matrixB);
double **multiply3Matrices(int n, double **matrixA, double **matrixB, double **matrixC);
double **subtractIbyMatrix(int n, double ** matrix);
double **createNormalizedGraphLaplacian(int n, double** diagonalDegreeMatrix, double ** weightedAdjacencyMatrix);
double **transposeMatrix(int n, int m, double **matrix);
double **multipleFromBothSides(int n, double** matrixA, double ** matrixP);
void copy_matrix_into_another(int n, double **from, double **into);
int *getIndicesOfLargestAbsoluteValueInOffDiagonal(int n, double **matrix);
double *obtainCAndT(int n, double **matrix, int pivotI, int pivotJ);
double **getUnitMatrix(int n);
int run_eigengap_heuristic(int n,double *eigenvalues);
double **get_t_matrix(double **eigen_matrix, int n, int k);
int compare_doubles_reversed(const void * first, const void * second);
int compare_doubles_vect_reversed(const void * first, const void * second);
#endif 
