#include "useful.h"
#include "math.h"
#include "stdlib.h"
#include "matrix_calculations.h"

/* creates the weighted adjacency matrix from dXn data_points matrix */
double **createWeightedAdjacencyMatrix(int n, int d, double** data_points) 
{
    double **matrix = allocate_memory_array_of_points(n, n);
    double val;
    int i;
    int j;
    for (i=0; i < n; i++) {
        for (j=i+1; j < n; j++) {
            val = exp(-euclidean_norm(d, data_points[i], data_points[j])/2);
            matrix[i][j] = val;
            matrix[j][i] = val;
        }
    }
    return matrix;
}

/* creates the diagonal degree matrix from nXn matrix weightedMatrix */
double **createDiagonalDegreeMatrix(int n, double ** weightedMatrix) 
{
    double **matrix = allocate_memory_array_of_points(n, n);
    double val;
    int i;
    int j;
    for (i=0; i < n; i++) {
        val = 0;
        for (j=0; j < n; j++) {
            val += weightedMatrix[i][j];
        }
        matrix[i][i] = val;
    }
    return matrix;
}

/* return matrix^(-1/2) for the diagonal nXn matrix - matrix */
double ** getDiagonalMatrixPoweredByMinusHalf(int n, double ** matrix) 
{
    double **resultMatrix = allocate_memory_array_of_points(n, n);
    int i;
    for (i=0; i < n; i++) {
        resultMatrix[i][i] = 1 / sqrt(matrix[i][i]);
    }
    return resultMatrix;
}

/* returns matrixA*matrixB for nXn matrices matrixA,matrixB */
double **multiplyMatrix(int n, double **matrixA, double **matrixB) 
{
    double **matrix = allocate_memory_array_of_points(n, n);
    int i;
    int j;
    int k;
    double val;
    for (i=0; i < n ; i++) {
        for (j=0; j < n ; j++) {
            val = 0;
            for (k=0; k < n; k++) {
                val += matrixA[i][k] * matrixB[k][j];
            }
            matrix[i][j] = val;
        }
    }
    return matrix;
}

/* returns matrixA*matrixB*matrixC for nXn matrices matrixA,matrixB,matrixC */
double **multiply3Matrices(int n, double **matrixA, double **matrixB, double **matrixC) 
{
    double ** tmp = multiplyMatrix(n, matrixA, matrixB);
    double ** result = multiplyMatrix(n, tmp, matrixC);
    free_matrix(tmp);
    return result;
}

/* returns I-matrix for nXn matrix - matrix */
double **subtractIbyMatrix(int n, double **matrix) 
{
    double **result = allocate_memory_array_of_points(n, n);
    int i;
    int j;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if (i==j) {
                result[i][j] = 1;
            }
            result[i][j] -= matrix[i][j];
        }
    }
    return result;
}

/* returns the nXn unit matrix - I */
double **getUnitMatrix(int n) 
{
    double **result = allocate_memory_array_of_points(n, n);
    int i;
    int j;
    for (i=0;i<n;i++)
    {
        for (j=0; j<n; j++)
        {
            if (i==j)
                result[i][j]=1;
            else
                result[i][j]=0;
        }
    }
    return result;
}

/* returns the normalized graph laplacian matrix from nXn matrices diagonalDegreeMatrix,weightedAdjacencyMatrix */
double **createNormalizedGraphLaplacian(int n, double** diagonalDegreeMatrix, double ** weightedAdjacencyMatrix) 
{
    double ** d = getDiagonalMatrixPoweredByMinusHalf(n, diagonalDegreeMatrix);
    double ** tmp = multiply3Matrices(n, d, weightedAdjacencyMatrix, d);
    double ** result = subtractIbyMatrix(n, tmp);
    free_matrix(tmp);
    free_matrix(d);
    return result;
}

/* returns the transpose of the nXm matrix - matrix */
double **transposeMatrix(int n, int m, double **matrix) 
{
    double **result = allocate_memory_array_of_points(n,m);
    int i;
    int j;
    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {
                result[i][j] = matrix[j][i];
        }
    }
    return result;
}

/* returns matrixP^t*matrixA*matrixP for nXn matrices matrixA,matrixP */
double **multipleFromBothSides(int n, double** matrixA, double ** matrixP) 
{
    double **matrixPTransposed = transposeMatrix(n,n, matrixP);
    double **result = multiply3Matrices(n, matrixPTransposed, matrixA, matrixP);
    free_matrix(matrixPTransposed);
    return result;
}

/* copies the from nXn matrix to the into nXn matrix */
void copy_matrix_into_another(int n, double **from, double **into)
{
    int i;
    int j;
    
    for (i=0;i<n;i++)
    {
        for (j=0; j<n; j++)
        {
            into[i][j] = from[i][j];
        }
    }

    return;
}

/* returns the indices of the element in nXn matrix with the largest absolute value off the diagonal */
int *getIndicesOfLargestAbsoluteValueInOffDiagonal(int n, double **matrix) 
{
    int *indices = calloc(2,sizeof(int));
    int i;
    int j;
    if (indices==NULL)
        error();
    indices[0]=0;
    indices[1]=1;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if ((fabs(matrix[indices[0]][indices[1]]) <= fabs(matrix[i][j]))&&(i!=j))
            {
                indices[0] = i;
                indices[1] = j;
            }
        }
    }
    return indices;
}


/* finds the best k based on the Eigengap Heuristic in project description */
int run_eigengap_heuristic(int n, double *eigenvalues)
{
    double *copy_of_eigenvalues;
    int i;
    int res = 1;
    double best_delta = 0;
    double curr_delta = 0;

    copy_of_eigenvalues = allocate_memory_array_of_doubles_of_size(n);
    for (i=0;i<n;i++)
    {
        copy_of_eigenvalues[i] = eigenvalues[i];
    }

    qsort(copy_of_eigenvalues,n,sizeof(double),compare_doubles_reversed); /*sorts eigen values by descending order */

    i=0;
    while (2*(i+1)<n)
    {
        curr_delta = fabs(copy_of_eigenvalues[i]-copy_of_eigenvalues[i+1]);
        if (curr_delta>best_delta)
        {
            best_delta = curr_delta;
            res = i+1;
        }
        i++;
    }

    free(copy_of_eigenvalues);
    
    return res;
}

/* finds the T matrix as in algo description from the result of the jacobi algorithm
   eigen matrix is (n+1)Xn first row is eigen values */
double **get_t_matrix(double **eigen_matrix, int n, int k)
{
    double **tmp_mat;
    double **sorted_eigen_matrix;
    double **k_eigen_matrix; 
    int i;
    int j;
    double row_norm;
    double *zero_vect;

    k_eigen_matrix = allocate_memory_array_of_points(k,n);

    tmp_mat = transposeMatrix(n+1,n,eigen_matrix);

    qsort(tmp_mat,n,sizeof(double*),compare_doubles_vect_reversed); /*sort by eigenvalue */

    sorted_eigen_matrix = transposeMatrix(n,n+1,tmp_mat);
    free_matrix(tmp_mat);

    for (i=0;i<n;i++)
    {
        for (j=0;j<k;j++)
        {
            k_eigen_matrix[i][j] = sorted_eigen_matrix[i+1][j]; /* make k_eigen_matrix be U in the algo description */
        }
    }


    free_matrix(sorted_eigen_matrix);


    zero_vect = allocate_memory_array_of_doubles_of_size(k);
    for (i=0;i<n;i++)
    {
        zero_vect[i]=0;
    }


    for (i=0;i<n;i++)
    {
        row_norm = euclidean_norm(k,k_eigen_matrix[i],zero_vect);

        if(row_norm==0) /* what if we have row of zeros? rami said in forums that this should not happen in tests and we can ignore this case */
            error();

        for (j=0;j<k;j++)
        {
            k_eigen_matrix[i][j]=k_eigen_matrix[i][j]/row_norm; /* normalize rows */
        }
    } 

    free(zero_vect);

    return k_eigen_matrix;
}

/* compares doubles for qsort in reversed order */
int compare_doubles_reversed(const void * first, const void * second)
{
    double f = *(double *)first;
    double s = *(double *)second;
    double r = f-s;
    if (r==0)
        return 0;
    if (r>0)
        return -1;
    else /* r<0 */
        return 1;
}

/* compares vectors for qsort in reverse order based on first entry of vector */
int compare_doubles_vect_reversed(const void * first, const void * second)
{
    double f = **(double **)first;
    double s = **(double **)second;
    double r = f-s;
    if (r==0)
        return 0;
    if (r>0)
        return -1;
    else /* r<0 */
        return 1;
}

