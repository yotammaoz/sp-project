#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "useful.h"
#include "matrix_calculations.h"
#include "jacobi.h"

#define MIN(X,Y) (((X)<(Y))?(X):(Y)) /* macros for min and max */
#define MAX(X,Y) (((X)<(Y))?(Y):(X))

const int MAX_ROTATION = 100;
const double EPSILON = 0.00001;


/* preforms the jacobi algorithm on the nXn real symmetric matrix - matrix */
double **jacobi_alg(int n, double **matrix)
{
    int i;
    int j;
    double PrevOff;
    double Off;
    int iter;
    int *indices;
    double *c_and_s;
    double c;
    double s;
    double **working_mat;
    double **tmp_working_mat;
    double **eigen_vectors_mat;
    double **tmp_eigen_vectors_mat;
    double **rotation_mat;
    double **res;


    working_mat = allocate_memory_array_of_points(n,n); /* copy matrix into matrix working_mat */
    if (working_mat==NULL)
        error();
    
    copy_matrix_into_another(n, matrix, working_mat);
    

    eigen_vectors_mat  = getUnitMatrix(n); /* set the matrix for the eigenvectors to I */
    
    iter=0;

    PrevOff = calcOff(working_mat,n); /* calculates the off-diagonal sum of squares of working_mat */
    Off = 0;

    while ((iter<MAX_ROTATION)&&(PrevOff-Off>=EPSILON)) 
    /* while the algorithm did not converge or do more than 100 iterations */
    {
        indices = getIndicesOfLargestAbsoluteValueInOffDiagonal(n,working_mat); /* find largest off-diagonal element of working_mat*/
        i = MIN(indices[0],indices[1]); /* choose largest element in top part of working_mat - symmetric so can be done */
        j = MAX(indices[0],indices[1]); 
        free(indices);


        c_and_s = get_c_and_s_at_ij(n,working_mat,i,j); /* compute c and s as in the algorithm description */
        c = c_and_s[0];
        s = c_and_s[1];
        free(c_and_s);


        rotation_mat = get_rotation_mat(n,i,j,c,s); /* get the rotation mat for the element and c&s */

        tmp_working_mat = multipleFromBothSides(n,working_mat,rotation_mat); /* working_mat = rot_mat^T*working_mat*rot_mat */

        copy_matrix_into_another(n,tmp_working_mat,working_mat);

        free_matrix(tmp_working_mat);

        tmp_eigen_vectors_mat = multiplyMatrix(n, eigen_vectors_mat, rotation_mat); /* eigen_vectors_mat=eigen_vectors_mat*rot_mat */

        copy_matrix_into_another(n,tmp_eigen_vectors_mat,eigen_vectors_mat);

        free_matrix(tmp_eigen_vectors_mat);

        free_matrix(rotation_mat); /* free the rotation matrix rotation_mat */

        if (iter==0) /* if in the first iteration */
        {
            Off = calcOff(working_mat,n);
        }
        else /* not in first iteration - compute Off() of new working_mat */
        {
            PrevOff = Off;
            Off = calcOff(working_mat,n);
        }

        iter++;
    }

    res = allocate_memory_array_of_points(n,n+1); /* make space for eigen values-vectors matrix */

    if (res==NULL)
        error();

    for (i=0; i<n+1; i++)
    {
        for (j=0; j<n; j++)
        {
            if (i==0)
                res[i][j] = working_mat[j][j]; /* first line is eigen values */
            else
                res[i][j] = eigen_vectors_mat[i-1][j]; /* second is eigen vectors */
        }
    }

    free_matrix(working_mat); /* free working matrix */
    free_matrix(eigen_vectors_mat); /* free eigen matrix */

    return res; /* return eigen values-vectors matrix */

}



/* calculate the sum of squares of the off-diagonal entries in nXn matrix - matrix */
double calcOff(double **matrix, int n) 
{
    double sum = 0;
    int i,j;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if (i!=j) {
                sum += pow(matrix[i][j], 2);
            }
        }
    }
    return sum;
}

/* calculate c and s for nXn matrix as in the algorithm description
   when choosing element matrix[i][j] for computations */
double *get_c_and_s_at_ij(int n,double **matrix, int i, int j)
{
    double theta;
    double t;
    double c;
    double s;
    double *res;

    if ((i>=n)||(j>=n))
        error(); /* matrix is nXn so we cant have i>=n or j>=n */

    theta = (matrix[j][j]-matrix[i][i])/(2*matrix[i][j]);
    t = sign(theta)/(fabs(theta) + sqrt(pow(theta,2)+1));
    c = 1/sqrt(pow(t,2)+1);
    s = t*c;

    res = calloc(2, sizeof(double));
    if (res==NULL)
        error();
    
    res[0]=c;
    res[1]=s;
    return res;
}

/* assuming j>i, calculate the nXn rotation matrix at entry i,j with c and s */
double **get_rotation_mat(int n, int i, int j, double c, double s) 
{
    double **res = getUnitMatrix(n);
    res[i][i] = c;
    res[j][j] = c;
    res[i][j] = s;
    res[j][i] = -s;
    return res;
}

