#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "useful.h"
#include "matrix_calculations.h"
#include "jacobi.h"
#include "spkmeans.h"

const char *WAM = "wam";
const char *DDG = "ddg";
const char *LNORM = "lnorm";
const char *JACOBI = "jacobi";
const int MAXDIM = 200;

int main(int argc, char *argv[]) 
{
    char *file_name;
    char *goal;

    if (argc!=3)
    {
        invalid_input(); /* in useful.c */
    }

    goal = argv[1]; 
    file_name = argv[2];

    if ((strcmp(goal,WAM)!=0)&&(strcmp(goal,DDG)!=0)&&(strcmp(goal,LNORM)!=0)&&((strcmp(goal,JACOBI)!=0))) {
        invalid_input(); /* if the goal is not one of the specified options, raise invalid_input */
    }

    run_goal(goal, file_name);
    return 0;
}

struct k_n_and_t_matrix *calc_k_n_and_t(char *file_name, int k)
{
    int dimension;
    int num_of_points;
    double **points_mat;
    double **wam;
    double **ddm;
    double **norm_lap;
    double **jacobi_res;
    struct k_n_and_t_matrix * res;
    dimension = get_dimension_from_file(file_name); /* get the dimension of the points i.e. the num of features */
    num_of_points = get_num_of_points_from_file(file_name); /* get the num of points */
    points_mat = allocate_memory_array_of_points(dimension,num_of_points);
    read_data_from_input_file_to_matrix(points_mat, file_name);
    wam = createWeightedAdjacencyMatrix(num_of_points,dimension,points_mat);
    ddm = createDiagonalDegreeMatrix(num_of_points,wam);
    norm_lap = createNormalizedGraphLaplacian(num_of_points,ddm,wam);
    jacobi_res = jacobi_alg(num_of_points,norm_lap); /* this runs the jacobi algorithm from jacobi.c */
    res =  malloc(sizeof(struct k_n_and_t_matrix));
    if (k == 0) {
        k = run_eigengap_heuristic(num_of_points,jacobi_res[0]);
    }
    res->k=k;
    res->n=num_of_points;
    res->t_matrix = get_t_matrix(jacobi_res, num_of_points, k);
    return res;
}

void run_goal(char *file_name, char *goal)
{
    int dimension;
    int num_of_points;
    double **points_mat;
    double **wam;
    double **ddm;
    double **norm_lap;
    double **jacobi_res;
    dimension = get_dimension_from_file(file_name); /* get the dimension of the points i.e. the num of features */
    num_of_points = get_num_of_points_from_file(file_name); /* get the num of points */

    /*if the goal was jacobi we should have dimension==num_of_points */

    points_mat = allocate_memory_array_of_points(dimension,num_of_points);

    read_data_from_input_file_to_matrix(points_mat, file_name);
    /* reads the input points into points_mat */


    if (strcmp(goal,WAM) == 0)  /* goal was wam */
    {
        wam = createWeightedAdjacencyMatrix(num_of_points,dimension,points_mat);
        print_mat(num_of_points,num_of_points,wam);
        free_matrix(wam);
    }

    if (strcmp(goal,DDG)==0) /* goal was ddg */
    {
        wam = createWeightedAdjacencyMatrix(num_of_points,dimension,points_mat);
        ddm = createDiagonalDegreeMatrix(num_of_points,wam);
        print_mat(num_of_points,num_of_points,ddm);
        free_matrix(wam);
        free_matrix(ddm);
    }

    if (strcmp(goal,LNORM) == 0) /* goal was lnorm */
    {
        wam = createWeightedAdjacencyMatrix(num_of_points,dimension,points_mat);
        ddm = createDiagonalDegreeMatrix(num_of_points,wam);
        norm_lap = createNormalizedGraphLaplacian(num_of_points,ddm,wam);
        print_mat(num_of_points,num_of_points,norm_lap);
        free_matrix(wam);
        free_matrix(ddm);
        free_matrix(norm_lap);

    }

    if (strcmp(goal,JACOBI) == 0) /* goal was jacobi */
    {
        jacobi_res = jacobi_alg(dimension,points_mat); /* this runs the jacobi algorithm from jacobi.c */
        print_mat(num_of_points,num_of_points+1,jacobi_res);
        free_matrix(jacobi_res);
    }

    free_matrix(points_mat);
}


int get_dimension_from_file(char *filename)
/* gets the dimension of the points in file with name==filename */
{
    int res = 1;
    char curr_char;
    FILE *input_file;

    input_file = fopen(filename,"r");

    if (input_file==NULL)
        error();

    while((curr_char=fgetc(input_file))!=EOF)
    {
        if (curr_char == ',')
            res++;
        if (curr_char == '\n')
            break;
    }
    fclose(input_file);
    return res;
    
}


int get_num_of_points_from_file(char *filename)
/* gets the number of points in file with name==filename */
{
    int res = 0;
    char curr_char;
    FILE *input_file;

    input_file = fopen(filename,"r");

    if (input_file==NULL)
        error();

    while((curr_char=fgetc(input_file))!=EOF)
    {
        if (curr_char == '\n')
            res++;
    }
    fclose(input_file);
    return res;
}


void read_data_from_input_file_to_matrix(double **mat, char *filename)
/* reads the data from file with name==filename into matrix mat.
has a buffer that fills up and writes to mat when finding a ',' or '\n' */
{
    char curr_char;
    char *buffer;
    int index = 0;
    int row = 0;
    int column = 0;
    FILE *input_file;

    input_file = fopen(filename,"r"); /* open file */

    if (input_file==NULL)
        error();

    buffer = calloc(50, sizeof(char)); /* set buffer */
    
    if (buffer==NULL)
        error();
     
    memset(buffer,0,50); /* set buffer to 0 */


    while ((curr_char = fgetc(input_file))!=EOF) /* as long as not at the end of the file */
    {
        if ((curr_char!=',')&&(curr_char!='\n')) /* we found a digit or '.' -> write to buffer */
        {
            buffer[index] = curr_char;
            index++;
        }

        else /* char is ',' or '\n' */
        {
            /*buffer[index] = '\0'; mabye non needed? */
            sscanf(buffer,"%lf",&mat[row][column]);
            memset(buffer,0,50); /* reset buffer */
            index = 0;

            if (curr_char==',') /* not at end of point so we need to move 1 column right */
                column++;
            
            if (curr_char == '\n') /* at end of point so we need to move to the start of next row */
            {
                row++;
                column = 0;
            }
        }
    }

    free(buffer);
    fclose(input_file);
    return;
}





