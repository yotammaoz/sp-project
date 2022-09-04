#ifndef FINAL_SPKMEANS_H
#define FINAL_SPKMEANS_H

int get_dimension_from_file(char *filename);
int get_num_of_points_from_file(char *filename);
void read_data_from_input_file_to_matrix(double **mat, char *filename);
void print_mat(int d, int n, double **mat);
void run_goal(char *, char *);
struct k_n_and_t_matrix *calc_k_n_and_t(char *file_name, int k);

struct k_n_and_t_matrix {
    int k;
    int n;
    double **t_matrix;
};
#endif
