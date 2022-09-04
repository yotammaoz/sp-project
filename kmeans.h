#ifndef FINAL_KMEANS_H
#define FINAL_KMEANS_H

void update_centroids(int d, int k, double **new_centroids, int* clusters_size);
int check_convergence(int d, int k, double epsilon, double **centroids, double **new_centroids);
void assign_closest_cluster(int d, int k, double *data_point, double **centroids, double **new_centroids, int *clusters_size);
int find_closest_cluster(int d, int k, double *data_point, double **centroids);
double** runAlg(int d, int k, int n, int max_iter, double epsilon, double **data_points, double **centroids,
                double **new_centroids, int* clusters_size);


#endif //FINAL_KMEANS_H
