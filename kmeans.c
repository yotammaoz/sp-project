#define PY_SSIZE_T_CLEAN
#include "useful.h"
#include "kmeans.h"


void clear_new_centroids(int d, int k, double **new_centroids, int *clusters_size) {
    int i, j;
    for (i=0 ; i < k ; i++) {
        clusters_size[i] = 0;
        for (j = 0; j < d; j++)
            new_centroids[i][j] = 0;
    }
}

double ** runAlg(int d, int k, int n, int max_iter, double epsilon, double **data_points, double **centroids,
            double **new_centroids, int* clusters_size){
    int iter, i;
    double **temp;
    for (iter=0 ; iter < max_iter ; iter++) {
        for (i=0; i < n ; i++) {
            assign_closest_cluster(d, k, data_points[i], centroids, new_centroids, clusters_size);
        }
        update_centroids(d, k, new_centroids, clusters_size);
        if (check_convergence(d, k, epsilon, centroids, new_centroids) == 1) {
            return new_centroids;
        }
        temp = new_centroids;
        new_centroids = centroids;
        centroids = temp;
        clear_new_centroids(d, k, new_centroids, clusters_size);
    }
    return centroids;
}


void assign_closest_cluster(int d, int k, double *data_point, double **centroids, double **new_centroids,
                            int *clusters_size) {
    int j, i = find_closest_cluster(d, k, data_point, centroids);
    clusters_size[i]++;
    for (j=0; j < d ; j++) {
        new_centroids[i][j] += data_point[j];
    }
}

int find_closest_cluster(int d, int k, double *data_point, double **centroids) {
    int i, closest_cluster = 0;
    double distance, min_distance = euclidean_norm_powered(d, data_point, centroids[0]);
    for (i=1 ; i < k ; i++) {
        distance = euclidean_norm_powered(d, data_point, centroids[i]);
        if (distance < min_distance) {
            min_distance = distance;
            closest_cluster = i;
        }
    }
    return closest_cluster;
}

void update_centroids(int d, int k, double **new_centroids, int* clusters_size) {
    int i, j;
    for (i=0; i < k ; i++) {
        for (j=0; j < d ; j++) {
            new_centroids[i][j] /= clusters_size[i];
        }
    }
}

int check_convergence(int d, int k, double epsilon, double **centroids, double **new_centroids) {
    int i;
    for (i=0; i < k ; i++) {
        double distance = euclidean_norm(d, centroids[i], new_centroids[i]);
        if (distance >= epsilon) {
            return 0;
        }
    }
    return 1;
}
