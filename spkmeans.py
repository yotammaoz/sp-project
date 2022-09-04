import math
import sys
import spkmeansc
import numpy as np
import pandas as pd

WEIGHTED_ADJACENCY_MATRIX = "wam"
DIAGONAL_DEGREE_MATRIX = "ddg"
L_NORM = "lnorm"
JACOBI = "jacobi"
FULL_SPECTRAL = "spk"
GOALS = [WEIGHTED_ADJACENCY_MATRIX, DIAGONAL_DEGREE_MATRIX, L_NORM, JACOBI]
MAX_ITER = 300
EPSILON = 0


# kmeans_pp algorithm start
def kmeans_pp(k, data_points):
    np.random.seed(0)
    indices = data_points[0].to_numpy()
    data_points = data_points.to_numpy()
    data = {indices[i]: data_points[i][1:] for i in range(indices.size)}
    choice = np.random.choice(indices)
    m = {choice: data[choice]}
    for i in range(k-1):
        d = [calc_min_dis(point, m) for point in data.values()]
        sum_d = sum(d)
        p = [x / sum_d for x in d]
        choice = np.random.choice(indices, p=p)
        m[choice] = data[choice]
    return m.keys(), [line.tolist() for line in m.values()]


def calc_min_dis(point, m):
    cur_min = math.inf
    for other_point in m.values():
        sum = 0
        for v, v2 in zip(point, other_point):
            sum += math.pow(v-v2, 2)
        if sum < cur_min:
            cur_min = sum
    return cur_min


def min_distance(data_point, m):
    diff = m.apply(lambda x: x - data_point, axis=1)
    diff.pop(0)
    power = diff.apply(lambda x: pow(x, 2), axis=1)
    distance = power.sum(axis=1)
    return distance.min()

# kmeans_pp algorithm end


def print_result(result):
    for centroid in result:
        print(','.join("%.4f" % data for data in centroid))


def parse_input_file(file_name, k):
    data_points = pd.read_csv(file_name, header=None)
    if k >= data_points.shape[0]:
        raise InvalidInput
    return data_points.values.tolist()


class InvalidInput(Exception):
    pass


def parse_args():
    if len(sys.argv) != 4:
        raise InvalidInput
    k, goal, file_name = sys.argv[1:]
    if not k.isnumeric():
        raise InvalidInput
    k = int(k)
    if k < 0:
        raise InvalidInput
    if not (file_name.endswith('.csv') or file_name.endswith('.txt')):
        raise InvalidInput
    if goal not in GOALS:
        raise InvalidInput
    return k, goal, k


def main():
    k, goal, file_name = parse_args()
    if goal == FULL_SPECTRAL:
        if k == 0:
            k = spkmeansc.get_k(file_name)
        data_points = spkmeansc.get_t_matrix(k, file_name)
        indices, init_centroids = kmeans_pp(k, data_points)
        print(','.join([str(int(i)) for i in indices]))
        data_points.pop(0)
        data_points = data_points.values.tolist()
        result = spkmeansc.fit(k, MAX_ITER, EPSILON, len(data_points[0]), init_centroids, len(data_points), data_points)
        print_result(result)
    else:
        spkmeansc.goal(goal, file_name)


if __name__ == '__main__':
    try:
        main()
    except InvalidInput:
        print('Invalid Input!')
    except:
        print('An Error Has Occurred')