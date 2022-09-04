#include <Python.h>
#include <ctype.h>
#include "useful.h"
#include "kmeans.h"
#include "stdlib.h"
#include "spkmeans.h"

void fill_data_list(int colNum, int rowNum, double **data, PyObject * data_python);
PyObject* save_to_output(int colNum, int rowNum, double **data);
static PyObject* fit(PyObject *self, PyObject *args);
static void goal(PyObject *self, PyObject *args);
static PyObject* get_k_and_t_matrix(PyObject *self, PyObject *args);


static PyObject* fit(PyObject *self, PyObject *args)
{
    int k, max_iter, n, d, *cluster_size;
    double epsilon;
    PyObject * init_centroids_python, * data_points_python, *result;
    double **data_points, **centroids, **new_centroids, **centroids_result;
    if (!PyArg_ParseTuple(args, "iidiOiO", &k, &max_iter, &epsilon, &d, &init_centroids_python, &n, &data_points_python))
        return NULL;
    data_points = allocate_memory_array_of_points(d, n);
    centroids = allocate_memory_array_of_points(d, k);
    new_centroids = allocate_memory_array_of_points(d, k);
    cluster_size = allocate_memory_array_of_size(k);
    fill_data_list(d, n, data_points, data_points_python);
    fill_data_list(d, k, centroids, init_centroids_python);
    save_to_output(d, k, centroids);
    centroids_result = runAlg(d, k, n, max_iter, epsilon, data_points, centroids, new_centroids, cluster_size);
    free_matrix(data_points);
    free(cluster_size);
    result = save_to_output(d, k, centroids_result);
    free_matrix(centroids);
    free_matrix(new_centroids);
    return result;
}

static PyObject * get_k_and_t_matrix(PyObject *self, PyObject *args)
{
    char *file_name;
    int k;
    struct k_n_and_t_matrix * k_n_t;
    PyObject * res, * t_python, * k_python;
    if (!PyArg_ParseTuple(args, "si", &file_name, &k)) {
        return NULL;
    }
    k_n_t = calc_k_n_and_t(file_name, k);
    t_python = save_to_output(k_n_t->k, k_n_t->n, k_n_t->t_matrix);
    res = PyList_New(2);
    k_python = Py_BuildValue("i", k);
    PyList_SetItem(res, 0, k_python);
    PyList_SetItem(res, 1, t_python);
    return res;
}

static void goal(PyObject *self, PyObject *args)
{
    char *file_name, *goal;
    if (!PyArg_ParseTuple(args, "ss", &goal, &file_name))
        return;
    run_goal(file_name, goal);
}

void fill_data_list(int colNum, int rowNum, double **data, PyObject * data_python) {
    int i, j;
    PyObject* inner_list, * py_float;
    for (i = 0 ; i < rowNum ; i++) {
        inner_list = PyList_GetItem(data_python, i);
        for (j = 0; j < colNum; j++) {
            py_float = PyList_GetItem(inner_list, j);
            data[i][j] = PyFloat_AsDouble(py_float);
        }
    }
}

PyObject* save_to_output(int colNum, int rowNum, double **data) {
    int i, j;
    PyObject *val, *inner_list, *result;
    result = PyList_New(rowNum);
    for (i = 0 ; i < rowNum ; i++) {
        inner_list = PyList_New(colNum);
        for (j = 0; j < colNum; j++) {
            val = PyFloat_FromDouble(data[i][j]);
            PyList_SetItem(inner_list, j, val);
        }
        PyList_SetItem(result, i, inner_list);
    }
    return result;
}

static PyMethodDef kmeansMethods[] = {
        {"fit",
                (PyCFunction) fit,
                     METH_VARARGS,
                        PyDoc_STR("kmeans alg")},
        {"goal",
                (PyCFunction) goal,
                METH_VARARGS,
                PyDoc_STR("run specific goal")},
        {"get_k_and_t_matrix",
                (PyCFunction) get_k_and_t_matrix,
                METH_VARARGS,
                PyDoc_STR("get k and T matrix")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansc",
        NULL,
        -1,
        kmeansMethods
};

PyMODINIT_FUNC
PyInit_spkmeansc(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


