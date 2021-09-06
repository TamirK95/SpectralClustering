#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject* preFit_capi(PyObject *self, PyObject *args)
{
    double **eigenVectorsMatrix, **DPs, dp;
    PyObject *DataPoints, *vectorsMatrix, *prevDataPoint, *prevCord, *T, *Trow;
    int k, numOfCords, N, t, j;
    char *goal;
    if(!PyArg_ParseTuple(args, "OiisiO", &DataPoints, &N, &numOfCords, &goal, &k, &vectorsMatrix))
    {
        return NULL;
    }
    DPs = (double **) malloc(N * sizeof(double *));
    checkMatrixIsntNULL(DPs);
    eigenVectorsMatrix = (double **) malloc((N+1) * sizeof(double *));
    checkMatrixIsntNULL(eigenVectorsMatrix);
    for (t = 0; t < N; t++)
    {
        DPs[t] = (double *) malloc(numOfCords * sizeof(double));
        checkArrayIsntNULL(DPs[t]);
        eigenVectorsMatrix[t] = (double *) malloc(N * sizeof(double));
        checkArrayIsntNULL(eigenVectorsMatrix[t]);
    }
    eigenVectorsMatrix[N] = (double *) malloc(N * sizeof(double));
    checkArrayIsntNULL(eigenVectorsMatrix[t]);
    /*initialize DPs:*/
    for (t = 0; t < N; t++)
    {
        prevDataPoint = PyList_GetItem(DataPoints, t);
        double *DataPoint = (double *) malloc(numOfCords * sizeof(double));
        checkArrayIsntNULL(DataPoint);
        for (j=0; j < numOfCords; j++)
        {
            prevCord = PyList_GetItem(prevDataPoint, j);
            DataPoint[j] = PyFloat_AsDouble(prevCord);
        }
        DPs[t] = DataPoint;
    }
    k = executeSPKmeans(DPs, N, numOfCords, goal, k, eigenVectorsMatrix);
    eigenVectorsMatrix[N][0] = k;
    T = PyList_New(N+1);
    for (t = 0; t < N+1; t++)
    {
        Trow = PyList_New(N);
        for (j=0; j < N; j++)
        {
            dp = eigenVectorsMatrix[t][j];
            PyList_SetItem(Trow, j, Py_BuildValue("d",dp));
        }
        PyList_SetItem(T, t, Trow);
    }
    return T;
}

static PyObject* fit_capi(PyObject *self, PyObject *args)
{
    PyObject *DataPoints, *initialCentroids, *prevDataPoint, *prevCord, *res, *resRow;
    double **centroids, **DPs, dp;
    int k, N, t, j;
    if(!PyArg_ParseTuple(args, "OOii", &DataPoints, &initialCentroids, &k, &N))
    {
        return NULL;
    }
    DPs = (double **) malloc(N * sizeof(double *));
    checkMatrixIsntNULL(DPs);
    centroids = (double **) malloc(k * sizeof(double *));
    checkMatrixIsntNULL(centroids);
    /*initialize DPs:*/
    for (t = 0; t < N; t++)
    {
        prevDataPoint = PyList_GetItem(DataPoints, t);
        double *DataPoint = (double *) malloc(k * sizeof(double));
        checkArrayIsntNULL(DataPoint);
        for (j=0; j < k; j++)
        {
            prevCord = PyList_GetItem(prevDataPoint, j);
            DataPoint[j] = PyFloat_AsDouble(prevCord);
        }
        DPs[t] = DataPoint;
    }
    /*initialize centroids:*/
    for (t = 0; t < k; t++)
    {
        prevDataPoint = PyList_GetItem(initialCentroids, t);
        double *DataPoint = (double *) malloc(k * sizeof(double));
        checkArrayIsntNULL(DataPoint);
        for (j=0; j < k; j++)
        {
            prevCord = PyList_GetItem(prevDataPoint, j);
            DataPoint[j] = PyFloat_AsDouble(prevCord);
        }
        centroids[t] = DataPoint;
    }
    centroids = computeKmeans(DPs, centroids, k, k, N);
    res = PyList_New(k);
    for (t = 0; t < k; t++)
    {
        resRow = PyList_New(N);
        for (j=0; j < N; j++)
        {
            dp = centroids[t][j];
            PyList_SetItem(resRow, j, Py_BuildValue("d",dp));
        }
        PyList_SetItem(res, t, resRow);
    }
    return res;
}

static PyMethodDef capiMethods[] = {
    {"fit", (PyCFunction) fit_capi, METH_VARARGS,
    PyDoc_STR("list of Data Points and k initial centroids for kmeans algorithm.")},
    {"preFit", (PyCFunction) preFit_capi, METH_VARARGS,
    PyDoc_STR("list of data points (DPs), number of data points(N), numOfCords(d), goal, k, eigenVectorsMatrix")},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    capiMethods
};


PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}
