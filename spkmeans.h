#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

double** computeAdjacencyMatrix(double **DPs, double **W, int d, int N);
double** computeDiagonalDegreeMatrix(double **W, double** D, int N);
double** computeNormalizedGraphLaplacian(double **W, double **D, double **Lnorm, int N);
double** computeKmeans(double **DPs, double **centroids, int k, int numOfCords, int numOfDPs);
double** jacobiAlgorithm(double **lnorm, int N, double **eigenVectorsMatrix, double **res);
double** getRotationMatrix(double **A, int N, double **P);
double** transpose(double **matrix, int N);
double** sortVectorsAndValues(double **matrix, double *vector, int N);
double** swap(double **matrix, double *vector, int *indexes, int N, int i, int j);
double** updateNewA(double **newA, double **A, int i, int j, int N, double s, double c);
double** normalizeT(double **eigenVectorsMatrix, int N, int k);
double** mallocDoubleMatrix(int N);
double* getBiggestNumberIndexes(double **A, int N, double *biggestNumberIndexes);
double* mallocDoubleArray(int N);
double computeFi(double **A, int i, int j);
double computeTeta(double fi);
int* turnArrayToZero(int* arr, int length);
int eigengapHeuristic(double* eigenvalues, int N);
int partition(double **matrix, double *vector, int*indexes, int N, int low, int high);
int countingCoordinates(char str[]);
int converged(double **A, double **B, int N);
int getN(char* file_name);
int getD(char* file_name);
int isDiagonal(double **matrix, int N);
double** readDPs(char* file_name, int N, double**DPs);
int isConverged(double **clusters1, double **clusters2, int k, int numOfCords);
void qSort(double **matrix, double *vector, int *indexes, int N, int low, int high);
void copyArr1ToArr2(double **arr1, double **arr2, int k, int numOfCords);
void makeAllZero(double **currClusters, int k , int numOfCords);
void printDiagonal(double **matrix, int N);
double** multiplyMatrices(double **A, double **B, double **C, int N);
void printFinalRes(double **res, int k, int numOfCords);
void printMatrix(double **matrix, int N);
void freeSpace (double **arr, int i);
int checkIfInteger(char n[]);
void doJacobi(double **DPs, int N);
int executeSPKmeans(double **DPs, int N, int numOfCords, char* goal, int k, double **eigenVectorsMatrix);
void callError();
void errorOccured();
void checkMatrixIsntNULL(double **matrix);
void checkArrayIsntNULL(double *array);
void assertNisGreaterThanK(int N, int k);

