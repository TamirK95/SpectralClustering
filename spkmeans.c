#include "spkmeans.h"

int main(int argc, char *argv[]){
    char *goal, *file_name;
    int t, N, k, numOfCords;
    double **DPs, **clusters, **eigenVectorsMatrix, **T;
    /* reading all arguments and check if they are valid: */
    if(argc != 4 || !checkIfInteger(argv[1])){
        callError();
    }
    k=atoi(argv[1]);
    if (k < 0){
        callError();
    }
    goal = argv[2];
    file_name = argv[3];
    N = getN(file_name);
    numOfCords = getD(file_name);
    DPs = mallocDoubleMatrix(N);
    eigenVectorsMatrix = mallocDoubleMatrix(N);
    for (t = 0; t < N; t++){
        DPs[t] = mallocDoubleArray(numOfCords);
        eigenVectorsMatrix[t] = mallocDoubleArray(N);
    }
    DPs = readDPs(file_name, N, DPs);
    /*now DPs contains all tha data points.*/
    k = executeSPKmeans(DPs, N, numOfCords, goal, k, eigenVectorsMatrix);
    /*now k is what it should be if it was 0, and eigenVectorsMatrix is updated.*/
    if (!strcmp(goal, "spk")){
        T = mallocDoubleMatrix(N);
        for (t = 0; t < N; t++){
            T[t] = mallocDoubleArray(k);
        }
        copyArr1ToArr2(eigenVectorsMatrix, T, N, k);
        clusters = mallocDoubleMatrix(k);
        for (t = 0; t < k; t++){
            clusters[t] = mallocDoubleArray(k);
        }
        /*take first k points as initial centroids:*/
        copyArr1ToArr2(T, clusters, k, k);
        /*compute clusters by kmeans algorithm:*/
        clusters = computeKmeans(T, clusters, k, k, N);
        printFinalRes(clusters, k, k);
        freeSpace(clusters, k);
        freeSpace(T, N);
    }
    freeSpace(DPs, N);
    freeSpace(eigenVectorsMatrix,N);
    return 1;
}

/*
* Function: getN
* -----------------------------
* Read the data points (=lines) from file and determine how many data points in file.
*
* file_name - name of the file to read from.
*/
int getN(char* file_name){
    int count = 1;
    char c;
    FILE* file = fopen(file_name, "r");
    if (file == NULL){
        errorOccured();
    }
    rewind(file);
    c = fgetc(file);
    while (!feof(file)){
        if (c == '\n'){
            count ++;
        }
        c = fgetc(file);
    }
    fclose(file);
    return count;
}

/*
* Function: getD
* -----------------------------
* Read the first data point from file and determine how many coordinates there are for each data point.
*
* file_name - name of the file to read from.
*/
int getD(char* file_name){
    char* data_point = (char*) malloc(1024 * sizeof(char));
    int count = 0;
    FILE* file = fopen(file_name, "r");
    if (file == NULL){
        errorOccured();
    }
    rewind(file);
    fscanf(file, "%s",data_point);
    count = countingCoordinates(data_point);
    fclose(file);
    free(data_point);
    return count;
}

/*
* Function: readDPs
* -----------------------------
* Read the data points from file and put them in the given array.
*
* file_name - name of the file to read from.
* N - number of datapoints to read.
* DPs - the array that the function will fill with data points.
*/
double** readDPs(char* file_name, int N, double**DPs){
    char* data_point = (char*) malloc(1024 * sizeof(char));
    const char s[2] = ",";
    int i, j;
    char *cord;
    double fcord;
    FILE* file = fopen(file_name, "r");
    if (file == NULL){
        errorOccured();
    }
    rewind(file);
    for(i=0; i < N; i++){
        fscanf(file, "%s", data_point);
        cord = strtok(data_point, s);
        j = 0;
        while (cord!=NULL){
            fcord = (double)atof(cord);
            DPs[i][j] = fcord;
            j++;
            cord = strtok(NULL,s);
        }
    }
    fclose(file);
    free(data_point);
    return DPs;
}

/*
* Function: computeKmeans
* -----------------------------
* Return final k centroids, according to values provided.
*
* DPs - all the data points.
* centroids - the first k centroids for the algorithm.
* k - number of centroids.
*/
double** computeKmeans(double **DPs, double **centroids, int k, 
int numOfCords, int numOfDPs){
    int t, a, b, iterCounter=0, maxIter = 300, *countersArr;
    double **currClusters, **oldClusters;
    currClusters = mallocDoubleMatrix(k);
    oldClusters = mallocDoubleMatrix(k);
    for (t = 0; t < k; t++){
        currClusters[t] = mallocDoubleArray(numOfCords);
        oldClusters[t] = mallocDoubleArray(numOfCords);
    }
    copyArr1ToArr2(centroids, currClusters, k, numOfCords);
    makeAllZero(oldClusters, k, numOfCords);
    /*now to compute the exact clusters, by the algo*/
    while((!isConverged(currClusters, oldClusters, k, numOfCords)) && 
    (iterCounter<maxIter)){
        copyArr1ToArr2(currClusters, oldClusters, k, numOfCords);
        iterCounter++;
        countersArr = (int*) malloc(k * sizeof(int));
        if (countersArr == NULL){
            errorOccured();
        }
        countersArr = turnArrayToZero(countersArr, k);
        makeAllZero(currClusters, k , numOfCords);
        for(t=0; t<numOfDPs; t++){
            double *DPoint = DPs[t];  /*the point ill search min for*/
            double minDistance=__FLT_MAX__;
            int closestClusterIndx = 0;
            for(a=0; a<k; a++){
                double *cluster = oldClusters[a];
                double minCandidate = 0;
                for(b=0; b<numOfCords; b++){
                    double x = DPoint[b]-cluster[b];
                    minCandidate+=x*x;
                }
                if(minCandidate<minDistance){ /*found new minimum*/
                    closestClusterIndx = a;
                    minDistance = minCandidate;
                }
            }
            countersArr[closestClusterIndx]++;
            for(a=0; a<numOfCords;a++){
                currClusters[closestClusterIndx][a]+=DPs[t][a];
            }
        }
        /*now to divide by |s_i|*/
        for(t=0;t<k;t++){
            for(a=0;a<numOfCords;a++){
                if(countersArr[t] == 0){ /*zero division*/
                    errorOccured();
                }
                currClusters[t][a] = (currClusters[t][a]/countersArr[t]);
            }
        }
        free(countersArr);
    }
    freeSpace(oldClusters,k);
    copyArr1ToArr2(currClusters, centroids, k, numOfCords);
    freeSpace(currClusters,k);
    return centroids;
}

/*
* Function: computeAdjacencyMatrix
* -----------------------------
* Returns the adjacency matrix according to the provided data points.
*
* DPs - the data points.
* W - some N*N matrix. will fill the correct values of W inside it.
* N - number of data points provided.
* d- number of coordinates in each data point.
*/
double** computeAdjacencyMatrix(double **DPs, double **W, int d, int N){
    int i,j,k;
    double w;
    for (i=0; i < N; i++){
        for(j=i; j < N; j++)
        {
            w = 0;
            for (k = 0; k < d; k++){ /*sum row elements*/
                w += pow(DPs[i][k] - DPs[j][k], 2);
            }
            w = sqrt(w);
            w *= (-0.5);
            w = exp(w);
            if(i == j){
                w = 0;
            }
            W[i][j] = w;
            W[j][i] = w;
        }
    }
    return W;
}

/*
* Function: printDiagonal
* -----------------------------
* Printing the diagonal values of the given matrix.
*
* matrix - some N*N matrix.
* N - dimensions of the matrix.
*/
void printDiagonal(double **matrix, int N)
{
    int i;
    double temp;
    for (i=0; i < N; i++){
        temp = matrix[i][i];
        if (temp < 0 && temp > -0.00005){
            temp = 0;
        }
        if (i < N-1){
            printf("%.4f,", temp);
        }
        else{
            printf("%.4f", temp);
        } 
    }
    printf("\n");
}

/*
* Function: computeDiagonalDegreeMatrix
* -----------------------------
* Returns the diagonal degree matrix according to the provided data points.
*
* W - the adjacency matrix.
* D - some N*N matrix. will fill the correct values of D inside it.
* N - number of data points provided.
*/
double** computeDiagonalDegreeMatrix(double **W, double** D, int N){
    int i, j, k;
    double sum;
    for (i=0; i < N; i++){
        for(j=0; j < N; j++){
            if (i != j){
                D[i][j] = 0;
            }
            else{
                sum = 0;
                for(k=0; k < N; k++){
                    sum += W[i][k];
                }
                D[i][j] = sum;
            }
        }
    }
    return D;
}

/*
* Function: computeNormalizedGraphLaplacian
* -----------------------------
* Returns the normalized graph laplacian according to the provided data points.
*
* W - the adjacency matrix.
* D - the diagonal degree matrix.
* Lnorm - some N*N matrix. will fill the correct values of Lnorm inside it.
* N - number of data points provided.
*/
double** computeNormalizedGraphLaplacian(double **W, double **D, 
double **Lnorm, int N){
    int i, j;
    double **DWD, **DWDcpy, **Dcpy;
    /* turning Dcpy to D^-0.5: */
    Dcpy = mallocDoubleMatrix(N);
    for (i = 0; i < N; i++){
        Dcpy[i] = mallocDoubleArray(N);
    }
    copyArr1ToArr2(D, Dcpy, N, N);
    for (i=0; i < N; i++){
        Dcpy[i][i] = sqrt(Dcpy[i][i]);
        Dcpy[i][i] = 1/Dcpy[i][i];
    }
    DWD = mallocDoubleMatrix(N);
    DWDcpy = mallocDoubleMatrix(N);
    for (i = 0; i < N; i++){
        DWD[i] = mallocDoubleArray(N);
        DWDcpy[i] = mallocDoubleArray(N);
    }
    DWDcpy = multiplyMatrices(Dcpy, W, DWDcpy, N);  /* DWDcpy = (D^-0.5)*W */
    DWD = multiplyMatrices(DWDcpy, Dcpy, DWD, N);   /* DWD = (D^-0.5)*W*(D^-0.5) */
    /* filling Lnorm: */
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                Lnorm[i][j] = 1 - DWD[i][j];
            }
            else{
                Lnorm[i][j] = -DWD[i][j];
            }
        }
    }
    freeSpace(DWD,N);
    freeSpace(DWDcpy,N);
    freeSpace(Dcpy,N);
    return Lnorm;
}

/*
* Function: MultiplyMatrices
* -----------------------------
* Returns the product of the two matrices provided.
*
* A - first matrix.
* B - second matrix.
* C - some N*N matrix. will fill the correct values inside it.
* N - Dimensions.
*/
double** multiplyMatrices(double **A, double **B, double **C, int N)
{
    int i, j, k;
    for(i=0; i < N; i++){
        for(j=0; j < N; j++){
            C[i][j] = 0;
            for(k=0; k < N; k++){
                C[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
    return C;
}

/*
* Function: CountingCoordinates
* -----------------------------
* Returns number of coordinates in given string.
*
* str - the given string.
*/
int countingCoordinates(char str[]){
    int numOfCords = 1;
    char *ch=strchr(str,',');
    while (ch!=NULL){
        numOfCords++;
        ch=strchr(ch+1,',');
    }
    return numOfCords;
}

/*
* Function: jacobiAlgorithm
* -----------------------------
* Returns A', which is a diagonal matrix and it's diagonal contains all eigenvalues of lnorm. it's values will be filled in res.
*
* lnorm - the lnorm matrix, which is a real symmetric matrix.
* N - dimentions of the matrix.
* eigenVectorsMatrix - an empty matrix which at the end of the function will contain all the eigenVectors of lnorm as columns.
*/
double** jacobiAlgorithm(double **lnorm, int N, 
double **eigenVectorsMatrix, double **res){
    int i, j, t, flag = 0, iterationCounter = 0;
    double c, s, fi, teta, **P, *biggestNumberIndexes,**A, **newA;
    A = mallocDoubleMatrix(N);
    for (t = 0; t < N; t++){
        A[t] = mallocDoubleArray(N);
    }
    copyArr1ToArr2(lnorm, A, N, N);
    while (!isDiagonal(A, N)){
        iterationCounter ++;
        P = mallocDoubleMatrix(N);
        for (t = 0; t < N; t++){
            P[t] = mallocDoubleArray(N);
        }
        P = getRotationMatrix(A, N, P);
        /*find index of biggest off-diagonal element:*/
        biggestNumberIndexes = mallocDoubleArray(2); 
        biggestNumberIndexes = getBiggestNumberIndexes(A, 
        N, biggestNumberIndexes);
        i = biggestNumberIndexes[0];
        j = biggestNumberIndexes[1];
        fi = computeFi(A, i, j);
        teta = computeTeta(fi);
        free(biggestNumberIndexes);
        c = 1/(sqrt(pow(teta, 2) + 1));
        s= teta*c;
        /* making A': */
        newA = mallocDoubleMatrix(N);
        for (t = 0; t < N; t++){
            newA[t] = mallocDoubleArray(N);
        }
        copyArr1ToArr2(A, newA, N, N);
        newA = updateNewA(newA, A, i, j, N, s, c);
        /* multiplying P: */
        if (!flag){
            copyArr1ToArr2(P, eigenVectorsMatrix, N, N);
            flag ++;
        }
        else{
            double** temp = mallocDoubleMatrix(N);
            for (i = 0; i < N; i++){
                temp[i] = mallocDoubleArray(N);
            }
            copyArr1ToArr2(eigenVectorsMatrix, temp, N, N);
            eigenVectorsMatrix = multiplyMatrices(temp, P, 
            eigenVectorsMatrix, N);
            freeSpace(temp, N);
        }
        freeSpace(P, N); 
        if(converged(A, newA, N) || (iterationCounter == 100)){ /*finished*/
            copyArr1ToArr2(newA, A, N, N);
            freeSpace(newA, N);
            break;
        }
        else{
            copyArr1ToArr2(newA, A, N, N);
            freeSpace(newA, N);
        }
    }
    copyArr1ToArr2(A, res, N, N);
    freeSpace(A, N);
    return res;
}

/*
* Function: converged
* -----------------------------
* Returns 1 if both matrixes are the same, as defined for jacobi algorithm, else 0.
*
* A - first matrix.
* B - second matrix.
* N - dimentions of the matrices.
*/
int converged(double **A, double **B, int N){
    int i, j;
    double offAsquared = 0, offBsquared = 0;
    for(i=0; i < N; i++){
        for(j=0; j < N; j++){
            if (i != j){
                offAsquared += pow(A[i][j], 2);
                offBsquared += pow(B[i][j], 2);
            }
        }
    }
    if ((offAsquared - offBsquared) <= pow(10, -15)){
        return 1;
    }
    return 0;
}

/*
* Function: isDiagonal
* -----------------------------
* Returns 1 if the given matrix is diagonal, else 0.
*
* matrix - the matrix that I will check if it's diagonal.
* N - dimentions of the matrix.
*/
int isDiagonal(double **matrix, int N){
    int i,j;
    for (i=0; i<N; i++){
        for(j=0; j<N; j++){
            if ((i != j) && (matrix[i][j] != 0)){
                return 0;
            }
        }
    }
    return 1;
}

/*
* Function: transpose
* -----------------------------
* Returns the transpose of the given matrix (and change it in-place).
*
* matrix - some N*N matrix.
* N - dimentions of the matrix.
*/
double** transpose(double **matrix, int N){
    int i, j;
    double temp;
    for(i=0; i<N; i++){
        for(j=i+1; j<N; j++){
            temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
    return matrix;
}

/*
* Function: sortVectorsAndValues
* -----------------------------
* Returns the sorted eigenvectors matrix and the sorted eigenvalues vector inplace.
*
* matrix - some N*N matrix in which each column is a eigen vector.
* vector - some vector of size N, representing the eigen values of the vectors in the matrix.
* N - dimentions of the vector.
*/
double** sortVectorsAndValues(double **matrix, double *vector, int N){
    int i;
    int* indexes = (int*) malloc(N*sizeof(int));
    if (indexes == NULL){
        errorOccured();
    }
    for(i=0; i < N; i++){
        indexes[i] = i;
    }
    qSort(matrix, vector, indexes, N, 0, N-1);
    free(indexes);
    return matrix;
}

/*
* Function: swap
* -----------------------------
* Swapping the i'th and j'th elements in the provided vector, 
* and swapping the i'th and j'th columns in the provided matrix.
*
* matrix - some N*N matrix.
* vector - some vector of size N.
* N - dimention of the vector.
* i - element location.
* j - element location.
*/
double** swap(double **matrix, double *vector, int *indexes, int N, 
int i, int j){
    int t, indx;
    double temp, *iColumnCpy;
    /* update indexes vector */
    indx = indexes[i];
    indexes[i] = indexes[j];
    indexes[j] = indx;
    /* copy i'th column in the matrix to iColumnCpy */
    iColumnCpy = mallocDoubleArray(N);
    for (t=0; t<N; t++){
        iColumnCpy[t] = matrix[t][i];
    }
    /* copy j'th column in the matrix to i'th column in the matrix */
    for (t=0; t<N; t++){
        matrix[t][i] = matrix[t][j];
    }
    /* copy the original i'th column in the matrix to j'th column in the matrix */
    for (t=0; t<N; t++){
        matrix[t][j] = iColumnCpy[t];
    }
    /* swapping the elements in the vector */
    temp = vector[i];
    vector[i] = vector[j];
    vector[j] = temp;
    free(iColumnCpy);
    return matrix;
}

/*
* Function: partition
* -----------------------------
* Doing a partiotion to the given array and keeping
* the matrix ordered accordingly (in-place).
*
* matrix - some N*N matrix.
* vector - some vector of size N.
* N - dimention of the vector.
* low - lower border.
* high - higher border.
*/
int partition(double **matrix, double *vector, int*indexes, int N, 
int low, int high){
    int i, j;
    double pivot = vector[high];
    i = low - 1;
    for (j=low; j<high; j++){
        if(vector[j] < pivot || (vector[j] == pivot && indexes[j] < indexes[high])){
            i ++;
            swap(matrix, vector, indexes, N, i, j);
        }
    }
    swap(matrix, vector, indexes, N, i+1, high);
    return i+1;
}

/*
* Function: qSort
* -----------------------------
* Using quick sort algorithm in order to sort the vector (and matrix's columns accordingly).
*
* matrix - some N*N matrix.
* vector - some vector of size N.
* N - dimention of the vector.
* low - lower border.
* high - higher border.
*/
void qSort(double **matrix, double *vector, int *indexes, int N, 
int low, int high){
    int pivot;
    if (low < high){
        pivot = partition(matrix, vector, indexes, N, low, high);
        qSort(matrix, vector, indexes, N, low, pivot-1);
        qSort(matrix, vector, indexes, N, pivot+1, high);
    }
}

/*
* Function: getRotationMatrix
* -----------------------------
* Returns the rotation matrix P, according to jacobi algorithm and the given matrix A.
*
* A - the matrix that I will return rotation matrix for.
* N - dimentions of the matrix.
* P - an empty N*N matrix, will be the rotation matrix.
*/
double** getRotationMatrix(double **A, int N, double **P){
    int i,j;
    double fi, teta, c, s;
    double* biggestNumberIndexes = mallocDoubleArray(2);
    biggestNumberIndexes = getBiggestNumberIndexes(A, N, biggestNumberIndexes);
    for(i=0; i< N; i++){
        for(j=0; j< N; j++){
            if(i==j){
                P[i][j] = 1;
            }
            else{
                P[i][j] = 0;
            }
        }
    }
    i = biggestNumberIndexes[0];
    j = biggestNumberIndexes[1];
    fi = computeFi(A, i, j);
    teta = computeTeta(fi);
    free(biggestNumberIndexes);
    c = 1/(sqrt(pow(teta, 2) + 1));
    s= teta*c;
    P[i][i] = c;
    P[i][j] = s;
    P[j][j] = c;
    P[j][i] = -s;
    return P;
}

/*
* Function: computeFi
* -----------------------------
* Returns fi, according to the values given.
*
* A - the given matrix.
* i - the row number of the biggest value of A.
* j - the column number of the biggest value of A.
*/
double computeFi(double **A, int i, int j){
    double fi = A[j][j] - A[i][i];
    fi = fi/A[i][j];
    fi *= 0.5;
    return fi;
}

/*
* Function: computeTeta
* -----------------------------
* Returns teta, according to the values given.
*
* fi - the fi value.
*/
double computeTeta(double fi){
    double teta = 1, temp, absFee;
    if(fi < 0){
        teta = -1;
    }
    absFee = fi;
    if (absFee < 0){
        absFee *= -1;
    }
    temp = absFee + sqrt(pow(fi, 2) + 1);
    teta = teta / temp;
    return teta;
}

/*
* Function: getBiggestNumberIndexes
* -----------------------------
* Returns indexes [i,j] of the biggest number in A's upper part as a list (fill them in biggestNumberIndexes).
*
* A - the given matrix.
* N - dimentions of the matrix.
* biggestNumberIndexes - an empty list of length 2. will contain the result indexes.
*/
double* getBiggestNumberIndexes(double **A, int N, 
double *biggestNumberIndexes){
    int i, j;
    double maximum, temp;
    maximum = -1.7 * pow(10, 50);
    for(i=0; i< N; i++){
        for(j=i+1; j< N; j++){
            temp = A[i][j];
            if (temp < 0){
                temp *= -1;
            }
            if (temp > maximum){
                maximum = temp;
                biggestNumberIndexes[0] = i;
                biggestNumberIndexes[1] = j;
            }
        }
    }
    return biggestNumberIndexes;
}

/*
* Function: isConverged
* -----------------------------
* Check if clusters1 equals to clusters2, return True if so, or False if not.
*
* clusters1 - first cluster.
* clusters1 - second cluster.
* k - number of data points in each cluster.
* numOfCords - number of coordinates in each data point.
*/
int isConverged(double **clusters1, double **clusters2, int k, int numOfCords){
    int i;
    int j;
    for(i=0; i<k; i++){
        for(j=0; j<numOfCords; j++){
            if(clusters1[i][j]!=clusters2[i][j]){
                return 0;
            }
        }
    }
    return 1;
}

/*
* Function: eigengapHeuristic
* -----------------------------
* Return k, according to the Eigengap Heuristic.
*
* eigenvalues - a list of all eigen values.
* N - number of eigen values.
*/
int eigengapHeuristic(double* eigenvalues, int N){
    int i, k;
    double temp, maxEigengap;
    k = 0;
    maxEigengap = -1.7 * pow(10, 50);
    for(i=0; i<(int)(floor(N/2)); i++){
        temp = eigenvalues[i+1] - eigenvalues[i];
        if (temp > maxEigengap){
            maxEigengap = temp;
            k = i;
        }
    }
    return k;
}

/*
* Function: copyArr1ToArr2
* -----------------------------
* Copy the contents of arr1 into arr2.
*
* arr1 - first array (the one which we will copy from).
* arr2 - second array (the one which we will copy to).
* k - number of data points in each array (cluster).
* numOfCords - number of coordinates in each data point.
*/
void copyArr1ToArr2(double **arr1, double **arr2, int k, int numOfCords){
    int i;
    int j;
    for(i=0; i<k; i++){
        for(j=0; j<numOfCords; j++){
            arr2[i][j] = arr1[i][j];
        }
    }
}

/*
* Function: makeAllZero
* -----------------------------
* Turn all slots of currClusters to zero.
*
* currClusters - the given cluster which will turn to zero-cluster.
* k - number of data points in currClusters.
* numOfCords - number of coordinates in each data point.
*/
void makeAllZero(double **currClusters, int k , int numOfCords){
    int i;
    int j;
    for(i=0; i<k; i++){
        for(j=0; j<numOfCords; j++){
            currClusters[i][j] = 0;
        }
    }
}

/*
* Function: printFinalRes
* -----------------------------
* Print the coordinates as needed for the exercise.
*
* res - the result cluster.
* k - number of data points in res.
* numOfCords - number of coordinates in each data point.
*/
void printFinalRes(double **res, int k, int numOfCords){
    int t, p;
    double temp;
    for(t=0;t<k;t++){
        for(p=0;p<numOfCords;p++){
            temp = res[t][p];
            if (temp < 0 && temp > -0.00005){
                temp = 0;
            }
            if(p!=numOfCords-1){
                printf("%.4f,",temp);
            }
            else{
                printf("%.4f",temp);
            }
        }
        printf("\n");
    }
}

/*
* Function: freeSpace
* -----------------------------
* Free a given 2D array that was created using malloc, in which wach line was also created using malloc.
*
* arr - the array which will be free.
* i - number of rows in arr.
*/
void freeSpace(double **arr, int i){
    int a;
    for(a=0; a<i; a++){
        free(arr[a]);
    }
    free(arr);
}

/*
* Function: checkIfInteger
* -----------------------------
* Check if n is non-negative integer.
*
* n - parameter that may be an integer
*/
int checkIfInteger(char n[]){
    int i;
    int len = strlen(n);
    for(i=0; i<len; i++ ){
        if((n[i] < '0') || (n[i] > '9')){
            return 0;
        }
    }
    return 1;
}

/*
* Function: callError
* -----------------------------
* Print error message and exit.
*/
void callError(){
    printf("Invalid Input!");
    exit(0);
}

/*
* Function: printMatrix
* -----------------------------
* print the matrix provided.
*
* matrix - the matrix to be printed.
* N - dimensions.
*/
void printMatrix(double **matrix, int N){
    int t, p;
    double temp;
    for(t=0;t<N;t++){
        for(p=0;p<N;p++){
            temp = matrix[t][p];
            if (temp < 0 && temp > -0.00005){
                temp = 0;
            }
            if(p!=N-1){
                printf("%.4f,",temp);
            }
            else{
                printf("%.4f",temp);
            }
        }
        printf("\n");
    }
}

/*
* Function: errorOccured
* -----------------------------
* Prints error statement and exit program.
*
*/
void errorOccured(){
    printf("An Error Has Occured");
    exit(0);
}

/*
* Function: turnArrayToZero
* -----------------------------
* Turn all elements of the array to zero.
*
* arr - the array i will turn to zeros.
* length - arr's length.
*/
int* turnArrayToZero(int* arr, int length){
    int t;
    for(t=0; t<length; t++){
        arr[t]=0;
    }
    return arr;
}

/*
* Function: updateNewA
* -----------------------------
* Update A' (the given matrix) as needed for the iteration.
*
* newA - the matrix to be updated.
* A - previous matrix.
* i - biggest off-diagonal absolute value of matrix element's row.
* j - i - biggest off-diagonal absolute value of matrix element's column.
* N - matrix's dimensions.
* s - the parameter s.
* c - the parameter c. 
*/
double** updateNewA(double **newA, double **A, int i, int j, int N, 
double s, double c){
    int t;
    for(t=0; t<N; t++){
        if (t != i && t != j){
            newA[t][i] = c*A[t][i] - s*A[t][j];
            newA[i][t] = c*A[t][i] - s*A[t][j];
            newA[t][j] = c*A[t][j] + s*A[t][i];
            newA[j][t] = c*A[t][j] + s*A[t][i];
        }
    }
    newA[i][i] = pow(c, 2) * A[i][i] + pow(s, 2) * A[j][j] - 2* s* c* A[i][j];
    newA[j][j] = pow(s, 2) * A[i][i] + pow(c, 2) * A[j][j] + 2* s* c* A[i][j];
    newA[i][j] = 0;
    newA[j][i] = 0;
    return newA;
}

/*
* Function: checkMatrixIsntNULL
* -----------------------------
* Check that a given double** matrix isnt NULL.
*
* matrix - the matrix to be verified. 
*/
void checkMatrixIsntNULL(double **matrix){
    if (matrix == NULL){
        errorOccured();
    }
}

/*
* Function: checkArrayIsntNULL
* -----------------------------
* Check that a given double* array isnt NULL.
*
* array - the array to be verified. 
*/
void checkArrayIsntNULL(double *array){
    if (array == NULL){
        errorOccured();
    }
}

/*
* Function: executeSPKmeans
* -----------------------------
* Make the needed executions according to goal and the data provided. Return the right k, if relevant, or the original k if not.
*
* DPs - all the data points.
* N - number of data points.
* numOfCords - number of coordinations in each data point.
* goal - execution order. 
* eigenVectorsMatrix - an N*N matrix that will eventually contain the N eigenvectors of Lnorm as columns.
*/
int executeSPKmeans(double **DPs, int N, int numOfCords, char *goal, 
int k, double **eigenVectorsMatrix){
    int i;
    double **W, **D, **LN, **eigenValuesMatrix, *eigenValues;
    if ((!strcmp(goal, "wam")) || (!strcmp(goal, "ddg")) || 
    (!strcmp(goal, "lnorm")) || (!strcmp(goal, "spk"))){
        W = mallocDoubleMatrix(N);
        D = mallocDoubleMatrix(N);
        LN = mallocDoubleMatrix(N);
        eigenValuesMatrix = mallocDoubleMatrix(N);
        eigenValues = mallocDoubleArray(N);
        for (i = 0; i < N; i++){
            W[i] = mallocDoubleArray(N);
            D[i] = mallocDoubleArray(N);
            LN[i] = mallocDoubleArray(N);
            eigenValuesMatrix[i] = mallocDoubleArray(N);
        }
        W = computeAdjacencyMatrix(DPs, W, numOfCords, N);
        D = computeDiagonalDegreeMatrix(W, D, N);
        LN = computeNormalizedGraphLaplacian(W, D, LN, N);
        if (!strcmp(goal, "wam")){
            printMatrix(W, N);
        }
        else if (!strcmp(goal, "ddg")){
            printMatrix(D, N);
        }
        else if (!strcmp(goal, "lnorm")){
            printMatrix(LN, N);
        }
        else if (!strcmp(goal, "spk") && k == 0){
            eigenValuesMatrix = jacobiAlgorithm(LN, N, eigenVectorsMatrix, 
            eigenValuesMatrix);
            for(i=0; i<N; i++){
                eigenValues[i] = eigenValuesMatrix[i][i];
            }
            /* sorting the eigenvalues and rearrange the vectors in the matrix accordingly */
            sortVectorsAndValues(eigenVectorsMatrix, eigenValues, N);
            k = eigengapHeuristic(eigenValues, N) + 1;
            assertNisGreaterThanK(N,k);
            /* forming T matrix: */
            eigenVectorsMatrix = normalizeT(eigenVectorsMatrix, N, k);
        }
        freeSpace(W,N);
        freeSpace(D,N);
        freeSpace(LN,N);
        freeSpace(eigenValuesMatrix,N);
        free(eigenValues);
    }
    else if (!strcmp(goal, "jacobi")){
        doJacobi(DPs, N);
    }
    else{
        callError();
    }
    return k;
}

/*
* Function: assertNisGreaterThanK
* -----------------------------
* Make sure that N is greater than k.
*
* N - integer.
* k - integer.
*/
void assertNisGreaterThanK(int N, int k){
    if(k >= N){
        errorOccured();
    }
}

/*
* Function: doJacobi()
* -----------------------------
* Execute jacobi and print the result.
*
* DPs - the matrix read.
* N - integer.
*/
void doJacobi(double **DPs, int N){
    int i;
    double **eigenVectorsMatrix, **eigenValuesMatrix;
    eigenVectorsMatrix =  mallocDoubleMatrix(N);
    eigenValuesMatrix =  mallocDoubleMatrix(N);
    for (i = 0; i < N; i++){
        eigenVectorsMatrix[i] =  mallocDoubleArray(N);
        eigenValuesMatrix[i] = mallocDoubleArray(N);
    }
    eigenValuesMatrix = jacobiAlgorithm(DPs, N, eigenVectorsMatrix, 
    eigenValuesMatrix);
    printDiagonal(eigenValuesMatrix, N);
    eigenVectorsMatrix = transpose(eigenVectorsMatrix, N);
    printMatrix(eigenVectorsMatrix, N);
    freeSpace(eigenValuesMatrix, N);
    freeSpace(eigenVectorsMatrix, N);
}

/*
* Function: normalizeT()
* -----------------------------
* Normalize matrix containing eigenvectors as columns.
*
* DPs - the matrix read.
* N - integer.
* k - integer.
*/
double** normalizeT(double **eigenVectorsMatrix, int N, int k){
    double** T, sum;
    int i, j, t;
    T = mallocDoubleMatrix(N);
    for (i = 0; i < N; i++){
        T[i] = mallocDoubleArray(k);
    }
    for (i = 0; i < N; i++){
        sum = 0;
        for(t=0; t < k; t++){
            sum += pow(eigenVectorsMatrix[i][t],2);
        }
        sum = sqrt(sum);
        for (j=0; j < k; j++){
            T[i][j] = eigenVectorsMatrix[i][j] / sum;
        }
    }
    copyArr1ToArr2(T, eigenVectorsMatrix, N, k);
    freeSpace(T, N);
    return eigenVectorsMatrix;
}

/*
* Function: mallocDoubleMatrix()
* -----------------------------
* Initialize space for a malloc matrix containing double* elements.
*
* N - number of rows in matrix.
*/
double** mallocDoubleMatrix(int N){
    double **matrix = (double **) malloc(N * sizeof(double *));
    checkMatrixIsntNULL(matrix);
    return matrix;
}

/*
* Function: mallocDoubleArray()
* -----------------------------
* Initialize space for a malloc array containing double elements.
*
* N - number of elements in array.
*/
double* mallocDoubleArray(int N){
    double *arr = (double *) malloc(N * sizeof(double ));
    checkArrayIsntNULL(arr);
    return arr;
}
