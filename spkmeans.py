import sys
import pandas as pd
import numpy as np
import spkmeansmodule

np.random.seed(0)

def input_error():
    """Print error message and exit in case of invalid input."""
    print("Invalid Input!")
    sys.exit()

def printFinalRes(res, k, numOfCords):
    """print final result to stdo."""
    for t in range(k):
        for p in range (numOfCords):
            temp = res[t][p]
            #temp = np.round(temp, 4)
            if (temp < 0 and temp > -0.00005):
                temp = 0
            if(p!=numOfCords-1):
                #print(temp, ",", sep="" ,end="")
                print("%.4f" % temp, ",", sep="" ,end="")
            else:
                #print(temp, end="\n")
                print("%.4f" % temp, end="\n")

def checkLength(length):
    """assert length of arguments is valid."""
    if length != 4:
        input_error()

def checkNbiggerThanK(N, k):
    """assert N is bigger than k."""
    if k >= N:
        input_error()

def copyMatrix(matrix1, matrix2, N, k):
    """copy contents of matrix2 into matirx1."""
    for i in range(N):
        for j in range(k):
            matrix1[i][j] = matrix2[i][j]

def getTemp(k, dp, mus, j):
    """get temp, which is the value for the current observation."""
    temp = 0
    for i in range(k):
        temp += (dp[i] - mus[j][i]) ** 2
    return temp

def printRandoms():
    """print every element in array randoms."""
    started = True
    for random in randoms:
        if started:
            print(random, end="")
            started = False
        else:
            print(",", random, sep="", end="")
    print("")

length = len(sys.argv)
checkLength(length)
k = int(sys.argv[1])        #read k
goal = sys.argv[2]          #read goal
file_name = sys.argv[3]     #read file's name
with open(file_name) as f:
    first_line = f.readline()
f.close()
if first_line == "" or first_line is None:
    input_error()
first_line_list = first_line.split(",")
d = len(first_line_list)    #determine d
arr = [0]*d
for i in range(d):
    arr[i] = i
DPs = np.genfromtxt(file_name, delimiter=",")
N = DPs.shape[0]            #determine N
checkNbiggerThanK(N, k)
eigenVectorsMatrix = [[0 for i in range(N)] for j in range(N+1)]
DPsList = DPs.tolist()
if d == 1:                  #take care of special case
    newDPs = [[DPs[i]] for i in range(N)]
    DPsList = newDPs
eigenVectorsMatrix = spkmeansmodule.preFit(DPsList, N,
                                           d, goal, k, eigenVectorsMatrix)
k = int(eigenVectorsMatrix[N][0])
eigenVectorsMatrix = eigenVectorsMatrix[:N]
if goal == "spk":   #determine k initial centroids by kmeans++
    mus = [[0 for i in range(k)] for j in range(k)]
    randoms = [0] * k
    rnd = np.random.choice(N)
    randoms[0] = rnd
    T = [[0 for i in range(k)] for j in range(N)]
    copyMatrix(T, eigenVectorsMatrix, N, k)
    mus[0] = T[rnd]
    distances = [-1] * N
    probs = [0] * N
    z = 1
    j = 0
    rndIndex = 1
    while z < k:
        sumOfMins = 0
        t = 0
        for dp in T:
            minimum = distances[t]
            temp = getTemp(k, dp, mus, j)
            if temp < minimum or minimum == -1:
                minimum = temp
                distances[t] = minimum
            sumOfMins += distances[t]
            t += 1
        j += 1
        for x in range(N):
            probs[x] = distances[x] / sumOfMins
        rnd = np.random.choice(N, 1, p=probs)
        rnd = rnd[0]
        randoms[rndIndex] = rnd
        rndIndex += 1
        mus[z] = T[rnd]
        z += 1
    printRandoms()
    mus = spkmeansmodule.fit(T, mus, k, N)
    printFinalRes(mus, k, k)

