import sys
import pandas as pd
import numpy as np
import mykmeanssp as mkm

if len(sys.argv) == 3:
    goal = argv[1]
    inputfile = argv[2]
else:
    k = argv[1]
    goal = argv[2]
    inputfile = argv[3]

X = defX(inputfile)

if goal == "symnmf":
    W = norm(X, N, dim)
    H0 = initial_H(W)
    H = symnmf(H0, W, N, dim)
    printMatrix(H)

elif goal == "sym":
    A = sym(X, N, dim)
    printMatrix(H)

elif goal == "ddg":
    D = ddg(X, N, dim)
    printMatrix(D)

else:  # goal == "norm"
    W = norm(X, N, dim)
    printMatrix(W)


defX(filename):
    return pd.read_csv(filename, header=None)

def initial_H(W):
    H0 =

    return H0