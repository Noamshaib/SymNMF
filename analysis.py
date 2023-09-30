import sys
import pandas as pd
import numpy as np
import mysymnmfsp as msm
np.random.seed(0)

def defX(filename):
    return pd.read_csv(filename, header=None).values.tolist()

def initial_H(W, n, k):
    W_np = np.array(W)
    m = np.mean(W_np)
    upper_bound = 2 * np.sqrt(m / k)
    H0_np = np.random.uniform(0, upper_bound, size=(n, k))
    H0 = H0_np.tolist()
    return H0

def silhouette_score(centroids, belongX):
    
    return 
    

def belong_point_to_centroid(centroidsH):
    max_column_indices = np.argmax(H, axis=1)
    max_column_indices += 1                         #add 1 because indecies are 0 - k-1 and centroids are 1 - k
    belong_lst = max_column_indices.tolist()
    return belong_lst

k = sys.argv[1]
inputfile = sys.argv[2]

X = defX(inputfile)
dim = len(X[0])
N = len(X)
W = msm.norm(X, N, dim)
H0 = initial_H(W, N, k)
H = msm.symnmf(W, H0, k, N, dim)
belong_lst = belong_point_to_centroid(H)

print("nmf: ", silhouette_score(H, belong_lst))


# H = np.array([[0.2, 0.4, 0.6],
#               [0.1, 0.8, 0.3],
#               [0.5, 0.2, 0.7]])
# print(belong_point_to_centroid(H))
