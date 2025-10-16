# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 16:54:48 2025

@author: subha
"""
import numpy as np
import matplotlib.pyplot as plt
from hadinit import Had_init3,face_split


n = 16
m = 20
r = 4
maxit = 10
err_new = np.zeros(maxit)
err_WVG = np.zeros(maxit)
epsilon = 1e-2

for j in range(maxit):
    np.random.seed(j+1)
    Xtrue = np.random.randint(1, 11, size=(n, r))
    Ytrue = np.random.randint(1, 11, size=(m, r))
    Utrue = np.random.randint(1, 11, size=(n, r))
    Vtrue = np.random.randint(1, 11, size=(m, r))
    E = np.random.randint(1, 11, size=(n, m))
    A = (Xtrue @ Ytrue.T) * (Utrue @ Vtrue.T) + epsilon * E / np.linalg.norm(E, 'fro')
    X, U, Y, V = Had_init3(A, r)
    err_new[j] = np.linalg.norm(A - face_split(X, U) @ face_split(Y, V).T, 'fro') / np.linalg.norm(A, 'fro')
    M = np.sqrt(np.abs(A))
    N = np.sign(A) * M
    U1, S1, V1t = np.linalg.svd(M, full_matrices=False)
    S1 = np.diag(S1)
    W1 = U1[:, :r] @ np.sqrt(S1[:r, :r])
    H1 = V1t.T[:, :r] @ np.sqrt(S1[:r, :r])
    U2, S2, V2t = np.linalg.svd(N, full_matrices=False)
    S2 = np.diag(S2)
    W2 = U2[:, :r] @ np.sqrt(S2[:r, :r])
    H2 = V2t.T[:, :r] @ np.sqrt(S2[:r, :r])
    err_WVG[j] = np.linalg.norm(A - (W1 @ H1.T) * (W2 @ H2.T), 'fro') / np.linalg.norm(A, 'fro')

plt.plot(err_new, 'r-o')
plt.plot(err_WVG, 'b-o')
plt.legend(['New approach', 'Wertz et al.'])
plt.show()

score = np.sum(err_WVG > err_new) / maxit