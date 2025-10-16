# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 16:18:50 2025

@author: subha
"""

import numpy as np
from numpy.linalg import svd, pinv

def Had_init(A, r):
    # Had_init:
    # It is an initialization for the Hadamard decomposition problem, that is:
    # given a matrix A of size nxm, compute two matrices W and H
    # such that:
    #   - A ~ W*H' where W is nx(r^2) and H is mx(r^2),
    #   - W ~ face_split(X,U) and H ~ face_split(Y,V).
    # The problem is usually not solvable exactly if rank(A) > r^2.
    # The function assumes that min(m,n) > r^2.

    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    S_sqrt = np.sqrt(np.diag(S[:r**2]))
    W = U[:, :r**2] @ S_sqrt
    H = Vt.T[:, :r**2] @ S_sqrt

    n, m = A.shape
    X = np.zeros((n, r))
    Y = np.zeros((m, r))
    U_mat = np.zeros((n, r))
    V = np.zeros((m, r))

    for i in range(n):
        Pi, Sigma, Qi_t = np.linalg.svd(W[i, :].reshape(r, r))
        s = np.sqrt(Sigma[0])
        X[i, :] = (s * Pi[:, 0]).T
        U_mat[i, :] = (s * Qi_t.T[:, 0]).T

    for i in range(m):
        Fi, Theta, Gi_t = np.linalg.svd(H[i, :].reshape(r, r))
        t = np.sqrt(Theta[0])
        Y[i, :] = (t * Fi[:, 0]).T
        V[i, :] = (t * Gi_t.T[:, 0]).T

    return X, U_mat, Y, V

def face_split(Y, V):
    # Assuming face_split stacks horizontally the kronecker products of rows of Y and V
    # This is a placeholder; the actual implementation depends on the original face_split definition
    r = Y.shape[1]
    m = Y.shape[0]
    result = np.zeros((m, r*r))
    for i in range(m):
        result[i, :] = np.kron(Y[i, :], V[i, :])
    return result

 
def Had_init2(A, r):
    # SVD decomposition
    U_svd, S_svd, Vh_svd = np.linalg.svd(A, full_matrices=False)
    S_sqrt = np.sqrt(np.diag(S_svd[:r*r]))
    H = Vh_svd[:r*r, :].T @ S_sqrt

    n, m = A.shape
    X = np.zeros((n, r))
    Y = np.zeros((m, r))
    U = np.zeros((n, r))
    V = np.zeros((m, r))

    for i in range(m):
        F, Theta, Gt = np.linalg.svd(H[i, :].reshape(r, r))
        t = np.sqrt(Theta[0])
        Y[i, :] = t * F[:, 0]
        V[i, :] = t * Gt.T[:, 0]

    M = np.linalg.pinv(face_split(Y, V)) @ A.T
    M = M.T

    for i in range(n):
        P, Sigma, Q_t = np.linalg.svd(M[i, :].reshape(r, r))
        s = np.sqrt(Sigma[0])
        X[i, :] = s * P[:, 0]
        U[i, :] = s * Q_t.T[:, 0]

    return X, U, Y, V



def Had_init3(A, r):
    W, S, _ = svd(A, full_matrices=False)
    W = W[:, :r**2] @ np.sqrt(np.diag(S[:r**2]))
    n, m = A.shape
    X = np.zeros((n, r))
    Y = np.zeros((m, r))
    U = np.zeros((n, r))
    V = np.zeros((m, r))
    for i in range(n):
        Pi, Sigma, Qi = svd(W[i, :].reshape(r, r))
        s = np.sqrt(Sigma[0])
        X[i, :] = (s * Pi[:, 0]).T
        U[i, :] = (s * Qi[:, 0]).T
    M = pinv(face_split(X, U)) @ A
    M = M.T
    for i in range(m):
        Fi, Theta, Gi = svd(M[i, :].reshape(r, r))
        t = np.sqrt(Theta[0])
        Y[i, :] = (t * Fi[:, 0]).T
        V[i, :] = (t * Gi[:, 0]).T
    return X, U, Y, V