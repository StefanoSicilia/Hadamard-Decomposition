# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 16:00:32 2025

@author: subha
"""

import numpy as np
from hadinit import Had_init,Had_init2,Had_init3
def Upd_manifold(X, U, G, h):
    n, r = X.shape
    for i in range(n):
        Gi = G[i, :].reshape(r, r)
        Gixi = Gi @ X[i, :].T
        X[i, :] = X[i, :] + h * (-Gi.T @ U[i, :].T + (U[i, :] @ Gixi) * X[i, :].T).T
        X[i, :] = X[i, :] / np.linalg.norm(X[i, :])
        U[i, :] = U[i, :] + h * (-Gixi).T
    return X, U


def face_split(A, B):
    n, p = A.shape
    m, q = B.shape
    if m != n:
        raise ValueError('Incompatible sizes of the matrices.')
    C = np.zeros((n, p * q))
    for j in range(n):
        C[j, :] = np.kron(A[j, :], B[j, :])
    return C


from numpy.linalg import norm, svd

def Had_manBCD(X, opts):
    """
    Computes a(n) (approximate) Hadamard decomposition X=WH', where 
    W=face_split(W1,W2) and H=face_split(H1,H2) have rank r=opts.rank. 
    It uses a 2 block coordinate descend algorithm for W and H and optimizes 
    each rank-r matrix on the manifold such that they admit a face-split 
    decomposition.
    """

    n, m = X.shape
    r = opts['rank']
    maxit = opts['maxit']
    tol = opts['tol']
    extrapar = opts['momentum']
    beta, betat, gamma, gammat, eta = extrapar

    if opts['init'] == 'Wertz':
        M = np.sqrt(np.abs(X))
        N = np.sign(X) * M
        U1, S1, V1t = svd(M, full_matrices=False)
        S1_sqrt = np.sqrt(np.diag(S1[:r]))
        W1 = U1[:, :r] @ S1_sqrt
        H1 = V1t[:r, :].T @ S1_sqrt
        U2, S2, V2t = svd(N, full_matrices=False)
        S2_sqrt = np.sqrt(np.diag(S2[:r]))
        W2 = U2[:, :r] @ S2_sqrt
        H2 = V2t[:r, :].T @ S2_sqrt
    elif opts['init'] == 'FS':
        W1, W2, H1, H2 = Had_init(X, r)
    elif opts['init'] == 'FS2':
        W1, W2, H1, H2 = Had_init2(X, r)
    elif opts['init'] == 'FS3':
        W1, W2, H1, H2 = Had_init3(X, r)
    elif opts['init'] == 'given':
        W1 = opts['X']
        H1 = opts['Y']
        W2 = opts['U']
        H2 = opts['V']
    else:
        raise ValueError('Initialization not available.')

    # normalization of the rows of W1, W2
    for i in range(n):
        normXi = norm(W1[i, :])
        W1[i, :] = W1[i, :] / normXi
        W2[i, :] = W2[i, :] * normXi

    # normalization of the rows of H1, H2
    for i in range(m):
        normYi = norm(H1[i, :])
        H1[i, :] = H1[i, :] / normYi
        H2[i, :] = H2[i, :] * normYi

    W = face_split(W1, W2)
    H = face_split(H1, H2)
    normX = norm(X, 'fro')

    Hblock = opts['Hblock']
    Wblock = opts['Wblock']
    err = np.zeros(maxit * (Hblock + Wblock) + 1)
    err[0] = 0.5 * norm(X - W @ H.T, 'fro') ** 2 / normX
    i = 0
    j = 1
    tau = opts['tau']

    while j <= maxit and err[i] > tol:
        # W-update
        for k in range(Wblock):
            W1old = W1.copy()
            W2old = W2.copy()
            G = (W @ H.T - X) @ H
            alpha = tau / norm(H.T @ H, 2)
            W1, W2 = Upd_manifold(W1, W2, G, alpha)
            W1 = W1 + beta * (W1 - W1old)
            W2 = W2 + beta * (W2 - W2old)
            W = face_split(W1, W2)
            i += 1
            err[i] = 0.5 * norm(X - W @ H.T, 'fro') ** 2 / normX
            if err[i] < err[i - 1]:
                beta = min(betat, gamma * beta)
                betat = min(1, gammat * betat)
            else:
                beta = beta / eta
                betat = beta

        # H-update
        for k in range(Hblock):
            H1old = H1.copy()
            H2old = H2.copy()
            G = (H @ W.T - X.T) @ W
            alpha = tau / norm(W.T @ W, 2)
            H1, H2 = Upd_manifold(H1, H2, G, alpha)
            H1 = H1 + beta * (H1 - H1old)
            H2 = H2 + beta * (H2 - H2old)
            H = face_split(H1, H2)
            i += 1
            err[i] = 0.5 * norm(X - W @ H.T, 'fro') ** 2 / normX
            if err[i] < err[i - 1]:
                beta = min(betat, gamma * beta)
                betat = min(1, gammat * betat)
            else:
                beta = beta / eta
                betat = beta
        j += 1

    return W1, H1, W2, H2, err


