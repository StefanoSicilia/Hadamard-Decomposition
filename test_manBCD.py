# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 16:16:48 2025

@author: subha
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from hadamard_code_manBCD import Had_manBCD
n = 4
m = n
r = int(np.sqrt(n))
eta = 1e-5
nmax = 4

# Example type
example = 'oneseye'
np.random.seed(11)

def eH1e(n, m):
    # Placeholder for the eH1e function, as it's not defined in the Matlab code
    return np.eye(n, m)

if example == 'eye':
    X = eH1e(n, m)
elif example == 'oneseye':
    X = np.ones((n, m)) - np.eye(n, m)
elif example == 'random':
    X = np.random.rand(n, m)
elif example == 'randexact':
    W1 = np.random.randint(1, nmax + 1, size=(n, r))
    H1 = np.random.randint(1, nmax + 1, size=(m, r))
    W2 = np.random.randint(1, nmax + 1, size=(n, r))
    H2 = np.random.randint(1, nmax + 1, size=(m, r))
    X = (W1 @ H1.T) * (W2 @ H2.T)
elif example == 'randexactpert':
    W1 = np.random.randint(1, nmax + 1, size=(n, r))
    H1 = np.random.randint(1, nmax + 1, size=(m, r))
    W2 = np.random.randint(1, nmax + 1, size=(n, r))
    H2 = np.random.randint(1, nmax + 1, size=(m, r))
    X = (W1 @ H1.T) * (W2 @ H2.T) + eta * np.random.randint(1, nmax + 1, size=(n, m))
else:
    raise ValueError('Example type not available.')

X = X / np.linalg.norm(X, 'fro')

block = 10
opts = {
    'rank': r,
    'maxit': 50,
    'init': 'FS',
    'tau': 1,
    'Hblock': block,
    'Wblock': block,
    'tol': 1e-7,
    'momentum': [0, 0, 0, 0, 1]
}

# Placeholder for Had_manBCD function, as it's not defined in the Matlab code

W1_FS, H1_FS, W2_FS, H2_FS, err_FS = Had_manBCD(X, opts)
fin_err_FS = err_FS[-1]

opts['init'] = 'Wertz'
W1_Wertz, H1_Wertz, W2_Wertz, H2_Wertz, err_Wertz = Had_manBCD(X, opts)
fin_err_Wertz = err_Wertz[-1]

opts['init'] = 'FS2'
start_time = time.time()
W1_FS2, H1_FS2, W2_FS2, H2_FS2, err_FS2 = Had_manBCD(X, opts)
t_manifold = time.time() - start_time
fin_err_FS2 = err_FS2[-1]

opts['init'] = 'FS3'
W1_FS3, H1_FS3, W2_FS3, H2_FS3, err_FS3 = Had_manBCD(X, opts)
fin_err_FS3 = err_FS3[-1]

plt.close('all')
lw = 1.3
legendlabel = []

if fin_err_FS < 1e5:
    plt.semilogy(err_FS, 'r-', linewidth=lw)
    legendlabel.append('FS')
if fin_err_Wertz < 1e5:
    plt.semilogy(err_Wertz, 'b-', linewidth=lw)
    legendlabel.append('Wertz')
if fin_err_FS2 < 1e5:
    plt.semilogy(err_FS2, 'g-', linewidth=lw)
    legendlabel.append('FS2')
if fin_err_FS3 < 1e5:
    plt.semilogy(err_FS3, 'm-', linewidth=lw)
    legendlabel.append('FS3')

plt.legend(legendlabel, loc='best')
plt.show()