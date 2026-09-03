Given a matrix $X \in \mathbb{R}^{m \times n}$ and two ranks $1<r_1,r_2 \ll \min(m,n) $, the Hadamard Decomosition problem aims to solve 

$$\min_{W_1,H_1,W_2,H_2} ‖X - (W_1H_1^\top) \circ (W_2H_2^\top)‖_F^2=\min_{W_1,H_1,W_2,H_2} ‖X - (W_1\textbullet W_2)(H_1\textbullet H_2)^\top‖_F^2,$$

where $W_1\in \mathbb{R}^{m\times r_1}$, $W_2\in \mathbb{R}^{m\times r_2}$, $H_1\in \mathbb{R}^{n\times r_1}$ and $H_2\in \mathbb{R}^{n\times r_2}$, $\circ$ denotes the entry-wise product and $\textbullet$ denotes the face-splitting product.

This code implements the algorithms from the paper

[GSSV26] N. Gillis, S. Saha, S. Sicilia and A. Vandaele, Manifold-based Algorithms for the Hadamard Decomposition, May 2026. [arXiv](https://arxiv.org/pdf/2605.28980) 

The function HadDec.m implements all the algorithms discussed in the paper and, depending on the option chosen, it selects the method for the decomposition among
1) Manopt: it uses a product manifold for $X_1=W_1H_1^\top$ and $X_2=W_2H_2^\top$ and it optimizes through $\texttt{Manopt}$ software to find an approximation $X\approx X_1\circ X_2$.
2) manBCD: it uses a 2 block coordinate descend algorithm for the rank- $(r_1 r_2)$ matrices $W=W_1\textbullet W_2$ and $H=H_1\textbullet H_2$ and it optimizes on $\mathcal{B}\_{m,r}$ and $\mathcal{B}\_{n,r}$ respectively, where $\mathcal{B}_{m,r}$ is the manifold of $m\times (r_1 r_2)$ matrices that admit a face-split decomposition of size $r=[r_1,r_2]$.
3) projBCD: projected block coordinate descent onto on $\mathcal{B}\_{m,r}$ and $\mathcal{B}\_{n,r}$.
4) BCD: it uses the 4 block coordinate descent algorithm described in the paper [WVG25] S. Wertz, A. Vandaele, N.Gillis, Efficient algorithms for the Hadamard decomposition, IEEE International Workshop on Machine Learning for Signal Processing, 2025. [arXiv](https://arxiv.org/abs/2504.13633) 

In the repository 'tests' you can find the codes for generating the numerical experiments in the paper. The script 'example_HadDec' is a simple script to try the codes for small examples.

Before running any code, please run Install.m to have all paths added. To run the method 'Manopt', the codes need $\texttt{Manopt}$ software by

[BMAS14] N. Boumal, B. Mishra, P.-A. Absil, and R. Sepulchre. Manopt, a Matlab toolbox for optimization on manifolds. The Journal of Machine Learning Research, 15(1):1455–1459, 2014.
