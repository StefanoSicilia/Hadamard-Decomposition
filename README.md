Given a matrix $X \in \mathbb{R}^{m \times n}$ and a rank $r \ll \min\{m,n\} $, the Hadamard Decomosition problem aims to solve
$ \min_{W_1,H_1,W_2,H_2} \|X - (W_1H_1^T) \circ (W_2H_2^T)\|_F^2, $
where $W_1,W_2\in \mathbb{R}^{m\times r}$ and $H_1,H_2\in \mathbb{R}^{n\times r}$.

This code implements the algorithms from the paper

[GSSV26] N. Gillis, S. Saha, S. Sicilia and A. Vandaele, Manifold-based approaches for the Hadamard Decomposition of a matrix, January 2026.

The function HadDec.m implements all the algorithms discussed in the paper and, depending on the option selected, it chooses the method for the decomposition.

In the repository 'tests' you can find the codes for generating the numerical experiments in the paper.

Before running any code, please run Install.m to have all paths added. The codes need Manopt directory to run. 
