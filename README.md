# Whittaker-Fourier-series-in-N-4-Super-Yang-Mills-Theory

The code searches for the Fourier expansion of the generalized Eisenstein series $\mathcal{E}(\alpha, \beta, r)$ in a specific form. For the definition of the latter, see (1.6) in https://arxiv.org/abs/2008.02713
 
Conjecturaly, it is possible to find solutions for any $\alpha, \beta \in \mathbb{Z}+1/2$, $\alpha+\beta+r \in 2 \mathbb{Z}$ and $|\alpha-\beta| < r$. 
The output coincides up to a sign with a particular solution for the ordinary differential equation on the $(n_1, n_2)$-component in the Fourier series of $\mathcal{E}(\alpha, \beta, r)$. 
The code as in the repository gives an output for $\alpha = 3/2, beta = 7/2$ and $r = 9$ in the Wolfram Mathematica compatible format -- for other values of parameter and the signs of $n_1$ and $n_2$, please edit the main file. At the moment, only the nodes with $n_1 n_2 \neq 0$, $n_1 + n_2 \neq 0$ are being supported.
