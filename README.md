# Whittaker-Fourier-series-in-N-4-Super-Yang-Mills-Theory

The code searches for the Fourier expansion of the generalized Eisenstein series $\mathcal{E}(r, \alpha, \beta)$ in a specific form. It is complementary to the preprint https://arxiv.org/pdf/2209.09319

The folder "examples" contains .nb files with pre-calculated modes of the generalized Eisenstein series for small values of parameters for $n_1>0$ and $n_2>0$. Using those does not require you to have anything installed rather than Mathematica 12; if you require more solutions to be calculated, please feel free to email us!
 
The code in the repository finds an $(n_1, n_2)$ node in $\mathcal{E}(r, \alpha, \beta)$ for $\alpha = 3/2, \beta = 7/2$, $r = 9$, $n_1>0$ and $n_2>0$ in the Wolfram Mathematica compatible format. This corresponds to the generalized Eisenstein series with the biggest values of parameters in (2.13) in the mentioned article. For other values of parameters and signs of $n_1$ and $n_2$, please edit the main file. At the moment, only the nodes with $n_1 n_2 \neq 0$, $n_1 + n_2 \neq 0$ are being supported.

Conjecturaly, it is possible to find solutions for any $\alpha, \beta \in \mathbb{Z}+1/2$, $\alpha+\beta+r \in 2 \mathbb{Z}$ and $|\alpha-\beta| < r$ (we tested the conjecture for $\alpha < 30$, $\beta < 30$ and $r < 15$).

The code only requires sympy. For the installation instructions, see
https://docs.sympy.org/latest/guides/getting_started/install.html
