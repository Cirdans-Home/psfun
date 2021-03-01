---
title: Algorithms
author: F. Durastante
---

We describe here the main algorithms implemented for the solution of
\begin{equation}\label{eq:eq_problem_to_solve}
\mathbf{y} = f(A) \mathbf{x}, A \in \mathbb{R}^{n \times n}, \; \operatorname{nnz}(A) = O(n), \; f:\mathbb{R}\rightarrow \mathbb{R},
\end{equation}
and point to the relevant modules and subroutines in the library.

Krylov Methods
--------------

Let \(V_{k}\) be an orthogonal matrix whose columns \(\mathbf{x}_{1},\dots, \mathbf{x}_{k}\)
span an arbitrary Krylov subspace \(\mathcal{W}_{k}(A,\mathbf{x})\) of
dimension \(k\). We obtain an approximation of \(f(A)\mathbf{x}\) by
\begin{equation}
f(A)\mathbf{x} = V_{k}f(V_{k}^{T}AV_{k})V_{k}^{T}\mathbf{x}.
\label{eq:krylov_fA_approx}
\end{equation}

Different methods for the approximation of matrix functions are obtained for
different choices of the projection spaces \(\mathcal{W}_{k}(A,\mathbf{x})\).

Given a set of scalars \(\left\{ \sigma _{1},\dots,\sigma _{k-1}\right\} \subset \overline{\mathbb{C}}\)
in the the extended complex plane \(\overline{\mathbb{C}}\), that are not
eigenvalues of \(A\), let
\begin{equation*}
q_{k-1}(z)=\prod\nolimits_{j=1}^{k-1}(\sigma _{j}-z).
\end{equation*}
The **rational Krylov** subspace of order \(k\) associated with \(A\),
\(\mathbf{x}\) and \(q_{k-1}\) is defined by
\begin{equation*}
\mathcal{Q}_{k}(A,\mathbf{x})=\left[ q_{k-1}(A)\right] ^{-1}\mathcal{K}_{k}(A,\mathbf{x}),
\end{equation*}
where
\begin{equation*}
\mathcal{K}_{k}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x},A\mathbf{x},\ldots ,A^{k-1}\mathbf{x}\}
\end{equation*}
is the standard polynomial Krylov space.

By defining the matrices
\begin{equation*}
C_{j}=\left( \mu _{j}\sigma _{j}A-I\right) (\sigma _{j}I-A)^{-1},
\end{equation*}
where \(\left\{ \mu _{1},\dots, \mu _{k-1}\right\} \subset \overline{\mathbb{C}}\)
are such that \(\sigma _{j}\neq \mu _{j}^{-2}\), it is known that the
rational Krylov space can also be written as follows [[1]](#guttel2013rational)
\begin{equation*}
\mathcal{Q}_{k}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x},C_{1}\mathbf{x},\ldots ,C_{k-1}\cdots C_2 C_{1}\mathbf{x}\}.
\end{equation*}
This general formulation allows to recast most of the classical Krylov
methods in terms of a rational Krylov method with a specific choice of \(\sigma _{j}\)
and \(\mu _{j}\). In particular,

  * the **polynomial Krylov** method in which \(\mathcal{W}_{k}(A,\mathbf{x})=\mathcal{K}_{k}(A,\mathbf{x})\)
    can be recovered by defining \(\mu _{j}=1\) and \(\sigma_{j}=\infty\) for each \(j\).
  * The **extended Krylov** method [[2]](#EKSM1)[[3]](#EKSM2), in which
    \begin{equation*}
    \mathcal{W}_{2k-1}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x}, A^{-1}\mathbf{x}, A\mathbf{x}, \ldots, A^{-(k-1)}\mathbf{x}, A^{k-1}\mathbf{x}\},
    \end{equation*}
    is obtained by setting
    \begin{equation*}
    (\mu _j,\sigma_j)= \left\lbrace \begin{array}{ll} (1,\infty), & \mbox{for  $j$ even},\\ (0,0),  & \mbox{for  $j$  odd.} \end{array}\right.
    \end{equation*}

  * The **shift-and-invert** rational Krylov [[4]](#moret2004rd)[[5]](#van2006preconditioning), where
   \begin{equation*}
   \mathcal{W}_{k}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x},(\sigma I-A)^{-1}\mathbf{x},\ldots ,(\sigma I-A)^{-(k-1)}\mathbf{x}\},
   \end{equation*}
    is defined by taking \(\mu _{j}=0\) and \(\sigma _{j}=\sigma\) for each \(j\).

The PSFUN library contains the implementation of several flavour of these methods
that can be used for the computation of \eqref{eq:krylov_fA_approx}, the field
in the [[psfun_d_krylov]] type represent the options neeeded to for setting
up and applying the different implemented method for a given matrix function
`fun` (represented by an object of type [[psfun_d_serial]]).

The following table has the info on the available methods.

<table style="width:100%">
<tr>
  <th>Method</th>
  <th>Class</th>
  <th>Matrix type</th>
  <th>"kname"</th>
  <th>Source</th>
</tr>
<tr>
<td>Arnoldi</td>
<td>Polynomial</td>
<td>General</td>
<td>"ARNOLDI"</td>
<td><a href="#MR1149094">[6]</a></td>
</tr>
<tr>
<td>Lanczos</td>
<td>Polynomial</td>
<td>Symmetric</td>
<td>"LANCZOS"</td>
<td><a href="#MR1149094">[6]</a></td>
</tr>
</table>

Quadrature Methods
------------------

This module makes use of the matrix function definition based on the *Cauchy integral*:
given a closed contour \(\Gamma\) lying in the region of
analyticity of the function \(f(x)\) and containing the spectrum of \(A\),
\(f(A)\) can be defined as
\begin{equation}
f(A) =\frac{1}{2\pi i} \int_{\Gamma}  f(z) (zI - A)^{-1} dz.
\label{eq:eq_ContourInt_quadrature}
\end{equation}
By applying a quadrature formula on \(N\) points to \eqref{eq:eq_ContourInt_quadrature},
with weights \(\{c_j\}_{j=1}^{N}\) and nodes \(\{\xi_j\}_{j=1}^{N}\),
it is possible to approximate \eqref{eq:eq_problem_to_solve} as
\begin{equation}
\mathbf{y} = f(A)\mathbf{x} \approx \sum_{j=1}^{N} c_j (A+\xi_j I)^{-1}\mathbf{x},
\label{eq:eq_QuadratureFormula}
\end{equation}
that is then computationally equivalent to the solution of \(N\) linear
systems with the same right-hand side.


The construction of the quadrature module is made of several interconnected modules.
The base module is the [[psfun_base_quadrature_mod]], it contains the
base module of which the different quadratures are extensions.

Then the functions for working with the quadrature formula having either
real (`psb_dpk_`) or complex quadrature nodes and weights for
\eqref{eq:eq_QuadratureFormula} are contained in the relative modules

- [[psfun_d_quadrature_mod]]
- [[psfun_z_quadrature_mod]]

These two modules make use of `abstract interface` for both the
[[dquadrule]]/[[zquadrule]] and the generic
function for which we compute the \(f(A)\). This is implemented this way
to permit the user to implement its own quadrature rule and functions. An example
of how this can be achieved is contained in the functions included in the
submodules of [[psfun_z_quadrature_mod]] that implements the three routines
from [7](#MR2421045).

## References
<a id="guttel2013rational">[1]</a> Güttel, Stefan.
Rational Krylov approximation of matrix functions: numerical methods and optimal
pole selection. GAMM-Mitt. 36 (2013), no. 1, 8--31.

<a id="EKSM1">[2]</a> Druskin, Vladimir; Knizhnerman, Leonid. Extended Krylov subspaces: approximation of the matrix square root and related functions. SIAM J. Matrix Anal. Appl. 19 (1998), no. 3, 755--771.

<a id="EKSM2">[3]</a> Knizhnerman, L.; Simoncini, V. A new investigation of the extended Krylov subspace method for matrix function evaluations. Numer. Linear Algebra Appl. 17 (2010), no. 4, 615--638.

<a id="moret2004rd">[4]</a> Moret, I.; Novati, P. RD-rational approximations of the matrix exponential. BIT 44 (2004), no. 3, 595--615.

<a id="van2006preconditioning">[5]</a> van den Eshof, Jasper; Hochbruck, Marlis. Preconditioning Lanczos approximations to the matrix exponential. SIAM J. Sci. Comput. 27 (2006), no. 4, 1438--1457.

<a id="MR1149094">[6]</a> Saad, Y. Analysis of some Krylov subspace approximations to the matrix exponential operator. SIAM J. Numer. Anal. 29 (1992), no. 1, 209--228.

<a id="MR2421045" >[7]</a> Hale, Nicholas; Higham, Nicholas J.; Trefethen, Lloyd N. Computing \(\mathbf{A}^\alpha,\ \log(\mathbf{A})\), and related matrix functions by contour integrals. SIAM J. Numer. Anal. 46 (2008), no. 5, 2505--2523.
