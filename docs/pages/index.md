---
title: Matrix Functions
author: F. Durastante
---

Matrix function vector products
-------------------------------

This library is focused on the computation of matrix-function [[1]](#MR2396439) vector products
\begin{equation}\label{eq:eq_problem_to_solve}
\mathbf{y} = f(A) \mathbf{x}, A \in \mathbb{R}^{n \times n}, \; \operatorname{nnz}(A) = O(n), \; f:\mathbb{R}\rightarrow \mathbb{R},
\end{equation}
for large and sparse matrices in a distributed setting. Matrix functions are
ubiquitous in models for applied sciences. They are involved in the solution of
ordinary, partial, and fractional differential equations, systems of coupled
differential equations, hybrid differential-algebraic problems,
equilibrium problems, measures of complex networks, and many others.

To perform the computation in \eqref{eq:eq_problem_to_solve}, we consider here two
main approaches, the first one makes use of a definition based on the *Cauchy integral*
for a matrix function: given a closed contour \(\Gamma\) lying in the region of
analyticity of the function \(f(x)\) and containing the spectrum of \(A\),
\(f(A)\) can be defined as
\begin{equation}\label{eq:eq_ContourInt}
f(A) =\frac{1}{2\pi i} \int_{\Gamma}  f(z) (zI - A)^{-1} dz.
\end{equation}
By applying a quadrature formula on \(N\) points to \eqref{eq:eq_ContourInt},
with weights \(\{c_j\}_{j=1}^{N}\) and nodes \(\{\xi_j\}_{j=1}^{N}\),
it is possible to approximate \eqref{eq:eq_problem_to_solve} as
\begin{equation*}
\mathbf{y} = f(A)\mathbf{x} \approx \sum_{j=1}^{N} c_j (A+\xi_j I)^{-1}\mathbf{x},
\end{equation*}
that is then computationally equivalent to the solution of \(N\) linear
systems with the same right-hand side.

The second approach to problem \eqref{eq:eq_problem_to_solve} resides instead on the
use of projection algorithm. Specifically, we suppose having two $k$-th
dimensional subspaces \(\mathcal{V}\) and \(\mathcal{W}\) spanned by the
column of the matrices \(V,W \in \mathbb{R}^{n \times k}\). Then, problem
\eqref{eq:eq_problem_to_solve} can be projected and approximated on the two subspaces
by doing
\begin{equation*}
\mathbf{y} = f(A)\mathbf{x} \approx W f(V^T A W) V^T \mathbf{x},
\end{equation*}
where now \(A_k = V^T A W\) is a small matrix of size \(k \times k\),
to which we can apply many specific algorithms for the particular choice of
\(f(x)\), [[1]](#MR2396439), or again a quadrature formula.

The PSFUN Library
-----------------

The recent developments on softwares for sparse linear algebra have been made
essential for a wide variety of scientific applications. Specifically, they have
been dedicated to the construction of of massively parallel sparse solvers for
a particular matrix function \(f(x) = x^{-1}\), i.e., for the solution of
large and sparse linear system. A computational framework that lies at the core
of pretty much all multi-physics and multi-scale simulations.

With this library, we try to face the analogous challenge of computing
matrix-function vector products for more general functions than the inverse.

The library described here is substantially based on the parallel BLAS feature
for sparse matrices made available by the [PSBLAS library](https://psctoolkit.github.io/),
and is geared towards the possibility of running on machines with thousands
of high-performance cores, and is divided in three main modules,

Serial module:
   this module implements (or interfaces) the computation of \(f(A)\),
   \(f(A)\mathbf{x}\) for matrices of small-size that can be
   handled in a sequential way,

Krylov module:
   this module implements distributed Krylov based methods for the
   reduction of problem \eqref{eq:eq_problem_to_solve} to the solution of problems of
   small dimensions,

Quadrature module:
   this module implements the approach in \eqref{eq:eq_ContourInt} by implementing
   different quadrature formulas.

How To Install
--------------

The first step to install the PSFUN is to obtain and install the PSBLAS library
from [PSCTOOLKIT](https://psctoolkit.github.io/). All the relevant information
can be found there.

The actual version of the library works with the development version of PSBLAS,
this can be done obtained via GitHub by doing

```bash
git clone https://github.com/sfilippone/psblas3.git
cd psblas3
./configure -with-<stuff>=... -prefix=/path/to/psblas
make -j
make install
```

in which the various `-with-<stuff>=...` options can be read from the
output of the `./configure -h`, again please refer to the original
documentation of PSBLAS for all the relevant information.

Auxiliary packages that can be used with the library are:

   * the package for the computation of \(\varphi\)-functions from , that can be
   obtained from the [ACM website](https://doi.org/10.1145/1499096.1499101).

To **build the documentation** you need to use the [FORD](https://github.com/Fortran-FOSS-Programmers/ford)
Python package, it can be installed simply by doing:
```bash
pip install ford
```
Building the documentation is *optional*, and can
be skipped during the configuration phase. In every case a copy of the docs
is included with the code.

After having installed all the dependencies, and the auxiliary packages the PSFUN
library can be installed via `ccmake` (Version \(\geq\) 3.15), by setting
the position of PSBLAS, and all the auxiliary packages.

```bash
git clone https://github.com/Cirdans-Home/psfun.git
mkdir build
cd build
ccmake ../psfun/
make
make install
```

## References

<a id="">[1]</a> Higham, Nicholas J. Functions of matrices. Theory and computation. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA, 2008. xx+425 pp. ISBN: 978-0-89871-646-7
