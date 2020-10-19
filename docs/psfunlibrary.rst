****************
Matrix Functions
****************

This library is focused on the computation of matrix-function :cite:`MR2396439` vector products

.. math:: \mathbf{y} = f(A) \mathbf{x}, A \in \mathbb{R}^{n \times n}, \; \operatorname{nnz}(A) = O(n), \; f:\mathbb{R}\rightarrow \mathbb{R},
   :label: eq_problem_to_solve

for large and sparse matrices in a distributed setting. Matrix functions are
ubiquitous in models for applied sciences. They are involved in the solution of
ordinary, partial, and fractional differential equations, systems of coupled
differential equations, hybrid differential-algebraic problems,
equilibrium problems, measures of complex networks, and many others.

To perform the computation in :eq:`eq_problem_to_solve`, we consider here two
main approaches, the first one makes use of a definition based on the *Cauchy integral*
for a matrix function: given a closed contour :math:`\Gamma` lying in the region of
analyticity of the function :math:`f(x)` and containing the spectrum of :math:`A`,
:math:`f(A)` can be defined as

.. math:: f(A) =\frac{1}{2\pi i} \int_{\Gamma}  f(z) (zI - A)^{-1} dz.
   :label: eq_ContourInt

By applying a quadrature formula on :math:`N` points to :eq:`eq_ContourInt`,
with weights :math:`\{c_j\}_{j=1}^{N}` and nodes :math:`\{\xi_j\}_{j=1}^{N}`,
it is possible to approximate :eq:`eq_problem_to_solve` as

.. math:: \mathbf{y} = f(A)\mathbf{x} \approx \sum_{j=1}^{N} c_j (A+\xi_j I)^{-1}\mathbf{x},

that is then computationally equivalent to the solution of :math:`N` linear
systems with the same right-hand side.

The second approach to problem :eq:`eq_problem_to_solve` resides instead on the
use of projection algorithm. Specifically, we suppose having two :math:`k`-th
dimensional subspaces :math:`\mathcal{V}` and :math:`\mathcal{W}` spanned by the
column of the matrices :math:`V,W \in \mathbb{R}^{n \times k}`. Then, problem
:eq:`eq_problem_to_solve` can be projected and approximated on the two subspaces
by doing

.. math:: \mathbf{y} = f(A)\mathbf{x} \approx W f(V^T A W) V^T \mathbf{x},

where now :math:`A_k = V^T A W` is a small matrix of size :math:`k \times k`,
to which we can apply many specific algorithms for the particular choice of
:math:`f(x)`, :cite:`MR2396439`, or again a quadrature formula.

The PSFUN Library
==================

The recent developments on softwares for sparse linear algebra have been made
essential for a wide variety of scientific applications. Specifically, they have
been dedicated to the construction of of massively parallel sparse solvers for
a particular matrix function :math:`f(x) = x^{-1}`, i.e., for the solution of
large and sparse linear system. A computational framework that lies at the core
of pretty much all multi-physics and multi-scale simulations.

With this library, we try to face the analogous challenge of computing
matrix-function vector products for more general functions than the inverse.

The library described here is substantially based on the parallel BLAS feature
for sparse matrices made available by the `PSBLAS library <https://psctoolkit.github.io/>`_,
and is geared towards the possibility of running on machines with thousands
of high-performance cores, and is divided in three main modules,

Serial module:
   this module implements (or interfaces) the computation of :math:`f(A)`,
   :math:`f(A)\mathbf{x}` for matrices of small-size that can be
   handled in a sequential way,

Krylov module:
   this module implements distributed Krylov based methods for the
   reduction of problem :eq:`eq_problem_to_solve` to the solution of problems of
   small dimensions,

Quadrature module:
   this module implements the approach in :eq:`eq_ContourInt` by implementing
   different quadrature formulas.

.. only:: html

      .. figure:: librarystructure.png
         :alt: Structure of the PSFUN library.
         :scale: 80%

         Structure of the PSFUN library.

.. only:: latex

      .. figure:: librarystructure.png
         :alt: Structure of the PSFUN library.
         :scale: 50%

         Structure of the PSFUN library.
