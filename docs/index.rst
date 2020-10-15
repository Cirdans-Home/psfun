Welcome to PSFUN's documentation!
=================================

Parallel Sparse Matrix Function library. This library contains routines for the
computation of matrix function-vector products

.. math:: y = f(A) x, A \in \mathbb{R}^{n \times n}, \; \operatorname{nnz}(A) = O(n), \; f:\mathbb{R}\rightarrow \mathbb{R},

for large and sparse matrices in a distributed setting.

The parallel environment is managed through the `PSBLAS library <https://psctoolkit.github.io/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   serial
   krylov



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
