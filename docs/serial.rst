**************
Serial Module
**************

.. highlight:: fortran
.. role:: fortran(code)

This module contains the routines needed for the computation of :math:`f(A)x`
for :math:`A` a matrix of small size. It interfaces external codes and algorithms
that usually work with matrix memorized in dense storage. The intended use of
the functions contained here is to use them at the lower level of a Krylov
subspace method. The library directly contains the EXPOKIT code :cite:`expokit`
for the computation of the matrix exponential, together with the scaling and
squaring and Taylor algorithms :cite:`MR508383,MR1981253` by `J. Burkardt <https://people.sc.fsu.edu/~jburkardt/f_src/matrix_exponential/matrix_exponential.html>`_.
For using the :math:`\varphi`-functions, the code from :cite:`10.1145/1499096.1499101` is needed.
It can be `downloaded <https://dl.acm.org/doi/10.1145/1499096.1499101>`_,
compiled and linked to the main library in the install phase.

The module is centered on the :code:`psfun_d_serial` type, this module contains
all the options needed to set a specific matrix function to be computed. Not all
the options are needed for every type of matrix-function, e.g., the field
:fortran:`integer(psb_ipk_) :: padedegree` is used only if a Padè type algorithm is
employed. All the keywords needed to load the implemented functions and algorithic variants are
given in :numref:`implemented_methods`.

.. list-table:: Implemented Methods
   :name: implemented_methods
   :widths: 10 10 10 10 10 10
   :header-rows: 1

   * - Function
     - Variant
     - Matrix
     - fname
     - variant
     - Source
   * - :math:`f(\alpha A)`
     - Diagonalization
     - Symmetric
     - :code:`"USERF"`
     - :code:`"SYM"`
     -
   * - :math:`\exp(\alpha A)`
     - Taylor
     - General
     - :code:`"EXP"`
     - :code:`"TAYLOR"`
     - :cite:`MR508383,MR1981253`
   * -
     - Scaling and Squaring
     - General
     - :code:`"EXP"`
     - :code:`"SASQ"`
     - :cite:`MR508383,MR1981253`
   * -
     - Generalized Padè
     - General
     - :code:`"EXP"`
     - :code:`"GENPADE"`
     - :cite:`expokit`
   * -
     - Chebyshev
     - Hessenberg
     - :code:`"EXP"`
     - :code:`"CHBHES"`
     - :cite:`expokit`
   * -
     - Chebyshev
     - General
     - :code:`"EXP"`
     - :code:`"CHBGEN"`
     - :cite:`expokit`
   * -
     - Chebyshev
     - Symmetric
     - :code:`"EXP"`
     - :code:`"CHBSYM"`
     - :cite:`expokit`
   * - :math:`\varphi_k(\alpha A)`
     - Scaling and Squaring
     - Symmetric
     - :code:`"PHI"`
     - :code:`"NONE"`
     - :cite:`10.1145/1499096.1499101`



Module
======

.. f:automodule:: psfun_d_serial_mod
