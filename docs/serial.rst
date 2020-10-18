Serial Module
=================================

This module contains the routines needed for the computation of :math:`f(A)x`
for :math:`A` a matrix of small size. It interfaces external codes and algorithms
that usually work with matrix memorized in dense storage. The intended use of
the functions contained here is to use them at the lower level of a Krylov
subspace method. The library directly contains the EXPOKIT code :cite:`expokit`
for the computation of the matrix exponential, together with the scaling and
squaring and Taylor algorithms :cite:`MR508383,MR1981253` by `J. Burkardt <https://people.sc.fsu.edu/~jburkardt/f_src/matrix_exponential/matrix_exponential.html>`_.
For using the phi functions, the code from :cite:`10.1145/1499096.1499101` is needed.
It can be `downloaded <https://dl.acm.org/doi/10.1145/1499096.1499101>`_,
compiled and linked to the main library in the install phase. bla bla

All the implemented functions and the keywords needed to load the are given in
Table :ref:`implemented_methods`.

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
     - "USERF"
     - "SYM"
     -
   * - :math:`\exp(\alpha A)`
     - Taylor
     - General
     - "EXP"
     - "TAYLOR"
     - :cite:`MR508383,MR1981253`
   * -
     - Scaling and Squaring
     - General
     - "EXP"
     - "SASQ"
     - :cite:`MR508383,MR1981253`
   * -
     - Generalized Padè
     - General
     - "EXP"
     - "GENPADE"
     - :cite:`expokit`
   * -
     - Chebyshev
     - Hessenberg
     - "EXP"
     - "CHBHES"
     - :cite:`expokit`
   * -
     - Chebyshev
     - General
     - "EXP"
     - "CHBGEN"
     - :cite:`expokit`
   * -
     - Chebyshev
     - Symmetric
     - "EXP"
     - "CHBSYM"
     - :cite:`expokit`
   * - :math:`\phi_k(\alpha A)`
     - Scaling and Squaring
     - Symmetric
     - "PHI"
     - "NONE"
     - :cite:`10.1145/1499096.1499101`



Module
------------------------------------
.. f:automodule:: psfun_d_serial_mod

Bibliography
------------

.. bibliography:: refserial.bib
   :cited:
   :style: plain
