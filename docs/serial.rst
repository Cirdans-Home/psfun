Serial Module
=================================

This modules contains the routines needed for the computation of :math:`f(A)x`
for :math:`A` a matrix of small size. It interaces external codes and algorithms
that usually work with matrix memorized in dense storage. The intended use of
the functions contained here is to use them at the lower level of a Krylov
subspace method.

.. list-table:: Implemented Methods
   :widths: 25 25 25 50 10 10
   :header-rows: 1

   * - Function
     - Variant
     - Matrix
     - fname
     - variant
     - Description
   * - Exponential
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


Module
------------------------------------
.. f:automodule:: psfun_d_serial_mod

Bibliography
------------

.. bibliography:: refserial.bib
   :cited:
