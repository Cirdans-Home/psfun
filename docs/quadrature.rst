*****************
Quadrature Module
*****************

.. highlight:: fortran
.. role:: fortran(code)

This module makes use of the matrix function definition based on the *Cauchy integral*:
given a closed contour :math:`\Gamma` lying in the region of
analyticity of the function :math:`f(x)` and containing the spectrum of :math:`A`,
:math:`f(A)` can be defined as

.. math:: f(A) =\frac{1}{2\pi i} \int_{\Gamma}  f(z) (zI - A)^{-1} dz.
   :label: eq_ContourInt_quadrature

By applying a quadrature formula on :math:`N` points to :eq:`eq_ContourInt_quadrature`,
with weights :math:`\{c_j\}_{j=1}^{N}` and nodes :math:`\{\xi_j\}_{j=1}^{N}`,
it is possible to approximate :eq:`eq_problem_to_solve` as

.. math:: \mathbf{y} = f(A)\mathbf{x} \approx \sum_{j=1}^{N} c_j (A+\xi_j I)^{-1}\mathbf{x},
   :label: eq_QuadratureFormula

that is then computationally equivalent to the solution of :math:`N` linear
systems with the same right-hand side.

Module
======

The construction of the quadrature module is made of several interconnected modules.
The base module is the :fortran:`psfun_base_quadrature_mod`, it contains the
base module of which the different quadratures are extensions.

.. f:automodule:: psfun_base_quadrature_mod

Then the functions for working with the quadrature formula having either
real (:fortran:`psb_dpk_`) or complex quadrature nodes and weights for
:eq:`eq_QuadratureFormula` are contained in the relative modules.

.. f:automodule:: psfun_d_quadrature_mod

.. f:automodule:: psfun_z_quadrature_mod

These two modules make use of :fortran:`abstract interface` for both the
:fortran:`subroutine dquadrule`/:fortran:`subroutine zquadrule` and the generic
function for which we compute the :math:`f(A)`. This is implemented this way
to permit the user to implement its own quadrature rule and functions. An example
of how this can be achieved is contained in the functions included in the
submodule :fortran:`psfun_z_quadrules_mod` that implements the three routines
from :cite:`MR2421045`.

.. f:autosubroutine:: psfun_z_quadrules_mod/hhtmethod1
