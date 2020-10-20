*************
Krylov Module
*************

.. highlight:: fortran
.. role:: fortran(code)

Let :math:`V_{k}` be an orthogonal matrix whose columns :math:`\mathbf{x}_{1},\dots, \mathbf{x}_{k}`
span an arbitrary Krylov subspace :math:`\mathcal{W}_{k}(A,\mathbf{x})` of
dimension :math:`k`. We obtain an approximation of :math:`f(A)\mathbf{x}` by

.. math:: f(A)\mathbf{x} = V_{k}f(V_{k}^{T}AV_{k})V_{k}^{T}\mathbf{x}.
  :label: krylov_fA_approx

Different methods for the approximation of matrix functions are obtained for
different choices of the projection spaces :math:`\mathcal{W}_{k}(A,\mathbf{x})`.

Given a set of scalars :math:`\left\{ \sigma _{1},\dots,\sigma _{k-1}\right\} \subset \overline{\mathbb{C}}`
in the the extended complex plane :math:`\overline{\mathbb{C}}`, that are not
eigenvalues of :math:`A`, let

.. math:: q_{k-1}(z)=\prod\nolimits_{j=1}^{k-1}(\sigma _{j}-z).

The **rational Krylov** subspace of order :math:`k` associated with :math:`A`,
:math:`\mathbf{x}` and :math:`q_{k-1}` is defined by

.. math:: \mathcal{Q}_{k}(A,\mathbf{x})=\left[ q_{k-1}(A)\right] ^{-1}\mathcal{K}_{k}(A,\mathbf{x}),

where

.. math:: \mathcal{K}_{k}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x},A\mathbf{x},\ldots ,A^{k-1}\mathbf{x}\}

is the standard polynomial Krylov space.

By defining the matrices

.. math:: C_{j}=\left( \mu _{j}\sigma _{j}A-I\right) (\sigma _{j}I-A)^{-1},

where :math:`\left\{ \mu _{1},\dots, \mu _{k-1}\right\} \subset \overline{\mathbb{C}}`
are such that :math:`\sigma _{j}\neq $ $\mu _{j}^{-2}`, it is known that the
rational Krylov space can also be written as follows :cite:`guttel2013rational`

.. math:: \mathcal{Q}_{k}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x},C_{1}\mathbf{x},\ldots ,C_{k-1}\cdots C_2 C_{1}\mathbf{x}\}.

This general formulation allows to recast most of the classical Krylov
methods in terms of a rational Krylov method with a specific choice of :math:`\sigma _{j}`
and :math:`\mu _{j}`. In particular,

  * the **polynomial Krylov** method in which :math:`\mathcal{W}_{k}(A,\mathbf{x})=\mathcal{K}_{k}(A,\mathbf{x})`
    can be recovered by defining :math:`\mu _{j}=1` and :math:`\sigma_{j}=\infty` for each `j`.
  * The **extended Krylov** method :cite:`EKSM1,EKSM2`, in which

    .. math:: \mathcal{W}_{2k-1}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x}, A^{-1}\mathbf{x}, A\mathbf{x}, \ldots, A^{-(k-1)}\mathbf{x}, A^{k-1}\mathbf{x}\},

    is obtained by setting

    .. math:: (\mu _j,\sigma_j)= \left\lbrace \begin{array}{ll} (1,\infty), & \mbox{for  $j$ even},\\ (0,0),  & \mbox{for  $j$  odd.} \end{array}\right.

  * The **shift-and-invert** rational Krylov :cite:`moret2004rd,van2006preconditioning`, where

    .. math:: \mathcal{W}_{k}(A,\mathbf{x})=\operatorname{Span}\{\mathbf{x},(\sigma I-A)^{-1}\mathbf{x},\ldots ,(\sigma I-A)^{-(k-1)}\mathbf{x}\},

    is defined by taking :math:`\mu _{j}=0` and :math:`\sigma _{j}=\sigma` for each :math:`j`.

The PSFUN library contains the implementation of several flavour of these methods
that can be used for the computation of :eq:`krylov_fA_approx`, the field
in the :fortran:`psfun_d_krylov` type represent the options neeeded to for setting
up and applying the different implemented method for a given matrix function
:fortran:`fun` (represented by an object of type :fortran:`psfun_d_serial`).

:numref:`implemented_krylov_methods` has the info on the method available.

.. list-table:: Implemented Krylov Methods
   :name: implemented_krylov_methods
   :widths: 10 10 10 10 10
   :header-rows: 1

   *  - Method
      - Class
      - Matrix Type
      - :fortran:`kname`
      - Source
   *  - Arnoldi
      - Polynomial
      - General
      - :fortran:`"ARNOLDI"`
      - :cite:`MR1149094`
   *  - Lanczos
      - Polynomial
      - Symmetric
      - :fortran:`"LANCZOS"`
      - :cite:`MR1149094`

Stopping Criterion
==================


Module
======
.. f:automodule:: psfun_d_krylov_mod
