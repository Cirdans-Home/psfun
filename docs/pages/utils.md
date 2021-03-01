---
title: Utils
author: F. Durastante
---

The utils module contains some functions and subroutines which are used in
various places in the library and which do not specifically belong to any of the
other modules. Routines for the computation of some special functions, e.g.,
elliptic integrals, Jacobi polynomials, etc., together with some internal service
routines.

External libraries
==================

To make some plots with the [Gnuplot](http://www.gnuplot.info/) software
directly from the Fortran code, we distribute a modified version of the
[ogpf](https://github.com/kookma/ogpf) library.
