![splpak](media/splpak.png)
============

This library contains routines for fitting (least squares) a multidimensional cubic spline to arbitrarily located data.  It also contains routines for evaluating this spline (or its partial derivatives) at any point.
This is a modernization of the double precision SPLPAK files from [NCL](https://github.com/NCAR/ncl).

## Status

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/splpak.svg)](https://github.com/jacobwilliams/splpak/releases/latest)
[![Build Status](https://github.com/jacobwilliams/splpak/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/splpak/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/splpak/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/splpak)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/splpak)](https://github.com/jacobwilliams/splpak/commits/master)

## Compiling

A `fmp.toml` file is provided for compiling splpak with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To run the unit tests:

```
fpm test --profile release
```

To use `splpak` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
splpak = { git="https://github.com/jacobwilliams/splpak.git" }
```

or, to use a specific version:
```toml
[dependencies]
splpak = { git="https://github.com/jacobwilliams/splpak.git", tag = "2.0.0"  }
```

## Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/splpak/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford) (i.e. by running `ford ford.md`).

## See also
 *  [The NCAR Command Language ](https://github.com/NCAR/ncl) (specificially, the [csagrid](https://github.com/NCAR/ncl/tree/develop/ngmath/src/lib/gridpack/csagrid) directory)
 * [bspline](https://github.com/NCAR/bspline) - Cubic B-Spline implementation in C++ templates. Also has a copy of [splpak.f](https://github.com/NCAR/bspline/tree/master/Tests/Fortran)
 * [Ngmath Library Overview](https://ngwww.ucar.edu/ngmath/)
 * [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran) Multidimensional B-Spline Interpolation of Data on a Regular Grid
 * [regridpack](https://github.com/jacobwilliams/regridpack) Linear or cubic interpolation for 1D-4D grids
 * [finterp](https://github.com/jacobwilliams/finterp) 1D-6D linear or nearest-neighbor interpolation