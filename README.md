Modernization of the double precision SPLPAK files from [NCL](https://github.com/NCAR/ncl). A work in progress.

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

## Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/splpak/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford) (i.e. by running `ford ford.md`).

## See also
 * [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran) Multidimensional B-Spline Interpolation of Data on a Regular Grid
 * [bspline](https://github.com/NCAR/bspline) - Cubic B-Spline implementation in C++ templates. Also has a copy of [splpak.f](https://github.com/NCAR/bspline/tree/master/Tests/Fortran)
