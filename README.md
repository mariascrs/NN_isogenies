# Isogenies on Kummer Surfaces

This repository contains the code which is attached to *Isogenies on Kummer Surfaces* by [Maria Corte-Real Santos](https://www.mariascrs.com/) and [E. Victor Flynn](https://people.maths.ox.ac.uk/flynn/). All of the code in this repository was written in [`Magma`](http://magma.maths.usyd.edu.au/magma/).

# Organisation

This directory contains the following files:
- `benchmarks.m` contains functions needed to (re-)run the benchmarks given in the accompanying article.
- `functions.m` contains functions needed to perform arithmetic on Fast Kummer surfaces and to set up (pseudo-)random superspecial Kummer surfaces defined over $\mathbb{F}_{p^2}$ for testing and benchmarks.
- `isoNN.m` is the main file in this repository. It contains the two the main algorithms: 
    - $\textsf{GetIsogeny}$: computes the formulae defining the (N,N)-isogeny given two points generating its kernel. These formulae can then be used, for example, to evaluate points through the isogeny.
    - $\textsf{GetImage}$: computes the image of the (N,N)-isogeny (assuming $\textsf{GetIsogeny}$ has been run without any final scaling).
- `partition.m` computes the partition of monomials of degree-$N$ into four parts (corresponding to which coordinate of the isogeny they appear), as described in the accompanying paper. It is only used for the function `Scaling_sqrt()` in the file `scalings.m`, and can be precomputed for any fixed $N$.
- `runner.m` gives an example of how to run the main algorithms in `isoNN.m`.
- `scaling.m` gives algorithms for all three methods to find the final scaling map, as described in the accompanying paper:
    - Method 1: Finding the scaling for $N = 5$, given by `Scaling_5()`.
    - Method 2: Finding the scaling using Gaussian elimination for $N \geq 7$, given by `Scaling_GE()`.
    - Method 3: Finding the scaling by computing 2 square roots for $N \geq 7$, given by `Scaling_sqrt()`.
