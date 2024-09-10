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
- MapleFiles: This folder contains Maple files used to obtain the re-derivation of the fast Kummer surface, as well as develop the algorithms that are optimised in the Magma files. 

# Reproducing benchmarks

We now indicate how to reproduce the benchmarks in the paper using the file `benchmarks.m`.

### Figure 1
Load `benchmarks.m` in terminal and run 
```
benchmark_fig1(17, 50);
```

### Figure 2
Load `benchmarks.m` in terminal and run 
```
primes_logp100 := [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
benchmark_fig2(100, 50, 50 : Ns := [5, 7, 11, 13, 17, 19] , primes := primes_logp100);
```
Here, `primes_logp100` is a list of primes that we have generated for our benchmarks, and `Ns` is a list of odd prime $N$ that we benchmark. The $i$-th entry of the array correspends to the $i$-th $N$, and the prime is chosen so that we have $\mathbb{F}_{p^2}$-rational $N$-torsion. 

### Table 1 
Load `benchmarks.m` in terminal and run 
```
primes_logp100 := [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
benchmark_table1(100, 50, 50 : Ns := [5, 7, 11, 13, 17, 19] , primes := primes_logp100);
```

### Table 2
Load `benchmarks.m` in terminal and run 
```
primes_logp100 := [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];
benchmark_getimage(100, 50, 50 : Ns := [7, 11, 13, 17, 19], primes := primes_logp100);
```
