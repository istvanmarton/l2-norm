
This repository contains an implementation of a branch-and-bound algorithm calculating the L2 norm of a matrix.

A modified problem is solved at [https://github.com/divipp/kmn-programming](https://github.com/divipp/kmn-programming).

## L2 norm

Let *A* be a set of vectors.
The goal is to partition the set $A$ into two disjoint substest $B$ and $C$ such
that the sum of the Manhattan norms of the sum of vectors in $B$ and the sum of vectors in $C$ is maximal.


## Implementation

The branch-and-bound solution of the problem is implemented in `L2.hs`.

An optimized algorithm is implemented in `L2Opt.hs`.

## Compilation

    $ ghc -O2 L2

## Usage

Example usage:

    $ ./L2 test.mat

should return `282`.
