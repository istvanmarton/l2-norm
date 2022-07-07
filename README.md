
This repository contains an implementation of a branch-and-bound algorithm calculating the L2 norm of a matrix.

A modified problem is solved at [https://github.com/divipp/kmn-programming](https://github.com/divipp/kmn-programming).

## L2 norm

Let *A* be a set of vectors.
The goal is to partition the set $A$ into two disjoint substest $B$ and $C$ such
that the sum of the Manhattan norms of the sum of vectors in $B$ and the sum of vectors in $C$ is maximal.


## Implementation

The essence of the branch-and-bound algorithm is implemented in `L2Slow.hs`.

`L2.hs` is the optimized version of `L2Slow.hs`.


## Installation

1.  Install the GHC Haskell compiler (for example with [ghcup](https://www.haskell.org/ghcup/)).
2.  Execute the following command: 

        $ cabal install

## Usage


### Basic usage

Example usage:

    $ L2 test.mat

returns

    Row numbers in the two partitions:
      Partition A:  1. 5. 7. 9. 10.
      Partition B:  2. 3. 4. 6. 8.
    L2 norm:
      282


### Running with guessed result

Giving the guessed result does not change the result, but it is faster:

    $ L2 --guessed 282 test.mat

Note that `L2` will never give a result which is lower than the guessed value.


### Running on multiple cores

The number of cores can be given by the `+RTS -N` option.
Example:

    $ L2 test.mat +RTS -N2
