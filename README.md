
This repository contains an implementation of a branch-and-bound algorithm calculating the L_n norm of a matrix.

A modified problem is solved at [https://github.com/divipp/kmn-programming](https://github.com/divipp/kmn-programming).

## L_n norm

Let *A* be a set of vectors.
The goal is to partition the set $A$ into $n$ disjoint subsets such
that the sum of the Manhattan norms of the sum of the vectors in the subsets is maximal.
The maximum is the L_n norm.


## Installation

1.  Install the GHC Haskell compiler with [ghcup](https://www.haskell.org/ghcup/).
2.  Execute the following command: 

        $ cabal install

## Usage


### Basic usage

Example usage:

    $ Ln --order 2 test.mat

returns

    Row numbers in the partitions:
      1 5 7 9 10
      2 3 4 6 8
    L2 norm:
      282

Other example:

    $ Ln --order 3 test.mat

returns

    Row numbers in the partitions:
      1 9 10
      2 4 6 8
      3 5 7
    L3 norm:
      350


### Running with guessed result

Giving the guessed result does not change the result, but it is faster:

    $ Ln --order 2 --guessed 282 test.mat

The guessed value may be lower than the norm.
If the guessed value is higher than the norm, the program halts with an error.


### Running on multiple cores

The number of cores can be given by the `+RTS -N` option.
For example, running on 64 cores:

    $ Ln --order 2 test.mat +RTS -N64

You can fine-tune the granularity of the tasks for the cores by the `--depth` option.
Bigger `--depth` means smaller tasks. Minimum depth value is 0, maximum is the number of matrix rows.
Default depth value is 1/4*rows.


