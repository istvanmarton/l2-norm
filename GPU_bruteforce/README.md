This directory contains a brute-force GPU implementation of the Lk matrix norms.

## Basic usage

The compilation (needs the nvcc Nvidia CUDA compiler) and can be compiled with the following command:

     nvcc -o L L.cu

The compiled code can be invoked with the

     ./L num filename_of_matrix

command, where 'num' is a positive integer indicating the order of the L norm we want the code to calculate. The 'filename_of_matrix' is the name of the file containing the matrix. The maximal order of the L norm the code can calculate is defined in the 'RANK_OF_NORM' variable. If one wants to calculate an L norm of order bigger than what is defined by 'RANK_OF_NORM', then the 'RANK_OF_NORM' variable should be increased, and the code should be recompiled. The entries of the matrix must be integers. The data separator containing the matrix can be any character apart from a numeral, and any of the following characters: 'e', 'E', '.', '+', '-'. The 'length' variable defines the maximal number of columns of a matrix, the program can deal with in case the order of the L norm is bigger than one. If the order of the L norm is 1, the number of rows or columns (whichever is bigger) of the matrix can not be more than the 'length' variable. In case it is bigger, the 'length' variable should be increased, and the code should be recompiled.

The program has 2 outputs. It writes out the L norm of the matrix to the screen. It also writes the strategy vector belonging to the L norm to a file. The filename has the following format: strategy_Lnum.txt, where num means the order of the L norm. For example, the file 'strategy_L2.txt' contains the strategy vector belonging to the $L_2$ norm of the matrix. These files consist of $\pm 1$ in case the order is not bigger than 2. In case it is bigger than 2, the entries of the strategy files are {0,1,...,num-1}, where num is the order of the L norm.
