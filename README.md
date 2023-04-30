# Numerical Methods utilites

Set of numerical utilities developed in C and Python for both Numerical Methods courses in the Mathematics degree at the University of Barcelona. All code is my own and anyone is free to use it in any way. It is NOT memory efficient as I added a layer of pointers on top of all the structures in order to simplify the logic of the functions and for ease of use. The main goal was for it to be very modular and easy to use as the exams consisted in modifying the developed code for new applications. However, to the best of my knowledge the algorithms are correct.

## Contents:

# linalg.c
I implement several matrix operations:
- Initialization (zero and identity)
- Transposition
- Product 
- Solution of systems of equations via Gauss-Jordan elimination
- LU decomposition
- Cholesky decomposition

# functions.c
I implement several utilities related to continuous functions including:
- Representation of functions in finite subsets of R as x, y coordinates
- Saving these to files for representation using GNUPlot
- Newton and Lagrange polynomial interpolation
- Txebitxev coordinate calculations
- Root-finding via the Newton-Rhapson and Secant methods
- Simple approximations for derivatives and integrals
- Refining of these approximations in an efficient way via the Richardson and Romberg tricks
