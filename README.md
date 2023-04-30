# Numerical Methods utilites

set of numerical utilities developed in C for the Numerical Methods course in the Mathematics degree at the University of Barcelona. All code is my own and anyone is free to use it in any way. It is NOT memory efficient as I added a layer of pointers on top of all the structures in order to simplify the logic of the functions and for ease of use. The main goal was for it to be very modular and easy to use as the exams consisted in modifying the developed code for new applications. However, to the best of my knowledge the algorithms are correct.

## Contents:

# linalg.c
I implement several matrix operations:
- Initialization (zero and identity)
- Transposition
- Product 
- Solution of systems of equations via Gauss-Jordan elimination
- LU decomposition
- Cholesky decomposition
