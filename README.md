# Ade-Fem_MATLAB_NOT

Solving the Advection-Diffusion Equation based on Bagus Muljadi's Fortran Code. By Emir Muhammad

## Table of Contents
* [Demo](#Demo) 
1.  Introduction
2.  Usage
3. Technical Descriptions 
    * Classes
    * Part 1: Initialization
    * Part 2: Matrix A
    * Part 3: RHS Matrix
    * Part 4: Solving
    * Part 5: Saving
    * Part 6: Displaying
## Demo

## Classes

### AIJ
### BasFunc
### Param
### Quad
The quad class contains the weights and the x points, used for Gaussian Quadrature to compute any integrals required. The Theory can be found in XXX. 

quad_x0: [Array], a one-dimensional array where the ith element is xi.

quad_w: [Array], a one-dimensional array where the ith element is wi.
