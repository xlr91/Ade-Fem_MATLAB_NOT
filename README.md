# Ade-Fem_MATLAB_NOT

Solving the Advection-Diffusion Equation based on Bagus Muljadi's Fortran Code. By Emir Muhammad

## Table of Contents
* [Demo](#Demo) 
1.  Introduction
2.  User Guide
3. Technical Descriptions 
    * Definitions 
        (these are like global node numbers, local node numbers, etc that terms will be used throughout this guide)
    * Classes
    * Part 1: Initialization
    * Part 2: Matrix A
    * Part 3: RHS Matrix
## Demo


## User Guide
- talk about the myparam.dat

## Technical Guide
## Definitions
- list down important things with figures, such as global node numbers, local node numbers, neighbor elements, and shit likt that
## Classes

### AIJ
The AIJ (A for Matrix A, Irow and Jcolumn indices) class holds all the sparse matrix information. This is what is fed into the matrix solver. 

GML(k,m,n): [l]
A 3D array where k is the element number, m and n are local node numbers, and l is a counting integer. This is used to transform the local answers into global answers by finding which local node pairs found in one element are also found in another element.  
EXAMPLE.

IRN(l): [n]
A 1D array, the I Row Node, returning the global node number from the counting integer l, originating from GML. Together IRN and JCN make up the row and column indices of the full matrix A

JCN(l): [n]
A 1D array, the J Column Node, returning the global node number from the counting integer l, originating from GML. Together IRN and JCN make up the row and column indices of the full matrix A.

A(l): [value]
A 1D array containing the numerical value of the full matrix element A with indices IRN(l) and JCN(l). Overall, the array A holds the numerical value, IRN and JCN are the row and column indices for the numerical value in a full matrix. 

nonzero= The number of nonzero elements in the matrix A, temporary number used by the code.

nbdof= The total number of nodes in the domain.
RHS= ?


### BasFunc
The BasFunc (Basis Function) class contains information on anything related to using the basis function and the results.

f, dxf, dyf: Temporary arrays used by the code to store values for the basis functions
rhsLoc= ?

Aloc(k,m,n) = [value of integral thing]
A 3d array where k is the element number, m and n are local node numbers. This array contains the values of the integrals for each local element. 

### Param
The Param (parameter) class is the largest file and contains information from the .dat file to create the domain. This class is basically the entire description of the simulated environment. 
Some definitions 
 
Nbe(k,n): [Neighbor Element]
 A 2d array where k is the element number, n is the neighbor number described in the picture below, and the result shows the element number of the neighbor

INSERT PICTURE HERE ABOUT THE NEIGHBOR NUMBERS

Lgm(K,n): [Global Node]
 A 2d array where k is the element number, n is the local node number, and the result is the global node number. 


NumCst(n): [Constant]
A 1d array taken directly from the dat file. Refer to the comments in the dat file for a description of what each Numerical Constant is.

tnode, lnode, rnode, and bnode(n): [Global Node]
A 1d array listing all the global nodes that make up a top, left, right or bottom boundary node. n is just a counting integer. 

qbc(n): [Global Node]
A 1d array listing all the dirichlet boundary nodes. n is a counting integer. Together qbc and qbcval shows the global node boundary condition and the value at that global node. 

nbNC = ?
nbPC = ?
Nbc = ?

neX and neY: Number of elements in the X or Y axis

Tne: Total number of elements in the domain
Tnp: Total number of points in the domain

leX and leY(k,n): [x/y coordinates] 
A 2d array where k is the element number, n is the local node number, and returns the x and y coordinates of the nodes. They are used to convert from local node numbers to x/y coordiantes. 

h(n): [distance]
A 2d array which returns the distance between 2 adjacent nodes in an axis. h(1) is for x, while h(2) is for y.

PhysCst(n): [Constant]
A 1d array taken directly from the dat file. Refer to the comments in the dat file for a description of what each Physical Constant is.

qbcval(n): [Value of Boundary condition]
A 1d array listing all the dirichlet boundary nodes. n is a counting integer. Together qbc and qbcval shows the global node boundary condition and the value at that global node. 

xmin,xmax,ymin,ymax:  Maximum or minimum possible value of the x or y coordinate. This is the extremes of the domain.

uerr = ?
x,y(n) = [coordinates]
A 1d array listing all the coordinates. n is a counting integer.

xg,yg(n): [coordinates]
A 1d array listing the coordinates of the nodes. n is the global node. together, xg and yg are used to convert from global nodes to coordinates

uex = ?
delta: a number that is determined by multiplying the smallest h with the penalization coefficient

bctop, bcbottom, bcleft, bcright: [letter]
A charater representing the type of boundary condition at the top, bottom, left, and right. 'd' for dirichlet and 'n' for neumann.

wx and wy: [function]
These represent the x and y components of the w velocity vector.

### Quad
The quad (quadrature) class contains the weights and the x points, used for Gaussian Quadrature to compute any integrals required. The Theory can be found in XXX. 

quad_x0: [Array], a one-dimensional array where the ith element is xi.

quad_w: [Array], a one-dimensional array where the ith element is wi.

## Part 1: Initialization
### Introduction
The Initialization part of the program reads the param.dat data file, loads it into the program, and processes the information found in the file to create the environment. It loads the 4 classes into the program and alters the Param class extensively. Each subroutine has its own docstring in the .m file and some comment blocks labelling each major section. As most of the code here is straightforward, a better resource would be to look at the class documentation HERE LINK IT

### extractdata
In order to read the data from the dat file properly, a custom function was defined. It converts the MATLAB data from the importdata function and returns the needed value, whether it is a number, string, or function.

## Part 2: Matrix A
### Introduction
The Matrix A part uses the parameter data to create the sparse matrix A. It does this by calculating the local element matrix (see theory guide), and then assembling them into a single matrix. Due to the sparsity of the generated matrix, the matrix is stored as a sparse matrix consisting of sp.A, sp.IRN, sp.JCN. Details can be found in the sparse class description. 
### calAloc
This is the bit that calculates the local matrix elements. a, b, c, and d are the obstacle parameters. The for loop logic goes like 
for all elements...
    check if the element is an 'obstacle' element or not
    for all local node pairs...
        for all quadrature points...
            calculate the basis funtions and sum it all up
It sequentially calculates the basis functions for each node pairs. Why the code does this can be found in the theory guide. 
An efficiency trick here used in this, compared to the original calAloc, was that the switch case was used, which allows the function to only evaluate the required basis functions. Furthermore, the bfdxf and bfdyf values are no longer stored in the BasFunc class but are instead stored within this function which increases efficiency. MATLAB now only needs to look inside this function instead of calling an external class, which takes time. 

### GlobalMap
What global map does is it generates a 'map' of the domain. EXAMPLE HERE EYET YEET AND EXPLAIN SOMETHING SOMETHING SUPPORT. 

### bcond

### assembly

### lagmul

## Part 3 Matrix RHS

