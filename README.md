# Ade-Fem_MATLAB_NOT
By Emir Muhammad
Solving the Advection-Diffusion Equation based on Dr Bagus Muljadi's Fortran Code. This code does not require any dependencies, so it can be immediately run after downloading. 

## User Guide
The main way the code works is it interprets the parameters from the myparam.dat file, solves the Advection Diffusion Equation, and saves the solution in the soltec.mat. The code also displays the solution as a 3d graph. 

### Using the code
To use the code, you need to change the parameters found in the myparam.dat file, which can be opened using any text editor. The values are already labelled and explain what they are for. The majority of the properties that can be changed are numerical, however some of them can only be integers or specific numbers which have been labelled. The velocity vector must be written in a format that MATLAB can understand, as it will be directly parsed as a function. The source terms also need to be formatted as labelled, as the code will take in the arrays as is. 

Once the parameters have been edited, the script main.m is can be loaded in and ran. This script automatically generates the graph and saves the solution. However, for domains with more elements (say about 1024x1024 elements), the script can take a long time to run. Therefore, it may be more desirable to use the script files initialize.m, solve.m, and Graph.m. Each of the scripts do a part of main.m, and achieve the exact same result. However, the load command in the beginning of each of these scripts can be edited so that other solutions can be loaded in quickly without having to recalculate every parameter. 

Furthermore, the script also contains graphingerror.m, which is a script that compares the generated MATLAB Graph with the solution from the original ade fem fortran code, if you have access to it. Simply copy the sol.tec file into the same directory as the graphingerror.m script, change the values of inum and jnum to the values number of nodes in x and y, and run it. 


### Documentations
The repository contains several documentations. 

The AnnotatedCode.pdf is my personal notes on Dr Muljadi's own codes for me to understand what is going on so I could translate it. 

The Theory Guide is a brief and very crash course introduction to what Finite Elements are and how they were used in this code. This is by no means a proper  introduction but it should be enough to get started on. I highly recommend Sayas' Gentle Introduction to Finite Elements as a proper introduction to it, as well as Bagus' own works. 

The (UNFINISHED) Technical Guide is meant to be a guide to the MATLAB code. As of now it contains a mostly finished description of the four classes used in the code, however the individual functions have not been properly documented yet. It has *some* useful information, but its still highly disorganized and unpolished so be warned.

The Report is somewhat scientific document that reports what has been done on this code, as well as the math used. 