Readme



The project is translating Bagus' Advection-Diffusion Equation using the 
Finite Elements Method from Fortran into MATLAB, so that he can use it 
as a teaching tool for his masters class. This requires reading and 
using his code to understand how it behaves. Reading the code can be done
using Visual Studios Code on the macOS, but because it was meant to run 
in ubuntu the mac cannot run his code.

The best (and so far only) workaround is to use Emir's personal mac to 
access his machine remotely via terminal. All the reading, editing,
and viewing the results of the code is done in the mac, but running
it requires Bagus' machine. Therefore, once all the necessary files
have been edited, the entire code will be securely copied over to 
his machine over the internet, compiled and executed, and shipped back
into Emir's machine to be analyzed. 

How to remotely access Bagus' machine 
1. Open terminal
2. Login to his machine using 
    ssh -X username@DUIP67638
3. Now you are in Bagus' machine on your own account. Transfer the entire
to his machine using the scp command (documentation link)
4. (i think) delete the makefile using make clean command
5. compile the entire code using the make command
6. Once the code finishes compiling, it will create a fem file. This 
fem file is the executable 
7. execite the fem file using ./fem and wait until finished
8. It will give a sol.tec file as the solution. using scp, move the solution 
into Emir's macbook
9. View the sol.tech using paraview 


File Index
- param.dat is the input file containing the input parameters (what you edit)
- gllc_bs.f90 is the main program, calling in other shits
- fem is the executable to be run in Bagus' machine
- sol.tec is the solution the program returns
- makefile is a byproduct (i think) of the make command