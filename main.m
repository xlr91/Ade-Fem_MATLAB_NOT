%Some introductionary lines here

%% Part 1: Initialization
clear; close; clc;
ActualParam = class.Param;
qd = class.quad;
bf = class.BasFunc;
sparse = class.AIJ;
filename = 'myparam.dat';

ActualParam = subroutines.init(ActualParam, filename); 
ActualParam = subroutines.LGM(ActualParam);
ActualParam = subroutines.NBE(ActualParam);
%% Part 2: Matrix A

%% Part 3: Matrix RHS
%% Part 4: Solving
%% Part 5: Graphing/Saving idk yet whatever floats the bill