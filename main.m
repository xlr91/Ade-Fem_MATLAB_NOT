%Some introductionary lines here
%% Part 1: Initialization
clear; close; clc;
ActualParam = class.Param;
qd = class.quad;
bf = class.BasFunc;
sp = class.AIJ;
filename = 'myparam.dat';

ActualParam = subroutines.init(ActualParam, filename); 
ActualParam = subroutines.LGM(ActualParam);
ActualParam = subroutines.NBE(ActualParam);
%% Part 2: Matrix A

qd = subroutines.quad_calc(ActualParam, qd);
bf = subroutines.calAloc(ActualParam,qd,bf);
sp = subroutines.GlobalMap(ActualParam, sp);
ActualParam = subroutines.bcond(ActualParam);
sp = subroutines.assembly(ActualParam, sp, bf);
sp = subroutines.lagmul(ActualParam, sp);

matrix = sparse(sp.IRN, sp.JCN, sp.A);
matrix = full(matrix);
%% Part 3: Matrix RHS
bf = subroutines.calRHSloc(ActualParam, qd, bf);
sp = subroutines.GRHS(sp, ActualParam, bf);
%% Part 4: Solving
A = sparse(sp.IRN, sp.JCN, sp.A);
RHS = sparse(sp.RHS);
RHS = RHS.';
x = A\RHS;
%% Part 5: Graphing/Saving idk yet whatever floats the bill
