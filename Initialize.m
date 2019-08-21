%Some introductionary lines here
%% Part 1: Initialization
tic

clear; close; clc;
ActualParam = class.Param;
qd = class.quad;
bf = class.BasFunc;
sp = class.AIJ;
filename = 'myparam.dat';

ActualParam = subroutines.init(ActualParam, filename); 
ActualParam = subroutines.LGM(ActualParam);
ActualParam = subroutines.NBE(ActualParam);

toc
%% Part 2: Matrix A
tic

qd = subroutines.quad_calc(ActualParam, qd);
bf = subroutines.calAloc(ActualParam,qd,bf);
sp = subroutines.GlobalMap(ActualParam, sp);
ActualParam = subroutines.bcond(ActualParam);
sp = subroutines.assembly(ActualParam, sp, bf);
sp = subroutines.lagmul(ActualParam, sp);

toc
%% Part 3: Matrix RHS
tic

bf = subroutines.calRHSloc(ActualParam, qd, bf);
sp = subroutines.GRHS(sp, ActualParam, bf);
toc
save('soltec.mat')
