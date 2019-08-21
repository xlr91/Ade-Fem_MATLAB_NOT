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
%% Part 4: Solving
tic

A = sparse(sp.IRN, sp.JCN, sp.A);
RHS = sparse(sp.RHS);
RHS = RHS.';
x = A\RHS;

toc
%% Part 5: Outputting
% call output
par = ActualParam;
[xg,yg] = meshgrid(par.xmin:par.h(1):par.xmax, par.ymin:par.h(2):par.ymax);
i = 1;
zg = zeros(ActualParam.NumCst(3),ActualParam.NumCst(2));
for yi = 1:ActualParam.NumCst(3)
    for xi = 1:ActualParam.NumCst(2)
        zg(yi,xi) = x(i);
        i = i+1;
    end
end
%xg = par.xg;
%yg = par.yg;
%zg = x(1:par.Tnp).';
surf(xg, yg, zg, 'FaceColor', 'interp', 'EdgeAlpha', '0.5')
title('Flow of Transport Obeying the ADE, simulated using Finite Elements')
xlabel('x axis')
ylabel('y axis')
zlabel('Concentration u')
colorbar
%% Part 6: Graphing/Saving idk yet whatever floats the bill
save('soltec.mat', 'xg', 'yg', 'zg')