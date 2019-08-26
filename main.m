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
%% Part 4: Solving
tic
A = sparse(sp.IRN, sp.JCN, sp.A);
RHS = sparse(sp.RHS);
RHS = RHS.';
x = A\RHS;
toc
%% Part 5: Output
tic
[xg,yg] = meshgrid(ActualParam.xmin:ActualParam.h(1):ActualParam.xmax,...
ActualParam.ymin:ActualParam.h(2):ActualParam.ymax);
i = 1;
zg = zeros(ActualParam.NumCst(3),ActualParam.NumCst(2));
for yi = 1:ActualParam.NumCst(3)
    for xi = 1:ActualParam.NumCst(2)
        zg(yi,xi) = x(i);
        i = i+1;
    end
end

surf(xg, yg, zg, 'FaceColor', 'interp', 'EdgeAlpha', '0.25')
title('Flow of Transport Obeying the ADE, simulated using Finite Elements')
xlabel('x axis')
ylabel('y axis')
zlabel('Concentration u')
colorbar
toc
%% Part 6: Saving
save('soltec.mat')