%This script solves the matrix after it has been initialized
%% Part 4: Solving
load('soltec.mat')
tic

A = sparse(sp.IRN, sp.JCN, sp.A);
RHS = sparse(sp.RHS);
RHS = RHS.';
x = A\RHS;

toc
%% Part 5: Output
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
save('soltec.mat')
