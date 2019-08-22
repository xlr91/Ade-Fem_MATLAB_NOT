% This script takes the solved matrix and graphs it
% It is useful to use this if the environment is large

%Put in the finished soltec.mat file here
load('soltec.mat')

surf(xg, yg, zg, 'FaceColor', 'interp', 'EdgeAlpha', '0.5')
title('Flow of Transport Obeying the ADE, simulated using Finite Elements')
xlabel('x axis')
ylabel('y axis')
zlabel('Concentration u')
colorbar
