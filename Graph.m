

%Put in the finished soltec.mat file here
load('soltec.mat')

surf(xg, yg, zg, 'FaceColor', 'interp', 'EdgeAlpha', '0.5')
title('Flow of Transport Obeying the ADE, simulated using Finite Elements')
xlabel('x axis')
ylabel('y axis')
zlabel('Concentration u')
colorbar
