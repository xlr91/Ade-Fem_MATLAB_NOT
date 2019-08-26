D = importdata('sol.tec');
%%
inum = 1025;
jnum = 1025;
n = 1;
testxg = zeros(inum,jnum);
testyg = testxg;
testzg = testxg;
for xi = 1:inum
    for yj = 1:jnum
        thexnumber = str2double(D.textdata{2+n,1});
        theynumber = str2double(D.textdata{2+n,2});
        theznumber = str2double(D.textdata{2+n,3});
        testxg(xi, yj) = thexnumber;
        testyg(xi, yj) = theynumber;
        testzg(xi,yj) = theznumber;
        n = n+1;
    end
end
errorzg = zg-testzg;

%%
surf(testxg, testyg, errorzg, 'FaceColor', 'interp', 'EdgeAlpha', '0.25')
title('Flow of Transport Obeying the ADE, simulated using Finite Elements')
xlabel('x axis')
ylabel('y axis')
zlabel('Concentration u')
colorbar


