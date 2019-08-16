function prm = init(prm, filename)
% INIT Initializes the Environment.
%   prm = An object belonging to the Param Class
%   filename = name of file to be loaded
%
%   Each section does a specific part in loading the parameters
%   from the file into the program.
%   See the init subroutine from Bagus' Code
    
    C = importdata(filename);

    %uses function that converts arguments from text file
    prm.xmin = extractstr(C, 4, 'num'); 
    prm.xmax = extractstr(C, 5, 'num');
    prm.ymin = extractstr(C, 6, 'num');
    prm.ymax = extractstr(C, 7, 'num');
    prm.nbNC = extractstr(C, 11, 'num');
    prm.NumCst = zeros(1,prm.nbNC);
    
    for i = 1:prm.nbNC
        prm.NumCst(i) = extractstr(C, i+11, 'num');
    end

    prm.nbPC = extractstr(C, 21, 'num');
    for i = 1:prm.nbPC
        prm.PhysCst(i) = extractstr(C, i+21, 'num');
    end

    prm.h(1) = (prm.xmax-prm.xmin) / (prm.NumCst(2)-1);
    prm.h(2) = (prm.ymax-prm.ymin) / (prm.NumCst(3)-1);

    %penalization coefficient
    prm.delta = prm.PhysCst(6) * min(prm.h);

    %Boundary Conditions
    prm.bctop = extractstr(C, 30, 'str');
    prm.bcbottom = extractstr(C, 31, 'str');
    prm.bcleft = extractstr(C, 32, 'str');
    prm.bcright = extractstr(C, 33, 'str');

    %creating the x and y coordinates
    prm.x = zeros(1,prm.NumCst(2));
    prm.y = zeros(1,prm.NumCst(3));
    prm.x(1) = prm.xmin;
    prm.y(1) = prm.ymin;

    for i = 2:prm.NumCst(2)
        prm.x(i) = prm.x(i-1) + prm.h(1);
    end

    for i = 2:prm.NumCst(3)
        prm.y(i) = prm.y(i-1) + prm.h(2);
    end

    %Creating Node Coordinates
    prm.neX = prm.NumCst(2) - 1;
    prm.neY = prm.NumCst(3) - 1;
    prm.Tne = prm.neX*prm.neY;
    prm.Tnp = prm.NumCst(2)*prm.NumCst(3);
    prm.xg = zeros(1, prm.NumCst(2)*prm.NumCst(3));
    prm.yg = prm.xg;

    m = 1;
    for i = 1:prm.NumCst(2)*prm.NumCst(3)
        prm.xg(i) = prm.x(m);
        m = m+1;
        if m>prm.NumCst(2)
            m=1;
        end
    end

    m = 1;
    for i = 1:prm.NumCst(2)*prm.NumCst(3)
        prm.yg(i) = prm.y(m);
        if i == m*prm.NumCst(2)
            m=m+1;
        end
    end

    %Creating local element coordinates
        %leXY(k,n) = x/y coordinates 
        %k = element number
        %n = local node number
    prm.leX = zeros(prm.neX*prm.neY, prm.NumCst(1));
    prm.leY = prm.leX;

    kcounter = 1;
    for k = 1:prm.Tne %for all elements
        prm.leX(k,1) = prm.x(kcounter);
        prm.leX(k,2) = prm.x(kcounter) + prm.h(1);
        prm.leX(k,3) = prm.x(kcounter) + prm.h(1);
        prm.leX(k,4) = prm.x(kcounter);
        kcounter = kcounter + 1;
        if kcounter > prm.neX
            kcounter = 1;
        end
    end

    kcounter = 1;
    j = 1;
    i = 1;
    for k = 1:prm.Tne
        prm.leY(k,1) = prm.y(j);
        prm.leY(k,2) = prm.y(j);  
        prm.leY(k,3) = prm.y(j) + prm.h(2);
        prm.leY(k,4) = prm.y(j) + prm.h(2);
        kcounter = kcounter + 1;
        if kcounter > i*prm.neX
            j = j+1;
            i = i+1;
        end
    end    
end
