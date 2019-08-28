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
    prm.xmin = subroutines.extractdata(C, 4, 'num'); 
    prm.xmax = subroutines.extractdata(C, 5, 'num');
    prm.ymin = subroutines.extractdata(C, 6, 'num');
    prm.ymax = subroutines.extractdata(C, 7, 'num');
    prm.nbNC = subroutines.extractdata(C, 11, 'num');
    prm.NumCst = zeros(1,prm.nbNC);
    
    for i = 1:prm.nbNC
        prm.NumCst(i) = subroutines.extractdata(C, i+11, 'num');
    end

    prm.nbPC = subroutines.extractdata(C, 21, 'num');
    for i = 1:prm.nbPC
        prm.PhysCst(i) = subroutines.extractdata(C, i+21, 'num');
    end

    prm.h(1) = (prm.xmax-prm.xmin) / (prm.NumCst(2)-1);
    prm.h(2) = (prm.ymax-prm.ymin) / (prm.NumCst(3)-1);

    %penalization coefficient
    prm.delta = prm.PhysCst(6) * min(prm.h);

    %Boundary Conditions
    prm.bctop = subroutines.extractdata(C, 30, 'str');
    prm.bcbottom = subroutines.extractdata(C, 31, 'str');
    prm.bcleft = subroutines.extractdata(C, 32, 'str');
    prm.bcright = subroutines.extractdata(C, 33, 'str');

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
    
    %Adding the w vector
    prm.wx = subroutines.extractdata(C,43,'func');
    prm.wy = subroutines.extractdata(C,44,'func');
    
    %Adding Source Terms
    prm.nSource = subroutines.extractdata(C, 48, 'num');
    prm.sourceterms = zeros(prm.nSource, 5);
    for k = 1:prm.nSource
        str = split(C{48+k});
        tempstr = string(str(1));
        newstr = extractAfter(tempstr, '[');
        prm.sourceterms(k,1) = str2double(newstr);
        for i = 2:4
            tempstr = string(str(i));
            prm.sourceterms(k,i) = str2double(tempstr);
        end
        tempstr = string(str(5));
        newstr = extractBefore(tempstr, ']');
        prm.sourceterms(k,5) = str2double(newstr);
    end
end
