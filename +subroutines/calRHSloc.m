function bf = calRHSloc(par, qd, bf)
% CALRHSLOC Generates the local RHS matrix for each element.
%   par = An object belonging to the Param Class
%   qd = An object belonging to the quad class
%   bf = An object belonging to the BasFunc class

    bf.rhsLoc = zeros(par.Tne, par.NumCst(1));

    %a = ((par.xmax - par.xmin )/3) + par.xmin;
    %b = a + ((par.xmax - par.xmin ) / 3);
    %c = ((par.ymax - par.ymin ) / 3) + par.ymin;
    %d = c + ((par.ymax - par.ymin)/3);
    
    bff = zeros(1,par.NumCst(1));
    st = par.sourceterms;
    qdquad_x0 = qd.quad_x0;
    qdquad_w = qd.quad_w;


    for k = 1:par.Tne
         %element parameters
         xmin = par.leX(k,1);    
         xmax = par.leX(k,2);
         ymin = par.leY(k,1);
         ymax = par.leY(k,4);
         xe = (xmin + xmax) / 2D0;
         ye = (ymin + ymax) / 2D0;
         dx = xmax - xmin;
         dy = ymax - ymin;
         
         for m = 1:par.NumCst(1)
             for i = 1:par.NumCst(4)
                    %Quadrature Coordinates
                    xr = (dx/2)* qdquad_x0(i) + dx/2;
                    
                 for j = 1:par.NumCst(4)
                    yr = (dy/2)* qdquad_x0(j) + dy/2;
                    
                    switch m
                        case 1
                            bff(1) = (1-xr/dx)*(1-yr/dy);
                        case 2
                            bff(2) = (xr/dx)*(1-yr/dy);
                        case 3
                            bff(3) = (xr/dx)*(yr/dy);
                        case 4
                            bff(4) = (1-xr/dx)*(yr/dy);
                    end
                    
                    %Source terms
                    sr = 0; 
                    for s = 1:par.nSource
                        if xe > st(s,1) && xe < st(s,2) 
                            if ye > st(s,3) && ye < st(s,4)
                                sr = st(s,5);
                            end
                        end
                    end  
                    
                    wsr = sr*qdquad_w(i)*qdquad_w(j)*bff(m);
                    bf.rhsLoc(k,m) = bf.rhsLoc(k,m) + ((dx*dy)/4)*wsr;
                 end
             end
         end
    end
end