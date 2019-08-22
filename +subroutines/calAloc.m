function bf = calAloc(par, qd, bf)
% CALALOC Creates the local solutions for Matrix A.
%   par = An object belonging to the Param Class
%   qd = An object belonging to the quad class
%   bf = An object belonging to the BasFunc class
%
%   This function also has a fancy wait bar due to how long it takes to 
%   generate the local matrix solutions
%
%   This is a slightly modified code of Bagus' original, in order to
%   increase the efficiency. The original code can be found in calAlocori.


    diam = par.PhysCst(7);
    nobs = par.NumCst(6) + 1;

    bf.Aloc = zeros(par.Tne, par.NumCst(1), par.NumCst(1));
    bff = zeros(1,par.NumCst(1));
    bfdxf = zeros(1,par.NumCst(1));
    bfdyf = zeros(1,par.NumCst(1));
    
    a = zeros(1,nobs);
    eps = (par.xmax - par.xmin) / (nobs-1);
    if nobs == 0
        eps = 0; 
    end
    for m = 1:nobs
        a(m) = par.xmin + m*eps-eps;
    end 
    
    b = a + 0.5D0*diam*eps; %xmax of boundaries
    a = a - 0.5D0*diam*eps; %xmin of boundaries
    c = a; %ymin of boundaries
    d = b; %ymax of boundaries
    
    wxf = par.wx;
    wyf = par.wy;
    qdquad_x0 = qd.quad_x0;
    qdquad_w = qd.quad_w;
    D = par.PhysCst(8);
    
    for k = 1:par.Tne %for all elements
         %Element parameters
         xmin = par.leX(k,1);    
         xmax = par.leX(k,2);
         ymin = par.leY(k,1);
         ymax = par.leY(k,4);
         xe = (xmin + xmax) / 2D0;
         ye = (ymin + ymax) / 2D0;
         dx = xmax - xmin;
         dy = ymax - ymin;
         wx = wxf(xe,ye);
         wy = wyf(xe,ye);
         sigma = 0D0;
         
         %Fancy Wait Bar
         if round(k/par.Tne*100) > round((k-1)/par.Tne*100)
            wbmsg = sprintf('calAloc In progress: %.1f%%', k/par.Tne*100);
            if ~exist('wb', 'var')
                wb = waitbar((k-1)/par.Tne, wbmsg);
            else
                waitbar((k-1)/par.Tne, wb, wbmsg);
            end
         end
         
         %Creating obstacle elements
         for m = 1:nobs
             for n = 1:nobs
                 if xe>a(m) && xe<b(m) && ye>c(n) && ye<d(n)
                     sigma = 1/par.delta^3;
                 end
             end
         end
         
         %Calculating the integrals
         for m = 1:par.NumCst(1)
             for n = 1:par.NumCst(1)
                 bf.Aloc(k,m,n) = 0;
                 for i = 1:par.NumCst(4)
                     xr = (dx/2)* qdquad_x0(i) + dx/2;
                     for j = 1:par.NumCst(4)
                         yr = (dy/2)* qdquad_x0(j) + dy/2;
                         
                         switch m
                             case 1
                                 bff(1) = (1-xr/dx)*(1-yr/dy);
                                 bfdxf(1) = -1/dx + (yr/(dx*dy));
                                 bfdyf(1) = -1/dy + (xr/(dx*dy));
                             case 2
                                 bff(2) = (xr/dx)*(1-yr/dy);
                                 bfdxf(2) = 1/dx - (yr/(dx*dy));
                                 bfdyf(2) = -xr/(dx*dy);
                             case 3
                                 bff(3) = (xr/dx)*(yr/dy);
                                 bfdxf(3) = yr/(dx*dy);
                                 bfdyf(3) = xr/(dx*dy);
                             case 4
                                 bff(4) = (1-xr/dx)*(yr/dy);
                                 bfdxf(4) = -yr/(dx*dy);
                                 bfdyf(4) = 1/dy - (xr/(dx*dy));
                         end
                               
                         switch n
                             case 1
                                 bff(1) = (1-xr/dx)*(1-yr/dy);
                                 bfdxf(1) = -1/dx + (yr/(dx*dy));
                                 bfdyf(1) = -1/dy + (xr/(dx*dy));
                             case 2
                                 bff(2) = (xr/dx)*(1-yr/dy);
                                 bfdxf(2) = 1/dx - (yr/(dx*dy));
                                 bfdyf(2) = -xr/(dx*dy);
                             case 3
                                 bff(3) = (xr/dx)*(yr/dy);
                                 bfdxf(3) = yr/(dx*dy);
                                 bfdyf(3) = xr/(dx*dy);
                             case 4
                                 bff(4) = (1-xr/dx)*(yr/dy);
                                 bfdxf(4) = -yr/(dx*dy);
                                 bfdyf(4) = 1/dy - (xr/(dx*dy));
                         end

                         %Adding up the integral
                         conv = (wx*bfdxf(m) + wy*bfdyf(m)) * bff(n);
                         diff = D*(bfdxf(m)*bfdxf(n) + ...
                                bfdyf(m)*bfdyf(n));
                         F_xy = diff + conv + sigma*bff(m)*bff(n);
                         WF_xy= qdquad_w(i)*qdquad_w(j)*F_xy*((dx*dy)/4);
                         bf.Aloc(k,m,n) = bf.Aloc(k,m,n) + WF_xy;
                     end
                 end
             end
         end
    end
close(wb)
end
