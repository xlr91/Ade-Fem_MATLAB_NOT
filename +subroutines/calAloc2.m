function bf = calAloc2(par, qd, bf)
% CALALOC Creates the local solutions for Matrix A.
%   par = An object belonging to the Param Class
%   qd = An object belonging to the quad class
%   bf = An object belonging to the BasFunc class
%
%   This function also has a fancy wait bar due to how long it takes to 
%   generate the local matrix solutions
%
%   This one switches out the repeating basis functions by defining some
%   equations in the for loop, so it only calculates the necessary things
%   instead of wasting computing resource on calculating BFs that are
%   ultimately not needed. 
    %call calAloc
    diam = par.PhysCst(7);
    nobs = par.NumCst(6) + 1;

    bf.Aloc = zeros(par.Tne, par.NumCst(1), par.NumCst(1));
    %bf.f = zeros(1,par.NumCst(1));
    %bf.dxf = zeros(1,par.NumCst(1));
    %bf.dyf = zeros(1,par.NumCst(1));
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

     %basis functions
     bff{1} = @(xr,yr, dx, dy) (1-xr/dx)*(1-yr/dy);
     bff{2} = @(xr,yr, dx, dy) (xr/dx)*(1-yr/dy);
     bff{3} = @(xr,yr, dx, dy) (xr/dx)*(yr/dy);
     bff{4} = @(xr,yr, dx, dy) (1-xr/dx)*(yr/dy);

     %derivatives
     bfdxf{1} = @(xr,yr, dx, dy) -1/dx + (yr/(dx*dy));
     bfdyf{1} = @(xr,yr, dx, dy) -1/dy + (xr/(dx*dy));
     bfdxf{2} = @(xr,yr, dx, dy) 1/dx - (yr/(dx*dy));
     bfdyf{2} = @(xr,yr, dx, dy) -xr/(dx*dy);
     bfdxf{3} = @(xr,yr, dx, dy) yr/(dx*dy);
     bfdyf{3} = @(xr,yr, dx, dy) xr/(dx*dy);
     bfdxf{4} = @(xr,yr, dx, dy) -yr/(dx*dy);
     bfdyf{4} = @(xr,yr, dx, dy) 1/dy - (xr/(dx*dy));
     
     %Quadrture stuff
     qdquad_x0 = qd.quad_x0;
     qdquad_w = qd.quad_w;
                         
    for k = 1:par.Tne
         xmin = par.leX(k,1);    
         xmax = par.leX(k,2);
         ymin = par.leY(k,1);
         ymax = par.leY(k,4);
         xe = (xmin + xmax) / 2D0;
         ye = (ymin + ymax) / 2D0;
         dx = xmax - xmin;
         dy = ymax - ymin;
         wx = par.wx(xe,ye);
         wy = par.wy(xe,ye);
         sigma = 0D0;
         
         
         %Fancy Wait bar
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
         %Calculating the Integral
         for m = 1:par.NumCst(1)
             for n = 1:par.NumCst(1)
                 bf.Aloc(k,m,n) = 0;
                 for i = 1:par.NumCst(4)
                     for j = 1:par.NumCst(4)
                         %quadrature coordinates
                         xr = (dx/2)* qdquad_x0(i) + dx/2;
                         yr = (dx/2)* qdquad_x0(j) + dx/2;
                         
                         bfdxfm = bfdxf{m}(xr,yr, dx, dy);
                         bfdyfm = bfdyf{m}(xr,yr, dx, dy);
                        
                         bfdxfn = bfdxf{n}(xr,yr, dx, dy);
                         bfdyfn = bfdyf{n}(xr,yr, dx, dy);
                         
                         bffn = bff{n}(xr,yr, dx, dy);
                         bffm = bff{m}(xr,yr, dx, dy);        
                         
                         %Adding up the integral
                         conv = (wx*bfdxfm + wy*bfdyfm) * bffn;
                         diff = par.PhysCst(8)*(bfdxfm*bfdxfn + ...
                                bfdyfm*bfdyfn);
                         F_xy = diff + conv + sigma*bffm*bffn;
                         WF_xy = qdquad_w(i)*qdquad_w(j)*F_xy*((dx*dy)/4);
                         bf.Aloc(k,m,n) = bf.Aloc(k,m,n) + WF_xy;
                     end
                 end
             end
         end        
    end
close(wb)
end
