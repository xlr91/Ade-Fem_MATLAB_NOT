function bf = calAloc(par, qd, bf)
    %call calAloc
    diam = par.PhysCst(7);
    nobs = par.NumCst(6) + 1;

    bf.Aloc = zeros(par.Tne, par.NumCst(1), par.NumCst(1));
    bf.f = zeros(1,par.NumCst(1));
    bf.dxf = zeros(1,par.NumCst(1));
    bf.dyf = zeros(1,par.NumCst(1));
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

tic
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

         %Creating obstacle elements
         for m = 1:nobs
             for n = 1:nobs
                 if xe>a(m) && xe<b(m) && ye>c(n) && ye<d(n)
                     sigma = 1/par.delta^3;
                 end
             end
         end
         for m = 1:par.NumCst(1)
             for n = 1:par.NumCst(1)
                 bf.Aloc(k,m,n) = 0;
                 for i = 1:par.NumCst(4)
                     for j = 1:par.NumCst(4)
                         %quadrature coordinates
                         xr = (dx/2)* qd.quad_x0(i) + dx/2;
                         yr = (dx/2)* qd.quad_x0(j) + dx/2;

                         %basis functions
                         bf.f(1) = (1-xr/dx)*(1-yr/dy);
                         bf.f(2) = (xr/dx)*(1-yr/dy);
                         bf.f(3) = (xr/dx)*(yr/dy);
                         bf.f(4) = (1-xr/dx)*(yr/dy);

                         %derivatives
                         bf.dxf(1) = -1/dx + (yr/(dx*dy));
                         bf.dyf(1) = -1/dy + (xr/(dx*dy));
                         bf.dxf(2) = 1/dx - (yr/(dx*dy));
                         bf.dyf(2) = -xr/(dx*dy);
                         bf.dxf(3) = yr/(dx*dy);
                         bf.dyf(3) = xr/(dx*dy);
                         bf.dxf(4) = -yr/(dx*dy);
                         bf.dyf(4) = 1/dy - (xr/(dx*dy));

                         %Adding up the integral
                         conv = (wx*bf.dxf(m) + wy*bf.dyf(m)) * bf.f(n);
                         diff = par.PhysCst(8)*(bf.dxf(m)*bf.dxf(n) + ...
                                bf.dyf(m)*bf.dyf(n));
                         F_xy = diff + conv + sigma*bf.f(m)*bf.f(n);
                         WF_xy = qd.quad_w(i)*qd.quad_w(j)*F_xy*((dx*dy)/4);
                         bf.Aloc(k,m,n) = bf.Aloc(k,m,n) + WF_xy;
                     end
                 end
             end
         end

         if round(k/par.Tne*100) > round((k-1)/par.Tne*100)
            wbmsg = sprintf('calAloc In progress: %.1f%%', k/par.Tne*100);
            if ~exist('wb', 'var')
                wb = waitbar((k-1)/par.Tne, wbmsg);
            else
                waitbar((k-1)/par.Tne, wb, wbmsg);
            end
        end
         
         
         
         
    end
toc
close(wb)
end
