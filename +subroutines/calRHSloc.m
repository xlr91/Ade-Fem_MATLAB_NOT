function bf = calRHSloc(par, qd, bf)

    bf.rhsLoc = zeros(par.Tne, par.NumCst(1));

    %actually idk what this is for 
    %a = ((par.xmax - par.xmin )/3) + par.xmin;
    %b = a + ((par.xmax - par.xmin ) / 3);
    %c = ((par.ymax - par.ymin ) / 3) + par.ymin;
    %d = c + ((par.ymax - par.ymin)/3);
    st = par.sourceterms;

    for k = 1:par.Tne
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
                 for j = 1:par.NumCst(4)
                    %quadrature coordinates
                    xr = (dx/2)* qd.quad_x0(i) + dx/2;
                    yr = (dx/2)* qd.quad_x0(j) + dx/2;

                    %basis functions
                    bf.f(1) = (1-xr/dx)*(1-yr/dy);
                    bf.f(2) = (xr/dx)*(1-yr/dy);
                    bf.f(3) = (xr/dx)*(yr/dy);
                    bf.f(4) = (1-xr/dx)*(yr/dy);

                    sr = 0;

                    for s = 1:par.nSource
                        if xe > st(s,1) && xe < st(s,2) 
                            if ye > st(s,3) && ye < st(s,4)
                                sr = st(s,5);
                            end
                        end
                    end     
                    wsr = sr*qd.quad_w(i)*qd.quad_w(j)*bf.f(m);
                    bf.rhsLoc(k,m) = bf.rhsLoc(k,m) + ((dx*dy)/4)*wsr;
                 end
             end
         end
    end
end