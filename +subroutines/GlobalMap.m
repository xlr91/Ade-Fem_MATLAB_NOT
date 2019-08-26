function sp = GlobalMap(par, sp)
% GLOBALMAP Generates the Global Map matrix.
%   par = An object belonging to the Param Class
%   sparse = An object belonging to the AIJ class

    sp.nonzero = 0;
    sp.GML = zeros(par.Tne, par.NumCst(1), par.NumCst(1));
    for k = 1:par.Tne
        
         %Fancy Wait Bar
         if round(k/par.Tne*100, -1) > round((k-1)/par.Tne*100, -1)
            wbmsg = sprintf('GlobalMap In progress: %.1f%%', k/par.Tne*100);
            if ~exist('wb', 'var')
                wb = waitbar((k-1)/par.Tne, wbmsg);
            else
                waitbar((k-1)/par.Tne, wb, wbmsg);
            end
         end
       
        for n = 1:par.NumCst(1)
            for m = 1:par.NumCst(1)
                if sp.GML(k,n,m) == 0 %checks if a value already exists
                    sp.nonzero = sp.nonzero + 1;
                    sp.GML(k,n,m) = sp.nonzero;
                    %checks all neighboring elements
                    for k1 = 1:8 
                        k2 = par.Nbe(k,k1);
                        
                        %justincase k2 returns a nonexistent element
                        if k2 > par.Tne
                            k2 = 0;
                        end
                       
                        if k2 ~= 0 %if neighbor exists
                            n1 = 0;
                            m1 = 0;
                            for n2 = 1:par.NumCst(1)
                                if par.Lgm(k,n) == par.Lgm(k2,n2)
                                    n1 = n2;
                                end
                                if par.Lgm(k,m) == par.Lgm(k2,n2)
                                    m1 = n2;
                                end
                            end
                            
                            %if n1 and m1 nodes in the neighbors
                            %and n and m nodes in the original
                            %return the same global nodes
                            if n1~=0 && m1 ~=0
                                sp.GML(k2,n1,m1) = sp.nonzero;
                            end
                        end
                    end
                end
            end
        end
    end
close(wb)
end