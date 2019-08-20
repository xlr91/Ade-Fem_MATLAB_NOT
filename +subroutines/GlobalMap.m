function sparse = GlobalMap(par, sparse)
% GLOBALMAP Generates the Global Map matrix.
%   par = An object belonging to the Param Class
%   sparse = An object belonging to the AIJ class

    sparse.nonzero = 0;
    sparse.GML = zeros(par.Tne, par.NumCst(1), par.NumCst(1));
    for k = 1:par.Tne
        for n = 1:par.NumCst(1)
            for m = 1:par.NumCst(1)
                if sparse.GML(k,n,m) == 0 %checks if a value already exists
                    sparse.nonzero = sparse.nonzero + 1;
                    sparse.GML(k,n,m) = sparse.nonzero;
                    %checks all neighboring elements
                    for k1 = 1:8 
                        k2 = par.Nbe(k,k1);
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
                                sparse.GML(k2,n1,m1) = sparse.nonzero;
                            end
                        end
                    end
                end
            end
        end
    end
end