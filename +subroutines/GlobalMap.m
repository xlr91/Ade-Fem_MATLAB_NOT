function sparse = GlobalMap(par, sparse)
    sparse.nonzero = 0;
    sparse.GML = zeros(par.Tne, par.NumCst(1), par.NumCst(1));
    for k = 1:par.Tne
        for n = 1:par.NumCst(1)
            for m = 1:par.NumCst(1)
                if sparse.GML(k,n,m) == 0
                    sparse.nonzero = sparse.nonzero + 1;
                    sparse.GML(k,n,m) = sparse.nonzero;
                    for k1 = 1:8
                        k2 = par.Nbe(k,k1);
                        if k2 ~= 0 
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