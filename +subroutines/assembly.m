function sparse = assembly(par, sparse, bf)
% Assembly creates the sparse matrix elements A, IRN, and JCN.
%   prm = An object belonging to the Param Class
    sparse.A = zeros(1, (sparse.nonzero + 2*par.Nbc));
    sparse.IRN = sparse.A;
    sparse.JCN = sparse.A;
    for k = 1:par.Tne
        for n = 1: par.NumCst(1)
            for m = 1:par.NumCst(1)
                l = sparse.GML(k,n,m);
                sparse.A(l) = sparse.A(l) + bf.Aloc(k,n,m);
                sparse.IRN(l) = par.Lgm(k,n);
                sparse.JCN(l) = par.Lgm(k,m);
            end
        end
    end
end
