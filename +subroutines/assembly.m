function sp = assembly(par, sp, bf)
% Assembly creates the sparse matrix elements A, IRN, and JCN.
%   prm = An object belonging to the Param Class
    sp.A = zeros(1, (sp.nonzero + 2*par.Nbc));
    sp.IRN = sp.A;
    sp.JCN = sp.A;
    for k = 1:par.Tne
        for n = 1: par.NumCst(1)
            for m = 1:par.NumCst(1)
                l = sp.GML(k,n,m);
                sp.A(l) = sp.A(l) + bf.Aloc(k,n,m);
                sp.IRN(l) = par.Lgm(k,n);
                sp.JCN(l) = par.Lgm(k,m);
            end
        end
    end
end
