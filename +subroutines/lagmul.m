function sp = lagmul(prm,sp)
    sp.nbdof = prm.Tnp;
    sp.nonzero = nnz(sp.A);

    for i = 1:prm.Nbc
        sp.A(sp.nonzero+i) = 1;
        sp.IRN(sp.nonzero+i) = prm.qbc(i);
        sp.JCN(sp.nonzero+i) = i + sp.nbdof;
    end

    sp.nonzero = nnz(sp.A);
    for i = 1:prm.Nbc
        sp.A(sp.nonzero+i) = 1;
        sp.IRN(sp.nonzero+i) = i + sp.nbdof;
        sp.JCN(sp.nonzero+i) = prm.qbc(i);
    end

end