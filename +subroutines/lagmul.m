function sp = lagmul(prm,sp)
% LAGMUL Generates the Lagrangian Multiplier Matrix.
%   prm = An object belonging to the Param Class
%   filename = name of file to be loaded
%
%   Each section does a specific part in loading the parameters
%   from the file into the program.
%   See the init subroutine from Bagus' Code
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