function sp = GRHS(sp, par, bf)
% GRHS Generates the Global Right Hand Side matrix.
%   par = An object belonging to the Param Class
%   bf = An object belonging to the BasFunc class
%   sp = An object belonging to the AIJ class
%
    sp.RHS = zeros(1, (par.Tnp+par.Nbc));

    %Add in the values for the nodes in the RHS
    for k = 1:par.Tne
        for n = 1:par.NumCst(1)
            m = par.Lgm(k,n);
            sp.RHS(m) = sp.RHS(m) + bf.rhsLoc(k,n);
        end
    end

    %Add in values for the lagrangan multipliers
    j = par.Tnp;
    j1 = par.Nbc;

    for i = (j+1):(j+j1)
        j2 = i-j;
        sp.RHS(i) = par.qbcval(j2);
    end
end