function prm = LGM(prm)
%call LGM
prm.Lgm = zeros(prm.neX*prm.neY,prm.NumCst(1));
m = 1;
for k = 1:prm.Tne
    if k == m*prm.neX
        prm.Lgm(k,1) = floor(k+ (k/prm.neX) - 1);
        m = m+1 ;
    else
        prm.Lgm(k,1) = floor(k + k/prm.neX);
    end
prm.Lgm(k,2) = prm.Lgm(k,1) + 1;
prm.Lgm(k,3) = prm.Lgm(k,2) + prm.NumCst(2);
prm.Lgm(k,4) = prm.Lgm(k,1) + prm.NumCst(2);
end
end
