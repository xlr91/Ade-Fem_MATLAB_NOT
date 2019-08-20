function prm = bcond(prm)
    prm.tnode = zeros(1, prm.NumCst(2));
    prm.bnode = prm.tnode;
    prm.lnode = zeros(1, prm.NumCst(3));
    prm.rnode = prm.lnode;

    j = 0;
    if string(prm.bctop) == 'd'
        j = j+prm.NumCst(2);
    end
    if string(prm.bcbottom) == 'd'
        j = j+prm.NumCst(2);
    end
    if string(prm.bcleft) == 'd'
        j = j+prm.NumCst(3);
    end
    if string(prm.bcright) == 'd'
        j = j+prm.NumCst(3);
    end

    prm.Nbc = j; %number of boundary conditions grossly overestimated
    prm.qbc = zeros(1,prm.Nbc);
    prm.qbcval = prm.qbc;


    tn = prm.Tnp; %total number of nodes
    lt = tn - prm.neX; %readys it for top node

    %Top node
    j = 0;
    for i = 1:prm.NumCst(2)
        prm.tnode(i) = lt + (i-1);
    end

    if string(prm.bctop) == 'd'
        for i = 1:prm.NumCst(2)
            prm.qbc(i) = prm.tnode(i);
            prm.qbcval(i) = prm.PhysCst(1);
        end
        j = j + prm.NumCst(2);
    end

    %Bottom node
    for i = 1:prm.NumCst(2)
        prm.bnode(i) = i;
    end

    if string(prm.bcbottom) == 'd'
        for i = 1:prm.NumCst(2)
            prm.qbc(i+j) = prm.bnode(i);
            prm.qbcval(i+j) = prm.PhysCst(2);
        end
        j = j + prm.NumCst(2);
    end

    %Right node
    for i = 1:prm.NumCst(3)
        k = i*prm.NumCst(2);
        prm.rnode(i) = k;
    end

    if string(prm.bcright) == 'd'
        if string(prm.bctop) == 'd'
            if string(prm.bcbottom) == 'd' 
               for i = 2:prm.NumCst(3)-1
                   prm.qbc(j+(i-1)) = prm.rnode(i);
                   prm.qbcval(j+(i-1)) = prm.PhysCst(4);
               end
               j = j + prm.NumCst(3) - 2;
            else
                for i = 1:prm.NumCst(3) - 1
                    prm.qbc(i+j) = prm.rnode(i);
                    prm.qbcval(i+j) = prm.PhysCst(4);
                end
                j = j+prm.NumCst(3) - 1;
            end
        elseif string(prm.bottom) == 'd'
            for i = 2, prm.NumCst(3)
                prm.qbc(j+(i-1)) = prm.rnode(i);
                prm.qbcval(j+(i-1)) = prm.PhysCst(4);
            end
            j = j+prm.NumCst(3)-1;
        end
    end

    %Left
    for i = 1:prm.NumCst(3)
        prm.lnode(i) = prm.rnode(i) - prm.neX;
    end

    if string(prm.bcleft) == 'd'
        if string(prm.bctop) == 'd'
            if string(prm.bcbottom) == 'd' 
               for i = 2:prm.NumCst(3)-1
                   prm.qbc(j+(i-1)) = prm.lnode(i);
                   prm.qbcval(j+(i-1)) = prm.PhysCst(3);
               end
               j = j + prm.NumCst(3) - 2;
            else
                for i = 1:prm.NumCst(3) - 1
                    prm.qbc(i+j) = prm.lnode(i);
                    prm.qbcval(i+j) = prm.PhysCst(3);
                end
                j = j+prm.NumCst(3) - 1;
            end
        elseif string(prm.bottom) == 'd'
            for i = 2, prm.NumCst(3)
                prm.qbc(j+(i-1)) = prm.lnode(i);
                prm.qbcval(j+(i-1)) = prm.PhysCst(3);
            end
            j = j+prm.NumCst(3)-1;
        end
    end
    prm.Nbc = j;
end