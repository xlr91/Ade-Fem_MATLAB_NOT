function prm = NBE(prm)
prm.Nbe = zeros(prm.Tne,8);
left = prm.xmin;
right = prm.xmax;
bottom = prm.ymin;
top = prm.ymax;
for k = 1:prm.Tne
    if prm.leX(k,1) > left && prm.leY(k,1) > bottom
        prm.Nbe(k,1) = k- (prm.neX +1);
    end

    if prm.leY(k,1) > bottom
        prm.Nbe(k,2) = k-prm.neX;
    end

    if prm.leX(k,2) < right && prm.leY(k,2) > bottom
        prm.Nbe(k,3) = k-(prm.neX -1);
    end
    
    if prm.leX(k,2) < right
        prm.Nbe(k,4) = k + 1;
    end
    
    if prm.leX(k,3) < right && prm.leY(k,3) < top
        prm.Nbe(k,5) = k+(prm.neX +1);
    end
    
    if  prm.leY(k,3) < top
        prm.Nbe(k,6) = k + (prm.neX);
    end
    
    if prm.leX(k,4) > left && prm.leY(k,4) < top
        prm.Nbe(k,7) = k + (prm.neX -1);
    end
    
    if prm.leX(k,1) > left
        prm.Nbe(k,8) = k-1 ;
    end
end

end
