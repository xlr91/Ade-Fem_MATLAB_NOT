classdef Param
    properties
    %Integer, 2d 
    Nbe;Lgm
    
    %Integer 1D
    NumCst;tnode;lnode;rnode;bnode;qbc
    
    %Integers
    nbNC; nbPC; Nbc;
    neX;neY;Tne;Tnp;
    
    %Double Precision 2D
    leX;leY
    
    %Double Precision 1D-2
    h
    
    %Double Precision 1D
    PhysCst; qbcval
    
    %Double Precision
    xmin;xmax;ymin;ymax;uerr
    x;y;xg;yg;uex
    delta
    
    %Characters
    bctop
    bcbottom
    bcleft
    bcright
    
    end
end
