function qd = quad_calc(par, qd)
    qd.quad_x0 = zeros(1, par.NumCst(4));
    qd.quad_w = qd.quad_x0;
    
    %selecting case
    switch par.NumCst(4)
        case 1
            qd.quad_w(1) = 1;
        case 2
            qd.quad_w(1) = 1;
            qd.quad_w(2) = 1;

            qd.quad_x0(1) = sqrt(1/3);
            qd.quad_x0(2) = -sqrt(1/3);
        case 3
            qd.quad_w(1) = 0.555555555555555555555555555556;
            qd.quad_w(2) = 0.888888888888888888888888888889;
            qd.quad_w(3) = 0.555555555555555555555555555556;

            qd.quad_x0(1) = 0.774596669241483377035853079956;
            qd.quad_x0(3) = -0.774596669241483377035853079956;
        case 4
            qd.quad_w(1) = 0.347854845137453857373063949222;
            qd.quad_w(2) = 0.652145154862546142626936050778;
            qd.quad_w(3) = 0.652145154862546142626936050778;
            qd.quad_w(4) = 0.347854845137453857373063949222;

            qd.quad_x0(1) = 0.861136311594052575223946488893;
            qd.quad_x0(2) = 0.339981043584856264802665759103;
            qd.quad_x0(3) = -0.339981043584856264802665759103;
            qd.quad_x0(4) = -0.861136311594052575223946488893;
        case 5
            qd.quad_w(1) = 0.236926885056189087514264040720;
            qd.quad_w(2) = 0.478628670499366468041291514836;
            qd.quad_w(3) = 0.568888888888888888888888888889;
            qd.quad_w(4) = 0.478628670499366468041291514836;
            qd.quad_w(5) = 0.236926885056189087514264040720;

            qd.quad_x0(1) = 0.906179845938663992797626878299D0;
            qd.quad_x0(2) = 0.538469310105683091036314420700;
            qd.quad_x0(3) = 0;
            qd.quad_x0(4) = -0.538469310105683091036314420700;
            qd.quad_x0(5) = -0.906179845938663992797626878299D0;
        case 7
            qd.quad_w(1) = 1;
            qd.quad_w(2) = 1;
            qd.quad_w(3) = 1;
            qd.quad_w(4) = 1;
            qd.quad_w(5) = 1;
            qd.quad_w(6) = 1;
            qd.quad_w(7) = 1;

            qd.quad_x0(1) = something;
            qd.quad_x0(2) = something;
            qd.quad_x0(3) = something;
            qd.quad_x0(4) = something;
            qd.quad_x0(5) = something;
            qd.quad_x0(6) = something;
            qd.quad_x0(7) = something;
        otherwise
            fprintf('Quadrature number not recognized\n')
            fprintf('Please double check\n')
    end
end