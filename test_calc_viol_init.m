G2dU = rand(50, 200);
w = rand(8e4,1);
M_r = rand(5000,100);
viol = calc_violations(G2dU, w, 0.1, M_r);