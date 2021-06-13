function dydt = chebyshev_reaction_diffusion_system(t, y, D1, D2, lap_Mat, n, beta)

    U = y(1:n^2);
    V = y(n^2+1:end);
    
    A2 = U.^2+V.^2;
    lambda = 1-A2;
    w = -beta*A2;
    
    dU = lambda.*U-w.*V+D1*lap_Mat*U;
    dV = w.*U+lambda.*V+D2*lap_Mat*V;
    
    dydt = [dU; dV];

end