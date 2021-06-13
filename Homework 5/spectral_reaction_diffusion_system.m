function dydt = spectral_reaction_diffusion_system(t, y, D1, D2, k2, n, beta)

    U_f = y(1:n^2);
    V_f = y(n^2+1:end);
    
    U_t = real(ifft2(reshape(U_f, n, n)));
    V_t = real(ifft2(reshape(V_f, n, n)));
    
    A2 = U_t.^2+V_t.^2;
    lambda = 1-A2;
    w = -beta*A2;
    
    dU_f = reshape(fft2(lambda.*U_t-w.*V_t), n^2, 1)-D1*k2.*U_f; 
    dV_f = reshape(fft2(w.*U_t+lambda.*V_t), n^2, 1)-D2*k2.*V_f;
    
    dydt = [dU_f; dV_f];
    
end