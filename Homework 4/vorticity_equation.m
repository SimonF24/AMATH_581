function dwdt = vorticity_equation(t, w, A, B, C, L, P, U, nu, method)
    if method == "direct"
        psi = A\w;
    elseif method == "LU"
        y = L\(P*w);
        psi = U\y;
    elseif method == "bicgstab"
        psi = bicgstab(A, w, 1e-6, 600);
    elseif method == "gmres"
        psi = gmres(A, w, 100, 1e-6);
    elseif method == "FFT"
        [N, ~] = size(w);
        n = sqrt(N);
        omega = reshape(w, n, n);
        
        L = 20;
        kx = (2*pi/L)*[0:(n/2-1) (-n/2):-1];
        kx(1) = 1e-6;
        [KX, KY] = meshgrid(kx, kx);
        
        transformed_omega = fft2(omega);
        transformed_psi = -transformed_omega./(KX.^2+KY.^2);
        matrix_psi = ifft2(transformed_psi);
        psi = real(reshape(matrix_psi, N, 1)); 
        % Real is only included as protection against errors, it shouldn't
        % be used
        
    end
    dwdt = -(B*psi).*(C*w)+(C*psi).*(B*w)+nu*(A*w);
end