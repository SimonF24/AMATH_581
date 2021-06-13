function dy = fitzhugh(t, y, a1, a2, b, c, I, d12, d21)
    % Simulates two linearly coupled fitzhugh neurons
    % t is only present to satisfy the form of matlab's ode solvers
    yCell = num2cell(y);
    [v1, v2, w1, w2] = yCell{:};
    dv1 = -v1^3 + (1 + a1)*v1^2 - a1*v1 - w1 + I + d12*v2;
    dw1 = b*v1 - c*w1;
    dv2 = -v2^3 + (1 + a2)*v2^2 - a2*v2 - w2 + I + d21*v1;
    dw2 = b*v2 - c*w2;
    dy = [dv1; dv2; dw1; dw2];
end