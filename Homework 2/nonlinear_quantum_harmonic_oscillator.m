function yout = nonlinear_quantum_harmonic_oscillator(x, yin, K, eps, gamma)
    % yout is a vector of the eigenfunction and its first derivative
    % yin is a vector of the current values of the eigenfunction and its
    % first derivative
    % x is the position at which the expression is evaluated
    % K and eps are parameters of the funciton
    y1_current = yin(1); y2_current = yin(2);
    yout1 = y2_current;
    yout2 = (gamma*abs(y1_current)^2+K*x^2-eps)*y1_current;
    yout = [yout1; yout2];
end