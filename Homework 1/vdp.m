function dy = vdp(t, y, eps)
    % Calculates the van der Pol oscillator equation and saves the result
    % as a two element vector where the first element is y' and the
    % second element is y''
    % the input y is expected to be a 2 element vector
    % t is present only for use with ode45, which expects that form
    dy = [y(2); eps*(1-y(1)^2)*y(2)-y(1)];
end