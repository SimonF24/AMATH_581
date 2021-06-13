function x = forward_euler(t, dx, x0)
    % t is expected to be a vector with a component for every timestep 
    % that the function should be evaluated at
    % dx is expected to be a callable function with form dx(t,x)
    % dx represents the differential equation being solved
    % x0 is expected to be a vector of initial conditions
    % x0 is included as the first row of x
    % x is a matrix where the rows correspond to timesteps and the columns
    % correspond to different variables as input in x0
    x_current = x0;
    x = zeros(length(t), length(x0));
    x(1,:) = x0;
    for i=2:length(t)
        current_time = t(i-1);
        tdelta = t(i)-t(i-1);
        x_next = x_current + tdelta * dx(current_time, x_current);
        x(i,:) = x_next;
        x_current = x_next;
    end
end