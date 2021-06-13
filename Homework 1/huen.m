function x = huen(t, dx, x0)
    % t is expected to be a vector with a component for every timestep 
    % that the function should be evaluated at
    % dx is expected to be a callable function with form dx(t,x)
    % dx represents the differential equation being solved
    % x0 is expected to be a vector of initial conditions
    % x0 is included as the first row of x
    % x is a matrix where the rows correspond to timesteps and the columns
    % correspond to different variables as input in x0
    % For the prediction, we assume that the next time step will be the 
    % same size as the previous one
    x_current = x0;
    x = zeros(length(t), length(x0));
    x(1,:) = x0;
    for i=2:length(t)
        current_time = t(i-1);
        tdelta = t(i)-t(i-1);
        est_next_time = current_time + tdelta;
        euler_t_vec = [current_time est_next_time];
        x_pred_full = forward_euler(euler_t_vec, dx, x_current);
        x_pred = x_pred_full(2,:);
        x_next = x_current + tdelta/2 * (dx(current_time, x_current) + dx(est_next_time, x_pred));
        x(i,:) = x_next;
        x_current = x_next;
    end
end