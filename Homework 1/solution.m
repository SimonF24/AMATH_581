% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.
% your extra functions, if you need them, can be in other files (don't forget to upload them too!)

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15] = solution()
 [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15] = student_solution(dummy_argument)
    % your solution code goes here
    % assign the variables you are asked to save here 
    
    %% Problem 1
    
    [A1, A2, A3] = eval_method("euler", false);
    [A4, A5, A6] = eval_method("huen", false);
    
    %% Problem 2
    
    tspan = [0:0.5:32];
    
    eps_array = [0.1 1 20];
    sol = zeros(length(tspan), 3);
    y0 = [sqrt(3) 1];
    for i=1:3
        eps = eps_array(i);
        [~, y] = ode45(@(t,y) vdp(t, y, eps), tspan, y0);
        sol(:,i) = y(:, 1);
    end
    A7 = sol;
    
    A8 = eval_solver("ode45", false);
    A9 = eval_solver("ode23", false);
    A10 = eval_solver("ode113", false);
    
    %% Problem 3
    
    A11 = eval_interaction_params(0,0, false);
    A12 = eval_interaction_params(0, 0.2, false);
    A13 = eval_interaction_params(-0.1, 0.2, false);
    A14 = eval_interaction_params(-0.3, 0.2, false);
    A15 = eval_interaction_params(-0.5, 0.2, false);
    
end

function [ans1, ans2, ans3] = eval_method(method, show_plots)
    % Evaluates a solution method in the manner expected by problem 1
        
    df1 = @(t, y) (-3*y*sin(t));
    y0 = pi/sqrt(2);
    
    sol = @(t) (pi*exp(3*(cos(t)-1))/sqrt(2));
    
    deltat_vals = zeros(7, 1);
    errors = zeros(1, 7);
    for i=2:8
        deltat_vals(i-1) = 2^-i;
        t = [0:2^-i:5];
        if method == "euler"
            x = forward_euler(t, df1, y0);
        elseif method == "huen"
            x = huen(t, df1, y0);
        end
        if i == 8
            ans1 = x;
        end
        y_true = sol(t);
        y_num = transpose(x);
        errors(i-1) = mean(abs(y_true-y_num));
        if show_plots
            % Checking solutions
            figure('Name', 'True Solution Comparison');
            plot(t(2:end), y_true)
            hold on
            plot(t, x)
            txt = sprintf('Solution Comparison for i=%d', i);
            title(txt);
            xlabel('t');
            ylabel('y');
            legend('True Values', 'Solved Values')
        end
    end
    
    ans2 = errors;
    
    p = polyfit(log(deltat_vals), log(errors), 1);
    
    if show_plots
        % Plots mentioned in the assignment
        figure('Name', 'Delta-t vs error');
        plot(deltat_vals, errors)
        xlabel('{\Delta}t'); ylabel('Error');
        figure('Name', 'Log(Delta-t) vs log(error)');
        plot(log(deltat_vals), log(errors));
        xlabel('log({\Delta}t)'); ylabel('log(error)');
        hold on
        plot(log(deltat_vals), polyval(p, log(deltat_vals)), 'r')
        legend('Original', 'Line Fit')
    end
    
    ans3 = p(1);
end

function answer = eval_solver(solver, show_plot)
    % Evaluates a solver in the manner required in problem 2
    av_steps = zeros(length(4:10), 1);
    counter = 1;
    eps = 1;
    tols = zeros(length(4:10), 1);
    tspan = [0 32];
    y0 = [2 pi^2];
    for i=4:10
        tol = 10^-i;
        tols(counter) = tol;
        options = odeset('AbsTol', tol, 'RelTol', tol);
        if solver == "ode45"
            [t, ~] = ode45(@(t, y) vdp(t, y, eps), tspan, y0, options);
        elseif solver == "ode23"
            [t, ~] = ode23(@(t, y) vdp(t, y, eps), tspan, y0, options);
        elseif solver == "ode113"
            [t, ~] = ode113(@(t, y) vdp(t, y, eps), tspan, y0, options);
        else
            error('Unsupported Solver')
        end
        av_step = mean(diff(t));
        av_steps(counter) = av_step;
        counter = counter + 1;
    end
    
    p = polyfit(log(av_steps), log(tols), 1);
    
    if show_plot
        txt = sprintf('Solver Evaluation for %s', solver);
        figure('Name', txt)
        plot(log(av_steps), log(tols))
        hold on
        plot(log(av_steps), polyval(p, log(av_steps)), 'r')
        legend('Original', 'Line Fit')
    end
    
    answer = p(1);
end

function answer = eval_interaction_params(d12, d21, show_plot)
    % Evaluates a pair of interaction parameters in the form required by
    % problem 3
    
    a1 = 0.05;
    a2 = 0.25;
    b = 0.01;
    c = 0.01;
    I = 0.1;
    tspan = [0:0.5:100];
    v10 = 0.1;
    v20 = 0.1;
    w10 = 0;
    w20 = 0;
    
    y0 = [v10; v20; w10; w20];
    [t, y] = ode15s(@(t, y) fitzhugh(t, y, a1, a2, b, c, I, d12, d21), tspan, y0);
    
    if show_plot
        
        hold on
        plot(t, y(:, 1))
        plot(t, y(:, 2))
        plot(t, y(:, 3))
        plot(t, y(:, 4))
        
        %{
        plot(y(:, 1), y(:, 3))
        plot(y(:, 2), y(:, 4))
        %}
    end
    
    answer = y;
    
end