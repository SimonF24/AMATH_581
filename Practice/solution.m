% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = solution()
 [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = student_solution(dummy_argument)

    A = [34 45; 17 6];

    A1 = A;
    
    A = [1 2; -1 1];
    B = [2 0; 0 2];
    C = [2 0 -3; 0  0 -1];
    D = [1 2; 2 3; -1 0];
    x = [1; 0];
    y = [0; 1];
    z = [1; 2; -1];
    
    A2 = A+B;
    A3 = 3*x-4*y;
    A4 = A*x;
    A5 = B*(x-y);
    A6 = D*x;
    A7 = D*y+z;
    A8 = A*B;
    A9 = B*C;
    A10 = C*D;
    
    f = @(x) (-x-cos(x));
    df = @(x) (-1+sin(x));
    
    err = 1;
    n_iterations = 0;
    previous_x = -3; % this is the initial guess
    tol = 10^-6;
    n_x_values = [previous_x];
    while abs(err) > tol
        next_x = previous_x - f(previous_x)/df(previous_x);
        n_x_values = [n_x_values; next_x];
        err = f(next_x);
        previous_x = next_x;
        n_iterations = n_iterations + 1;
    end
    
    A11 = n_x_values;
    
    b_iterations = 0;
    b_x_values = [];
    err = 1;
    left = -3;
    right = 1;
    tol = 10^-6;
    while abs(err) > tol
        mid = (left + right)/2;
        mid_val = f(mid);
        if mid_val > 0
            left = mid;
        else % mid_val < 0 
            % don't need to deal with equality since loop will have broken
            % on previous loop unless first loop
            right = mid;
        end
        err = f(mid);
        b_x_values = [b_x_values; mid];
        b_iterations = b_iterations + 1;
    end
    A12 = b_x_values;
    
    A13 = [n_iterations b_iterations];
    
    %{
    A1
    A2
    A3
    A4
    A5
    A6
    A7
    A8
    A9
    A10
    A11
    A12
    A13
    %}
    
end

% your extra functions, if you need them, can be in another files (don't forget to upload them too!)