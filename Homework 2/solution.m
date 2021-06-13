% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = solution()
 [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = student_solution(dummy_argument)
    % your solution code goes here
    % assign the variables you are asked to save here
    
    %% Problem 1
    
    A = 1; 
    eigenvalues = zeros(5, 1);
    K = 1;
    L = 4;
    last_sign = 1;
    show_plot = false;
    tol = 10^-4;
    xspan = -L:0.1:L;
    eps = 1;
    
    answers = zeros(length(xspan), 5);
    
    % The first guess is below the target
    for i=1:5 
        deps = 1;
        while true
            y0 = [A A*sqrt(K*L^2-eps)];
            [~, y] = ode45(@(t, y) quantum_harmonic_oscillator(t, y, K, eps), xspan, y0);
            
            if abs(y(end, 2) + sqrt(K*L^2-eps)*y(end,1)) < tol
                eigenvalues(i) = eps;
                answers(:,i) = abs(y(:,1)/sqrt(trapz(xspan,y(:,1).*y(:,1))));
                last_sign = -last_sign;
                break
            elseif last_sign*(y(end, 2) + sqrt(K*L^2-eps)*y(end,1)) > 0
                eps = eps + deps;
            else % last_sign*(y(end, 2) + sqrt(K*L^2-eps)*y(end,1)) < 0
                eps = eps - deps/2;
                deps = deps/2;
            end
        end
        eps = eps + 0.1;
    end
    
    A1 = answers(:, 1);
    A2 = answers(:, 2);
    A3 = answers(:, 3);
    A4 = answers(:, 4);
    A5 = answers(:, 5);
    A6 = eigenvalues;
    
    if show_plot
        figure()
        plot(xspan, A1)
        hold on
        plot(xspan, A2)
        plot(xspan, A3)
        plot(xspan, A4)
        plot(xspan, A5)
        xlabel('x')
        ylabel('|\phi(x)|')
        title('Shooting Method')
        legend('A1', 'A2', 'A3', 'A4', 'A5')
    end
    
    %% Problem 2
    
    K = 1;
    show_plot = true;
    deltax = 0.1;
    
    % Note that the below will break if xspan has less than 4 components
    %{
    A = zeros(length(xspan)-2);
    for i=1:length(A(:,1))
        if i==1
            A(1,1:2) = [2/3+deltax^2*K*xspan(i)^2 -2/3];
        elseif i==length(A(:,1))
            A(i,end-1:end) = [-2/3 2/3+deltax^2*K*xspan(end)^2];
        else
            A(i, i-1:i+1) = [-1 2+deltax^2*K*xspan(i)^2 -1];
        end
    end
    %}
    %{
    A = zeros(length(xspan)-2);
    for i=2:length(xspan)-1
        j = i-1;
        if i==2
            A(1,1:2) = [2/3+deltax^2*K*xspan(i)^2 -2/3];
        elseif i==length(xspan)-1
            A(j,end-1:end) = [-2/3 2/3+deltax^2*K*xspan(i)^2];
        else
            A(j, j-1:j+1) = [-1 2+deltax^2*K*xspan(i)^2 -1];
        end
    end
    
    [V,D] = eig(A);
    [d, I] = sort(diag(D/deltax^2));
    
    % These were the old answers with the old matrix A
    ans7 = abs(V(:, I(1))/sqrt(trapz(xspan(2:end-1), V(:, I(1)).*V(:, I(1)))));
    ans8 = abs(V(:, I(2))/sqrt(trapz(xspan(2:end-1), V(:, I(2)).*V(:, I(2)))));
    ans9 = abs(V(:, I(3))/sqrt(trapz(xspan(2:end-1), V(:, I(3)).*V(:, I(3)))));
    ans10 = abs(V(:, I(4))/sqrt(trapz(xspan(2:end-1), V(:, I(4)).*V(:, I(4)))));
    ans11 = abs(V(:, I(5))/sqrt(trapz(xspan(2:end-1), V(:, I(5)).*V(:, I(5)))));
    %}
    
    % A going from 2-N-1 approach
    
    A = diag(2+deltax^2*K*xspan(2:end-1).^2)+diag(-ones(length(xspan(3:end-1)),1),-1)+diag(-ones(length(xspan(2:end-2)),1),1);
    
    % Original/TA Slides
    A(1,1:2) = [2/3+deltax^2*K*xspan(2)^2 -2/3];
    A(end,end-1:end) = [-2/3 2/3+deltax^2*K*xspan(end-1)^2];
    
    % Alternate Assumptions

    % Full 2nd order forward and backwards
%     A(1,1:3) = [deltax^2*K*xspan(2)^2-1 2 -1];
%     A(end,end-2:end) = [-1 2 deltax^2*K*xspan(end-1)^2-1];

    % 2nd order forwards and backwards with assumption
%     A(1,1:3) = [0 4/3 -1/3];
%     A(end, end-2:end) = [-1/3 4/3 0];
    
    [V,D] = eig(A);
    [d, I] = sort(diag(D));
    
%     ans7 = abs(V(:, I(1))/sqrt(trapz(xspan(2:end-1), V(:, I(1)).*V(:, I(1)))));
%     ans8 = abs(V(:, I(2))/sqrt(trapz(xspan(2:end-1), V(:, I(2)).*V(:, I(2)))));
%     ans9 = abs(V(:, I(3))/sqrt(trapz(xspan(2:end-1), V(:, I(3)).*V(:, I(3)))));
%     ans10 = abs(V(:, I(4))/sqrt(trapz(xspan(2:end-1), V(:, I(4)).*V(:, I(4)))));
%     ans11 = abs(V(:, I(5))/sqrt(trapz(xspan(2:end-1), V(:, I(5)).*V(:, I(5)))));
    
    % With correction from slides
%     A7 = [1/3*(4*ans7(1)-ans7(2)); ans7; 1/3*(4*ans7(end)-ans7(end-1))];
%     A8 = [1/3*(4*ans8(1)-ans8(2)); ans8; 1/3*(4*ans8(end)-ans8(end-1))];
%     A9 = [1/3*(4*ans9(1)-ans9(2)); ans9; 1/3*(4*ans9(end)-ans9(end-1))];
%     A10 = [1/3*(4*ans10(1)-ans10(2)); ans10; 1/3*(4*ans10(end)-ans10(end-1))];
%     A11 = [1/3*(4*ans11(1)-ans11(2)); ans11; 1/3*(4*ans11(end)-ans11(end-1))];
%     A12 = d(1:5)/deltax^2;
    
    % Normalize again after bootstrapping
    
    % With full correction
%     A7 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(1))+3)*(4*ans7(2)-ans7(3)); ans7; 1/(2*deltax*sqrt(K*xspan(end)^2-d(1))+3)*(4*ans7(end-1)-ans7(end-2))];
%     A8 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(2))+3)*(4*ans8(2)-ans8(3)); ans8; 1/(2*deltax*sqrt(K*xspan(end)^2-d(2))+3)*(4*ans8(end-1)-ans8(end-2))];
%     A9 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(3))+3)*(4*ans9(2)-ans9(3)); ans9; 1/(2*deltax*sqrt(K*xspan(end)^2-d(3))+3)*(4*ans9(end-1)-ans9(end-2))];
%     A10 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(4))+3)*(4*ans10(2)-ans10(3)); ans10; 1/(2*deltax*sqrt(K*xspan(end)^2-d(4))+3)*(4*ans10(end-1)-ans10(end-2))];
%     A11 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(5))+3)*(4*ans11(2)-ans11(3)); ans11; 1/(2*deltax*sqrt(K*xspan(end)^2-d(5))+3)*(4*ans11(end-1)-ans11(end-2))];
    
    % Set answers to normalized eigenfunctions
%     ans7 = abs(V(:, I(1)))/sqrt(trapz(xspan(2:end-1), V(:, I(1)).*V(:, I(1))));
%     ans8 = abs(V(:, I(2)))/sqrt(trapz(xspan(2:end-1), V(:, I(2)).*V(:, I(2))));
%     ans9 = abs(V(:, I(3)))/sqrt(trapz(xspan(2:end-1), V(:, I(3)).*V(:, I(3))));
%     ans10 = abs(V(:, I(4)))/sqrt(trapz(xspan(2:end-1), V(:, I(4)).*V(:, I(4))));
%     ans11 = abs(V(:, I(5)))/sqrt(trapz(xspan(2:end-1), V(:, I(5)).*V(:, I(5))));
    
    % Set answers to eigenfunctions (previous approach)

    ans7 = V(:, I(1));
    ans8 = V(:, I(2));
    ans9 = V(:, I(3));
    ans10 = V(:, I(4));
    ans11 = V(:, I(5));
    
    % Full Forward and Backward Difference based on Boundary Condition
    % Use this one
  
    ans7 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(1))+3)*(4*ans7(1)-ans7(2)); ans7; 1/(2*deltax*sqrt(K*xspan(end)^2-d(1))+3)*(4*ans7(end)-ans7(end-1))];
    ans8 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(2))+3)*(4*ans8(1)-ans8(2)); ans8; 1/(2*deltax*sqrt(K*xspan(end)^2-d(2))+3)*(4*ans8(end)-ans8(end-1))];
    ans9 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(3))+3)*(4*ans9(1)-ans9(2)); ans9; 1/(2*deltax*sqrt(K*xspan(end)^2-d(3))+3)*(4*ans9(end)-ans9(end-1))];
    ans10 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(4))+3)*(4*ans10(1)-ans10(2)); ans10; 1/(2*deltax*sqrt(K*xspan(end)^2-d(4))+3)*(4*ans10(end)-ans10(end-1))];
    ans11 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(5))+3)*(4*ans11(1)-ans11(2)); ans11; 1/(2*deltax*sqrt(K*xspan(end)^2-d(5))+3)*(4*ans11(end)-ans11(end-1))];
    
    % Forward and Backward Difference not based on boundary condition
%     ans7 = [1/(d(1)*deltax^2-K*xspan(1)^2*deltax^2+1)*(2*ans7(1)-ans7(2)); ans7; 1/(d(1)*deltax^2-K*xspan(end)^2*deltax^2+1)*(2*ans7(end)-ans7(end-1))];
%     ans8 = [1/(d(2)*deltax^2-K*xspan(1)^2*deltax^2+1)*(2*ans8(1)-ans8(2)); ans8; 1/(d(2)*deltax^2-K*xspan(end)^2*deltax^2+1)*(2*ans8(end)-ans8(end-1))];
%     ans9 = [1/(d(3)*deltax^2-K*xspan(1)^2*deltax^2+1)*(2*ans9(1)-ans9(2)); ans9; 1/(d(3)*deltax^2-K*xspan(end)^2*deltax^2+1)*(2*ans9(end)-ans9(end-1))];
%     ans10 = [1/(d(4)*deltax^2-K*xspan(1)^2*deltax^2+1)*(2*ans10(1)-ans10(2)); ans10; 1/(d(4)*deltax^2-K*xspan(end)^2*deltax^2+1)*(2*ans10(end)-ans10(end-1))];
%     ans11 = [1/(d(5)*deltax^2-K*xspan(1)^2*deltax^2+1)*(2*ans11(1)-ans11(2)); ans11; 1/(d(5)*deltax^2-K*xspan(end)^2*deltax^2+1)*(2*ans11(end)-ans11(end-1))];

    % Calculate answers with normalization
    
    A7 = abs(ans7)/sqrt(trapz(xspan, ans7.*ans7));
    A8 = abs(ans8)/sqrt(trapz(xspan, ans8.*ans8));
    A9 = abs(ans9)/sqrt(trapz(xspan, ans9.*ans9));
    A10 = abs(ans10)/sqrt(trapz(xspan, ans10.*ans10));
    A11 = abs(ans11)/sqrt(trapz(xspan, ans11.*ans11));
    A12 = d(1:5)/deltax^2;
    
    % A going all the way approach
    %{
    A = diag(2+deltax^2*K*xspan.^2)+diag(-ones(length(xspan(2:end)),1),-1)+diag(-ones(length(xspan(1:end-1)),1),1);
    A(1,1:3) = [deltax^2*K*xspan(1)^2-1 2 -1];
    A(end,end-2:end) = [-1 2 deltax^2*K*xspan(end)^2];
    
    [V,D] = eig(A);
    [d, I] = sort(diag(D));
    
    ans7 = abs(V(:, I(1))/sqrt(trapz(xspan, V(:, I(1)).*V(:, I(1)))));
    ans8 = abs(V(:, I(2))/sqrt(trapz(xspan, V(:, I(2)).*V(:, I(2)))));
    ans9 = abs(V(:, I(3))/sqrt(trapz(xspan, V(:, I(3)).*V(:, I(3)))));
    ans10 = abs(V(:, I(4))/sqrt(trapz(xspan, V(:, I(4)).*V(:, I(4)))));
    ans11 = abs(V(:, I(5))/sqrt(trapz(xspan, V(:, I(5)).*V(:, I(5)))));
    
    % To see results without correction
%     A7 = abs(V(:, I(1))/sqrt(trapz(xspan, V(:, I(1)).*V(:, I(1)))));
%     A8 = abs(V(:, I(2))/sqrt(trapz(xspan, V(:, I(2)).*V(:, I(2)))));
%     A9 = abs(V(:, I(3))/sqrt(trapz(xspan, V(:, I(3)).*V(:, I(3)))));
%     A10 = abs(V(:, I(4))/sqrt(trapz(xspan, V(:, I(4)).*V(:, I(4)))));
%     A11 = abs(V(:, I(5))/sqrt(trapz(xspan, V(:, I(5)).*V(:, I(5)))));
    
    % Correcting first and last point
    A7 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(1))+3)*(4*ans7(2)-ans7(3)); ans7(2:end-1); 1/(2*deltax*sqrt(K*xspan(end)^2-d(1))+3)*(4*ans7(end-1)-ans7(end-2))];
    A8 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(2))+3)*(4*ans8(2)-ans8(3)); ans8(2:end-1); 1/(2*deltax*sqrt(K*xspan(end)^2-d(2))+3)*(4*ans8(end-1)-ans8(end-2))];
    A9 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(3))+3)*(4*ans9(2)-ans9(3)); ans9(2:end-1); 1/(2*deltax*sqrt(K*xspan(end)^2-d(3))+3)*(4*ans9(end-1)-ans9(end-2))];
    A10 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(4))+3)*(4*ans10(2)-ans10(3)); ans10(2:end-1); 1/(2*deltax*sqrt(K*xspan(end)^2-d(4))+3)*(4*ans10(end-1)-ans10(end-2))];
    A11 = [1/(2*deltax*sqrt(K*xspan(1)^2-d(5))+3)*(4*ans11(2)-ans11(3)); ans11(2:end-1); 1/(2*deltax*sqrt(K*xspan(end)^2-d(5))+3)*(4*ans11(end-1)-ans11(end-2))];
    A12 = d(1:5)/deltax^2;
    
    % Renormalizing
    A7 = abs(ans7)/sqrt(trapz(xspan, ans7.*ans7));
    A8 = abs(ans8)/sqrt(trapz(xspan, ans8.*ans8));
    A9 = abs(ans9)/sqrt(trapz(xspan, ans9.*ans9));
    A10 = abs(ans10)/sqrt(trapz(xspan, ans10.*ans10));
    A11 = abs(ans11)/sqrt(trapz(xspan, ans11.*ans11));
    A12 = d(1:5)/deltax^2;
    %}
    
    % Highest scoring submission
    K = 1;
    deltax = 0.1;
    old_A = zeros(length(xspan)-2);

    % Note that the below will break if xspan has less than 4 components
    %{
    for i=1:length(A(:,1))
        if i==1
            A(1,1:2) = [2/3+deltax^2*K*xspan(i)^2 -2/3];
        elseif i==length(A(:,1))
            A(i,end-1:end) = [-2/3 2/3+deltax^2*K*xspan(end)^2];
        else
            A(i, i-1:i+1) = [-1 2+deltax^2*K*xspan(i)^2 -1];
        end
    end
    %}

    for i=2:length(xspan)-1
        j = i-1;
        if i==2
            old_A(1,1:2) = [2/3+deltax^2*K*xspan(i)^2 -2/3];
        elseif i==length(xspan)-1
            old_A(j,end-1:end) = [-2/3 2/3+deltax^2*K*xspan(end-1)^2];
        else
            old_A(j, j-1:j+1) = [-1 2+deltax^2*K*xspan(i)^2 -1];
        end
    end

    [V,D] = eig(old_A);
    [d, I] = sort(diag(D/deltax^2));

    % These were the old answers with the old matrix A
    ans7 = abs(V(:, I(1))/sqrt(trapz(xspan(2:end-1), V(:, I(1)).*V(:, I(1)))));
    ans8 = abs(V(:, I(2))/sqrt(trapz(xspan(2:end-1), V(:, I(2)).*V(:, I(2)))));
    ans9 = abs(V(:, I(3))/sqrt(trapz(xspan(2:end-1), V(:, I(3)).*V(:, I(3)))));
    ans10 = abs(V(:, I(4))/sqrt(trapz(xspan(2:end-1), V(:, I(4)).*V(:, I(4)))));
    ans11 = abs(V(:, I(5))/sqrt(trapz(xspan(2:end-1), V(:, I(5)).*V(:, I(5)))));

    old_A7 = [1/3*(4*ans7(1)-ans7(2)); ans7; 1/3*(4*ans7(end)-ans7(end-1))];
    old_A8 = [1/3*(4*ans8(1)-ans8(2)); ans8; 1/3*(4*ans8(end)-ans8(end-1))];
    old_A9 = [1/3*(4*ans9(1)-ans9(2)); ans9; 1/3*(4*ans9(end)-ans9(end-1))];
    old_A10 = [1/3*(4*ans10(1)-ans10(2)); ans10; 1/3*(4*ans10(end)-ans10(end-1))];
    old_A11 = [1/3*(4*ans11(1)-ans11(2)); ans11; 1/3*(4*ans11(end)-ans11(end-1))];
    old_A12 = d(1:5);
    % End highest scoring submission
    
    show_old_v_new_plot = false;
    if show_old_v_new_plot
        figure()
        plot(xspan, A7)
        hold on
        plot(xspan, A8)
        plot(xspan, A9)
        plot(xspan, A10)
        plot(xspan, A11)
        plot(xspan, old_A7)
        plot(xspan, old_A8)
        plot(xspan, old_A9)
        plot(xspan, old_A10)
        plot(xspan, old_A11)
        legend('A7', 'A8', 'A9', 'A10', 'A11', 'Old A7', 'Old A8', 'Old A9', 'Old A10', 'Old A11')
    end
    
    if show_plot
        figure()
        plot(xspan, A7)
        hold on
        plot(xspan, A8)
        plot(xspan, A9)
        plot(xspan, A10)
        plot(xspan, A11)
        xlabel('x')
        ylabel('|\phi(x)|')
        title('Direct Method')
        legend('A7', 'A8', 'A9', 'A10', 'A11')
    end
    
    %% Comparison Between Problem 1 and 2
    
    show_comparison_plot = false;
    
    if show_comparison_plot
        figure()
        plot(xspan, A1)
        hold on
        plot(xspan, A2)
        plot(xspan, A3)
        plot(xspan, A4)
        plot(xspan, A5)
        plot(xspan, A7)
        plot(xspan, A8)
        plot(xspan, A9)
        plot(xspan, A10)
        plot(xspan, A11)
        legend('A1', 'A2', 'A3', 'A4', 'A5', 'A7', 'A8', 'A9', 'A10', 'A11')
    end
    
    show_old_comparison_plot = false;
    
    if show_old_comparison_plot
        figure()
        plot(xspan, A1)
        hold on
        plot(xspan, A2)
        plot(xspan, A3)
        plot(xspan, A4)
        plot(xspan, A5)
        plot(xspan, old_A7)
        plot(xspan, old_A8)
        plot(xspan, old_A9)
        plot(xspan, old_A10)
        plot(xspan, old_A11)
        legend('A1', 'A2', 'A3', 'A4', 'A5', 'Old A7', 'Old A8', 'Old A9', 'Old A10', 'Old A11')
    end
    
    %% Problem 3
    
    eigenvalues = zeros(2, 1);
    K = 1;
    L = 2;
    show_plot = false;
    tol = 10^-4;
    xspan = -L:0.1:L;
    
    answers = zeros(length(xspan), 2);
    
    for gamma=[-0.05 0.05]
        eps = 1;
        if gamma==-0.05
            last_sign = 1;
        elseif gamma==0.05
            last_sign = 1;
        end
        for i=1:2
            A = 1;
            deps = 1;
            while true
                y0 = [A A*sqrt(K*L^2-eps)];
                [x, y] = ode45(@(t,y) nonlinear_quantum_harmonic_oscillator(t, y, K, eps, gamma), xspan, y0);
                norm = sqrt(trapz(x, y(:,1).*y(:,1)));
                if abs(norm-1)<tol
                    
                else
                    A = A/norm;
                end
                
                y0 = [A A*sqrt(K*L^2-eps)];
                [x, y] = ode45(@(t,y) nonlinear_quantum_harmonic_oscillator(t, y, K, eps, gamma), xspan, y0);

                if abs(y(end, 2) + sqrt(K*L^2-eps)*y(end,1)) < tol
                    eigenvalues(i) = eps;
                    answers(:,i) = abs(y(:,1)/sqrt(trapz(x,y(:,1).*y(:,1))));
                    last_sign = -last_sign;
                    break
                elseif last_sign*(y(end, 2) + sqrt(K*L^2-eps)*y(end,1)) > 0
                    eps = eps + deps;
                else % last_sign*(y(end, 2) + sqrt(K*L^2-eps)*y(end,1)) < 0
                    eps = eps - deps/2;
                    deps = deps/2;
                end 
            end
            eps = eps + 0.1;

        end 
        
        if gamma==0.05
            ans13 = answers(:,1);
            ans14 = answers(:,2);
            ans15 = eigenvalues;
        elseif gamma==-0.05
            ans16 = answers(:,1);
            ans17 = answers(:,2);
            ans18 = eigenvalues;
        end
    end
    
    A13 = ans13;
    A14 = ans14;
    A15 = ans15;
    A16 = ans16;
    A17 = ans17;
    A18 = ans18;
    
    if show_plot
        figure()
        plot(xspan, A13)
        hold on
        plot(xspan, A14)
        plot(xspan, A16)
        plot(xspan, A17)
        legend('A13', 'A14', 'A16', 'A17')
    end
end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)