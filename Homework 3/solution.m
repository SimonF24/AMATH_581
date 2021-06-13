% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5] = solution()
 [consoleout, A1, A2, A3, A4, A5] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5] = student_solution(dummy_argument)
    
    % keep these lines to load Fmat and permvec
    % DO NOT submit the mat files to Gradescope
    load Fmat.mat 
    load permvec.mat
    
    %% Problem 1
    
    max_val = 10;
    min_val = -10;
    n = 8;
    
    delta = (abs(max_val) + abs(min_val))/n;
    
    % Generate a sample vector of the correct dimensions
    x = linspace(min_val, max_val, n);
    vec = repmat(x', n, 1);
    
    laplacian = generate_2d_laplacian(vec, delta);
    partial_x_derivative = generate_partial_x_derivative(vec, delta);
    partial_y_derivative = generate_partial_y_derivative(vec, delta);
     
    A1 = full(laplacian); 
    A2 = full(partial_x_derivative);
    A3 = full(partial_y_derivative);
    
    %% Compare to Notes Version 
    
    run_comparison_section = false;
    
    % Notes Version
    if run_comparison_section
        m=8;    % N value in x and y directions
        n=m*m;  % total size of matrix

        e0=zeros(n,1);  % vector of zeros
        e1=ones(n,1);   % vector of ones

        e2=e1;    % copy the one vector
        e4=e0;    % copy the zero vector

        for j=1:m
            e2(m*j)=0;  % overwrite every m^th value with zero
            e4(m*j)=1;  % overwrite every m^th value with one
        end

        e3(2:n,1)=e2(1:n-1,1); e3(1,1)=e2(n,1); % shift to correct
        e5(2:n,1)=e4(1:n-1,1); e5(1,1)=e4(n,1); % positions

        % place diagonal elements

        matA=spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1], ...
            [-(n-m) -m -m+1 -1 0 1 m-1 m (n-m)],n,n);

        % Compare notes matrix A to my matrix A
    
        figure()
        spy(matA) % view the matrix structure
        figure()
        spy(A1)
    end

    % disp(all(laplacian == matA)) 
    % Doesn't work in script but saved for reference in command line
    
    %% Problem 2
    
    show_unencrypted_image = false;
    
    interior = Fmat(161:240,161:240);
    interior_array = mat2cell(interior, 20*ones(4,1), 20*ones(4,1));
    working_array = reshape(interior_array, [16 1]);
    permuted_linear_array = working_array(permvec);
    permuted_array = reshape(permuted_linear_array, [4 4]);
    Fmat(161:240, 161:240) = cell2mat(permuted_array);
    
    A4 = abs(Fmat);
    
    shifted_matrix = ifftshift(Fmat);
    transformed_matrix = ifft2(shifted_matrix);
    
    A5 = abs(transformed_matrix);
    
    if show_unencrypted_image
        imshow(uint8(A5))
    end
end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)