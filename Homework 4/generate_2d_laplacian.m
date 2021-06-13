function laplacian_matrix = generate_2d_laplacian(vec, delta)
    % Returns the 2D Laplacian matrix for the given vector vec as a sparse
    % matrix
    % delta is expected to be the (uniform) discretization of space
    % Applying this matrix to a vector with the same form as vec will
    % approximate the Laplcian 
    % vec is expected to be a vector sampling the values of a function over a
    % square domain with equal discretization in the x and y directions
    % where the values for given x values have been concatenated with each 
    % other
    % The Laplcian is calculated using a central difference approximation
    % and assuming periodic boundary conditions
    
    N = sqrt(length(vec)); % This is the number of elements of the discretization of space in the x and y directions
    
    A_D = spdiags(-4*ones(N^2,1), 0, N^2, N^2);
    A_x = spdiags(ones(N^2,4), [-N*(N-1) -N N N*(N-1)], N^2, N^2);
    
    empty_pattern = zeros(N,1);
    empty_pattern(1) = 1;
    full_pattern = ones(N,1);
    full_pattern(N) = 0;
    empty_diag = repmat(empty_pattern, N, 1);
    full_diag = repmat(full_pattern, N, 1);
    A_y = spdiags([empty_diag full_diag flipud(full_diag) flipud(empty_diag)], [-(N-1) -1 1 N-1], N^2, N^2);
        
    laplacian_matrix = A_D + A_x + A_y;
    laplacian_matrix = laplacian_matrix/delta^2;
end