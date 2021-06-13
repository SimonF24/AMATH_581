function derivative_matrix = generate_partial_y_derivative(vec, delta)
    % Returns a matrix that takes the partial derivative with respect to x
    % of the vector vec as a sparse matrix
    % delta is expected to be the (uniform) discretization of space
    % vec is expected to be a vector sampling the values of a function over a
    % square domain with equal discretization in the x and y directions
    % where the values for given x values have been concatenated with each 
    % other 
    % Specifically the difference between the first two values of vec is
    % assumed to be the discretization of the space
    
    N = sqrt(length(vec));
    
    empty_pattern = zeros(N,1);
    empty_pattern(1) = 1;
    full_pattern = ones(N,1);
    full_pattern(N) = 0;
    empty_diag = repmat(empty_pattern, N, 1);
    full_diag = repmat(full_pattern, N, 1);
    derivative_matrix = spdiags([empty_diag -full_diag flipud(full_diag) flipud(-empty_diag)],[-(N-1) -1 1 N-1], N^2, N^2);
    
    derivative_matrix = derivative_matrix/(2*delta);
end