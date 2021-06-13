function derivative_matrix = generate_partial_x_derivative(vec, delta)
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
    
    derivative_matrix = spdiags(...
        [ones(N^2,1) -ones(N^2,1) ones(N^2,1) -ones(N^2,1)], ... 
        [-N*(N-1) -N N N*(N-1)], ...
        N^2, ...
        N^2); 
    derivative_matrix = derivative_matrix/(2*delta);
end