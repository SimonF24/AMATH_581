% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3] = solution()
 [consoleout, A1, A2, A3] = evalc('student_solution(0)'); 
end

function [A1, A2, A3] = student_solution(dummy_argument)

    %% Part a
    
    % Throughout this part quantities with a "_f" suffix are in Fourier
    % space and quantities with a "_t" suffix are in time (normal) space
    % Also uppercase quantities are matrices while lowercase quantities are
    % vectors (with the exception of L, D1, and D2)
    
    show_movie = false;
    
    D1 = 0.1;
    D2 = 0.1;
    n = 64;
    m = 1;
    tspan = 0:0.5:4;
    beta = 1;
    
    x = linspace(-10,10,n+1);
    x = x(1:n);

    [X,Y] = meshgrid(x,x);
    
    L = 20;
    kx = (2*pi/L)*[0:(n/2-1) (-n/2):-1];
    kx(1) = 1e-6;
    [KX, KY] = meshgrid(kx, kx);
    K2 = KX.^2+KY.^2;
    k2 = reshape(K2, n^2, 1);
    
    U_t = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+1i*Y)-sqrt(X.^2+Y.^2));
    V_t = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+1i*Y)-sqrt(X.^2+Y.^2));
   
    U_f = fft2(U_t);
    V_f = fft2(V_t);
    
    u_f = reshape(U_f, n^2, 1);
    v_f = reshape(V_f, n^2, 1);
    
    y0 = [u_f; v_f];
    
    [~, y] = ode45(@(t, y) spectral_reaction_diffusion_system(t, y, D1, D2, k2, n, beta), tspan, y0);
    
    A1 = real(y); 
    A2 = imag(y);
    
    if show_movie
        fig = figure;
        fig.Visible = 'off';
        frames = length(tspan);
        M(frames) = struct('cdata',[],'colormap',[]);
        for t=1:frames
            u = y(t, 1:n^2);
            frame = real(ifft2(reshape(u, n, n)));
            pcolor(X, Y, frame)
            shading interp
            drawnow
            M(t) = getframe;
        end
        fig.Visible = 'on';
        movie(M)
    end
    
    %% Part b
    
    show_movie = true;
    
    n = 31;
    
    [Dcheb, z] = cheb(n-1);
    Dcheb2 = Dcheb^2;
    Dcheb2(1,:) = zeros(1, n);
    Dcheb2(n,:) = zeros(1, n); % This applies the Dirichlet boundary conditions
    
    I = eye(length(Dcheb2));
    
    lap_Mat = (4/L^2)*(kron(Dcheb2, I)+kron(I,Dcheb2)); 
    % The factor of 4/L^2 comes from chain ruling the scaling below through
    % the Laplacian
    
    x = L*z/2;
    
    [X, Y] = meshgrid(x, x);
    
    U = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+1i*Y)-sqrt(X.^2+Y.^2));
    V = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+1i*Y)-sqrt(X.^2+Y.^2));
    
    u = reshape(U, n^2, 1);
    v = reshape(V, n^2, 1);
    
    y0 = [u; v];
    
    [~, y] = ode45(@(t, y) chebyshev_reaction_diffusion_system(t, y, D1, D2, lap_Mat, n, beta), tspan, y0);
    
    A3 = y;
    
    if show_movie
        fig = figure;
        fig.Visible = 'off';
        frames = length(tspan);
        M(frames) = struct('cdata',[],'colormap',[]);
        for t=1:frames
            pcolor(X, Y, reshape(y(t,1:n^2), n, n))
            shading interp
            drawnow
            M(t) = getframe;
        end
        fig.Visible = 'on';
        movie(M)
    end
    
    
end
