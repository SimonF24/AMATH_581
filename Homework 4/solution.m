% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5] = solution()
 [consoleout, A1, A2, A3, A4, A5] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5] = student_solution(dummy_argument)

    run_main_code = true;

    if run_main_code
        %% Section a

        max_val = 10;
        min_val = -10;
        n = 64;
        tspan = 0:0.5:4;
        nu = 0.001;

        delta = (abs(max_val) + abs(min_val))/n;

        % Generate a sample vector of the correct dimensions
        x = linspace(min_val, max_val, n+1);
        x = x(1:n);
        vec = repmat(x', n, 1);

        A = generate_2d_laplacian(vec, delta);
        B = generate_partial_x_derivative(vec, delta);
        C = generate_partial_y_derivative(vec, delta);
        A(1, 1) = 2;
        [L, U, P] = lu(A);

        % Generating initial vector w0
    %     x = repelem(x', n); % Old method, generates same results as below
    %     y = vec;
        [X, Y] = meshgrid(x, x);
        x = reshape(X, n^2, 1);
        y = reshape(Y, n^2, 1);
        w0 = exp(-x.^2-(y.^2/20));

        tic;
        [~, w] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, "direct"), tspan, w0);
        A1_time = toc;
        A1 = w; 
        tic;
        [~, w] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, "LU"), tspan, w0);
        A2_time = toc;
        A2 = w;

        %% Section b

        tic;
        [~, w] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, "bicgstab"), tspan, w0);
        A3_time = toc;
        A3 = w;
        tic;
        [~, w] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, "gmres"), tspan, w0);
        A4_time = toc;
        A4 = w;

        %% Section c

        tic;
        [~, w] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, "FFT"), tspan, w0);
        A5_time = toc;
        A5 = w;
    
    end
    
    %% As Cool as it Gets
    
    run_movie_code = true;
    show_only_custom_movie = true;
    
    if run_movie_code
        
        if run_main_code
            [~, fastest_method_num] = min([A1_time A2_time A3_time A4_time A5_time]);
            if fastest_method_num == 1
                fastest_method = "direct";
            elseif fastest_method_num == 2
                fastest_method = "LU";
            elseif fastest_method_num == 3
                fastest_method = "bicgstab";
            elseif fastest_method_num == 4
                fastest_method = "gmres";
            elseif fastest_method_num == 5
                fastest_method = "FFT";
            end
            calculation_method = fastest_method; 
        else
            calculation_method = "FFT";
            % FFT Method is fastest by factor of ~6 over second fastest (LU) on
            % my machine
            A1 = 0;
            A2 = 0;
            A3 = 0;
            A4 = 0;
            A5 = 0;
            % Not setting A1-A5 causes an error message to pop up
        end
        
        max_val = 10;
        min_val = -10;
        n = 64;
        movie_tspan = 0:0.5:15;
        nu = 0.001;

        delta = (abs(max_val) + abs(min_val))/n;

        % Generate a sample vector of the correct dimensions
        side = linspace(min_val, max_val, n+1);
        side = side(1:n);
        vec = repmat(side', n, 1);

        A = generate_2d_laplacian(vec, delta);
        B = generate_partial_x_derivative(vec, delta);
        C = generate_partial_y_derivative(vec, delta);
        A(1, 1) = 2;
        [L, U, P] = lu(A);
        
        [X, Y] = meshgrid(side, side);
        x = reshape(X, n^2, 1);
        y = reshape(Y, n^2, 1);
    
        random_gaussians = 10; % Number of Gaussians to randomly place in case 4
        spacing = 1;

        % Suggested Movie initial vectors
        w01 = gaussian(x-spacing, y)-gaussian(x+spacing, y); % Two oppositely charged Gaussians next to each other
        w02 = gaussian(x-spacing, y)+gaussian(x+spacing, y); % Two same charged Gaussians next to each other
        w03 = gaussian(x-spacing, y)-gaussian(x+spacing, y); % Two pairs of oppositely charged Gaussians that collide
        w03 = w03+gaussian(x-spacing, y+3*spacing)-gaussian(x+spacing, y+3*spacing);
        w04 = zeros(length(x),1); % random_gaussians random Gaussians
        for index=1:random_gaussians
            w04 = w04 + gaussian(x-10+20*rand(1,1), y-10+20*rand(1,1));
        end
        
        % Custom movie intial vectors
        
        custom_indices_1 = [];
        % A
        custom_indices_1 = append_point_line(custom_indices_1, [-9 -5], [-9 5]);
        custom_indices_1 = append_point_line(custom_indices_1, [-7 -5], [-7 5]);
        custom_indices_1 = [custom_indices_1 [-8; 0] [-8; 5]]; 
        % M
        custom_indices_1 = append_point_line(custom_indices_1, [-5 -5], [-5 5]);
        custom_indices_1 = [custom_indices_1 [-4; 4] [-3; 3] [-2; 4]];
        custom_indices_1 = [custom_indices_1 [-3; 3]]; % Included second time for emphasis
        custom_indices_1 = append_point_line(custom_indices_1, [-1 -5], [-1 5]);
        % A
        custom_indices_1 = append_point_line(custom_indices_1, [1 -5], [1 5]);
        custom_indices_1 = append_point_line(custom_indices_1, [3 -5], [3 5]);
        custom_indices_1 = [custom_indices_1 [2; 0] [2; 5]];
        % T
        custom_indices_1 = append_point_line(custom_indices_1, [5 -5], [5 5]);
        custom_indices_1 = [custom_indices_1 [4; 5] [6; 5]];
        % H
        custom_indices_1 = append_point_line(custom_indices_1, [7 -5], [7 5]);
        custom_indices_1 = append_point_line(custom_indices_1, [9 -5], [9 5]);
        custom_indices_1 = [custom_indices_1 [8; 0]];
        w05 = zeros(length(x),1);
        for index=1:size(custom_indices_1, 2)
            w05 = w05 + gaussian(x-custom_indices_1(1,index), y-custom_indices_1(2,index));
        end
        
        custom_indices_2 = [];
        % 5
        custom_indices_2 = append_point_line(custom_indices_2, [-7 0], [-7 5]);
        custom_indices_2 = append_point_line(custom_indices_2, [-6 5], [-4 5]);
        custom_indices_2 = append_point_line(custom_indices_2, [-6 0], [-4 0]);
        custom_indices_2 = append_point_line(custom_indices_2, [-4 -1], [-4 -5]);
        custom_indices_2 = append_point_line(custom_indices_2, [-7 -5], [-5 -5]);
        % 8
        custom_indices_2 = append_point_line(custom_indices_2, [-2 -5], [-2 5]);
        custom_indices_2 = append_point_line(custom_indices_2, [2 -5], [2 5]);
        custom_indices_2 = append_point_line(custom_indices_2, [-1 0], [1 0]);
        custom_indices_2 = append_point_line(custom_indices_2, [-1 -5], [1 -5]);
        custom_indices_2 = append_point_line(custom_indices_2, [-1 5], [1 5]);
        % 1
        custom_indices_2 = append_point_line(custom_indices_2, [4 -5], [8 -5]);
        custom_indices_2 = append_point_line(custom_indices_2, [6 -4], [6 5]);
        custom_indices_2 = append_point_line(custom_indices_2, [4 3], [5 4]);
        w06 = zeros(length(x),1);
        for index=1:size(custom_indices_2, 2)
            w06 = w06 + gaussian(x-custom_indices_2(1,index), y-custom_indices_2(2,index));
        end
        
        custom_indices_3 = [];
        % X (I think it looks cool)
        custom_indices_3 = append_point_line(custom_indices_3, [-5 -5], [5 5]);
        custom_indices_3 = append_point_line(custom_indices_3, [-5 5], [5 -5]); 
        w07 = zeros(length(x),1);
        for index=1:size(custom_indices_3, 2)
            w07 = w07 + gaussian(x-custom_indices_3(1,index), y-custom_indices_3(2,index));
        end
        
        if ~show_only_custom_movie
            [~, w1] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w01);
            [~, w2] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w02);
            [~, w3] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w03);
            [~, w4] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w04);
        end
        [~, w5] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w05);
        [~, w6] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w06);
        [~, w7] = ode45(@(t,w) vorticity_equation(t, w, A, B, C, L, P, U, nu, calculation_method), movie_tspan, w07);
        
        if ~show_only_custom_movie
            % Suggested Movies
            for w_num=1:4
                if w_num == 1
                    w = w1;
                elseif w_num == 2
                    w = w2;
                elseif w_num == 3
                    w = w3;
                elseif w_num == 4
                    w = w4;
                end
                fig = figure;
                fig.Visible = 'off';
                frames = length(movie_tspan);
                M(frames) = struct('cdata',[],'colormap',[]);
                for t=1:frames
                    pcolor(X, Y, reshape(w(t,:), n, n))
                    drawnow
                    M(t) = getframe;
                end
                fig.Visible = 'on';
                movie(M)
            end
        end
        
        % Custom movie
        start_delay = 24; % Number of frames to pause on the first frame of each section (2 seconds)
        fig = figure('Name', 'Custom Movie', 'NumberTitle', 'off');
        fig.Visible = 'off';
        frames = length(movie_tspan);
        M(frames) = struct('cdata',[],'colormap',[]);
        for w_num=1:3
            if w_num == 1
                w = w5;
                colormap jet
            elseif w_num == 2
                w = w6;
                colormap bone
            elseif w_num == 3
                w = w7;
                colormap copper
            end
            for i=1:start_delay
                frame_num = i+(w_num-1)*start_delay+(w_num-1)*frames;
                pcolor(X, Y, reshape(w(1,:), n, n))
                drawnow
                M(frame_num) = getframe;
            end
            for t=1:frames
                pcolor(X, Y, reshape(w(t,:), n, n))
                drawnow
                frame_num = t+w_num*start_delay+(w_num-1)*frames;
                M(frame_num) = getframe;
            end
        end
        fig.Visible = 'on';
        movie(M)
        v = VideoWriter('video'); 
        % The compression results in noticeable loss but we leave it
        % compressed for space reasons
        v.FrameRate = 12; % Match MATLAB Frame Rate
        open(v)
        writeVideo(v, M)
    end
    
end

function vec_out = append_point_line(vec_in, start_point, end_point)
    % Works only for horizontal lines, vertical lines, and diagonal lines
    if start_point(1) == end_point(1) && start_point(2) ~= end_point(2)
        if end_point(2)-start_point(2) > 0
            step = 1;
        else
            step = -1;
        end
        new_points = [repelem(start_point(1), length(start_point(2):step:end_point(2)));start_point(2):step:end_point(2)];
    elseif start_point(1) ~= end_point(1) && start_point(2) == end_point(2)
        if end_point(1)-start_point(1) > 0
            step = 1;
        else
            step = -1;
        end
        new_points = [start_point(1):step:end_point(1);repelem(start_point(2),length(start_point(1):step:end_point(1)))];
    elseif abs(end_point(1)-start_point(1)) == abs(end_point(2)-start_point(2))
        if end_point(1)-start_point(1) > 0
            x_step = 1;
        else 
            x_step = -1;
        end
        if end_point(2)-start_point(2) > 0
            y_step = 1;
        else 
            y_step = -1;
        end
        new_points = [start_point(1):x_step:end_point(1);start_point(2):y_step:end_point(2)];
    else
        error_message = ['This function works only when start_point and end_point ', ... 
            'are along a horizontal line, vertical line, or diagonal line.'];
        ME = MException('MyComponent:IncompatiblePoints', error_message);
        throw(ME)
    end
    vec_out = [vec_in new_points];
end
function answer = gaussian(x, y)
    answer = exp(-x.^2-y.^2);
end