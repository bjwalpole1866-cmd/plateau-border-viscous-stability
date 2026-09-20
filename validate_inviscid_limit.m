clear
clc
close all


%% ================================================================
%  INVISCID-LIMIT ERROR VALIDATION
%  ================================================================


%% Plot style

set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')


%% Parameters

nu = 2/sqrt(3);

Oh_vec = [1e-4 1e-3 1e-2 1e-1];

k_vec = 0.05:0.05:1.20;

Nk = length(k_vec);

NOh = length(Oh_vec);


%% fsolve options

options = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-11, ...
    'OptimalityTolerance',1e-11, ...
    'StepTolerance',1e-11, ...
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',3000);


%% Storage

s_viscous = zeros(NOh,Nk);

s_inviscid = zeros(NOh,Nk);

absolute_error = zeros(NOh,Nk);

relative_error = zeros(NOh,Nk);


%% ================================================================
%  Calculate solutions
%  ================================================================

for j = 1:NOh

    Oh = Oh_vec(j);

    previous_root = [];

    fprintf('\nOh = %.4g\n',Oh)


    for n = 1:Nk

        k = k_vec(n);


        % Inviscid root

        s_inv = inviscid_root(k,nu);


        % Full viscous root

        s_vis = solve_viscous_mode( ...
            k,Oh,nu,s_inv,previous_root,options);


        previous_root = s_vis;


        % Store solutions

        s_inviscid(j,n) = s_inv;

        s_viscous(j,n) = s_vis;


        % Errors

        absolute_error(j,n) = abs(s_vis-s_inv);

        relative_error(j,n) = ...
            abs(s_vis-s_inv)/abs(s_inv);

    end

end


%% ================================================================
%  Print error summary
%  ================================================================

fprintf('\n')
fprintf('===============================================================\n')
fprintf('INVISCID-LIMIT ERROR\n')
fprintf('===============================================================\n\n')

fprintf('      Oh        Maximum absolute error    Maximum relative error\n')
fprintf('----------------------------------------------------------------\n')


for j = 1:NOh

    fprintf('%10.4g        %12.6e             %10.6f %%\n', ...
        Oh_vec(j), ...
        max(absolute_error(j,:)), ...
        100*max(relative_error(j,:)));

end


%% ================================================================
%  Estimate error dependence on Oh
%  ================================================================

maximum_relative_error = max(relative_error,[],2);


fprintf('\n')
fprintf('Observed dependence of error on Oh\n')
fprintf('----------------------------------\n')


for j = 1:NOh-1

    p = log(maximum_relative_error(j+1) ...
        /maximum_relative_error(j)) ...
        /log(Oh_vec(j+1)/Oh_vec(j));


    fprintf('Oh = %.4g -> %.4g:   p = %.4f\n', ...
        Oh_vec(j), ...
        Oh_vec(j+1), ...
        p);

end


%% ================================================================
%  Plot
%  ================================================================

figure('Units','inches', ...
       'Position',[1 1 9.4 6.6], ...
       'Color','w')


tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact')


panel = {'(a)','(b)','(c)','(d)'};

Oh_text = {'0.0001','0.001','0.01','0.1'};


for j = 1:NOh

    nexttile

    error_percent = 100*relative_error(j,:);


    plot(k_vec,error_percent, ...
        'b-o', ...
        'LineWidth',1.3, ...
        'MarkerSize',4.5, ...
        'MarkerFaceColor','b', ...
        'MarkerEdgeColor','b')


    xlabel('$k$ (Wavenumber)', ...
        'Interpreter','latex', ...
        'FontSize',11)


    ylabel('Relative error (\%)', ...
        'Interpreter','latex', ...
        'FontSize',11)


    title(sprintf('%s Inviscid error, $Oh = %s$', ...
        panel{j},Oh_text{j}), ...
        'Interpreter','latex', ...
        'FontSize',11, ...
        'FontWeight','normal')


    xlim([0 1.2])

    xticks(0:0.2:1.2)


    ymax = max(error_percent);


    if ymax > 0

        ylim([0 1.08*ymax])

    end


    format_axes(gca)

end


%% ================================================================
%  Local functions
%  ================================================================


function s = inviscid_root(k,nu)

    A = besseli(1,k*nu)/besselk(1,k*nu);

    B = besseli(1,k) ...
        - A*besselk(1,k);

    C = besseli(0,k) ...
        + A*besselk(0,k);

    Q = k ...
        *(1/(1-nu)-k^2) ...
        *(B/C);

    s = sqrt(complex(-Q));

    if imag(s) < 0

        s = -s;

    end

end


function s = solve_viscous_mode( ...
    k,Oh,nu,s_reference,previous_root,options)

    guesses = s_reference;


    if ~isempty(previous_root)

        guesses(end+1) = previous_root;

    end


    guesses(end+1) = s_reference + 0.01;

    guesses(end+1) = s_reference - 0.01;

    guesses(end+1) = s_reference + 0.01i;

    guesses(end+1) = s_reference - 0.01i;


    candidates = [];


    for q = 1:length(guesses)

        s0 = guesses(q);

        x0 = [real(s0); imag(s0)];

        F = @(x) complex_system(x,k,Oh,nu);

        [x,~,exitflag] = fsolve(F,x0,options);


        if exitflag > 0

            s_trial = x(1) + 1i*x(2);

            residual = ...
                abs(dispersion_relation(s_trial,k,Oh,nu));


            if residual < 1e-8

                if imag(s_trial) < 0

                    s_conjugate = conj(s_trial);

                    residual_conjugate = ...
                        abs(dispersion_relation( ...
                        s_conjugate,k,Oh,nu));


                    if residual_conjugate < 1e-8

                        s_trial = s_conjugate;

                    end

                end


                candidates(end+1) = s_trial;

            end

        end

    end


    if isempty(candidates)

        error('No root found for Oh = %g and k = %g.',Oh,k)

    end


    [~,index] = min(abs(candidates-s_reference));

    s = candidates(index);

end


function F = complex_system(x,k,Oh,nu)

    s = x(1) + 1i*x(2);

    D = dispersion_relation(s,k,Oh,nu);


    if isfinite(real(D)) && isfinite(imag(D))

        F = [real(D); imag(D)];

    else

        F = [1e10; 1e10];

    end

end


function D = dispersion_relation(s,k,Oh,nu)

    lambda = sqrt(k^2+s/Oh);


    A_k = ...
        besseli(1,k*nu)/besselk(1,k*nu);


    B_k = ...
        besseli(1,k) ...
        - A_k*besselk(1,k);


    C_k = ...
        besseli(0,k) ...
        + A_k*besselk(0,k);


    E_k = ...
        besseli(2,k) ...
        + A_k*besselk(2,k);


    A_lambda = ...
        besseli(1,lambda*nu) ...
        /besselk(1,lambda*nu);


    B_lambda = ...
        besseli(1,lambda) ...
        - A_lambda*besselk(1,lambda);


    D_lambda = ...
        besseli(0,lambda) ...
        + A_lambda*besselk(0,lambda) ...
        + besseli(2,lambda) ...
        + A_lambda*besselk(2,lambda);


    viscous_term = ...
        1 ...
        + E_k/C_k ...
        - (2*k*lambda/(lambda^2+k^2)) ...
        *(B_k/B_lambda) ...
        *(D_lambda/C_k);


    capillary_term = ...
        k ...
        *(1/(1-nu)-k^2) ...
        *(B_k/C_k) ...
        *((lambda^2-k^2)/(lambda^2+k^2));


    D = ...
        s^2 ...
        + Oh*s*k^2*viscous_term ...
        + capillary_term;

end


function format_axes(ax)

    ax.FontSize = 10;

    ax.FontWeight = 'normal';

    ax.LineWidth = 0.8;

    ax.TickDir = 'in';

    ax.TickLength = [0.012 0.012];

    ax.Box = 'on';

    ax.XGrid = 'on';

    ax.YGrid = 'on';

    ax.GridAlpha = 0.18;

    ax.Layer = 'top';

    ax.TickLabelInterpreter = 'latex';

end
