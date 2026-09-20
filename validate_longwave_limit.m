clear
clc
close all


%% ================================================================
%  LONG-WAVE-LIMIT ERROR VALIDATION
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

s_longwave = zeros(NOh,Nk);

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


        % Long-wave analytical root

        s_LW = longwave_root(k,Oh,nu);


        % Full viscous solution

        s_vis = solve_viscous_mode( ...
            k,Oh,nu,s_LW,previous_root,options);


        previous_root = s_vis;


        % Store

        s_longwave(j,n) = s_LW;

        s_viscous(j,n) = s_vis;


        % Errors

        absolute_error(j,n) = abs(s_vis-s_LW);

        relative_error(j,n) = ...
            abs(s_vis-s_LW)/abs(s_vis);

    end

end


%% ================================================================
%  Print error summary
%  ================================================================

fprintf('\n')
fprintf('===============================================================\n')
fprintf('LONG-WAVE-LIMIT ERROR\n')
fprintf('===============================================================\n\n')

fprintf('      Oh       Max relative error, k <= 0.2     Max relative error, all k\n')
fprintf('--------------------------------------------------------------------------\n')


long_wave_region = k_vec <= 0.20;


for j = 1:NOh

    error_small_k = ...
        max(relative_error(j,long_wave_region));


    error_all_k = ...
        max(relative_error(j,:));


    fprintf('%10.4g              %10.6f %%                  %10.6f %%\n', ...
        Oh_vec(j), ...
        100*error_small_k, ...
        100*error_all_k);

end


%% ================================================================
%  Estimate convergence order
%  ================================================================

fprintf('\n')
fprintf('Observed small-k error dependence\n')
fprintf('---------------------------------\n')


for j = 1:NOh

    k_small = ...
        k_vec(long_wave_region);


    E_small = ...
        relative_error(j,long_wave_region);


    valid = E_small > 0;


    k_fit = k_small(valid);

    E_fit = E_small(valid);


    if length(k_fit) >= 2

        coefficients = ...
            polyfit(log(k_fit),log(E_fit),1);


        p = coefficients(1);


        fprintf('Oh = %.4g:    E_LW ~ k^{%.4f}\n', ...
            Oh_vec(j),p);

    end

end


%% ================================================================
%  Linear-scale figure
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


    error_percent = ...
        100*relative_error(j,:);


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


    title(sprintf('%s Long-wave error, $Oh = %s$', ...
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
%  Logarithmic-scale figure
%  ================================================================

figure('Units','inches', ...
       'Position',[1 1 9.4 6.6], ...
       'Color','w')


tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact')


for j = 1:NOh

    nexttile


    error_percent = ...
        100*relative_error(j,:);


    semilogy(k_vec,error_percent, ...
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


    title(sprintf('%s Long-wave error, $Oh = %s$', ...
        panel{j},Oh_text{j}), ...
        'Interpreter','latex', ...
        'FontSize',11, ...
        'FontWeight','normal')


    xlim([0 1.2])

    xticks(0:0.2:1.2)


    format_axes(gca)

end


%% ================================================================
%  Local functions
%  ================================================================


function s = longwave_root(k,Oh,nu)

    b = Oh*(nu^2+1)*k^2;

    c = 0.5*(nu+1)*k^2;

    discriminant = b^2-4*c;


    s = ...
        (-b+sqrt(complex(discriminant)))/2;


    if imag(s) < 0

        s = conj(s);

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


        [x,~,exitflag] = ...
            fsolve(F,x0,options);


        if exitflag > 0

            s_trial = ...
                x(1)+1i*x(2);


            residual = ...
                abs(dispersion_relation( ...
                s_trial,k,Oh,nu));


            if residual < 1e-8

                if imag(s_trial) < 0

                    s_conjugate = ...
                        conj(s_trial);


                    residual_conjugate = ...
                        abs(dispersion_relation( ...
                        s_conjugate,k,Oh,nu));


                    if residual_conjugate < 1e-8

                        s_trial = ...
                            s_conjugate;

                    end

                end


                candidates(end+1) = ...
                    s_trial;

            end

        end

    end


    if isempty(candidates)

        error('No viscous root found for Oh = %g and k = %g.',Oh,k)

    end


    [~,index] = ...
        min(abs(candidates-s_reference));


    s = candidates(index);

end


function F = complex_system(x,k,Oh,nu)

    s = x(1)+1i*x(2);


    D = ...
        dispersion_relation(s,k,Oh,nu);


    if isfinite(real(D)) && isfinite(imag(D))

        F = ...
            [real(D); imag(D)];

    else

        F = ...
            [1e10; 1e10];

    end

end


function D = dispersion_relation(s,k,Oh,nu)

    lambda = ...
        sqrt(k^2+s/Oh);


    A_k = ...
        besseli(1,k*nu) ...
        /besselk(1,k*nu);


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
