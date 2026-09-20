clear
clc
close all

%% ================================================================
%  PLOT STYLE
%  ================================================================

% Use LaTeX-style fonts throughout
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

% Font and line sizes chosen to match the attached figure
axesFontSize  = 11;
labelFontSize = 12;
titleFontSize = 12;
legendFontSize = 10;
curveLineWidth = 1.4;
markerSize     = 5;

% Figure dimensions
figureWidth  = 7.2;
figureHeight = 5.2;

%% ================================================================
%  PHYSICAL PARAMETERS
%  ================================================================

rho   = 1030;          % Density [kg/m^3]
gamma = 0.038;         % Surface tension [N/m]
mu    = 1.04e-3;       % Dynamic viscosity [Pa s]
a     = 5e-4;          % Plateau-border radius [m]

nu = 2/sqrt(3);

Oh = mu/sqrt(rho*gamma*a);

fprintf('Experimental Ohnesorge number = %.8f\n\n',Oh)

%% ================================================================
%  WAVENUMBERS
%  ================================================================

k_vec = 0.05:0.05:1.20;

%% ================================================================
%  SEARCH REGION IN COMPLEX s-PLANE
%  ================================================================

sigma_min = -8;
sigma_max =  2;
omega_max = 6;

N_sigma = 21;
N_omega = 31;

sigma_guess = linspace(sigma_min,sigma_max,N_sigma);
omega_guess = linspace(-omega_max,omega_max,N_omega);

%% ================================================================
%  ROOT ACCEPTANCE PARAMETERS
%  ================================================================

residual_tolerance = 1e-8;
root_tolerance     = 1e-5;
zero_tolerance     = 1e-7;

%% ================================================================
%  FSOLVE OPTIONS
%  ================================================================

options = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-11, ...
    'OptimalityTolerance',1e-11, ...
    'StepTolerance',1e-11, ...
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',3000);

%% ================================================================
%  STORAGE
%  ================================================================

modes = cell(length(k_vec),1);
number_modes = zeros(length(k_vec),1);
maximum_growth = NaN(length(k_vec),1);
frequency_of_maximum = NaN(length(k_vec),1);

previous_roots = [];

%% ================================================================
%  SYSTEMATIC SEARCH FOR ROOTS
%  ================================================================

for j = 1:length(k_vec)

    k = k_vec(j);

    fprintf('Searching k = %.3f ...\n',k)

    roots_found = [];

    %% Long-wave roots

    b = Oh*(nu^2 + 1)*k^2;

    c = 0.5*(nu + 1)*k^2;

    discriminant = b^2 - 4*c;

    s_LW_1 = (-b + sqrt(discriminant + 0i))/2;

    s_LW_2 = (-b - sqrt(discriminant + 0i))/2;

    %% Systematic grid of initial guesses

    initial_guesses = [];

    for p = 1:length(sigma_guess)

        for q = 1:length(omega_guess)

            initial_guesses(end+1) = ...
                sigma_guess(p) + 1i*omega_guess(q);

        end

    end

    %% Add long-wave roots

    initial_guesses(end+1) = s_LW_1;

    initial_guesses(end+1) = s_LW_2;

    %% Add previous roots for continuation in k

    for p = 1:length(previous_roots)

        initial_guesses(end+1) = previous_roots(p);

    end

    %% Solve from every initial guess

    for p = 1:length(initial_guesses)

        s0 = initial_guesses(p);

        x0 = [real(s0); imag(s0)];

        F = @(x) complex_system(x,k,Oh,nu);

        [x,~,exitflag] = fsolve(F,x0,options);

        if exitflag > 0

            s = x(1) + 1i*x(2);

            residual = ...
                abs(dispersion_relation(s,k,Oh,nu));

            inside_search_region = ...
                real(s) >= sigma_min && ...
                real(s) <= sigma_max && ...
                abs(imag(s)) <= omega_max;

            non_trivial_root = ...
                abs(s) > zero_tolerance;

            good_residual = ...
                residual < residual_tolerance;

            if inside_search_region && ...
               non_trivial_root && ...
               good_residual

                is_new_root = true;

                for m = 1:length(roots_found)

                    if abs(s - roots_found(m)) < root_tolerance

                        is_new_root = false;

                        break

                    end

                end

                if is_new_root

                    roots_found(end+1) = s;

                end

            end

        end

    end

    %% Sort roots from largest to smallest Re(s)

    if ~isempty(roots_found)

        [~,order] = ...
            sort(real(roots_found),'descend');

        roots_found = ...
            roots_found(order);

    end

    %% Store roots

    modes{j} = roots_found;

    number_modes(j) = ...
        length(roots_found);

    previous_roots = ...
        roots_found;

    %% Largest growth rate

    if ~isempty(roots_found)

        [maximum_growth(j),index_max] = ...
            max(real(roots_found));

        frequency_of_maximum(j) = ...
            imag(roots_found(index_max));

    end

    fprintf( ...
        '    Number of unique non-trivial roots = %d\n', ...
        number_modes(j))

    if ~isempty(roots_found)

        fprintf( ...
            '    Largest Re(s) = %.10e\n', ...
            maximum_growth(j))

    end

    fprintf('\n')

end

%% ================================================================
%  PRINT MODES
%  ================================================================

fprintf('\n')

fprintf('============================================================\n')

fprintf('MODES FOUND AT EACH WAVENUMBER\n')

fprintf('============================================================\n\n')

for j = 1:length(k_vec)

    fprintf('k = %.3f\n',k_vec(j))

    roots_found = modes{j};

    if isempty(roots_found)

        fprintf( ...
            '    No roots found in search region.\n')

    else

        for m = 1:length(roots_found)

            fprintf( ...
                ['    Mode %2d:  Re(s) = % .10e' ...
                 '    Im(s) = % .10e\n'], ...
                m, ...
                real(roots_found(m)), ...
                imag(roots_found(m)))

        end

    end

    fprintf('\n')

end

%% ================================================================
%  STABILITY SUMMARY
%  ================================================================

all_roots = [];

for j = 1:length(modes)

    all_roots = ...
        [all_roots modes{j}];

end

fprintf('\n')

fprintf('============================================================\n')

fprintf('STABILITY SUMMARY\n')

fprintf('============================================================\n\n')

if isempty(all_roots)

    fprintf( ...
        'No non-trivial roots were found.\n')

else

    largest_real_part = ...
        max(real(all_roots));

    fprintf( ...
        'Largest real part found = %.12e\n\n', ...
        largest_real_part)

    if largest_real_part < 0

        fprintf( ...
            'All non-trivial roots found have Re(s) < 0.\n')

        fprintf( ...
            'No unstable eigenvalue was detected in the searched region.\n')

    elseif largest_real_part > 0

        fprintf( ...
            'At least one root with Re(s) > 0 was found.\n')

        fprintf( ...
            'An unstable mode has therefore been detected.\n')

    else

        fprintf( ...
            'A neutrally stable root was found.\n')

    end

end

%% ================================================================
%  FIGURE 1
%
%  LARGEST GROWTH RATE
%  ================================================================

figure( ...
    'Units','inches', ...
    'Position',[1 1 figureWidth figureHeight], ...
    'Color','w')

% ---------------------------------------------------------------
% BLUE CURVE
% ---------------------------------------------------------------

plot( ...
    k_vec, ...
    maximum_growth, ...
    'Color',[0 0 1], ...
    'LineStyle','-', ...
    'LineWidth',curveLineWidth)

hold on

% ---------------------------------------------------------------
% BLUE CIRCULAR MARKERS
% ---------------------------------------------------------------

plot( ...
    k_vec, ...
    maximum_growth, ...
    'LineStyle','none', ...
    'Marker','o', ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor',[0 0 1], ...
    'MarkerEdgeColor',[0 0 1])

xlabel( ...
    '$k$ (Wavenumber)', ...
    'FontSize',labelFontSize)

ylabel( ...
    '$\max\,\Re(s)$', ...
    'FontSize',labelFontSize)

title( ...
    sprintf('Largest growth rate, $Oh = %.5f$',Oh), ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal')

xlim([0 1.2])

format_axes(gca,axesFontSize)

%% ================================================================
%  FIGURE 2
%
%  EIGENVALUES IN COMPLEX s-PLANE
%  ================================================================

figure( ...
    'Units','inches', ...
    'Position',[1 1 figureWidth figureHeight], ...
    'Color','w')

hold on

for j = 1:length(k_vec)

    roots_found = ...
        modes{j};

    if ~isempty(roots_found)

        plot( ...
            real(roots_found), ...
            imag(roots_found), ...
            'LineStyle','none', ...
            'Marker','o', ...
            'MarkerSize',markerSize, ...
            'MarkerFaceColor',[0 0 1], ...
            'MarkerEdgeColor',[0 0 1])

    end

end

xlabel( ...
    '$\Re(s)$', ...
    'FontSize',labelFontSize)

ylabel( ...
    '$\Im(s)$', ...
    'FontSize',labelFontSize)

title( ...
    sprintf('Modes found at $Oh = %.5f$',Oh), ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal')

format_axes(gca,axesFontSize)

%% ================================================================
%  FIGURE 3
%
%  REAL PART OF THE MODES
%  ================================================================

figure( ...
    'Units','inches', ...
    'Position',[1 1 figureWidth figureHeight], ...
    'Color','w')

hold on

% Store positive and negative frequency branches separately

k_positive = [];

real_positive = [];

k_negative = [];

real_negative = [];

for j = 1:length(k_vec)

    roots_found = ...
        modes{j};

    for m = 1:length(roots_found)

        if imag(roots_found(m)) >= 0

            k_positive(end+1) = ...
                k_vec(j);

            real_positive(end+1) = ...
                real(roots_found(m));

        else

            k_negative(end+1) = ...
                k_vec(j);

            real_negative(end+1) = ...
                real(roots_found(m));

        end

    end

end

% ---------------------------------------------------------------
% BLUE CURVE THROUGH THE ROOTS
% ---------------------------------------------------------------

plot( ...
    k_positive, ...
    real_positive, ...
    'Color',[0 0 1], ...
    'LineStyle','-', ...
    'LineWidth',curveLineWidth)

% ---------------------------------------------------------------
% BLUE NUMERICAL POINTS
% ---------------------------------------------------------------

plot( ...
    k_positive, ...
    real_positive, ...
    'LineStyle','none', ...
    'Marker','o', ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor',[0 0 1], ...
    'MarkerEdgeColor',[0 0 1])

xlabel( ...
    '$k$ (Wavenumber)', ...
    'FontSize',labelFontSize)

ylabel( ...
    '$\Re(s)$', ...
    'FontSize',labelFontSize)

title( ...
    sprintf('Real part of the growth rate, $Oh = %.5f$',Oh), ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal')

xlim([0 1.2])

format_axes(gca,axesFontSize)

%% ================================================================
%  FIGURE 4
%
%  IMAGINARY PART OF THE MODES
%  ================================================================

figure( ...
    'Units','inches', ...
    'Position',[1 1 figureWidth figureHeight], ...
    'Color','w')

hold on

k_upper = [];

omega_upper = [];

k_lower = [];

omega_lower = [];

for j = 1:length(k_vec)

    roots_found = ...
        modes{j};

    for m = 1:length(roots_found)

        if imag(roots_found(m)) >= 0

            k_upper(end+1) = ...
                k_vec(j);

            omega_upper(end+1) = ...
                imag(roots_found(m));

        else

            k_lower(end+1) = ...
                k_vec(j);

            omega_lower(end+1) = ...
                imag(roots_found(m));

        end

    end

end

% Sort the branches in increasing k

[k_upper,order] = ...
    sort(k_upper);

omega_upper = ...
    omega_upper(order);

[k_lower,order] = ...
    sort(k_lower);

omega_lower = ...
    omega_lower(order);

% ---------------------------------------------------------------
% BLUE CURVES
% ---------------------------------------------------------------

plot( ...
    k_upper, ...
    omega_upper, ...
    'Color',[0 0 1], ...
    'LineStyle','-', ...
    'LineWidth',curveLineWidth)

plot( ...
    k_lower, ...
    omega_lower, ...
    'Color',[0 0 1], ...
    'LineStyle','-', ...
    'LineWidth',curveLineWidth)

% ---------------------------------------------------------------
% BLUE NUMERICAL POINTS
% ---------------------------------------------------------------

plot( ...
    k_upper, ...
    omega_upper, ...
    'LineStyle','none', ...
    'Marker','o', ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor',[0 0 1], ...
    'MarkerEdgeColor',[0 0 1])

plot( ...
    k_lower, ...
    omega_lower, ...
    'LineStyle','none', ...
    'Marker','o', ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor',[0 0 1], ...
    'MarkerEdgeColor',[0 0 1])

xlabel( ...
    '$k$ (Wavenumber)', ...
    'FontSize',labelFontSize)

ylabel( ...
    '$\Im(s)$', ...
    'FontSize',labelFontSize)

title( ...
    sprintf('Imaginary part of the growth rate, $Oh = %.5f$',Oh), ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal')

xlim([0 1.2])

format_axes(gca,axesFontSize)

%% ================================================================
%  EXPERIMENTAL COMPARISON AT k = 0.4
%  ================================================================

k_experiment = 0.4;

[~,index_k] = ...
    min(abs(k_vec-k_experiment));

k_used = ...
    k_vec(index_k);

roots_experiment = ...
    modes{index_k};

oscillatory_modes = ...
    roots_experiment(abs(imag(roots_experiment)) > 1e-6);

if ~isempty(oscillatory_modes)

    [~,index_mode] = ...
        max(real(oscillatory_modes));

    s_experiment = ...
        oscillatory_modes(index_mode);

    T_scale = ...
        sqrt(rho*a^3/gamma);

    sigma = ...
        real(s_experiment)/T_scale;

    omega = ...
        abs(imag(s_experiment))/T_scale;

    t_decay = ...
        1/abs(sigma);

    c_wave = ...
        0.2;

    L_experiment = ...
        0.1;

    t_propagation = ...
        L_experiment/c_wave;

    L_decay = ...
        c_wave*t_decay;

    amplitude_ratio = ...
        exp(-L_experiment/L_decay);

    period = ...
        2*pi/omega;

    number_oscillations = ...
        t_decay/period;

    fprintf('\n')

    fprintf('============================================================\n')

    fprintf('EXPERIMENTAL COMPARISON\n')

    fprintf('============================================================\n\n')

    fprintf( ...
        'Oh                         = %.8f\n', ...
        Oh)

    fprintf( ...
        'k                          = %.3f\n', ...
        k_used)

    fprintf( ...
        's                          = %.10e %+.10ei\n', ...
        real(s_experiment), ...
        imag(s_experiment))

    fprintf( ...
        'Re(s)                      = %.10e\n', ...
        real(s_experiment))

    fprintf( ...
        '|Im(s)|                    = %.10e\n', ...
        abs(imag(s_experiment)))

    fprintf( ...
        'Capillary time T_scale     = %.6e s\n', ...
        T_scale)

    fprintf( ...
        'Dimensional damping sigma  = %.6e 1/s\n', ...
        sigma)

    fprintf( ...
        'Dimensional omega          = %.6e rad/s\n', ...
        omega)

    fprintf( ...
        'Oscillation period         = %.6e s\n', ...
        period)

    fprintf( ...
        'Decay time                 = %.6f s\n', ...
        t_decay)

    fprintf( ...
        'Propagation time over 0.1m = %.6f s\n', ...
        t_propagation)

    fprintf( ...
        'Attenuation length         = %.6f m\n', ...
        L_decay)

    fprintf( ...
        'A(0.1m)/A(0)               = %.6f\n', ...
        amplitude_ratio)

    fprintf( ...
        'Oscillations per decay time= %.6f\n', ...
        number_oscillations)

end

%% ================================================================
%  FUNCTIONS
%  ================================================================

function F = complex_system(x,k,Oh,nu)

    s = ...
        x(1) + 1i*x(2);

    D = ...
        dispersion_relation(s,k,Oh,nu);

    if isfinite(real(D)) && ...
       isfinite(imag(D))

        F = ...
            [real(D); imag(D)];

    else

        F = ...
            [1e10; 1e10];

    end

end


function D = dispersion_relation(s,k,Oh,nu)

    lambda = ...
        sqrt(k^2 + s/Oh);

    %% Bessel combinations depending on k

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

    %% Bessel combinations depending on lambda

    A_lambda = ...
        besseli(1,lambda*nu)/besselk(1,lambda*nu);

    B_lambda = ...
        besseli(1,lambda) ...
        - A_lambda*besselk(1,lambda);

    D_lambda = ...
        besseli(0,lambda) ...
        + A_lambda*besselk(0,lambda) ...
        + besseli(2,lambda) ...
        + A_lambda*besselk(2,lambda);

    %% Viscous contribution

    viscous_term = ...
        1 ...
        + E_k/C_k ...
        - (2*k*lambda/(lambda^2+k^2)) ...
        *(B_k/B_lambda) ...
        *(D_lambda/C_k);

    %% Capillary contribution

    capillary_term = ...
        k ...
        *(1/(1-nu) - k^2) ...
        *(B_k/C_k) ...
        *((lambda^2-k^2)/(lambda^2+k^2));

    %% Full viscous dispersion relation

    D = ...
        s^2 ...
        + Oh*s*k^2*viscous_term ...
        + capillary_term;

end


function format_axes(ax,fontSize)

    % Consistent style based on the attached figure

    ax.FontSize = fontSize;

    ax.FontWeight = 'normal';

    ax.LineWidth = 0.8;

    ax.TickDir = 'in';

    ax.TickLength = [0.012 0.012];

    ax.Box = 'on';

    ax.XGrid = 'on';

    ax.YGrid = 'on';

    ax.GridAlpha = 0.18;

    ax.MinorGridAlpha = 0.10;

    ax.Layer = 'top';

end
