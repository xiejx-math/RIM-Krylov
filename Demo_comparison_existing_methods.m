%% Figure 3: Comparison with existing randomized iterative methods
% This demo reproduces Figure 3 in the paper.
%
% It compares RABK, SCGP, RBK, Kaczmarz++, IS-Krylov-CS, and
% IS-Krylov-PS on synthetic consistent linear systems.
%
% The plots show the evolution of the relative solution error (RSE)
% with respect to the number of iterations and the CPU time.
% Three synthetic settings with different ranks and condition-number
% bounds can be selected below.
%
% RABK, SCGP, and IS-Krylov-PS use partition sampling, 
% whereas RBK, Kaczmarz++, and IS-Krylov-CS use uniform row-block sampling.
% For the partition-based methods, the partition-construction time is
% included in the reported CPU time.

clear
close all
clc

%% problem setup
m = 4096;          % number of rows
n = 1024;          % number of columns

% Choose one experimental setting.
% kappa = 40; r = 1024; Max_length = 80000; CPU_xmax = 12;
% kappa = 10; r = 1024; Max_length = 6500;  CPU_xmax = 2;
kappa = 10; r = 512;  Max_length = 3500;  CPU_xmax = 0.5;

run_time = 20;     % number of repeated runs
q = 64;            % row-block size

plot_tol = 1e-12;  % lower limit shown in the figures
stop_tol = 1e-14;  % stopping tolerance used by the solvers

% Kaczmarz++ parameters. In chol mode, the inner SRHT/LSQR parameters are
% ignored. KPP_use_RHT controls the outer RHT preprocessing.
KPP_reg = 1e-8;
% KPP_use_RHT = false;
KPP_use_RHT = true;
KPP_accelerated = true;
KPP_memoization = true;
KPP_inner_solver = 'chol';    % 'lsqr' or 'chol'
KPP_lsqr_maxit = 8;
KPP_lsqr_tol = 1e-8;
KPP_inner_tau = 2*q;
KPP_count_RHT_time = true;

%% arrays for storing numerical results
num_grid = Max_length + 1;

RABK_CPU = nan(run_time, num_grid);
RABK_error = nan(run_time, num_grid);

SCGP_CPU = nan(run_time, num_grid);
SCGP_error = nan(run_time, num_grid);

RBK_CPU = nan(run_time, num_grid);
RBK_error = nan(run_time, num_grid);

KPP_CPU = nan(run_time, num_grid);
KPP_error = nan(run_time, num_grid);

IS_Krylov_PS_CPU = nan(run_time, num_grid);
IS_Krylov_PS_error = nan(run_time, num_grid);

IS_Krylov_CS_CPU = nan(run_time, num_grid);
IS_Krylov_CS_error = nan(run_time, num_grid);

%% execute the algorithms run_time times
for ii = 1:run_time
    %% generate the matrix A
    [U, ~] = qr(randn(m, r), 0);
    [V, ~] = qr(randn(n, r), 0);
    D = diag(1 + (kappa - 1).*rand(r, 1));
    A = U*D*V';
    clear U V D

    %% generate the right-hand vector b
    [m, n] = size(A);
    x_true = randn(n, 1);
    b = A*x_true;
    xLS = lsqminnorm(A, b);

    %% parameter setup
    opts = struct;
    opts.xstar = xLS;
    opts.TOL1 = eps^2;
    opts.TOL = stop_tol;
    opts.Pre_iter = 50;
    opts.Max_iter = Max_length;

    % Parameters used only by Kaczmarz++.
    opts.KPP_reg = KPP_reg;
    opts.KPP_use_RHT = KPP_use_RHT;
    opts.KPP_accelerated = KPP_accelerated;
    opts.KPP_memoization = KPP_memoization;
    opts.KPP_inner_solver = KPP_inner_solver;
    opts.KPP_lsqr_maxit = KPP_lsqr_maxit;
    opts.KPP_lsqr_tol = KPP_lsqr_tol;
    opts.KPP_inner_tau = KPP_inner_tau;
    opts.KPP_count_RHT_time = KPP_count_RHT_time;

    % This controls only the KPP memoization budget. RBK and IS-Krylov-CS
    % are not forced to use the same number of sampled blocks; the intended
    % consistency is uniform sampling, not identical sample counts.
    % opts.KPP_max_memo_blocks = ceil(m/q);

    %% construct the partition used by RABK, SCGP, and IS-Krylov-PS
    tic
    tau = floor(m/q);
    opts.permS = randperm(m);

    % A row permutation does not change the linear system. We apply the
    % same permuted system to all six methods so that every method uses the
    % same coefficient matrix and right-hand side in each trial.
    A = A(opts.permS, :);
    b = b(opts.permS);

    Aarrs = cell(tau, 1);
    barrs = cell(tau, 1);
    blockAnormfro = zeros(tau, 1);
    normAfro = norm(A, 'fro')^2;

    for i = 1:tau
        if i == tau
            ps = ((i - 1)*q + 1):m;
        else
            ps = ((i - 1)*q + 1):(i*q);
        end

        Aarrs{i} = A(ps, :);
        barrs{i} = b(ps);
        blockAnormfro(i) = norm(Aarrs{i}, 'fro')^2;
    end

    prob = blockAnormfro/normAfro;
    opts.Aarrs = Aarrs;
    opts.barrs = barrs;
    opts.cumsumpro = cumsum(prob);
    opts.blockAnormfro = blockAnormfro;
    opts.probset = 1;
    Partition_CPU = toc;

    %% run algorithms
    [xRABK, OutRABK] = My_RABK(A, b, q, opts);
    [xSCGP, OutSCGP] = My_AmRABK(A, b, q, opts);
    [xRBK, OutRBK] = My_RBK(A, b, q, opts);
    [xKPP, OutKPP] = My_Kaczmarz_plus_plus(A, b, q, opts);
    [xIS_Krylov_PS, OutIS_Krylov_PS] = My_IS_Krylov_PS(A, b, q, opts);
    [xIS_Krylov_CS, OutIS_Krylov_CS] = My_IS_Krylov_CS(A, b, q, opts);

    %% store numerical results and pad the tails after early stopping
    % The partition time is included only for the three partition-based
    % methods. The uniform-sampling methods do not use this preprocessing.
    [RABK_error(ii, :), RABK_CPU(ii, :)] = ...
        pad_history(OutRABK, num_grid, Partition_CPU);
    [SCGP_error(ii, :), SCGP_CPU(ii, :)] = ...
        pad_history(OutSCGP, num_grid, Partition_CPU);
    [RBK_error(ii, :), RBK_CPU(ii, :)] = ...
        pad_history(OutRBK, num_grid, 0);
    [KPP_error(ii, :), KPP_CPU(ii, :)] = ...
        pad_history(OutKPP, num_grid, 0);
    [IS_Krylov_PS_error(ii, :), IS_Krylov_PS_CPU(ii, :)] = ...
        pad_history(OutIS_Krylov_PS, num_grid, Partition_CPU);
    [IS_Krylov_CS_error(ii, :), IS_Krylov_CS_CPU(ii, :)] = ...
        pad_history(OutIS_Krylov_CS, num_grid, 0);

    fprintf(['Run %d/%d, iterations: RABK = %d, SCGP = %d, ', ...
        'RBK = %d, Kaczmarz++ = %d, IS-Krylov-CS = %d, ', ...
        'IS-Krylov-PS = %d\n'], ...
        ii, run_time, OutRABK.iter, OutSCGP.iter, OutRBK.iter, ...
        OutKPP.iter, OutIS_Krylov_CS.iter, OutIS_Krylov_PS.iter);
end

%% compute curve statistics
iter_axis = 0:Max_length;

[min_RABK, max_RABK, q25_RABK, q75_RABK, med_RABK] = ...
    curve_stats(RABK_error);
[min_SCGP, max_SCGP, q25_SCGP, q75_SCGP, med_SCGP] = ...
    curve_stats(SCGP_error);
[min_RBK, max_RBK, q25_RBK, q75_RBK, med_RBK] = ...
    curve_stats(RBK_error);
[min_KPP, max_KPP, q25_KPP, q75_KPP, med_KPP] = ...
    curve_stats(KPP_error);
[min_ISCS, max_ISCS, q25_ISCS, q75_ISCS, med_ISCS] = ...
    curve_stats(IS_Krylov_CS_error);
[min_ISPS, max_ISPS, q25_ISPS, q75_ISPS, med_ISPS] = ...
    curve_stats(IS_Krylov_PS_error);

med_RABK_CPU = median(RABK_CPU, 1);
med_SCGP_CPU = median(SCGP_CPU, 1);
med_RBK_CPU = median(RBK_CPU, 1);
med_KPP_CPU = median(KPP_CPU, 1);
med_ISCS_CPU = median(IS_Krylov_CS_CPU, 1);
med_ISPS_CPU = median(IS_Krylov_PS_CPU, 1);

%% plot style
fig_pos = [100, 100, 560, 420];
axis_font = 14;
axis_lw = 1.2;
curve_lw = 1.2;
label_font = 16;
title_font = 15;
legend_font = 10;

%% plot: error versus iterations
figure('Position', fig_pos)
hold on
box on
set(gcf, 'Color', 'w')
set(gca, 'FontSize', axis_font, ...
    'LineWidth', axis_lw, ...
    'TickLabelInterpreter', 'latex', ...
    'Layer', 'top')

% Min-max and interquartile bands.
plot_band(iter_axis, min_RABK, max_RABK, q25_RABK, q75_RABK, 'k', .05, .10);
plot_band(iter_axis, min_SCGP, max_SCGP, q25_SCGP, q75_SCGP, 'c', .05, .10);
plot_band(iter_axis, min_RBK, max_RBK, q25_RBK, q75_RBK, 'm', .05, .10);
plot_band(iter_axis, min_KPP, max_KPP, q25_KPP, q75_KPP, 'g', .05, .10);
plot_band(iter_axis, min_ISCS, max_ISCS, q25_ISCS, q75_ISCS, 'b', .05, .10);
plot_band(iter_axis, min_ISPS, max_ISPS, q25_ISPS, q75_ISPS, 'r', .05, .10);

p1 = semilogy(iter_axis, med_RABK, 'k-', 'LineWidth', curve_lw, ...
    'DisplayName', 'RABK');
p2 = semilogy(iter_axis, med_SCGP, 'c-', 'LineWidth', curve_lw, ...
    'DisplayName', 'SCGP');
p3 = semilogy(iter_axis, med_RBK, 'm-', 'LineWidth', curve_lw, ...
    'DisplayName', 'RBK');
p4 = semilogy(iter_axis, med_KPP, 'g-', 'LineWidth', curve_lw, ...
    'DisplayName', 'Kaczmarz++');
p5 = semilogy(iter_axis, med_ISCS, 'b-', 'LineWidth', curve_lw, ...
    'DisplayName', 'IS-Krylov-CS');
p6 = semilogy(iter_axis, med_ISPS, 'r-', 'LineWidth', curve_lw, ...
    'DisplayName', 'IS-Krylov-PS');

set(gca, 'YScale', 'log')
ylim([plot_tol, 1])
xlim([0, Max_length])
ylabel('RSE', 'Interpreter', 'latex', 'FontSize', label_font)
xlabel('Number of iterations', 'Interpreter', 'latex', 'FontSize', label_font)

legend([p1 p2 p6 p3 p4 p5], ...
    {'RABK', 'SCGP', 'IS-Krylov-PS', ...
     'RBK', 'Kaczmarz$++$', 'IS-Krylov-CS'}, ...
    'Interpreter', 'latex', ...
    'NumColumns', 2, ...
    'Location', 'best', ...
    'FontSize', legend_font);

title(['$\kappa=$ ', num2str(kappa), ', $r=$ ', num2str(r)], ...
    'Interpreter', 'latex', 'FontSize', title_font)

%% plot: error versus CPU time
figure('Position', fig_pos)
hold on
box on
set(gcf, 'Color', 'w')
set(gca, 'FontSize', axis_font, ...
    'LineWidth', axis_lw, ...
    'TickLabelInterpreter', 'latex', ...
    'Layer', 'top')

% Min-max and interquartile bands.
plot_band(med_RABK_CPU, min_RABK, max_RABK, q25_RABK, q75_RABK, 'k', .05, .10);
plot_band(med_SCGP_CPU, min_SCGP, max_SCGP, q25_SCGP, q75_SCGP, 'c', .05, .10);
plot_band(med_RBK_CPU, min_RBK, max_RBK, q25_RBK, q75_RBK, 'm', .05, .10);
plot_band(med_KPP_CPU, min_KPP, max_KPP, q25_KPP, q75_KPP, 'g', .05, .10);
plot_band(med_ISCS_CPU, min_ISCS, max_ISCS, q25_ISCS, q75_ISCS, 'b', .05, .10);
plot_band(med_ISPS_CPU, min_ISPS, max_ISPS, q25_ISPS, q75_ISPS, 'r', .05, .10);

p1 = semilogy(med_RABK_CPU, med_RABK, 'k-', 'LineWidth', curve_lw, ...
    'DisplayName', 'RABK');
p2 = semilogy(med_SCGP_CPU, med_SCGP, 'c-', 'LineWidth', curve_lw, ...
    'DisplayName', 'SCGP');
p3 = semilogy(med_RBK_CPU, med_RBK, 'm-', 'LineWidth', curve_lw, ...
    'DisplayName', 'RBK');
p4 = semilogy(med_KPP_CPU, med_KPP, 'g-', 'LineWidth', curve_lw, ...
    'DisplayName', 'Kaczmarz++');
p5 = semilogy(med_ISCS_CPU, med_ISCS, 'b-', 'LineWidth', curve_lw, ...
    'DisplayName', 'IS-Krylov-CS');
p6 = semilogy(med_ISPS_CPU, med_ISPS, 'r-', 'LineWidth', curve_lw, ...
    'DisplayName', 'IS-Krylov-PS');

set(gca, 'YScale', 'log')
xlim([0, CPU_xmax])
ylim([plot_tol, 1])
ylabel('RSE', 'Interpreter', 'latex', 'FontSize', label_font)
xlabel('CPU time', 'Interpreter', 'latex', 'FontSize', label_font)

legend([p1 p2 p6 p3 p4 p5], ...
    {'RABK', 'SCGP', 'IS-Krylov-PS', ...
     'RBK', 'Kaczmarz$++$', 'IS-Krylov-CS'}, ...
    'Interpreter', 'latex', ...
    'NumColumns', 2, ...
    'Location', 'best', ...
    'FontSize', legend_font);

title(['$\kappa=$ ', num2str(kappa), ', $r=$ ', num2str(r)], ...
    'Interpreter', 'latex', 'FontSize', title_font)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err_row, time_row] = pad_history(Out, num_grid, time_offset)
%PAD_HISTORY Store an output history and pad the tail by the last value.

    if nargin < 3
        time_offset = 0;
    end

    err = Out.error(:)';
    times = time_offset + Out.times(:)';

    if isempty(err) || isempty(times)
        error('The output histories Out.error and Out.times must be nonempty.');
    end

    err_row = zeros(1, num_grid);
    time_row = zeros(1, num_grid);

    L = min([numel(err), numel(times), num_grid]);
    err_row(1:L) = err(1:L);
    time_row(1:L) = times(1:L);

    if L < num_grid
        err_row(L + 1:end) = err_row(L);
        time_row(L + 1:end) = time_row(L);
    end
end

function [ymin, ymax, yq25, yq75, ymed] = curve_stats(Y)
%CURVE_STATS Statistics over repeated runs. Rows are runs.

    ymin = min(Y, [], 1);
    ymax = max(Y, [], 1);
    yq25 = quantile(Y, 0.25, 1);
    yq75 = quantile(Y, 0.75, 1);
    ymed = median(Y, 1);
end

function plot_band(x, ymin, ymax, yq25, yq75, color_name, alpha_outer, alpha_inner)
%PLOT_BAND Plot min-max and interquartile uncertainty bands.

    h = fill([x fliplr(x)], [ymin fliplr(ymax)], ...
        color_name, 'EdgeColor', 'none');
    set(h, 'FaceAlpha', alpha_outer)

    h = fill([x fliplr(x)], [yq25 fliplr(yq75)], ...
        color_name, 'EdgeColor', 'none');
    set(h, 'FaceAlpha', alpha_inner)
end
