%% Table 1: Comparison on SuiteSparse matrices
% This demo reproduces one row of Table 1 in the paper.
%
% It compares RABK, SCGP, RBK, Kaczmarz++, IS-Krylov-CS, and
% IS-Krylov-PS on one selected matrix from the SuiteSparse Matrix
% Collection.
%
% Uncomment exactly one dataset below and comment out all the others.
% The variable tab reports the average number of iterations and average
% CPU time over 20 independent trials.
%
% For RABK, SCGP, and IS-Krylov-PS, the reported CPU time includes the partition-construction time.
%
% Required solver files on the MATLAB path:
%   My_RABK.m
%   My_AmRABK.m
%   My_RBK.m
%   My_Kaczmarz_plus_plus.m
%   My_IS_Krylov_CS.m
%   My_IS_Krylov_PS.m
%
% Required SuiteSparse data files:
%   abtaha2.mat, model1.mat, crew1.mat, WorldCities.mat,
%   lp_sctap1.mat, well1033.mat, cr42.mat, Franz1.mat,
%   GL7d11.mat, D_6.mat, rel6.mat, lp_ship04s.mat

clear
close all
clc

%% experiment setup
q = 30;                  % block size used in Table 1
ell = 50;                % number of previous iterates for IS-Krylov
run_time = 20;           % number of independent trials
stop_tol = 1e-12;        % stopping tolerance for RSE
Max_length = 1500000;    % maximum number of iterations

% Uncomment the following line to use a fixed random seed.
% rng(1, 'twister');

%% Kaczmarz++ parameters
KPP_reg = 1e-8;
KPP_use_RHT = true;
KPP_accelerated = true;
KPP_memoization = true;
KPP_inner_solver = 'chol';
KPP_lsqr_maxit = 8;
KPP_lsqr_tol = 1e-8;
KPP_inner_tau = 2*q;
KPP_count_RHT_time = true;

%% select one dataset
% Uncomment exactly one line.

data_file = 'abtaha2';

% data_file = 'model1';
% data_file = 'crew1';
% data_file = 'WorldCities';
% data_file = 'lp_sctap1';
% data_file = 'well1033';
% data_file = 'cr42';
% data_file = 'Franz1';
% data_file = 'GL7d11';
% data_file = 'D_6';
% data_file = 'rel6';
% data_file = 'lp_ship04s';

data_name = data_file;

%% load the selected dataset
fprintf('\n============================================================\n');
fprintf('Dataset: %s\n', data_name);
fprintf('============================================================\n');

data = load(data_file);

if ~isfield(data, 'Problem') || ~isfield(data.Problem, 'A')
    error('File %s.mat does not contain Problem.A.', data_file);
end

A_origin = data.Problem.A;
[m, n] = size(A_origin);

if q > m
    error('The block size q=%d exceeds the number of rows m=%d.', q, m);
end

if nnz(A_origin) == 0
    error('The coefficient matrix in %s.mat is zero.', data_file);
end

fprintf('Matrix size: %d-by-%d\n', m, n);

%% trial-level storage
Iter_RABK = nan(run_time, 1);
Iter_SCGP = nan(run_time, 1);
Iter_RBK = nan(run_time, 1);
Iter_KPP = nan(run_time, 1);
Iter_ISCS = nan(run_time, 1);
Iter_ISPS = nan(run_time, 1);

CPU_RABK = nan(run_time, 1);
CPU_SCGP = nan(run_time, 1);
CPU_RBK = nan(run_time, 1);
CPU_KPP = nan(run_time, 1);
CPU_ISCS = nan(run_time, 1);
CPU_ISPS = nan(run_time, 1);

%% run the experiments
for ii = 1:run_time
    %% construct a consistent linear system
    x_true = randn(n, 1);
    b_origin = A_origin*x_true;
    xLS = lsqminnorm(A_origin, b_origin);

    %% common options
    opts = struct;
    opts.xstar = xLS;
    opts.TOL = stop_tol;
    opts.TOL1 = eps^2;
    opts.Pre_iter = ell;
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

    % Uncomment to impose an explicit memoization cap.
    % opts.KPP_max_memo_blocks = ceil(m/q);

    %% construct a common partition for RABK, SCGP, and IS-Krylov-PS
    partition_timer = tic;

    permS = randperm(m);
    A_part = A_origin(permS, :);
    b_part = b_origin(permS);

    num_blocks = ceil(m/q);
    Aarrs = cell(num_blocks, 1);
    barrs = cell(num_blocks, 1);
    blockAnormfro = zeros(num_blocks, 1);

    for jj = 1:num_blocks
        first_row = (jj - 1)*q + 1;
        last_row = min(jj*q, m);
        rows = first_row:last_row;

        Aarrs{jj} = A_part(rows, :);
        barrs{jj} = b_part(rows);
        blockAnormfro(jj) = norm(Aarrs{jj}, 'fro')^2;
    end

    normAfro = sum(blockAnormfro);
    if normAfro == 0
        error('The coefficient matrix has zero Frobenius norm.');
    end

    prob = blockAnormfro/normAfro;
    Partition_CPU = toc(partition_timer);

    opts_part = opts;
    opts_part.permS = permS;
    opts_part.Aarrs = Aarrs;
    opts_part.barrs = barrs;
    opts_part.cumsumpro = cumsum(prob);
    opts_part.cumsumpro(end) = 1;
    opts_part.blockAnormfro = blockAnormfro;
    opts_part.probset = 1;

    %% RABK
    [~, Out_RABK] = My_RABK(A_part, b_part, q, opts_part);
    Iter_RABK(ii) = Out_RABK.iter;
    CPU_RABK(ii) = Partition_CPU + Out_RABK.times(end);

    %% SCGP
    [~, Out_SCGP] = My_AmRABK(A_part, b_part, q, opts_part);
    Iter_SCGP(ii) = Out_SCGP.iter;
    CPU_SCGP(ii) = Partition_CPU + Out_SCGP.times(end);

    %% RBK
    [~, Out_RBK] = My_RBK(A_origin, b_origin, q, opts);
    Iter_RBK(ii) = Out_RBK.iter;
    CPU_RBK(ii) = Out_RBK.times(end);

    %% Kaczmarz++
    [~, Out_KPP] = My_Kaczmarz_plus_plus(A_origin, b_origin, q, opts);
    Iter_KPP(ii) = Out_KPP.iter;
    CPU_KPP(ii) = Out_KPP.times(end);

    %% IS-Krylov-CS
    [~, Out_ISCS] = My_IS_Krylov_CS(A_origin, b_origin, q, opts);
    Iter_ISCS(ii) = Out_ISCS.iter;
    CPU_ISCS(ii) = Out_ISCS.times(end);

    %% IS-Krylov-PS
    [~, Out_ISPS] = My_IS_Krylov_PS(A_part, b_part, q, opts_part);
    Iter_ISPS(ii) = Out_ISPS.iter;
    CPU_ISPS(ii) = Partition_CPU + Out_ISPS.times(end);

    fprintf(['Run %2d/%2d: RABK=%g, SCGP=%g, RBK=%g, K++=%g, ', ...
        'IS-Krylov-CS=%g, IS-Krylov-PS=%g\n'], ...
        ii, run_time, ...
        Iter_RABK(ii), Iter_SCGP(ii), Iter_RBK(ii), ...
        Iter_KPP(ii), Iter_ISCS(ii), Iter_ISPS(ii));
end

%% compute average results
% The method order agrees with Table 1 of the paper.
trial_iter = [ ...
    Iter_RABK, ...
    Iter_SCGP, ...
    Iter_RBK, ...
    Iter_KPP, ...
    Iter_ISCS, ...
    Iter_ISPS];

trial_cpu = [ ...
    CPU_RABK, ...
    CPU_SCGP, ...
    CPU_RBK, ...
    CPU_KPP, ...
    CPU_ISCS, ...
    CPU_ISPS];

mean_iter = mean(trial_iter, 1);
mean_cpu = mean(trial_cpu, 1);

%% construct the final table
Method = [ ...
    "RABK"; ...
    "SCGP"; ...
    "RBK"; ...
    "Kaczmarz++"; ...
    "IS-Krylov-CS"; ...
    "IS-Krylov-PS"];

AverageIterations = mean_iter';
AverageCPU = mean_cpu';

tab = table(Method, AverageIterations, AverageCPU);

% Keep only the final result table in the MATLAB workspace.
clearvars -except tab

disp(tab)
