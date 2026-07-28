function [x, Out] = My_Kaczmarz_plus_plus(A, b, q, opts)
%MY_KACZMARZ_PLUS_PLUS_STORE  Practical Kaczmarz++ for Ax=b.
%
%   [x, Out] = My_Kaczmarz_plus_plus_store(A,b,q,opts)
%
% This MATLAB implementation is adapted from the public Python
% implementation of Kaczmarz++:
%
%   https://github.com/EdwinYang7/kaczmarz-plusplus
%
% The original Python implementation is distributed under the MIT License.
% This MATLAB version includes modifications for the numerical experiments
% in the accompanying paper.

% This implementation is intended for MATLAB comparison experiments.  It
% implements the practical Kaczmarz++ ingredients discussed in the paper:
%   1. optional row RHT preprocessing, A <- Q*A and b <- Q*b;
%   2. uniform row-block sampling from the transformed system;
%   3. Tikhonov-regularized block projections;
%   4. adaptive acceleration;
%   5. block memoization;
%   6. cached sampled blocks SA and Sb;
%   7. either exact Cholesky projection or SRHT-preconditioned LSQR.
%
% Default setting:
%   opts.KPP_inner_solver = 'lsqr', i.e., K++(LSQR-8)-type update.
%
% Common opts fields:
%   opts.xstar                       reference solution for error computation
%   opts.Max_iter                    maximum number of outer iterations
%   opts.TOL                         stopping tolerance
%   opts.x0                          initial point; default zeros(n,1)
%   opts.KPP_reg                     regularization parameter; default 1e-8
%   opts.KPP_use_RHT                 true/false; default true
%   opts.KPP_count_RHT_time          true/false; default true
%   opts.KPP_accelerated             true/false; default true
%   opts.KPP_memoization             true/false; default true
%   opts.KPP_cache_sampled_blocks    true/false; default true
%   opts.KPP_inner_solver            'lsqr' or 'chol'; default 'lsqr'
%   opts.KPP_lsqr_maxit              LSQR inner iterations; default 8
%   opts.KPP_lsqr_tol                LSQR tolerance; default 1e-8
%   opts.KPP_inner_tau               inner SRHT sketch size; default 2*q
%   opts.KPP_memoize_inner_precond   true/false; default true
%   opts.KPP_max_memo_blocks         optional cap on memoized blocks; default inf
%
% Output fields:
%   Out.error                        relative squared solution error history
%   Out.times                        CPU time history
%   Out.iter                         number of performed iterations
%   Out.num_memo_blocks              number of stored sampled blocks
%   Out.num_new_blocks               number of freshly sampled blocks
%   Out.num_reused_blocks            number of reused blocks
%   Out.KPP_*                        parameters actually used
%
% Notes:
%   In LSQR mode, the stored factor F is an SRHT-based approximate lower
%   Cholesky factor satisfying F*F' approx SA*SA' + lambda*I.  This is the
%   MATLAB counterpart of \tilde R[S] in Algorithm 1.1.  
%   If opts.KPP_memoize_inner_precond=false, the factor is recomputed whenever
%   a block is reused, which is closer to the current public Python code's
%   exact=false branch. 

    if nargin < 4
        opts = struct();
    end

    [m0, n] = size(A);
    b = b(:);
    q = min(q, m0);

    % General parameters.
    Max_iter = get_option(opts, 'Max_iter', 1000);
    TOL      = get_option(opts, 'TOL', 1e-12);
    xstar    = get_option(opts, 'xstar', []);

    % Kaczmarz++ parameters.
    lambda        = get_option(opts, 'KPP_reg', 1e-8);
    use_RHT       = get_option(opts, 'KPP_use_RHT', true);
    count_RHT     = get_option(opts, 'KPP_count_RHT_time', true);
    accelerated   = get_option(opts, 'KPP_accelerated', true);
    memoization   = get_option(opts, 'KPP_memoization', true);
    cache_blocks  = get_option(opts, 'KPP_cache_sampled_blocks', true);
    inner_solver  = lower(get_option(opts, 'KPP_inner_solver', 'lsqr'));
    lsqr_maxit    = get_option(opts, 'KPP_lsqr_maxit', 8);
    lsqr_tol      = get_option(opts, 'KPP_lsqr_tol', 1e-8);
    inner_tau     = get_option(opts, 'KPP_inner_tau', 2*q);
    memoize_precond = get_option(opts, 'KPP_memoize_inner_precond', true);
    max_memo_blocks = get_option(opts, 'KPP_max_memo_blocks', inf);

    if ~strcmp(inner_solver, 'lsqr') && ~strcmp(inner_solver, 'chol')
        error('opts.KPP_inner_solver must be either ''lsqr'' or ''chol''.');
    end

    % Initial point and reference solution.
    if isfield(opts, 'x0') && ~isempty(opts.x0)
        x = opts.x0(:);
    else
        x = zeros(n, 1);
    end
    x0 = x;

    if isempty(xstar)
        % Fallback for standalone testing.  In numerical experiments, pass
        % opts.xstar to avoid this extra least-squares solve.
        xstar = lsqminnorm(A, b);
    else
        xstar = xstar(:);
    end

    % Timer.  If count_RHT=true, Out.times includes the RHT preprocessing
    % cost.  Otherwise, the clock starts after RHT preprocessing.
    t_total = tic;

    % Optional row RHT preprocessing.
    if use_RHT
        % [Awork, bwork, rht_info] = apply_row_rht_with_padding(A, b);
        [Awork, bwork, rht_info] = My_RHT_sketch_C(A, b);
    else
        Awork = A;
        bwork = b;
        rht_info = struct('used', false, 'm_original', m0, 'm_work', m0, ...
                          'num_padded_rows', 0);
    end

    preprocess_time = toc(t_total);
    if count_RHT
        t_start = t_total;
        initial_time = preprocess_time;
    else
        t_start = tic;
        initial_time = 0;
    end

    [m, ncheck] = size(Awork);
    if ncheck ~= n
        error('Unexpected column dimension after RHT preprocessing.');
    end
    q = min(q, m);
    inner_tau = min(inner_tau, nextpow2_length(n));

    % Error/time histories.
    Out.error = zeros(1, Max_iter + 1);
    Out.times = zeros(1, Max_iter + 1);
    Out.error(1) = relative_squared_error(x, x0, xstar);
    Out.times(1) = initial_time;

    % Adaptive acceleration parameters, following the public K++ structure.
    z_param = double(accelerated);
    eta = 0.0;
    nu = 2*n/q;
    skip = max(1, round(n/q + 1));
    cnt = 1.0;
    ratio = 0.0;
    dist_old = 0.0;
    dist_new = 0.0;
    z = zeros(n, 1);

    % Memoization storage.
    block_indices = cell(0, 1);
    SA_blocks = cell(0, 1);
    Sb_blocks = cell(0, 1);
    factor_list = cell(0, 1);

    num_new_blocks = 0;
    num_reused_blocks = 0;
    last_lsqr_flag = NaN;
    last_lsqr_relres = NaN;
    last_lsqr_iter = NaN;

    iter = Max_iter;

    for k = 1:Max_iter
        % Online block memoization.  Draw a fresh block with decreasing
        % probability; otherwise reuse a previously sampled block.
        if ~memoization || isempty(block_indices)
            draw_new = true;
        elseif numel(block_indices) >= max_memo_blocks
            draw_new = false;
        else
            p_new = min(1, n*log(max(n, 2))/(q*k));
            draw_new = (rand < p_new);
        end

        if draw_new
            num_new_blocks = num_new_blocks + 1;

            idx = randperm(m, q);
            SA = Awork(idx, :);
            Sb = bwork(idx);

            if strcmp(inner_solver, 'chol')
                G = SA*SA' + lambda*eye(q);
                F = robust_chol_lower(G);
            else
                if memoize_precond
                    F = inner_srht_preconditioner(SA, inner_tau, lambda);
                else
                    F = [];
                end
            end

            if memoization
                block_indices{end + 1} = idx; %#ok<AGROW>
                if cache_blocks
                    SA_blocks{end + 1} = SA; %#ok<AGROW>
                    Sb_blocks{end + 1} = Sb; %#ok<AGROW>
                else
                    SA_blocks{end + 1} = []; %#ok<AGROW>
                    Sb_blocks{end + 1} = []; %#ok<AGROW>
                end
                factor_list{end + 1} = F; %#ok<AGROW>
            end
        else
            num_reused_blocks = num_reused_blocks + 1;

            memo_id = randi(numel(block_indices));
            idx = block_indices{memo_id};

            if cache_blocks && ~isempty(SA_blocks{memo_id})
                SA = SA_blocks{memo_id};
                Sb = Sb_blocks{memo_id};
            else
                SA = Awork(idx, :);
                Sb = bwork(idx);
            end

            F = factor_list{memo_id};
            if strcmp(inner_solver, 'lsqr') && (~memoize_precond || isempty(F))
                % GitHub-aligned option: recompute the inner SRHT
                % preconditioner rather than memoizing \tilde R[S].
                F = inner_srht_preconditioner(SA, inner_tau, lambda);
            end
        end

        rS = SA*x - Sb;

        if strcmp(inner_solver, 'chol')
            % Exact regularized projection.
            u = F' \ (F \ rS);
            w = SA' * u;
        else
            % SRHT-preconditioned LSQR projection.
            [w, last_lsqr_flag, last_lsqr_relres, last_lsqr_iter] = ...
                inner_regularized_lsqr(SA, rS, lambda, F, lsqr_tol, lsqr_maxit);
        end

        % Adaptive momentum update.
        z = z_param*(z + w);
        x = x - w - eta*z;

        Out.error(k + 1) = relative_squared_error(x, x0, xstar);
        Out.times(k + 1) = toc(t_start);

        if mod(k, 2*skip) <= skip
            dist_old = dist_old + norm(rS)^2;
        else
            dist_new = dist_new + norm(rS)^2;
        end

        if Out.error(k + 1) <= TOL
            iter = k;
            break
        end

        if accelerated && mod(k, 2*skip) == 0
            a_old = cnt^log(cnt);
            a_new = (cnt + 1)^log(cnt + 1);

            if dist_old <= eps
                local_ratio = 1;
            else
                local_ratio = min(1, dist_new/dist_old);
            end

            ratio = ratio*(a_old/a_new) + local_ratio*(1 - a_old/a_new);
            rho = max(0, 1 - ratio^(1/skip));
            z_param = (1 - rho)/(1 + rho);
            eta = 1/nu;

            cnt = cnt + 1;
            dist_old = 0.0;
            dist_new = 0.0;
        end
    end

    Out.error = Out.error(1:iter + 1);
    Out.times = Out.times(1:iter + 1);
    Out.iter = iter;
    Out.num_memo_blocks = numel(block_indices);
    Out.num_new_blocks = num_new_blocks;
    Out.num_reused_blocks = num_reused_blocks;
    Out.KPP_reg = lambda;
    Out.KPP_use_RHT = use_RHT;
    Out.KPP_count_RHT_time = count_RHT;
    Out.KPP_accelerated = accelerated;
    Out.KPP_memoization = memoization;
    Out.KPP_cache_sampled_blocks = cache_blocks;
    Out.KPP_inner_solver = inner_solver;
    Out.KPP_lsqr_maxit = lsqr_maxit;
    Out.KPP_lsqr_tol = lsqr_tol;
    Out.KPP_inner_tau = inner_tau;
    Out.KPP_memoize_inner_precond = memoize_precond;
    Out.KPP_max_memo_blocks = max_memo_blocks;
    Out.last_lsqr_flag = last_lsqr_flag;
    Out.last_lsqr_relres = last_lsqr_relres;
    Out.last_lsqr_iter = last_lsqr_iter;
    Out.RHT_info = rht_info;
    Out.preprocess_time = preprocess_time;
end

function [w, flag, relres, iter] = inner_regularized_lsqr(SA, rS, lambda, F, tol, maxit)
%INNER_REGULARIZED_LSQR  SRHT-preconditioned LSQR projection step.
%
% F is a lower factor satisfying approximately
%     F*F' approx SA*SA' + lambda*I.
%
% The routine approximates
%     w = SA'*(SA*SA' + lambda*I)^{-1}*rS
% by solving the augmented least-squares problem
%     [SA, sqrt(lambda)*I] y = rS
% after left preconditioning by F^{-1}, and returning y(1:n).

    [s, n] = size(SA);
    B = [SA, sqrt(max(lambda, 0))*eye(s)];
    M = F \ B;
    rhs = F \ rS;

    [y, flag, relres, iter] = lsqr(M, rhs, tol, maxit);
    w = y(1:n);
end

function F = inner_srht_preconditioner(SA, tau, lambda)
%INNER_SRHT_PRECONDITIONER  SRHT approximate lower Cholesky factor.
%
% Given SA in R^{s x n}, construct Ahat = sqrt(N/tau)*(SA*D*H)_J,
% where N is the padded power-of-two column dimension.  Return a lower
% factor F such that F*F' approx Ahat*Ahat' + lambda*I.

    [s, n] = size(SA);
    N = nextpow2_length(n);
    tau = min(tau, N);

    Z = zeros(s, N);
    Z(:, 1:n) = SA;

    signs = 2*(rand(1, N) > 0.5) - 1;
    Z = Z .* signs;
    Z = fwht_dim2(Z) / sqrt(N);

    cols = randperm(N, tau);
    Ahat = sqrt(N/tau) * Z(:, cols);

    G = Ahat*Ahat' + lambda*eye(s);
    F = robust_chol_lower(G);
end

% function [Aout, bout, info] = apply_row_rht_with_padding(A, b)
% %APPLY_ROW_RHT_WITH_PADDING  Apply normalized row RHT to [A,b].
% %
% % If m is not a power of two, zero rows are appended to reach the next
% % power-of-two length.  The transformation is orthogonal on the padded
% % system and therefore preserves the solution set of Ax=b.
% 
%     [m, n] = size(A);
%     M = nextpow2_length(m);
% 
%     Aout = zeros(M, n);
%     bout = zeros(M, 1);
%     Aout(1:m, :) = A;
%     bout(1:m) = b;
% 
%     signs = 2*(rand(M, 1) > 0.5) - 1;
%     Aout = Aout .* signs;
%     bout = bout .* signs;
% 
%     Aout = fwht_dim1(Aout) / sqrt(M);
%     bout = fwht_dim1(bout) / sqrt(M);
% 
%     info = struct('used', true, 'm_original', m, 'm_work', M, ...
%                   'num_padded_rows', M - m);
% end

function X = fwht_dim1(X)
%FWHT_DIM1  Unnormalized Walsh-Hadamard transform along rows.
    N = size(X, 1);
    h = 1;
    while h < N
        for i = 1:2*h:N
            idx1 = i:(i+h-1);
            idx2 = (i+h):(i+2*h-1);
            U = X(idx1, :);
            V = X(idx2, :);
            X(idx1, :) = U + V;
            X(idx2, :) = U - V;
        end
        h = 2*h;
    end
end

function X = fwht_dim2(X)
%FWHT_DIM2  Unnormalized Walsh-Hadamard transform along columns.
    N = size(X, 2);
    h = 1;
    while h < N
        for j = 1:2*h:N
            idx1 = j:(j+h-1);
            idx2 = (j+h):(j+2*h-1);
            U = X(:, idx1);
            V = X(:, idx2);
            X(:, idx1) = U + V;
            X(:, idx2) = U - V;
        end
        h = 2*h;
    end
end

function L = robust_chol_lower(G)
%ROBUST_CHOL_LOWER  Lower Cholesky factor with diagonal jitter if needed.
    G = (G + G')/2;
    s = size(G, 1);
    jitter = 0;
    scale = max(1, norm(G, 'fro'));

    for attempt = 1:10 
        [L, p] = chol(G + jitter*eye(s), 'lower');
        if p == 0
            return
        end
        if jitter == 0
            jitter = 1e-14*scale;
        else
            jitter = 10*jitter;
        end
    end

    % Last-resort fallback: project eigenvalues to a tiny positive floor.
    [V, D] = eig(G);
    d = real(diag(D));
    d = max(d, 1e-14*scale);
    Gpd = V*diag(d)*V';
    Gpd = real((Gpd + Gpd')/2);

    [L, p] = chol(Gpd, 'lower');
    if p ~= 0
        L = chol(Gpd + 1e-12*scale*eye(s), 'lower');
    end
end

function err = relative_squared_error(x, x0, xstar)
%RELATIVE_SQUARED_ERROR  ||x-xstar||^2 / ||x0-xstar||^2.
    denom = norm(x0 - xstar)^2;
    if denom <= eps
        denom = max(norm(xstar)^2, 1);
    end
    err = norm(x - xstar)^2 / denom;
end

function N = nextpow2_length(n)
%NEXTPOW2_LENGTH  Smallest power of two not smaller than n.
    N = 2^nextpow2(n);
end

function val = get_option(opts, field_name, default_val)
    if isfield(opts, field_name) && ~isempty(opts.(field_name))
        val = opts.(field_name);
    else
        val = default_val;
    end
end
