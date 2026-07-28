function [x, Out] = My_RBK(A, b, q, opts)
%MY_RANDOMIZED_BLOCK_KACZMARZ_UNIFORM Randomized block Kaczmarz with uniform row-block sampling.
%
%   [x, Out] = My_randomized_block_kaczmarz_uniform(A,b,q,opts)
%   samples a size-q row subset uniformly without replacement at each
%   iteration, i.e., tau = randperm(m,q).  This matches the uniform
%   row-subset sampling used by Kaczmarz++ for generating new blocks.
%
%   Required fields in opts:
%       opts.xstar    - reference solution for error computation
%       opts.Max_iter - maximum number of iterations
%   Optional fields:
%       opts.TOL      - stopping tolerance. Default: 1e-12
%       opts.x0       - initial point. Default: zeros(n,1)
%       opts.initial  - alternative name for the initial point

    [m, n] = size(A);

    if isfield(opts, 'x0') && ~isempty(opts.x0)
        x = opts.x0;
    elseif isfield(opts, 'initial') && ~isempty(opts.initial)
        x = opts.initial;
    else
        x = zeros(n, 1);
    end

    if isfield(opts, 'xstar') && ~isempty(opts.xstar)
        xstar = opts.xstar;
    else
        xstar = lsqminnorm(A, b);
    end

    if isfield(opts, 'Max_iter')
        Max_iter = opts.Max_iter;
    else
        Max_iter = 10000;
    end

    if isfield(opts, 'TOL')
        TOL = opts.TOL;
    else
        TOL = 1e-12;
    end

    q = min(q, m);
    normxstar2 = max(norm(xstar)^2, eps);

    Out.error = zeros(1, Max_iter + 1);
    Out.times = zeros(1, Max_iter + 1);
    Out.error(1) = norm(x - xstar)^2 / normxstar2;
    Out.times(1) = 0;

    tic;
    iter = Max_iter;

    for k = 1:Max_iter
        % Uniform row-block sampling, consistent with KPP new-block sampling.
        tau = randperm(m, q);
        AS = A(tau, :);
        rS = AS*x - b(tau);

        % Exact block Kaczmarz correction:
        %   d = A_S^dagger (A_S x - b_S).
        % For the usual case q <= rank(A_S), use the cheaper normal-equation
        % form d = A_S' (A_S A_S')^{-1} rS; otherwise fall back to lsqminnorm.
        G = AS*AS';
        G = (G + G')/2;
        [L, p] = chol(G, 'lower');
        if p == 0
            y = L' \ (L \ rS);
            d = AS' * y;
        else
            d = lsqminnorm(full(AS), rS);
        end

        x = x - d;

        Out.error(k + 1) = norm(x - xstar)^2 / normxstar2;
        Out.times(k + 1) = toc;

        if Out.error(k + 1) < TOL
            iter = k;
            Out.error = Out.error(1:k + 1);
            Out.times = Out.times(1:k + 1);
            break;
        end
    end

    Out.iter = iter;
    Out.sampling = 'uniform row subset without replacement';
    Out.block_size = q;
end
