function [SA, Sb] = My_SRHT_sketch_C(X,m,n, ell,srht_scale)
    %MY_SRHT_SKETCH_C Apply an SRHT sketch using srht.mexw64.
    %
    % This MATLAB wrapper uses the SRHT MEX routine provided in the
    % arLMM repository:
    %
    %   https://github.com/ztanml/arLMM
    %
    % The core SRHT source files `tran_srht.c` and `tran_srht.h` in the
    % original repository include an MIT License notice:
    %
    %   Copyright (C) 2017 Zilong Tan
    %
    % This MATLAB wrapper was written or modified for the numerical
    % experiments in the accompanying paper.

    % [m, n] = size(A);

    % if issparse(A) || issparse(b)
    %     error('This C-MEX SRHT implementation expects a full dense matrix. Do not pass sparse A directly.');
    % end

    % 
    % if ~isa(X, 'double')
    %     X = double(X);
    % end

    % MATLAB controls the Rademacher signs.
    % Use D = diag(d)/sqrt(m), so srht applies the normalized Hadamard transform.
    d = 2 * randi([0, 1], m, 1) - 1;
    D = d / sqrt(m);

    % Use the MEX routine only for the full Hadamard transform.
    % The third argument is only used internally for the unused random index.
    [HX, ~, ~] = srht(X, D, ell);

    % MATLAB controls the sampled rows.
    idx = randperm(m, ell);

    % Standard SRHT scaling:
    % sqrt(m/ell) * R * (H_un/sqrt(m)) * D_sign * X

    SA = srht_scale * HX(idx, 1:n);
    Sb = srht_scale * HX(idx, n+1);
end