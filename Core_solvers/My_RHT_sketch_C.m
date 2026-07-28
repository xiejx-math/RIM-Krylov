function [Aout, bout, info] = My_RHT_sketch_C(A, b)
%MY_RHT_SKETCH_C Apply a normalized randomized Hadamard transform to [A,b].
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
%
% The function applies a randomized Hadamard transform to [A,b].
% If the number of rows is not a power of two, zero rows are appended
% before applying the transform.

    [m, n] = size(A);
    b = b(:);

    if length(b) ~= m
        error('The length of b must be equal to the number of rows of A.');
    end

    if exist('srht', 'file') ~= 3 && exist('srht', 'file') ~= 2
        error('srht.mexw64 or srht.m must be on the MATLAB path.');
    end

    M = 2^nextpow2(m);

    % The RHT output is dense, so full storage is unavoidable here.
    X = zeros(M, n + 1);
    X(1:m, 1:n) = full(A);
    X(1:m, n + 1) = full(b);

    % MATLAB controls the Rademacher signs.
    % If srht applies the unnormalized Hadamard transform, then
    % D = d/sqrt(M) gives the normalized RHT.
    d = 2*randi([0, 1], M, 1) - 1;
    D = d / sqrt(M);

    % Full RHT. Do not sample rows here.
    % The third input is kept as M for compatibility with srht.mexw64.
    [HX, ~, ~] = srht(X, D, M);

    Aout = HX(:, 1:n);
    bout = HX(:, n + 1);

    info = struct('used', true, ...
                  'm_original', m, ...
                  'm_work', M, ...
                  'num_padded_rows', M - m, ...
                  'implementation', 'mex');
end