function [SA, Sb] = My_Gaussian_sketch(A, b, ell)
    % Gaussian Sketch
    % Input:
    %   A: m x n matrix
    %   b: m x 1 vector
    %   ell: target dimension (ell < m)
    % Output:
    %   SA: ell x n compressed matrix
    %   Sb: ell x 1 compressed vector

    % Get the size of matrix A
    [m, n] = size(A);  

    % Check if ell is less than m
    if ell >= m
        error('Target dimension ell must be less than m.');
    end

    % Step 1: Generate a Gaussian random matrix S
    % S is an ell x m matrix with entries drawn from N(0, 1)
    S = randn(ell, m) / sqrt(ell);  

    % Step 2: Compute the Gaussian sketch for A, and b
    % Compute the sketched matrices and vector
    SA = S * A;  
    Sb = S * b;  
end
