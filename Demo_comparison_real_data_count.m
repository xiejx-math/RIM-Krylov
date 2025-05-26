% This Matlab file is used to compare RABK, SCGP, and IS-Krylov-PS
% based on the data from SuiteSparse Matrix Collection

clear
close all

%% generated the matrix A using the data from SuiteSparse Matrix Collection
%load abtaha2
%load model1
%load crew1
%load WorldCities
%load well1033
%load cr42
%load Franz1
load GL7d11
%load D_6
%load rel6
%load lp_ship04s

A=Problem.A;
Max_length=1500000; % max iterations
q=30; % size of the block
run_time=20; % average times


%% executing "run_time" times of the algorithms
for ii=1:run_time
    
    %% generated the right-hand vector b
    [m,n]=size(A);
    x=randn(n,1);
    b=A*x;
    xLS=lsqminnorm(A,b);

    %% parameter setup
    opts.xstar=xLS;
    opts.TOL=10^(-12);
    opts.TOL1=eps^2;
    opts.Pre_iter=50;  % number of previous iterations
    opts.Max_iter=Max_length;

    %% run RABK
    [xRABK,OutRABK]=My_RABK(A,b,q,opts);

    %% run SCGP
    [xSCGP,OutSCGP]=My_AmRABK(A,b,q,opts);

    %% run IS-Krylov-PS
    [xIS_Krylov_PS,OutIS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);

    %% store the numerical results
    Iter_RABK(ii)=OutRABK.iter;
    Iter_SCGP(ii)=OutSCGP.iter;
    Iter_IS_Krylov_PS(ii)=OutIS_Krylov_PS.iter;

    CPU_RABK(ii)=OutRABK.times(end);
    CPU_SCGP(ii)=OutSCGP.times(end);
    CPU_IS_Krylov_PS(ii)=OutIS_Krylov_PS.times(end);

    fprintf('Done, ii=%d\n',ii)
    fprintf('Number of iterations: %d,%d,%d\n',OutRABK.iter,OutSCGP.iter,OutIS_Krylov_PS.iter)
end

%% store the result
tab=[mean(Iter_RABK),mean(CPU_RABK),mean(Iter_SCGP),mean(CPU_SCGP),mean(Iter_IS_Krylov_PS),mean(CPU_IS_Krylov_PS)];