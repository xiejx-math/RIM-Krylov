function [x,Out]=My_IS_Krylov_GS(A,b,q,opts)

% The iterative-sketching-based Krylov subspace method with Gaussian sketch (IS-Krylov-GS) for solving linear systems
%              Ax=b
%Input: the coefficent matrix A, the vector b, block size q and opts
% opts.initial: the initial vector x^0
% opts.TOL: the stopping rule
% opts.Pre_iter: the number of previous iterations
% opts.Max_iter: the number of maximum iterations
%.....
%
%Output: the approximate solution x and Out
% Out.error: the relative iterative error \|x^k-x^*\|^2/\|x^k\|^2
% Out.iter: the total number of iteration
% Out.times: the CPU time
%
% Based on the manuscript:
%
%
%
% Coded by Yonghan Sun, Beihang University, sunyonghan@buaa.edu.cn

tic

%% setting some parameters
[m,n]=size(A);

flag=exist('opts');

%%%% setting the max iteration
if (flag && isfield(opts,'Max_iter'))
    Max_iter=opts.Max_iter;
else
    Max_iter=200000;
end

%%%% setting the tolerance
if (flag && isfield(opts,'TOL'))
    TOL=opts.TOL;
else
    TOL=10^(-12);
end

%%%% setting the initial point
if (flag && isfield(opts,'initial'))
    initialx=opts.initial;
else
    initialx=zeros(n,1);
end
x=initialx;

%%%% setting the tolerance for checking the values of S_k(Ax_k-b)
if (flag && isfield(opts,'TOL1'))
    TOL1=opts.TOL1;
else
    %TOL1=10^-20;
    TOL1=eps^2;
end

%%%% determining what to use as the stopping rule
if (flag && isfield(opts,'xstar'))
    xstar=opts.xstar;
    if m>=n
        normxstar=sum(xstar.^2);
        error1=sum((xstar-x).^2)/normxstar;
        strategy=1;
    else
        strategy=0;
    end
else
    strategy=0;
end

if (flag && isfield(opts,'strategy'))
    strategy=opts.strategy;
    normxstar=sum(xstar.^2);
end


if (flag && isfield(opts,'Pre_iter'))
    l=opts.Pre_iter;
else
    l=10;
end

if ~strategy
    normb=sum(b.^2)+1;
    error1=sum((A*x-b).^2)/normb;
end

RSE(1)=error1;


%% executing the IS-Krylov-GS method
stopc=0;
iter=-1;
times(1)=toc;
current_idx = 1;
%parpool("local",4);
while ~stopc
    tic
    iter=iter+1;
    %% Randomly select a sampling matrix S_k until $S_k(A*x_k − b) \neq 0$
    stopsampling=0;
    while ~stopsampling
        %indexR=randperm(m,q);
        %AindexR=A(indexR,:);
        %bindexR=b(indexR);
        [AindexR,bindexR]=My_Gaussian_sketch(A, b, q);

        Axb=AindexR*x-bindexR;
        gamma=sum(Axb.^2);

        % If normAxb is greater than TOL1, we consider it to be equal to 0
        if gamma>TOL1
            stopsampling=1;
            dk=AindexR'*Axb;
            dk=-dk;
            %norm_dk=norm(dk)^2;
        end
    end

   %% compute the x
   if l==1
       norm_dk=sum(dk.^2);
       delta=gamma/norm_dk;
       x=x+delta*dk;
   elseif l==2
        if iter==0
           pk=dk;
           norm_pk=sum(pk.^2);
           delta=gamma/norm_pk;
           x=x+delta*pk;
           %norm_Pk(iter+1)=norm_pk;
           %Pk(:,iter+1)=pk;
        else
           %Pk=Pk(:,2:l-1);
           %norm_Pk=norm_Pk(:,2:l-1);
           eta=pk'*dk;
           eta=eta/norm_pk;
           pk=dk-eta*pk;
           norm_pk=norm(pk)^2;
           delta=gamma/norm_pk;
           x=x+delta*pk;
        end
    else
       if iter==0
           pk=dk;
           norm_pk=sum(pk.^2);
           delta=gamma/norm_pk;
           x=x+delta*pk;
           norm_Pk(iter+1)=norm_pk;
           Pk(:,iter+1)=pk;
       elseif iter<=l-2&&iter>=1
           %for ii=1:iter
           %    eta(ii)=dk'*Pk(:,ii)/norm_Pk(ii);
           %    pk=dk-eta(ii)*Pk(:,ii);
           %end
           eta=Pk'*dk;
           eta=eta./norm_Pk';
           pk=dk-Pk*eta;
           Pk(:,iter+1)=pk;
           norm_pk=sum(pk.^2);
           norm_Pk(iter+1)=norm_pk;
           delta=gamma/norm_pk;
           x=x+delta*pk;
       else
           eta=Pk'*dk;
           eta=eta./norm_Pk';
           pk=dk-Pk*eta;
           Pk(:,current_idx)=pk;
           norm_pk=sum(pk.^2);
           norm_Pk(current_idx)=norm_pk;
           delta=gamma/norm_pk;
           x=x+delta*pk;
           current_idx = mod(current_idx, l) + 1;
       end
   end
    %% stopping rule
    if strategy
        error1=sum((x-xstar).^2)/normxstar;
        RSE(iter+2)=error1;% RSE
        if error1<TOL  || iter>=Max_iter-1
            stopc=1;
        end
    else
        %%%% Note that we do not use this stopping rule during our test
        error1=sum((A*x-b).^2)/normb;
        RSE(iter+2)=error1;
        if  error1<TOL || iter>=Max_iter-1
            stopc=1;
        end
    end
    %% store the CPU time
    times(iter+2)=times(iter+1)+toc;
end
%% setting Output
Out.error=RSE;
Out.iter=iter+2;
Out.times=times;
end

