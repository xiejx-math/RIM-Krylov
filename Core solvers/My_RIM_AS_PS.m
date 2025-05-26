function [x,Out]=My_RIM_AS_PS(A,b,q,opts)

% The efficient randomized iterative methods with affine subspace search and parition sampling (RIM-AS-PS) for solving linear systems
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
    TOL1=10^-20;
end

%%%% determining what to use as the stopping rule
if (flag && isfield(opts,'xstar'))
    xstar=opts.xstar;
    if m>=n
        normxstar=norm(xstar)^2;
        error1=norm(xstar-x)^2/normxstar;
        strategy=1;
    else
        strategy=0;
    end
else
    strategy=0;
end

if (flag && isfield(opts,'strategy'))
    strategy=opts.strategy;
    normxstar=norm(xstar)^2;
end


if (flag && isfield(opts,'Pre_iter'))
    l=opts.Pre_iter;
else
    l=10;
end

if ~strategy
    normb=norm(b)^2+1;
    error1=norm(A*x-b)^2/normb;
end

RSE(1)=error1;

%% a uniform random permutation for both A and b
if (flag && isfield(opts,'permS'))
    S=opts.permS;
    A=A(S,:);
    b=b(S);
else
    S=randperm(m);
    A=A(S,:);
    b=b(S);
end

%% setting the probability
if (flag && isfield(opts,'probset'))
    probset=opts.probset;
else
    probset=0;
end

if probset
    Aarrs=opts.Aarrs;
    barrs=opts.barrs;
    cumsumpro=opts.cumsumpro;
else

    normAfro=norm(A,'fro')^2;
    tau=floor(m/q);
    blockAnormfro=zeros(tau,1);
    %prob=zeros(tau,1);
    for i=1:tau
        if i==tau
            ps=((i-1)*q+1):1:m;
        else
            ps=((i-1)*q+1):1:(i*q);
        end
        Aps=A(ps,:);
        blockAnormfro(i)=norm(A(ps,:),'fro')^2;
        Aarrs{i}=Aps;
        barrs{i}=b(ps);
    end
    prob=blockAnormfro/normAfro;
    cumsumpro=cumsum(prob);
end

%% executing the RIM-AS-PS method
stopc=0;
iter=-1;
times(1)=toc;
V0(:,1)=x;
%parpool("local",4);
while ~stopc
    tic
    iter=iter+1;
    %% Randomly select a sampling matrix S_k until $S_k(A*x_k − b) \neq 0$
    stopsampling=0;
    while ~stopsampling
        w=sum(cumsumpro<rand)+1;
        AindexR=Aarrs{w};
        bindexR=barrs{w};

        Axb=AindexR*x-bindexR;
        gamma=norm(Axb)^2;

        % If normAxb is greater than TOL1, we consider it to be equal to 0
        if gamma>TOL1
            stopsampling=1;
            dk=AindexR'*Axb;
            dk=-dk;
            norm_dk=norm(dk)^2;
        end
    end

   %% compute the x
    if l==1
        alpha=gamma/norm_dk;
        %%%% update x for k=1
        x=x+alpha*dk;
    elseif l==2
        if iter==0
            alpha=gamma/norm_dk;
            %xold=x;
            %%%% update x for k=1
            x=x+alpha*dk;
            xold=x;
            xoold=initialx;
            under_old=alpha;
            normAxb_old=gamma;
        else
            wk=(xoold-xold)'*dk;
            hk=wk/(normAxb_old*under_old);
            wh=wk'*hk;
            under_s=gamma/(norm_dk-wh);
            over_s=-under_s*hk;
            x=x+over_s*(xoold-xold)+under_s*dk;
            xoold=xold;
            xold=x;
            under_old=under_s;
            normAxb_old=gamma;
        end
    else
        if iter==0
            V0(:,iter+1)=x;
            under_s=gamma/norm_dk;
            %xold=x;
            %%%% update x for k=1
            x=x+under_s*dk;
            % Update the matrix C
            C0(iter+1)=1/(gamma*under_s);
            C(iter+1,iter+1)=C0(iter+1);
        elseif iter<=l-2&&iter>=1
            V=V0-x;
            %for i=1:iter
            %    V(:,i)=V0(:,i)-x;
            %end
            V0(:,iter+1)=x;
            wk=V'*dk;
            hk=C*wk;
            wh=wk'*hk;
            under_s=gamma/(norm_dk-wh);
            over_s=-under_s*hk;
            x=x+V*over_s+under_s*dk;
            % Update the matrix C
            C0(iter+1)=1/(gamma*under_s);
            C(iter+1,iter+1)=C0(iter)+C0(iter+1);
            C(iter+1,iter)=-C0(iter);
            C(iter,iter+1)=-C0(iter);
            
        elseif iter==l-1
            V=V0-x;
            %for i=1:iter
            %    V(:,i)=V0(:,i)-x;
            %end
            V0(:,iter+1)=x;
            wk=V'*dk;
            hk=C*wk;
            wh=wk'*hk;
            under_s=gamma/(norm_dk-wh);
            over_s=-under_s*hk;
            x=x+V*over_s+under_s*dk;
            
            % Update the matrix C
            %C0(1:l-2)=C0(2:l-1);
            C0(l)=1/(gamma*under_s);
            C(1:l-2,1:l-2)=C(2:l-1,2:l-1);
            C(1,1)=-C(1,2);
            C(l-1,l-2)=-C0(l-1);
            C(l-2,l-1)=C(l-1,l-2);
            C(l-1,l-1)=-C(l-1,l-2)+C0(l);
        else
            V0(:,l+1)=x;
            V0(:,1:l)=V0(:,2:l+1);
            V=V0-x;
            V=V(:,1:l-1);
            %for i=1:l-1
            %    V(:,i)=V0(:,i)-x;
            %end
            wk=V'*dk;
            hk=C*wk;
            wh=wk'*hk;
            under_s=gamma/(norm_dk-wh);
            over_s=-under_s*hk;
            x=x+V*over_s+under_s*dk;
            %x=x+initialx+V*over_s-under_s*dk;
            % Update the matrix C
            C0(1:l-1)=C0(2:l);
            C0(l)=1/(gamma*under_s);
            C(1:l-2,1:l-2)=C(2:l-1,2:l-1);
            C(1,1)=-C(1,2);
            C(l-1,l-2)=-C0(l-1);
            C(l-2,l-1)=C(l-1,l-2);
            C(l-1,l-1)=-C(l-1,l-2)+C0(l);
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
        %%%% Note that we do not us this stopping rule during our test
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
Out.iter=iter+1;
Out.times=times;
end

