%% Local function: matrix-free CGNE
function [x,Out]=My_CGNE(A,b,opts)

% Matrix-free conjugate-gradient method for the normal equations.
% This version implements the CGNE recurrence without explicitly forming
% A*A'.  It updates x by directions of the form A'*p.

tic

[~,n]=size(A);

if nargin < 3 || isempty(opts)
    opts=struct();
end

if isfield(opts,'initial')
    x=opts.initial;
else
    x=zeros(n,1);
end

if isfield(opts,'Max_iter')
    Max_iter=opts.Max_iter;
else
    Max_iter=200000;
end

if isfield(opts,'TOL')
    TOL=opts.TOL;
else
    TOL=1e-12;
end

if isfield(opts,'xstar')
    xstar=opts.xstar;
    normxstar=norm(xstar)^2;
    if normxstar==0
        normxstar=1;
    end
    error1=norm(x-xstar)^2/normxstar;
    strategy=1;
else
    normb=norm(b)^2+1;
    error1=norm(A*x-b)^2/normb;
    strategy=0;
end

RSE(1)=error1;
times(1)=toc;

% CGNE starts from the residual r=b-A*x and applies CG to the
% normal equations A*A'*y = r, with x updated by A'*y.
r=b-A*x;
p=r;
rho=r'*r;

iter=0;
stopc=0;

% if error1 < TOL || rho <= eps
%     stopc=1;
% end
if error1 < TOL 
    stopc=1;
end

while ~stopc

    tic

    z=A'*p;
    denom=z'*z;

    % if denom <= eps
    %     stopc=1;
    % else
        alpha=rho/denom;
        x=x+alpha*z;
        r=r-alpha*(A*z);
        rho_new=r'*r;

        % if rho <= eps
        %     beta=0;
        % else
        %     beta=rho_new/rho;
        % end
        
        beta=rho_new/rho;
        p=r+beta*p;
        rho=rho_new;
    % end

    iter=iter+1;

    if strategy
        error1=norm(x-xstar)^2/normxstar;
    else
        error1=norm(A*x-b)^2/normb;
    end

    RSE(iter+1)=error1;
    times(iter+1)=times(iter)+toc;

    % if error1 < TOL || iter >= Max_iter || rho <= eps
    %     stopc=1;
    % end
    if error1 < TOL || iter >= Max_iter 
        stopc=1;
    end
end

Out.error=RSE;
Out.iter=iter+1;
Out.times=times;

end