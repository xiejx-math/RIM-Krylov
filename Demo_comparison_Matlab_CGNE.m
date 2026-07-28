%% Figure 5: Comparison with MATLAB solvers and CGNE
% This demo reproduces Figure 5 in the paper.
%
% It compares pinv, lsqminnorm, CGNE, and IS-Krylov-PS with their
% Gaussian-sketched variants on synthetic consistent linear
% systems with an increasing number of rows.
%
% The figure reports the CPU time versus the number of rows.
% For the Gaussian-sketched variants, the preprocessing time required
% to form the sketched system and the solution time on the sketched
% system are reported separately.
%
% Three experimental settings corresponding to the panels in Figure 5
% can be selected below.

clear
close all
clc

n=500;       % the number of columns of the coefficient matrix
mi=10000;    % the initial setting for the number of rows
kappa=10;    % desired condition number
opts.TOL=5.2*10^(-29);  % change the tolerance if the size of problem changes
Pre_iter=10;
% 
% n=1000;       % the number of columns of the coefficient matrix
% mi=20000;    % the initial setting for the number of rows
% kappa=10;    % desired condition number
% opts.TOL=8*10^(-29);  % change the tolerance if the size of problem changes
% Pre_iter=10;

% n=1000;       % the number of columns of the coefficient matrix
% mi=10000;    % the initial setting for the number of rows
% kappa=30;    % desired condition number
% opts.TOL=9.5*10^(-28);  % change the tolerance if the size of problem changes
% Pre_iter=50;

r=n;         % the rank of the coefficient matrix
testsize=1:1:10; % used for updating the number of rows
run_time=20;  % number of repeated runs; change it to 20 for the paper figure
q=50;        % size of the block

%opts.TOL=5*10^(-28); % change the tolerance if the size of problem changes


% One-shot Gaussian sketching parameter. We take s = sketch_ratio*n rows.
sketch_ratio=10;

%% Arrays for storing CPU time
CPU_LS=zeros(run_time,length(testsize));
CPU_pinv=zeros(run_time,length(testsize));
CPU_CGNE=zeros(run_time,length(testsize));
CPU_IS_Krylov_PS=zeros(run_time,length(testsize));

CPU_sketched_LS=zeros(run_time,length(testsize));
CPU_sketched_pinv=zeros(run_time,length(testsize));
CPU_sketched_CGNE=zeros(run_time,length(testsize));
CPU_sketched_IS_Krylov_PS=zeros(run_time,length(testsize));

CPU_sketch=zeros(run_time,length(testsize));

%% Arrays for storing iteration numbers and success/failure flags
ITER_CGNE=zeros(run_time,length(testsize));
ITER_IS_Krylov_PS=zeros(run_time,length(testsize));
ITER_sketched_CGNE=zeros(run_time,length(testsize));
ITER_sketched_IS_Krylov_PS=zeros(run_time,length(testsize));

HITMAX_CGNE=false(run_time,length(testsize));
HITMAX_IS_Krylov_PS=false(run_time,length(testsize));
HITMAX_sketched_CGNE=false(run_time,length(testsize));
HITMAX_sketched_IS_Krylov_PS=false(run_time,length(testsize));

SUCCESS_CGNE=false(run_time,length(testsize));
SUCCESS_IS_Krylov_PS=false(run_time,length(testsize));
SUCCESS_sketched_CGNE=false(run_time,length(testsize));
SUCCESS_sketched_IS_Krylov_PS=false(run_time,length(testsize));
%% Main loop
for ii=1:length(testsize)
    m=mi*testsize(ii); % the number of rows of the coefficient matrix
    s=min(m,ceil(sketch_ratio*n)); % the number of rows after sketching
    
    for jj=1:run_time
        %% Generate the matrix A
        [U,~]=qr(randn(m,r),0);
        [V,~]=qr(randn(n,r),0);
        D=diag(1+(kappa-1).*rand(r,1));
        A=U*D*V';
        clear U V D
        
        %% Generate the right-hand vector b
        x=randn(n,1);
        b=A*x;
        
        %% Parameter setup for iterative solvers
        opts.xstar=x;
        opts.TOL1=eps^2;
        %opts.TOL1=10^(-30);
        opts.Max_iter=100000;
        opts.Pre_iter=Pre_iter;
        opts.initial=zeros(n,1);
        
        %% Original Matlab function solvers on (A,b)
        tic
        xLS=lsqminnorm(A,b);
        MTls_CPU=toc;
        
        tic
        xpinv=pinv(A)*b;
        MTpinv_CPU=toc;
        
        %% Original CGNE on (A,b)
        [xCGNE,OutCGNE]=My_CGNE(A,b,opts);
        
        %% Original IS-Krylov-PS on (A,b)
        [xIS_Krylov_PS,OutIS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);
        
        %% One-shot Gaussian sketch
        % We form As = S*A and bs = S*b, where S is a Gaussian sketching
        % matrix.  Since this is a synthetic consistent problem with
        % b = A*x, we then set bs = As*x to enforce floating-point
        % consistency of the sketched system.  This is algebraically
        % equivalent to bs = S*b in exact arithmetic.
        tic
        % [As,bs]=My_Gaussian_sketch(A,b,s);
        S = randn(s, m) / sqrt(s);  
        As=S*A;
        bs=As*x;
        Sketch_CPU=toc;
        
        %% Sketched Matlab function solvers on (As,bs)
        tic
        x_sketched_LS=lsqminnorm(As,bs);
        % MT_sketched_ls_CPU=Sketch_CPU+toc;
        MT_sketched_ls_CPU=toc;
        
        tic
        x_sketched_pinv=pinv(As)*bs;
        % MT_sketched_pinv_CPU=Sketch_CPU+toc;
        MT_sketched_pinv_CPU=toc;
        
        %% Sketched CGNE on (As,bs)
        [x_sketched_CGNE,Out_sketched_CGNE]=My_CGNE(As,bs,opts);
        
        %% Sketched IS-Krylov-PS on (As,bs)
        [x_sketched_IS_Krylov_PS,Out_sketched_IS_Krylov_PS]=My_IS_Krylov_PS(As,bs,q,opts);
        
        %% Store CPU time
        CPU_LS(jj,ii)=MTls_CPU;
        CPU_pinv(jj,ii)=MTpinv_CPU;
        CPU_CGNE(jj,ii)=OutCGNE.times(end);
        CPU_IS_Krylov_PS(jj,ii)=OutIS_Krylov_PS.times(end);
        
        CPU_sketch(jj,ii)=Sketch_CPU;
        CPU_sketched_LS(jj,ii)=MT_sketched_ls_CPU;
        CPU_sketched_pinv(jj,ii)=MT_sketched_pinv_CPU;
        % CPU_sketched_CGNE(jj,ii)=Sketch_CPU+Out_sketched_CGNE.times(end);
        % CPU_sketched_IS_Krylov_PS(jj,ii)=Sketch_CPU+Out_sketched_IS_Krylov_PS.times(end);
        % 
        CPU_sketched_CGNE(jj,ii)=Out_sketched_CGNE.times(end);
        CPU_sketched_IS_Krylov_PS(jj,ii)=Out_sketched_IS_Krylov_PS.times(end);
        
        %% Store iteration numbers and success/failure flags
        iter_CGNE=get_iter_from_out(OutCGNE);
        iter_IS_Krylov_PS=get_iter_from_out(OutIS_Krylov_PS);
        iter_sketched_CGNE=get_iter_from_out(Out_sketched_CGNE);
        iter_sketched_IS_Krylov_PS=get_iter_from_out(Out_sketched_IS_Krylov_PS);
        
        ITER_CGNE(jj,ii)=iter_CGNE;
        ITER_IS_Krylov_PS(jj,ii)=iter_IS_Krylov_PS;
        ITER_sketched_CGNE(jj,ii)=iter_sketched_CGNE;
        ITER_sketched_IS_Krylov_PS(jj,ii)=iter_sketched_IS_Krylov_PS;
        
        HITMAX_CGNE(jj,ii)=(iter_CGNE>=opts.Max_iter);
        HITMAX_IS_Krylov_PS(jj,ii)=(iter_IS_Krylov_PS>=opts.Max_iter);
        HITMAX_sketched_CGNE(jj,ii)=(iter_sketched_CGNE>=opts.Max_iter);
        HITMAX_sketched_IS_Krylov_PS(jj,ii)=(iter_sketched_IS_Krylov_PS>=opts.Max_iter);
        
        SUCCESS_CGNE(jj,ii)=~HITMAX_CGNE(jj,ii);
        SUCCESS_IS_Krylov_PS(jj,ii)=~HITMAX_IS_Krylov_PS(jj,ii);
        SUCCESS_sketched_CGNE(jj,ii)=~HITMAX_sketched_CGNE(jj,ii);
        SUCCESS_sketched_IS_Krylov_PS(jj,ii)=~HITMAX_sketched_IS_Krylov_PS(jj,ii);
        %% Print errors
        fprintf(['m=%d, s=%d, run=%d: ', ...
            'pinv: %8e, lsqminnorm: %8e, IS-Krylov-PS: %8e, CGNE: %8e, ', ...
            'sketched pinv: %8e, sketched lsqminnorm: %8e, ', ...
            'sketched IS-Krylov-PS: %8e, sketched CGNE: %8e, sketch time: %8e\n'], ...
            m,s,jj, ...
            norm(x-xpinv),norm(x-xLS),norm(x-xIS_Krylov_PS),norm(x-xCGNE), ...
            norm(x-x_sketched_pinv),norm(x-x_sketched_LS), ...
            norm(x-x_sketched_IS_Krylov_PS),norm(x-x_sketched_CGNE),Sketch_CPU);
        fprintf(['Iteration status: ', ...
            'CGNE iter=%d, success=%d; ', ...
            'IS-Krylov-PS iter=%d, success=%d; ', ...
            'sketched CGNE iter=%d, success=%d; ', ...
            'sketched IS-Krylov-PS iter=%d, success=%d\n'], ...
            iter_CGNE,SUCCESS_CGNE(jj,ii), ...
            iter_IS_Krylov_PS,SUCCESS_IS_Krylov_PS(jj,ii), ...
            iter_sketched_CGNE,SUCCESS_sketched_CGNE(jj,ii), ...
            iter_sketched_IS_Krylov_PS,SUCCESS_sketched_IS_Krylov_PS(jj,ii));
        clear A b As bs
    end
    fprintf('Done m=%d, s=%d\n',m,s)
end

%% Plot CPU time with min/max and quartile bands for all nine curves
num_iter_array=mi*testsize;

% Each matrix is stored as: rows = different m values, columns = repeated runs.
Y_pinv=CPU_pinv';
Y_LS=CPU_LS';
Y_IS=CPU_IS_Krylov_PS';
Y_CGNE=CPU_CGNE';

Y_sketched_pinv=CPU_sketched_pinv';
Y_sketched_LS=CPU_sketched_LS';
Y_sketched_IS=CPU_sketched_IS_Krylov_PS';
Y_sketched_CGNE=CPU_sketched_CGNE';
Y_sketch=CPU_sketch';

% Plot the original and sketched variants of the same algorithm as a pair:
% original = solid line, sketched = dashed line.  The color and marker are
% kept unchanged within each pair.  The sketching time is plotted separately.
arrsIter={ ...
    Y_pinv, Y_sketched_pinv, ...
    Y_LS, Y_sketched_LS, ...
    Y_IS, Y_sketched_IS, ...
    Y_CGNE, Y_sketched_CGNE, ...
    Y_sketch};

display_names={ ...
    'pinv', ...
    'sketched pinv', ...
    'lsqminnorm', ...
    'sketched lsqminnorm', ...
    'IS-Krylov-PS', ...
    'sketched IS-Krylov-PS', ...
    'CGNE', ...
    'sketched CGNE', ...
    'Gaussian sketch'};
% display_names={ ...
%     'pinv', ...
%     'Sk-pinv', ...
%     'lsqminnorm', ...
%     'Sk-lsqminnorm', ...
%     'IS-Krylov-PS', ...
%     'Sk-IS-Krylov-PS', ...
%     'CGNE', ...
%     'Sk-CGNE', ...
%     'Gaussian sketch'};
% display_names={ ...
%     'pinv', ...
%     'G-pinv', ...
%     'lsqminnorm', ...
%     'G-lsqminnorm', ...
%     'IS-Krylov-PS', ...
%     'G-IS-Krylov-PS', ...
%     'CGNE', ...
%     'G-CGNE', ...
%     'Gaussian sketch'};
% Use MATLAB built-in color specifications:
% pinv: black, lsqminnorm: green, IS-Krylov-PS: red, CGNE: blue,
% Gaussian sketching time: magenta.
line_colors={ ...
    'k', 'k', ...
    'g', 'g', ...
    'r', 'r', ...
    'b', 'b', ...
    'm'};

line_styles={ ...
    '-', '--', ...
    '-', '--', ...
    '-', '--', ...
    '-', '--', ...
    ':'};

markers={ ...
    's', 's', ...
    '^', '^', ...
    'o', 'o', ...
    'd', 'd', ...
    'x'};

num_methods=length(arrsIter);

figure('Position',[100,100,560,420])
hold on
box on
set(gcf,'Color','w')
set(gca,'FontSize',14, ...
    'LineWidth',1.2, ...
    'TickLabelInterpreter','latex')

%% Plot min/max bands and 25%-75% quantile bands
for method_id=1:num_methods

    y=arrsIter{method_id};

    miny=min(y,[],2)';
    maxy=max(y,[],2)';
    yq25=quantile(y,0.25,2)';
    yq75=quantile(y,0.75,2)';

    % Min-max shaded region.
    h=fill([num_iter_array fliplr(num_iter_array)], ...
        [miny fliplr(maxy)], ...
        line_colors{method_id}, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');
    set(h,'facealpha',0.05)

    % 25%-75% quantile shaded region.
    h=fill([num_iter_array fliplr(num_iter_array)], ...
        [yq25 fliplr(yq75)], ...
        line_colors{method_id}, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');
    set(h,'facealpha',0.10)

end

%% Plot median CPU curves
legend_handles=gobjects(num_methods,1);

for method_id=1:num_methods

    y=arrsIter{method_id};

    legend_handles(method_id)=semilogy(num_iter_array,median(y,2)', ...
        'Color',line_colors{method_id}, ...
        'LineWidth',1.5, ...
        'LineStyle',line_styles{method_id}, ...
        'Marker',markers{method_id}, ...
        'MarkerSize',6, ...
        'DisplayName',display_names{method_id});

end

xlim([num_iter_array(1),num_iter_array(end)]);
% ylim([1e-3, 4])
xlabel('Number of rows $(m)$','Interpreter','latex','FontSize',16)
ylabel('CPU time','Interpreter','latex','FontSize',16)

txt=title(['$n=',num2str(n), ...
    ',\ m_{\rm sk}=',num2str(ceil(sketch_ratio*n)), ...
    ',\ \kappa=',num2str(kappa), ...
    ',\ \ell=',num2str(opts.Pre_iter),'$'], ...
    'Interpreter','latex', ...
    'FontSize',15);
set(txt,'Interpreter','latex');

legend(legend_handles,display_names, ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',9)

% Use a positive lower bound for the logarithmic y-axis.
all_CPU=[CPU_pinv(:);CPU_sketched_pinv(:); ...
    CPU_LS(:);CPU_sketched_LS(:); ...
    CPU_IS_Krylov_PS(:);CPU_sketched_IS_Krylov_PS(:); ...
    CPU_CGNE(:);CPU_sketched_CGNE(:); ...
    CPU_sketch(:)];
positive_CPU=all_CPU(isfinite(all_CPU) & all_CPU>0);

if ~isempty(positive_CPU)
    ylim([0.01*min(positive_CPU),1*max(positive_CPU)]);
end

%% Save numerical data and figure
save('Demo_comparison_Matlab_PS_sketched_CGNE_with_sketch_time_results.mat', ...
    'n','r','kappa','mi','testsize','run_time','q','sketch_ratio', ...
    'CPU_pinv','CPU_LS','CPU_CGNE','CPU_IS_Krylov_PS', ...
    'CPU_sketched_pinv','CPU_sketched_LS','CPU_sketched_CGNE', ...
    'CPU_sketched_IS_Krylov_PS','CPU_sketch');

savefig('Demo_comparison_Matlab_PS_sketched_CGNE_with_sketch_time_CPU.fig');
exportgraphics(gcf,'Demo_comparison_Matlab_PS_sketched_CGNE_with_sketch_time_CPU.pdf','ContentType','vector');
function iter=get_iter_from_out(Out)
    if isfield(Out,'iter')
        iter=Out.iter;
    elseif isfield(Out,'iters')
        iter=Out.iters;
    elseif isfield(Out,'iteration')
        iter=Out.iteration;
    elseif isfield(Out,'times')
        iter=length(Out.times)-1;
    elseif isfield(Out,'error')
        iter=length(Out.error)-1;
    else
        error('Cannot find the iteration number from the output structure.');
    end
end