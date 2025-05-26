% This Matlab file is used to compare pinv, lsqminnorm, and IS-Krylov-PS

clear
close all

n=1000;% the column of the coefficient matrix
r=n; % the rank of the coefficient matrix
kappa=20; % desired condition number
mi=10000; % the initial setting for the number of rows
testsize=1:1:10; % used for updating the number of rows
run_time=20;  % average times
q=50;% size of the block % here we set q=50,  we can set it as other values,
%opts.TOL=5*10^(-28);% change the tolence if the size of problem change
opts.TOL=5*10^(-29);% change the tolence if the size of problem change

%% some vectors are used to store the desired numerical results
lstime=zeros(length(testsize),1);
pinvtime=zeros(length(testsize),1);
IS_Krylov_PStime=zeros(length(testsize),1);

CPU_LS=zeros(run_time,length(testsize));
CPU_pinv=zeros(run_time,length(testsize));
CPU_IS_Krylov_PS=zeros(run_time,length(testsize));

for ii=1:length(testsize)
    m=mi*testsize(ii); % the row of the coefficient matrix
    
    for jj=1:run_time 
         %% generated the matrix A
        [U,~]=qr(randn(m, r), 0);
        [V,~]=qr(randn(n, r), 0);
        D = diag(1+(kappa-1).*rand(r, 1));
        A=U*D*V';
        clear U V D
        %% generated the right-hand vector b
        
        x=randn(n,1);
        b=A*x;

        %% Matlab function solvers
        tic
        xLS=lsqminnorm(A,b);
        MTls_CPU=toc;

        tic
        xpinv=pinv(A)*b;
        MTpinv_CPU=toc;

        %% parameter setup
        opts.xstar=x;
        opts.TOL1=eps^2;
        %opts.TOL1=10^(-30);
        opts.Max_iter=500000;
        opts.Pre_iter=10;
        
         %% IS-Krylov-PS
        [xIS_Krylov_PS,OutIS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);
        %%
        CPU_IS_Krylov_PS(jj,ii)=OutIS_Krylov_PS.times(end);
        CPU_LS(jj,ii)=MTls_CPU;
        CPU_pinv(jj,ii)=MTpinv_CPU;

        %%
        fprintf('lsqminnorm: %8e, pinv: %8e,IS-Krylov-PS: %8e\n',norm(x-xLS),norm(x-xpinv),norm(x-xIS_Krylov_PS))
        
    end
    fprintf('Done %d\n',ii)
end

%% plot errors
lightgray =   [0.8 0.8 0.8];
mediumgray =  [0.6 0.6 0.6];
lightred =    [1 0.9 0.9];
mediumred =   [1 0.6 0.6];
lightgreen =  [0.9 1 0.9];
mediumgreen = [0.6 1 0.6];
lightblue =   [0.9 0.9 1];
mediumblue =  [0.6 0.6 1];
lightmagenta =   [1 0.9 1];
mediummagenta =  [1 0.6 1];

%%
display_names = {'pinv','lsqminnorm','IS-Krylov-PS'};
arrsIter = {CPU_pinv',CPU_LS',CPU_IS_Krylov_PS'};

num_methods = length(arrsIter);

%line_colors = {'blue','green','magenta'};
%minmax_colors = { lightblue, lightgreen,lightmagenta};
%quant_colors = { mediumblue,mediumgreen,mediummagenta};

line_colors = {'black','green','red','blue'};
minmax_colors = { lightgray, lightgreen,lightred,lightblue};
quant_colors = { mediumgray,mediumgreen,mediumred,mediumblue};

%line_colors = {'black', 'green', 'red', 'blue'};
%minmax_colors = {lightgray, lightgreen, lightred, lightblue};
%quant_colors = {mediumgray, mediumgreen, mediumred, mediumblue};

display_legend = true;
max_val_in_plot = 1e5;

%% plot CPU
%%
num_iter_array=mi*testsize;
%num_iter_array=num_iter_array';

%%%
y1=CPU_IS_Krylov_PS';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2)';
y1q75=quantile(y1,0.75,2)';

y2=CPU_pinv';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2)';
y2q75=quantile(y2,0.75,2)';

y3=CPU_LS';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2)';
y3q75=quantile(y3,0.75,2)';

figure
h = fill([num_iter_array  fliplr(num_iter_array)], [miny1 fliplr(maxy1)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([num_iter_array  fliplr(num_iter_array)], [y1q25 fliplr(y1q75)],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([num_iter_array  fliplr(num_iter_array)], [miny2 fliplr(maxy2)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([num_iter_array  fliplr(num_iter_array)], [y2q25 fliplr(y2q75)],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([num_iter_array  fliplr(num_iter_array)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([num_iter_array  fliplr(num_iter_array)], [y3q25 fliplr(y3q75)],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on

hold on
p1=semilogy( num_iter_array, median(y1'), 'black', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','o','DisplayName', 'IS-Krylov-PS');
p2=semilogy( num_iter_array, median(y2'), 'red', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','s','DisplayName', 'pinv');
p3=semilogy( num_iter_array, median(y3'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','^','DisplayName', 'lsqminnorm');

xlim([num_iter_array(1),num_iter_array(end)]);
ylim([0,0.01]);
ylabel('CPU')
xlabel('Number of rows (m)')
hold on
txt=title( ['$n=$ ',num2str(n),', $\kappa=$ ',num2str(kappa),', $\ell=$ ',num2str(opts.Pre_iter)]);
set(txt, 'Interpreter', 'latex');
legend([p2 p3 p1],{'pinv','lsqminnorm','IS-Krylov-PS'},'Interpreter', 'latex','location', 'best')
