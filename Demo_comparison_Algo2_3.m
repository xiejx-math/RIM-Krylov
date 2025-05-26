% Comparison of RIM-AS-PS and IS-Krylov-PS

clear all;
close all;

%% setup
m=5000; % number of rows
n=800; % number of columns
q=30; % size of the block
l=10; % number of previous iterations
kappa=10; % upper bound for the condition number
Max_length=10000; % max iterations
run_time=20; % average times

%% some vectors are used to store the computed results
Times_RIM_AS_PS=zeros(run_time,Max_length+1);
Error_RIM_AS_PS=zeros(run_time,Max_length+1);

Times_IS_Krylov_PS=zeros(run_time,Max_length+1);
Error_IS_Krylov_PS=zeros(run_time,Max_length+1);

%% run and store the numerical results
for ii=1:run_time
    %% generated the matrix A
    r=n;
    [U,~]=qr(randn(m, r), 0);
    [V,~]=qr(randn(n, r), 0);
    D = diag(1+(kappa-1).*rand(r, 1));
    A=U*D*V';
    clear U V D
    %% generated the right-hand vector b
    x=rand(n,1);
    b=A*x;
    xLS=lsqminnorm(A,b);
    %% parameter setup
    opts.xstar=xLS;
    opts.TOL1=eps^2;
    opts.TOL=10^(-20);
    opts.Pre_iter=l; % number of previous iterations
    opts.Max_iter=Max_length;

    %% run RIM-AS-PS
    [x_RIM_AS_PS,Out_RIM_AS_PS]=My_RIM_AS_PS(A,b,q,opts);

    %% run IS-Krylov-PS
    [x_IS_Krylov_PS,Out_IS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);
    
    CPU_RIM_AS_PS(ii,1:length(Out_RIM_AS_PS.times))=Out_RIM_AS_PS.times;
    Error_RIM_AS_PS(ii,1:length(Out_RIM_AS_PS.error))=Out_RIM_AS_PS.error;
    
    CPU_IS_Krylov_PS(ii,1:length(Out_IS_Krylov_PS.times))=Out_IS_Krylov_PS.times;
    Error_IS_Krylov_PS(ii,1:length(Out_IS_Krylov_PS.error))=Out_IS_Krylov_PS.error;
    
    fprintf('Done, ii=%d\n',ii)
end

%% plot the iterations
xlable=1:Max_length+1;

e1=Error_RIM_AS_PS;
mine1=min(e1);
maxe1=max(e1);
e1q25=quantile(e1',0.25,2)';
e1q75=quantile(e1',0.75,2)';

e2=Error_IS_Krylov_PS;
mine2=min(e2);
maxe2=max(e2);
e2q25=quantile(e2',0.25,2)';
e2q75=quantile(e2',0.75,2)';

%%%
figure
h = fill([xlable  fliplr(xlable)], [mine2 fliplr(maxe2)],'blue','EdgeColor', 'none');
hold on
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [e2q25 fliplr(e2q75)],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [mine1 fliplr(maxe1)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [e1q25 fliplr(e1q75)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)

%%%
hold on
p1=semilogy( xlable, median(e1), 'magenta', 'LineWidth', 1,...
    'LineStyle', '--','DisplayName', 'RIM-PS-AS');
p2=semilogy( xlable, median(e2), 'blue', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'IS-Krylov-PS');
hold on
set(gca, 'YScale', 'log')

xlim([1,Max_length]);
ylim([10^(-12),1]);
ylabel('RSE','Interpreter', 'latex')
xlabel('Number of iterations','Interpreter', 'latex')
hold on
txt=title(['$\kappa=$ ',num2str(kappa), ', $r=$ ', num2str(r)]);
set(txt, 'Interpreter', 'latex');
legend([p1 p2],{'RIM-AS-PS','IS-Krylov-PS'},'Interpreter', 'latex','location', 'best')

%% plot CPU
xlable=1:(Max_length+1);

y1=CPU_RIM_AS_PS;
miny1=min(y1);
maxy1=max(y1);
y1q25=quantile(y1',0.25,2)';
y1q75=quantile(y1',0.75,2)';

y2=CPU_IS_Krylov_PS;
miny2=min(y2);
maxy2=max(y2);
y2q25=quantile(y2',0.25,2)';
y2q75=quantile(y2',0.75,2)';

figure
h = fill([median(y2)  fliplr(median(y2))], [mine2 fliplr(maxe2)],'blue','EdgeColor', 'none');
hold on
set(h,'facealpha', .05)
hold on
h = fill([median(y2)  fliplr(median(y2))], [e2q25 fliplr(e2q75)],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([median(y1)  fliplr(median(y1))], [mine1 fliplr(maxe1)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(y1)  fliplr(median(y1))], [e1q25 fliplr(e1q75)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)

hold on
p1=semilogy( median(y1), median(e1), 'magenta', 'LineWidth', 1,...
    'LineStyle', '--','DisplayName', 'RIM-AS-PS');
p2=semilogy( median(y2), median(e2), 'blue', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'IS-Krylov-PS');
hold on
set(gca, 'YScale', 'log')

xlim([0,11]);
ylim([10^(-12),1]);
ylabel('RSE','Interpreter', 'latex')
xlabel('CPU','Interpreter', 'latex')
legend([p1 p2],{'RIM-AS-PS','IS-Krylov-PS'},'Interpreter', 'latex','location', 'best')

hold on
txt=title(['$\kappa=$ ',num2str(kappa), ', $r=$ ', num2str(r)]);
set(txt, 'Interpreter', 'latex');
