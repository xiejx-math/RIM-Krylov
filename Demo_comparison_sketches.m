% This Matlab file is used to compare IS_Krylov_PS, IS_Krylov_CS, IS_Krylov_GS, and IS_Krylov_SRHT

clear
close all

%% setup
m=256; % number of rows
n=128; % number of columns
r=128; % number of rank
kappa=100; % upper bound for the condition number
q=30; % size of the block
run_time=20;  % average times
Max_length=3000; % max iterations

%% some vectors are used to store the computed results
IS_Krylov_PS_CPU=zeros(run_time,Max_length+1);
IS_Krylov_PS_error=zeros(run_time,Max_length+1);

IS_Krylov_CS_CPU=zeros(run_time,Max_length+1);
IS_Krylov_CS_error=zeros(run_time,Max_length+1);

IS_Krylov_GS_CPU=zeros(run_time,Max_length+1);
IS_Krylov_GS_error=zeros(run_time,Max_length+1);

IS_Krylov_SRHT_CPU=zeros(run_time,Max_length+1);
IS_Krylov_SRHT_error=zeros(run_time,Max_length+1);

%% run and store the numerical results
for ii=1:run_time
    %% generated the matrix A
    [U,~]=qr(randn(m, r), 0);
    [V,~]=qr(randn(n, r), 0);
    D = diag(1+(kappa-1).*rand(r, 1));
    A=U*D*V';
    clear U V D
    %% generated the right-hand vector b
    x=randn(n,1);
    b=A*x;
    xLS=lsqminnorm(A,b);
    %% parameter setup
    opts.xstar=xLS;
    opts.Max_iter=Max_length;
    opts.TOL1=eps^2;
    opts.TOL=10^(-30);
    opts.Pre_iter=10;

    %% IS-Krylov method with parition sampling
    [xIS_Krylov_PS,OutIS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);

    %% IS-Krylov method with Countsketch
    [xIS_Krylov_CS,OutIS_Krylov_CS]=My_IS_Krylov_CS(A,b,q,opts);

    %% IS-Krylov method with Guassian sketch
    [xIS_Krylov_GS,OutIS_Krylov_GS]=My_IS_Krylov_GS(A,b,q,opts);

    %% IS-Krylov method with SRHT
    [xIS_Krylov_SRHT,OutIS_Krylov_SRHT]=My_IS_Krylov_SRHT(A,b,q,opts);

    IS_Krylov_PS_error(ii,1:length(OutIS_Krylov_PS.error))=OutIS_Krylov_PS.error;
    IS_Krylov_PS_CPU(ii,1:length(OutIS_Krylov_PS.times))=OutIS_Krylov_PS.times;

    IS_Krylov_CS_error(ii,1:length(OutIS_Krylov_CS.error))=OutIS_Krylov_CS.error;
    IS_Krylov_CS_CPU(ii,1:length(OutIS_Krylov_CS.times))=OutIS_Krylov_CS.times;

    IS_Krylov_GS_error(ii,1:length(OutIS_Krylov_GS.error))=OutIS_Krylov_GS.error;
    IS_Krylov_GS_CPU(ii,1:length(OutIS_Krylov_GS.times))=OutIS_Krylov_GS.times;
    
    IS_Krylov_SRHT_error(ii,1:length(OutIS_Krylov_SRHT.error))=OutIS_Krylov_SRHT.error;
    IS_Krylov_SRHT_CPU(ii,1:length(OutIS_Krylov_SRHT.times))=OutIS_Krylov_SRHT.times;

    fprintf('Done, ii=%d\n',ii)
end

%% setting the colour
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
lightyellow =[0.99 0.99 0.62];
mediumyellow =[0.91 0.91 .02];

%% plot
close all
%%%
xlable=1:(Max_length+1);

%%%
y1=IS_Krylov_PS_error';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%%
y2=IS_Krylov_CS_error';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%%
y3=IS_Krylov_GS_error';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%%
y4=IS_Krylov_SRHT_error';
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2);
y4q75=quantile(y4,0.75,2);

%% plot the iterations
figure
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y4q25' fliplr(y4q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
%%%
p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'IS_Krylov_CS');
p3=semilogy( xlable, median(y3'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'IS_Krylov_GS');
p4=semilogy( xlable, median(y4'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'DisplayName', 'IS_Krylov_SRHT');
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'IS_Krylov_PS');
set(gca, 'YScale', 'log')
ylim([10^(-12), 1])
xlim([0, 1200])
ylabel('RSE','Interpreter', 'latex')
xlabel('Number of iterations','Interpreter', 'latex')
legend([p1 p2 p3 p4],{'IS-Krylov-PS','IS-Krylov-CS','IS-Krylov-GS','IS-Krylov-SRHT'},'Interpreter', 'latex','location', 'best')
txt=title(['$\kappa=$ ',num2str(kappa),', $r=$ ',num2str(r)]);
set(txt, 'Interpreter', 'latex');
%% plot the CPU
figure
h = fill([median(IS_Krylov_CS_CPU)  fliplr(median(IS_Krylov_CS_CPU))], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(IS_Krylov_CS_CPU)  fliplr(median(IS_Krylov_CS_CPU))], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([median(IS_Krylov_GS_CPU)  fliplr(median(IS_Krylov_GS_CPU) )], [miny3 fliplr(maxy3)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(IS_Krylov_GS_CPU)  fliplr(median(IS_Krylov_GS_CPU) )], [y3q25' fliplr(y3q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([median(IS_Krylov_SRHT_CPU)  fliplr(median(IS_Krylov_SRHT_CPU))], [miny4 fliplr(maxy4)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(IS_Krylov_SRHT_CPU)  fliplr(median(IS_Krylov_SRHT_CPU))], [y4q25' fliplr(y4q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([median(IS_Krylov_PS_CPU)  fliplr(median(IS_Krylov_PS_CPU) )], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(IS_Krylov_PS_CPU)  fliplr(median(IS_Krylov_PS_CPU) )], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
%%%
p2=semilogy( median(IS_Krylov_CS_CPU), median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'IS_Krylov_CS');
hold on
p3=semilogy( median(IS_Krylov_GS_CPU), median(y3'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'IS_Krylov_GS');
hold on
p4=semilogy( median(IS_Krylov_SRHT_CPU), median(y4'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'DisplayName', 'IS_Krylov_SRHT');
hold on
p1=semilogy( median(IS_Krylov_PS_CPU), median(y1'), 'red', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'IS_Krylov_PS');
hold on
set(gca, 'YScale', 'log')

xlim([-0.01,0.7]);
ylim([10^(-12),10^(0)]);
ylabel('RSE','Interpreter', 'latex')
xlabel('CPU','Interpreter', 'latex')

txt=title(['$\kappa=$ ',num2str(kappa),', $r=$ ',num2str(r)]);
set(txt, 'Interpreter', 'latex');
legend([p1 p2 p3 p4],{'IS-Krylov-PS','IS-Krylov-CS','IS-Krylov-GS','IS-Krylov-SRHT'},'Interpreter', 'latex','location', 'best')
set(gca, 'YScale', 'log')
