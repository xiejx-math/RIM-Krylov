% This Matlab file is used to compare RABK, SCGP, and IS-Krylov-PS
% based on the data from LIBSVM

close all;
clear all;

%%
run_time=20; % average times

%% generated the matrix A using the data from LIBSVM
%load aloi;A=aloi_inst;Max_length=2000;opts.Pre_iter=50;

%load a9a;A=a9a_inst;Max_length=600;opts.Pre_iter=50;

%load cod-rna;A=cod_rna_inst;Max_length=80;opts.Pre_iter=5;

load protein.mat;A=protein_inst;Max_length=3000;opts.Pre_iter=50;

%% some vectors are used to store the desired numerical results

RABK_CPU=zeros(run_time,Max_length+1);
RABK_error=zeros(run_time,Max_length+1);

SCGP_CPU=zeros(run_time,Max_length+1);
SCGP_error=zeros(run_time,Max_length+1);

IS_Krylov_PS_CPU=zeros(run_time,Max_length+1);
IS_Krylov_PS_error=zeros(run_time,Max_length+1);

%% executing "run_time" times of the algorithms
for ii=1:run_time
     
    %% generated the right-hand vector b
    [m,n]=size(A);
    x=randn(n,1);
    b=A*x;
    xLS=lsqminnorm(A,b);
    q=300; % size of the block
    %% parameter setup
    opts.xstar=xLS;
    opts.TOL1=eps^2;
    opts.TOL=10^(-29);
    opts.Max_iter=Max_length;

    opts.permS=randperm(m);
    S=opts.permS;
    
    %%
    tau=floor(m/q);
    A=A(S,:);
    b=b(S);
    blockAnormfro=zeros(tau,1);
    normAfro=norm(A,'fro')^2;
    beta_max=0;
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
    %%
    opts.Aarrs=Aarrs;
    opts.barrs=barrs;
    opts.cumsumpro=cumsumpro;
    opts.blockAnormfro=blockAnormfro;
    opts.probset=1;
   
    %% run 
    [xRABK,OutRABK]=My_RABK(A,b,q,opts);
    [xSCGP,OutSCGP]=My_AmRABK(A,b,q,opts);
    [xIS_Krylov_PS,OutIS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);
    %% store the numerical results

    RABK_error(ii,1:length(OutRABK.error))=OutRABK.error;
    RABK_CPU(ii,1:length(OutRABK.times))=OutRABK.times;

    SCGP_error(ii,1:length(OutSCGP.error))=OutSCGP.error;
    SCGP_CPU(ii,1:length(OutSCGP.times))=OutSCGP.times;

    IS_Krylov_PS_error(ii,1:length(OutIS_Krylov_PS.error))=OutIS_Krylov_PS.error;
    IS_Krylov_PS_CPU(ii,1:length(OutIS_Krylov_PS.times))=OutIS_Krylov_PS.times;

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

%% plot RSE
close all
%%%
xlable=1:(Max_length+1);

%%%
y1=RABK_error';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%%
y2=SCGP_error';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%%
y3=IS_Krylov_PS_error';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
%%%
hold on
p1=semilogy( xlable, median(y1'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'DisplayName', 'RABK');
set(gca, 'YScale', 'log')
hold on
p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'SCGP');
hold on
p3=semilogy( xlable, median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'IS-Krylov-PS');
hold on
ylim([10^(-20), 1])
xlim([1, Max_length])
ylabel('RSE','Interpreter', 'latex')
xlabel('Number of iterations','Interpreter', 'latex')
legend([ p1 p2 p3],{'RABK','SCGP','IS-Krylov-PS, $\ell=50$'},'Interpreter', 'latex','location', 'best')
txt=title(['\texttt{protein}, ','$m=$ ',num2str(m),', $n=$ ',num2str(n)]);
set(txt, 'Interpreter', 'latex');
%% plot CPU
figure
h = fill([median(RABK_CPU)  fliplr(median(RABK_CPU) )], [miny1 fliplr(maxy1)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(RABK_CPU)  fliplr(median(RABK_CPU) )], [y1q25' fliplr(y1q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([median(SCGP_CPU)  fliplr(median(SCGP_CPU) )], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(SCGP_CPU)  fliplr(median(SCGP_CPU) )], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([median(IS_Krylov_PS_CPU)  fliplr(median(IS_Krylov_PS_CPU))], [miny3 fliplr(maxy3)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(IS_Krylov_PS_CPU)  fliplr(median(IS_Krylov_PS_CPU))], [y3q25' fliplr(y3q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)

hold on
set(gca, 'YScale', 'log')
%%%
hold on
p1=semilogy( median(RABK_CPU), median(y1'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'DisplayName', 'RABK');
hold on
p2=semilogy( median(SCGP_CPU), median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'SCGP');
hold on
p3=semilogy( median(IS_Krylov_PS_CPU), median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'IS-Krylov-PS');
set(gca, 'YScale', 'log')

xlim([0,0.13]);
ylim([10^(-20),10^(0)]);
ylabel('RSE','Interpreter', 'latex')
xlabel('CPU','Interpreter', 'latex')

txt=title(['\texttt{protein}, ','$m=$ ',num2str(m),', $n=$ ',num2str(n)]);
legend([ p1 p2 p3],{'RABK','SCGP','IS-Krylov-PS, $\ell=50$'},'Interpreter', 'latex','location', 'best')
set(txt, 'Interpreter', 'latex');
