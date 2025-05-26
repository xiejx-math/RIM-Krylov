%%
clear
close all

m=1024; % number of rows
n=128; % number of columns
r=128; % number of rank
kappa=10; % upper bound for the condition number
for i=1:10
    q_set(i)=2^i; % size of the block
end
l_set=[2,5,10,50,128]; % number of previous iterations
run_time=20; % average times
%% some vectors are used to store the computed results
Times_IS_Krylov_PS_set=zeros(10,run_time);
Iter_IS_Krylov_PS_set=zeros(10,run_time);

Times_RIMCS_AS_Krylov_set=zeros(10,run_time);
Iter_RIMCS_AS_Krylov_set=zeros(10,run_time);

%% execute the iteration
for ii=1:run_time
    %% generated the matrix A
    [U,~]=qr(randn(m, r), 0);
    [V,~]=qr(randn(n, r), 0);
    D = diag(1+(kappa-1).*rand(r, 1));
    A=U*D*V';
    clear U V D
    
    %% generated the right-hand vector b
    x=rand(n,1);
    b=A*x;
    xLS=lsqminnorm(A,b);
    
    for jj=1:10
        q=q_set(jj); % size of the block
        for kk=1:5
            %% parameter setup
            opts.Pre_iter=l_set(kk); % number of previous iterations
            %opts.TOL1=eps^2;
            opts.TOL1=10^(-20);
            opts.TOL=10^(-12);
            opts.xstar=xLS;
            
            %% run IS-Krylov-PS
            [x_IS_Krylov_PS,Out_IS_Krylov_PS]=My_IS_Krylov_PS(A,b,q,opts);

            %% store the compute results
            Times_IS_Krylov_PS_set(jj,ii,kk)=Out_IS_Krylov_PS.times(end);
            Iter_IS_Krylov_PS_set(jj,ii,kk)=Out_IS_Krylov_PS.iter*q/m;
        end
    end
    fprintf('Done, ii=%d\n',ii)
end
%% plot the result

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


%% number of iterations
num_iter_array=log2(q_set)';
close all
xlable=1:(size(q_set,2));

%%%
y1=Iter_IS_Krylov_PS_set(:,:,1);
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2)';
y1q75=quantile(y1,0.75,2)';

y2=Iter_IS_Krylov_PS_set(:,:,2);
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2)';
y2q75=quantile(y2,0.75,2)';

y3=Iter_IS_Krylov_PS_set(:,:,3);
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2)';
y3q75=quantile(y3,0.75,2)';

y4=Iter_IS_Krylov_PS_set(:,:,4);
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2)';
y4q75=quantile(y4,0.75,2)';

y5=Iter_IS_Krylov_PS_set(:,:,5);
miny5=min(y5');
maxy5=max(y5');
y5q25=quantile(y5,0.25,2)';
y5q75=quantile(y5,0.75,2)';

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25 fliplr(y1q75)],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y2q25 fliplr(y2q75)],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y3q25 fliplr(y3q75)],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y4q25 fliplr(y4q75)],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny5 fliplr(maxy5)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y5q25 fliplr(y5q75)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)

hold on
p1=semilogy( xlable, median(y1'), 'black', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','o');
p2=semilogy( xlable, median(y2'), 'red', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','s');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','^');
p4=semilogy( xlable, median(y4'), 'blue', 'LineWidth', 1,...
    'LineStyle', '-','Marker','v');
p5=semilogy( xlable, median(y5'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','p');
xlim([1,10]);
ylim([0,100]);
ylabel('$k\cdot\frac{q}{m}$','Interpreter', 'latex')
xlabel('$\log_2(q)$','Interpreter', 'latex')
hold on
txt=title(['$\kappa=$ ',num2str(kappa),', $r=$ ',num2str(r)]);
set(txt, 'Interpreter', 'latex');
legend([p1 p2 p3 p4 p5],{'$\ell=2$', '$\ell=5$', '$\ell=10$', '$\ell=50$', '$\ell=128$'},'Interpreter', 'latex','location', 'best')

%% CPU
num_iter_array=log2(q_set)';
xlable=1:(size(q_set,2));

%%%
y1=Times_IS_Krylov_PS_set(:,:,1);
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2)';
y1q75=quantile(y1,0.75,2)';

y2=Times_IS_Krylov_PS_set(:,:,2);
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2)';
y2q75=quantile(y2,0.75,2)';

y3=Times_IS_Krylov_PS_set(:,:,3);
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2)';
y3q75=quantile(y3,0.75,2)';

y4=Times_IS_Krylov_PS_set(:,:,4);
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2)';
y4q75=quantile(y4,0.75,2)';

y5=Times_IS_Krylov_PS_set(:,:,5);
miny5=min(y5');
maxy5=max(y5');
y5q25=quantile(y5,0.25,2)';
y5q75=quantile(y5,0.75,2)';

figure

h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25 fliplr(y1q75)],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y2q25 fliplr(y2q75)],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y3q25 fliplr(y3q75)],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y4q25 fliplr(y4q75)],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [miny5 fliplr(maxy5)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y5q25 fliplr(y5q75)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)

hold on
p1=semilogy( xlable, median(y1'), 'black', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','o');
p2=semilogy( xlable, median(y2'), 'red', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','s');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','^');
p4=semilogy( xlable, median(y4'), 'blue', 'LineWidth', 1,...
    'LineStyle', '-','Marker','v');
p5=semilogy( xlable, median(y5'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker','p');
xlim([1,10]);
ylim([0,0.2]);
ylabel('CPU','Interpreter', 'latex')
xlabel('$\log_2(q)$','Interpreter', 'latex')
hold on
txt=title(['$\kappa=$ ',num2str(kappa),', $r=$ ',num2str(r)]);
set(txt, 'Interpreter', 'latex');
legend([p1 p2 p3 p4 p5],{'$\ell=2$', '$\ell=5$', '$\ell=10$', '$\ell=50$', '$\ell=128$'},'Interpreter', 'latex','location', 'best')

