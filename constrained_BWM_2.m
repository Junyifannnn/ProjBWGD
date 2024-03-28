% comparison of RFW and BWGD for Wasserstein mean %
% RFW is Riemannian Frank-wolfe method from "Riemannian Frank-Wolfe algorithm for optimization on manifolds".

clc
clear 
rng(100);
addpath(genpath('utilities'));

%% Parameter
n = 60; % dim
m = 20; % #size
iter_max = 50; % #max iterations

run_count = 1;  

obj_BWGD = zeros(iter_max, run_count);
obj_RFW = zeros(iter_max, run_count);
obj_BWGD_noproj = zeros(iter_max, run_count);

%% Generate collection of PSD
A = genPosdef2(n,m);

%% Main Iteration
for r = 1 : run_count
    
    % Riemannian Frank-Wolfe
    % Initialization
    fprintf('Runing Riemannian Frank-Wolfe \n');
    am = arithmeticMean(A);
    lmin = smallest_eval(A, n, m);
    lb = lmin * eye(size(am));
%     ub = mean(eig(am)) * eye(size(am));
    ub = am;
    x  = lb;  % only use the lower bound for initialization
    % Riemannian Frank-Wolfe main iteration
    for i = 1 : iter_max
        obj_RFW(i, r) = wmobj(x,A);
        s = gradient_WM(x,A);
        zk = FWR_dir(x,real(s),lb,ub);
        al = 0.7 / i; % alpha_RFW = 2 / (i + 1);
        xh = x^-0.5;
        xh2 = x^0.5;
        x = xh2 * (xh * zk * xh)^al * xh2;
    end
    
    % BWGD
    % Initialization
    fprintf('Runing BWGD \n');
    eta_BWGD = 0.2;
    x = lb;
    % BWGD main iteration
    
    for i = 1 : iter_max
        obj_BWGD(i, r) = wmobj(x,A);     
        grad_BWGD = 2 * gradient_WM(x,A);
        x_tmp = (eye(n) - eta_BWGD * grad_BWGD) * x * (eye(n) - eta_BWGD * grad_BWGD);
        x = SDP_constrained_BWM(x_tmp, lb, ub);
    end
    
    fprintf('Runing BWGD_noproj \n');
    eta_BWGD = 0.2;
    x = lb;
    % BWGD main iteration
    for i = 1 : iter_max
        obj_BWGD_noproj(i, r) = wmobj(x,A);     
        grad_BWGD = 2 * gradient_WM(x,A);
        x = (eye(n) - eta_BWGD * grad_BWGD) * x * (eye(n) - eta_BWGD * grad_BWGD);
    end
end

prc = 0;
alphaa = 0.1;
font_size = 20;
colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.45, 0.6, 0.25; 0.2, 0.75, 0.35];
fig = figure;
set(fig, 'Units', 'normalized', 'Position', [0.35, 0.25, 0.4, 0.55])
hold on
p1 = plot_with_shade(1:iter_max, obj_RFW, prc, alphaa, colors(1,:));
p2 = plot_with_shade(1:iter_max, obj_BWGD, prc, alphaa, colors(2,:));
set(gca, 'XScale', 'linear', 'YScale', 'linear');
set(gca, 'FontSize', font_size - 2);
xlabel('iterations', 'FontSize', font_size);
ylabel('objective function','FontSize', font_size)
grid on
lgd = legend([p1, p2], 'RFW','BWGD', 'Location', 'northeast');
lgd.FontSize = font_size;
remove_border()