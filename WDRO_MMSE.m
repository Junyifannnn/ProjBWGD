clc
clear
rng(100);
addpath(genpath('utilities'));

d = 200;
run_count = 5; % 10 as default
tol = 1e-9;
verbose = true;
opts.iter_max = 500; % 1000 as default
opts.delta = 0.95;
opts.tol = tol;
opts.oracle = 'svd';
opts.verbose = verbose;
opts_PGD = opts;
opts_APGD = opts;
opts_EPGD = opts;

res_FW = NaN(opts.iter_max, run_count);
res_FAFW = NaN(opts.iter_max, run_count);
res_PGD = NaN(opts.iter_max, run_count);
res_APGD = NaN(opts.iter_max, run_count);
res_EGD = NaN(opts.iter_max, run_count);

obj_FW = NaN(opts.iter_max, run_count);
obj_FAFW = NaN(opts.iter_max, run_count);
obj_PGD = NaN(opts.iter_max, run_count);
obj_APGD = NaN(opts.iter_max, run_count);
obj_EGD = NaN(opts.iter_max, run_count);

for r = 1 : run_count

    rho_x = sqrt(d);
    rho_w = sqrt(d);
    fprintf('Running Iteration %d for d = %d \n', r, d);

    A = randn(d);
    [R_A, ~] = eig(A + A');
    lambda_x = 1 + 4 * rand(d,1);
    cov_x = R_A * diag(lambda_x) * R_A';

    B = randn(d);
    [R_B, ~] = eig(B + B');
    lambda_w = 1 + rand(d,1);
    cov_w = R_B * diag(lambda_w) * R_B';
    
    H = randn(d) / sqrt(d);

    opts.step_size = 'vanilla';
    [~, ~, obj_FW(:,r), res_FW(:,r)] = FrankWolfe(zeros(d,1), cov_x, rho_x, zeros(d,1), cov_w, rho_w, H, opts);
   
    opts_PGD.step_size = 'vanilla';  
    opts_PGD.strategy = 'BWGD';
    [~, ~, obj_PGD(:,r), res_PGD(:,r)] = PGD(zeros(d,1), cov_x, rho_x, zeros(d,1), cov_w, rho_w, H, opts_PGD);

    opts.step_size = 'full_adaptive';
    [~, ~, obj_FAFW(:,r), res_FAFW(:,r)] = FrankWolfe(zeros(d,1), cov_x, rho_x, zeros(d,1), cov_w, rho_w, H, opts);
    
    opts_APGD.step_size = 'Armijo';
    opts_APGD.strategy = 'BWGD';
    [~, ~, obj_APGD(:,r), res_APGD(:,r)] = PGD(zeros(d,1), cov_x, rho_x, zeros(d,1), cov_w, rho_w, H, opts_APGD);
    
    opts_EGD.step_size = 'vanilla';  
    opts_EGD.strategy = 'EGD';
    [~, ~, obj_EGD(:,r), res_EGD(:,r)] = PGD(zeros(d,1), cov_x, rho_x, zeros(d,1), cov_w, rho_w, H, opts_EGD);

end
%%
prc = 0;
alphaa = 0.1;
font_size = 20;
iter_print = 500;
colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.45, 0.6, 0.25; 0.2, 0.75, 0.35];
fig = figure;
set(fig, 'Units', 'normalized', 'Position', [0.35, 0.25, 0.4, 0.55])
hold on
p1 = plot_with_shade(1:iter_print, res_FW(1:iter_print,:), prc, alphaa, colors(1,:));
p2 = plot_with_shade(1:iter_print, res_PGD(1:iter_print,:), prc, alphaa, colors(2,:));
p3 = plot_with_shade(1:iter_print, res_FAFW(1:iter_print,:), prc, alphaa, colors(3,:));
p4 = plot_with_shade(1:iter_print, res_APGD(1:iter_print,:), prc, alphaa, colors(4,:));
p5 = plot_with_shade(1:iter_print, res_EGD(1:iter_print,:), prc, alphaa, colors(5,:));
set(gca, 'XScale', 'linear', 'YScale', 'log');
set(gca, 'FontSize', font_size - 2);
xlabel('iterations', 'FontSize', font_size);
ylabel('Surrogate duality gap','FontSize', font_size)
ylim([1e-8, 1e3]);
grid on
lgd = legend([p1, p2, p3, p4, p5], 'FW', 'BWGD', 'FAFW', 'Armijo BWGD','EGD', 'Location', 'east');
lgd.FontSize = font_size;
remove_border()
saveas(gcf, 'convergence', 'svg')