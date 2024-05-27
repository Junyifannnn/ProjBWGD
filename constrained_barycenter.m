clc
clear
rng(100);
addpath(genpath('utilities'));

d = 30;
n = 50; % number of matrices
run_count = 1; 
tol = 1e-7;
verbose = 1;
lambda_min = 0.1; lambda_max = 100;
iter_max = 50; 

res_FW = zeros(iter_max, run_count);
res_PGD = zeros(iter_max, run_count);
res_EGD = zeros(iter_max, run_count);
res_FAFW = zeros(iter_max, run_count);
res_APGD = zeros(iter_max, run_count);

obj_FW = zeros(iter_max, run_count);
obj_PGD = zeros(iter_max, run_count);
obj_EGD = zeros(iter_max, run_count);
obj_FAFW = zeros(iter_max, run_count);
obj_APGD = zeros(iter_max, run_count);

for r = 1 : run_count

    rho = 4 * sqrt(d); % 4 sqrt(d)
    fprintf('Running constrained barycenter problem in %d-th Iteration for rho = %d \n', r, d);

    noise_hat = randn(d);
    [R_hat, ~] = eig(noise_hat + noise_hat');
    lambda_hat = 1 + 1 * rand(d,1);
    Sigma_hat = R_hat * diag(lambda_hat) * R_hat';
    Sigma_hat_half = R_hat * diag(sqrt(lambda_hat)) * R_hat';
    
    noise_ini = randn(d);
    
    Sigma_i = zeros(d, d, n);
    Sigma_i_half = zeros(d, d, n);
    beta = rand(n, 1);
    beta = beta / sum(beta);
    for num = 1 : n
        noise_i = randn(d);
        [R_i, ~] = eig(noise_i + noise_i');
        lambda_i = lambda_min + lambda_max * rand(d,1);
        Sigma_i(:, :, num) = R_i * diag(lambda_i) * R_i';
        Sigma_i_half(:, :, num) = R_i * diag(sqrt(lambda_i)) * R_i';
    end
        
    % BWGD
    % Initialization
    fprintf('Runing BWGD \n');
    eta_PGD = 0.5;
    Sigma_PGD = Sigma_hat;
    grad_tmp_PGD = zeros(d, d, n);
    % BWGD main iteration
    for i = 1 : iter_max
        res_PGD(i, r) = 0;
        obj_PGD(i, r) = 0;
        grad_BW = eye(d);
        for num = 1 : n
            grad_tmp_PGD(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_PGD * Sigma_i_half(:, :, num));
            obj_PGD(i, r) = obj_PGD(i, r) + beta(num) * trace(Sigma_PGD + Sigma_i(:, :, num) - 2 * grad_tmp_PGD(:, :, num));
            grad_BW = grad_BW - beta(num) * (Sigma_i_half(:, :, num) * (grad_tmp_PGD(:, :, num) \ Sigma_i_half(:, :, num)));
        end
        L_PGD = SDP_CBW(grad_BW, Sigma_hat, Sigma_hat_half, rho);
        res_PGD(i, r) = abs(trace((L_PGD - Sigma_PGD)' * grad_BW));
        if res_PGD(i ,r) <= tol
            break
        end
        Sigma_PGD_1_2 = (eye(d) - eta_PGD * grad_BW) * Sigma_PGD * (eye(d) - eta_PGD * grad_BW);
        [Sigma_PGD_next, ~] = projected_wasserstein(Sigma_PGD_1_2, Sigma_hat, rho);
        Sigma_PGD = Sigma_PGD_next;
    end
%     
    % Frank-Wolfe
    % Initialization
    fprintf('Runing Frank-Wolfe \n');
    Sigma_FW = Sigma_hat;
    grad_tmp_FW = zeros(d, d, n);
    % Frank-Wolfe main iteration
    for i = 1 : iter_max
        eta_FW = 2 / (i + 1);
        res_FW(i, r) = 0;
        obj_FW(i, r) = 0;
        grad_E = eye(d);
        for num = 1 : n
            grad_tmp_FW(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_FW * Sigma_i_half(:, :, num));
            obj_FW(i, r) = obj_FW(i, r) + beta(num) * trace(Sigma_FW + Sigma_i(:, :, num) - 2 * grad_tmp_FW(:, :, num));
            grad_E = grad_E - beta(num) * (Sigma_i_half(:, :, num) * (grad_tmp_FW(:, :, num) \ Sigma_i_half(:, :, num)));
        end
        L_FW = SDP_CBW(grad_E, Sigma_hat, Sigma_hat_half, rho);
        res_FW(i, r) = abs(trace((L_FW - Sigma_FW)' * grad_E));
        if res_FW(i, r) <= tol
            break
        end
        Sigma_FW_next = (1 - eta_FW) * Sigma_FW + eta_FW * L_FW;
        Sigma_FW = Sigma_FW_next;
    end

    % Fully adaptive Frank-Wolfe
    fprintf('Runing Fully adaptive Frank-Wolfe \n');
    % Initialization
    mu = 0.6; % 0.8
    tau = 1.8;  % 1.5
    beta_FAFW = 1e3; %1e2
    beta_k = beta_FAFW / tau^10;
    Sigma_FAFW = Sigma_hat;
    grad_tmp_FAFW = zeros(d, d, n);
    grad_tmp2_FAFW = zeros(d, d, n);
    % main iteration
    for i = 1 : iter_max
        res_FAFW(i, r) = 0;
        obj_FAFW(i, r) = 0;
        grad_E = eye(d);
        for num = 1 : n
            grad_tmp_FAFW(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_FAFW * Sigma_i_half(:, :, num));
            obj_FAFW(i, r) = obj_FAFW(i, r) + beta(num) * trace(Sigma_FAFW + Sigma_i(:, :, num) - 2 * grad_tmp_FAFW(:, :, num));
            grad_E = grad_E - beta(num) * (Sigma_i_half(:, :, num) * (grad_tmp_FAFW(:, :, num) \ Sigma_i_half(:, :, num)));
        end
        L_FAFW = SDP_CBW(grad_E, Sigma_hat, Sigma_hat_half, rho);
        g_k = abs(trace((L_FAFW - Sigma_FAFW)' * grad_E));
        P_k = L_FAFW - Sigma_FAFW;
        res_FAFW(i, r) = abs(g_k);
        if res_FAFW(i ,r) <= tol
            break
        end
        P_norm = P_k(:)' * P_k(:);
        beta_k = beta_k * mu;
            while true
                eta_k = min(g_k/(beta_k * P_norm), 1);
                obj_tmp = 0; % calculate f(Sigma_FAFW + eta_k * P_k)
                Sigma_tmp_FAFW = Sigma_FAFW + eta_k * P_k;
                for num = 1 : n
                    grad_tmp2_FAFW(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_tmp_FAFW * Sigma_i_half(:, :, num));
                    obj_tmp = obj_tmp + beta(num) * trace(Sigma_tmp_FAFW + Sigma_i(:, :, num) - 2 * grad_tmp2_FAFW(:, :, num));
                end
                if obj_tmp <= obj_FAFW(i, r) - eta_k * g_k + 0.5 * eta_k^2 * beta_k * P_norm
                    break
                end
                beta_k = beta_k * tau;
            end              
        Sigma_FAFW_next = Sigma_FAFW + eta_k * P_k;
        Sigma_FAFW = Sigma_FAFW_next;
    end

    % Armoji BWGD
    fprintf('Runing Armoji BWGD \n');    
    % Initialization 
    eta_APGD = 0.5;
    tau = 2;
    Sigma_APGD = Sigma_hat;
    grad_tmp_APGD = zeros(d, d, n);
    grad_tmp2_APGD = zeros(d, d, n);
    % main iteration
    for i = 1 : iter_max
        res_APGD(i, r) = 0;
        obj_APGD(i, r) = 0;
        grad_BW = eye(d);
        for num = 1 : n
            grad_tmp_APGD(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_APGD * Sigma_i_half(:, :, num));
            obj_APGD(i, r) = obj_APGD(i, r) + beta(num) * trace(Sigma_APGD + Sigma_i(:, :, num) - 2 * grad_tmp_APGD(:, :, num));
            grad_BW = grad_BW - beta(num) * (Sigma_i_half(:, :, num) * (grad_tmp_APGD(:, :, num) \ Sigma_i_half(:, :, num)));
        end
        beta_k = 1;
        while true
            eta_APGD = beta_k;
            obj_tmp = 0;
            Sigma_tmp_APGD_1_2 = (eye(d) - eta_APGD * grad_BW) * Sigma_APGD * (eye(d) - eta_APGD * grad_BW);
            [Sigma_tmp_APGD, ~] = projected_wasserstein(Sigma_tmp_APGD_1_2, Sigma_hat, rho);
            cur_half = sqrtm(Sigma_APGD);
            for num = 1 : n
                grad_tmp2_APGD(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_tmp_APGD * Sigma_i_half(:, :, num));
                obj_tmp = obj_tmp + beta(num) * trace(Sigma_tmp_APGD + Sigma_i(:, :, num) - 2 * grad_tmp2_APGD(:, :, num));
            end
            tmp_val = trace(inv(cur_half)*grad_BW*cur_half*sqrtm(cur_half * Sigma_APGD * cur_half))-trace(grad_BW*Sigma_APGD);
            if obj_tmp <= obj_APGD(i, r) + 0.5 * tmp_val
                break
            end
            beta_k = beta_k / tau;
        end 
        L_APGD = SDP_CBW(grad_BW, Sigma_hat, Sigma_hat_half, rho);
        res_APGD(i, r) = abs(trace((L_APGD - Sigma_APGD)' * grad_BW));
        if i > 1
            if res_APGD(i, r) <= tol || (abs(res_APGD(i, r) - res_APGD(i - 1, r)) < tol)
                break
            end
        end
        Sigma_APGD_1_2 = (eye(d) - eta_APGD * grad_BW) * Sigma_APGD * (eye(d) - eta_APGD * grad_BW);
        [Sigma_APGD_next, ~] = projected_wasserstein(Sigma_APGD_1_2, Sigma_hat, rho);
        Sigma_APGD = Sigma_APGD_next;
    end
    
    % EGD
    fprintf('Runing EGD \n');
    % Initialization
    eta_EGD = 0.5;
    Sigma_EGD = Sigma_hat;
    grad_tmp_EGD = zeros(d, d, n);
    % BWGD main iteration
    for i = 1 : iter_max
        res_EGD(i, r) = 0;
        obj_EGD(i, r) = 0;
        grad_EGD = eye(d);
        for num = 1 : n
            grad_tmp_EGD(:, :, num) = sqrtm(Sigma_i_half(:, :, num) * Sigma_EGD * Sigma_i_half(:, :, num));
            obj_EGD(i, r) = obj_EGD(i, r) + beta(num) * trace(Sigma_EGD + Sigma_i(:, :, num) - 2 * grad_tmp_EGD(:, :, num));
            grad_EGD = grad_EGD - beta(num) * (Sigma_i_half(:, :, num) * (grad_tmp_EGD(:, :, num) \ Sigma_i_half(:, :, num)));
        end
        L_EGD = SDP_CBW(grad_EGD, Sigma_hat, Sigma_hat_half, rho);
        res_EGD(i, r) = abs(trace((L_EGD - Sigma_EGD)' * grad_EGD));
        if res_EGD(i ,r) <= tol
            break
        end
        Sigma_EGD_1_2 = Sigma_EGD - eta_k * grad_EGD;
        [Sigma_EGD_next, ~] = projected_wasserstein(Sigma_EGD_1_2, Sigma_hat, rho);
        Sigma_EGD = Sigma_EGD_next;
    end
end

prc = 0;
alphaa = 0.1;
font_size = 20;
colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.45, 0.6, 0.25; 0.2, 0.75, 0.35];
fig = figure;
set(fig, 'Units', 'normalized', 'Position', [0.35, 0.25, 0.4, 0.55])
hold on
p1 = plot_with_shade(1:iter_max, res_FW, prc, alphaa, colors(1,:));
p2 = plot_with_shade(1:iter_max, res_PGD, prc, alphaa, colors(2,:));
p3 = plot_with_shade(1:iter_max, res_FAFW, prc, alphaa, colors(3,:));
p4 = plot_with_shade(1:iter_max, res_APGD, prc, alphaa, colors(4,:));
p5 = plot_with_shade(1:iter_max, res_EGD, prc, alphaa, colors(5,:));
set(gca, 'XScale', 'linear', 'YScale', 'log');
set(gca, 'FontSize', font_size - 2);
xlabel('iterations', 'FontSize', font_size);
ylabel('Surrogate duality gap','FontSize', font_size)
ylim([1e-7, 1e4]);
grid on
lgd = legend([p1, p2, p3, p4, p5], 'FW', 'BWGD', 'FAFW', 'Armijo BWGD','EGD', 'Location', 'southeast');
lgd.FontSize = font_size;
remove_border()
