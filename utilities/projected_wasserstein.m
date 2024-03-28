function [S, flag] = projected_wasserstein(X, Sigma_hat, rho)
% Calculate the solution of projected wasserstein distance: 
%          min_S  W(S,X) s.t. W(S,Sigma_hat) < rho
%
% Syntax: res = projected_wasserstein(X, Sigma_hat, rho)

X_Sigma_hat_1_2 = sqrtm(X * Sigma_hat);
Sigma_hat_X_1_2 = sqrtm(Sigma_hat * X);
W_D_X_Sigma_hat = sqrt(trace(X) + trace(Sigma_hat) - 2 * trace(X_Sigma_hat_1_2));
if W_D_X_Sigma_hat <= rho
    S = X;
    flag = 0;
else
    gamma = W_D_X_Sigma_hat / rho - 1;
    S = 1 / (1 + gamma)^2 * X + (gamma / (1 + gamma))^2 * Sigma_hat + gamma / (1 + gamma)^2 * (X_Sigma_hat_1_2 + Sigma_hat_X_1_2);
    flag = 1;
end

end