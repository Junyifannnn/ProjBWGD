function L_star = SDP_linear_oracle(D, Sigma_hat, Sigma_hat_half, lambda_min, rho)
% SDP_CBW - Solving problem 
%        max_{L} <L-Sigma, D> s.t. W(L,Sigma_hat) <= rho and L >= lambda_min*I 
%
% Syntax: S = SDP_linear_oracle(X, Sigma_hat, rho)
%
% Calculate the solution of projection into a BW-ball via SDP method 

    % Define auxiliary variable
    n = size(D, 1);
    
    % Initialization
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);

    % Define decision variables
    L = sdpvar(n, n);
    E = sdpvar(n, n);

    % Declare objective function
    objective = trace(L * D);

    % Define constraints
    constraints = [L >= 0, E >= 0, trace(L + Sigma_hat - 2*E) <= rho^2, ... 
                    [Sigma_hat_half * L * Sigma_hat_half, E; E, eye(n)] >= 0];

    % Solving the Optimization Problem           
    diagnosis = optimize(constraints, objective, ops);
%     A_star = value(A);
%     b_star = mu_x - A_star * (H * mu_x + mu_w);
    L_star = value(L);
end