function Z_star = SDP_CBW(D, Sigma_hat, Sigma_hat_half, rho)
% SDP_CBW - Solving problem 
%        min_{\Sigma} <D, \Sigma> s.t. W(\Sigma,\Sigma_hat) <= rho 
%
% Syntax: S = SDP_CBW(X, Sigma_hat, rho)
%
% Calculate the solution of projection into a BW-ball via SDP method 

    % Define auxiliary variable
    n = size(D, 1);
    
    % Initialization
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);

    % Define decision variables
    Z = sdpvar(n, n);
    E = sdpvar(n, n);

    % Declare objective function
    objective = trace(Z * D);

    % Define constraints
    constraints = [Z >= 0, E >= 0, trace(Z + Sigma_hat - 2*E) <= rho^2, ... 
                    [Sigma_hat_half * Z * Sigma_hat_half, E; E, eye(n)] >= 0];

    % Solving the Optimization Problem           
    diagnosis = optimize(constraints, objective, ops);
%     A_star = value(A);
%     b_star = mu_x - A_star * (H * mu_x + mu_w);
    Z_star = value(Z);
end