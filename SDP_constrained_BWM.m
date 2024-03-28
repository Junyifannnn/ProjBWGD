function Z_star = SDP_constrained_BWM(D, lb, ub)
% SDP_CBW - Solving problem 
%        min_{\Sigma} d(D, \Sigma) s.t. W(\Sigma,\Sigma_hat) <= rho 
%
% Syntax: S = SDP_CBW(X, Sigma_hat, rho)
%
% Calculate the solution of projection into a BW-ball via SDP method 

    % Define auxiliary variable
    D = (D + D') / 2;
    n = size(D, 1);
    
    % Initialization
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);

    % Define decision variables
    Z = sdpvar(n, n);
    C = sdpvar(n, n, 'full');

    % Declare objective function
    objective = trace(Z + D - 2 * C);

    % Define constraints
    constraints = [Z >= 0, ub - Z >= 0, Z >= lb,... 
                    [D, C; C', Z] >= 0];

    % Solving the Optimization Problem           
    diagnosis = optimize(constraints, objective, ops);
%     A_star = value(A);
%     b_star = mu_x - A_star * (H * mu_x + mu_w);
    Z_star = value(Z);
end