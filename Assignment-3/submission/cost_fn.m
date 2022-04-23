function [fval] = cost_fn(sigma,A,z)
    [m, n_samples] = size(A);
    % Construct RHS (A*sigma_e*A')
    RHS = A*diag(sigma.^2)*A';
    % Since matrix is symmetric we just take the unique elements
    RHS_vec = [];
    for i = 1:m
        RHS_vec = [RHS_vec; RHS(i:end,i)];
    end
    % Estimate Sr from our fit
    n_samples = size(z,2);
    r_mat = A*z;
    Srhat = r_mat*r_mat'/n_samples;
    % Once again, since matrix is symmetric we just take the unique elements
    LHS_vec = [];
    for i = 1:m
        LHS_vec = [LHS_vec; Srhat(i:end,i)];
    end
    error = LHS_vec - RHS_vec;
    fval = norm(error);
end

