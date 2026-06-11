function bm = benchmarks(ret_w, F_w, basis_win, w_prev, p)
%benchmark portfolio weights - same signal logic as layer2

cvx_clear;
M = 4; K = p.K;

Sigma = cov(ret_w);
Sigma = lw_shrink(ret_w, Sigma);

B = ((F_w'*F_w) + p.ridge_lam*eye(K)) \ (F_w'*ret_w);
B = B';

%factor fallback
n_smooth  = min(4, size(F_w,1));
f_smooth  = mean(F_w(end-n_smooth+1:end,:))';
mu_factor = B * f_smooth;
mu_mv     = mean(ret_w)';

%per-maturity basis z-score
mu_hat    = mu_factor;
basis_now = basis_win(end,:);
for m = 1:M
    if ~isnan(basis_now(m)) && size(basis_win,1) >= 8
        bvals = basis_win(:,m);
        bvals = bvals(~isnan(bvals));
        if length(bvals) >= 8
            mu_b  = mean(bvals);
            sig_b = max(std(bvals), 1e-8);
            z_b   = (basis_now(m) - mu_b) / sig_b;
            z_b   = max(-2, min(2, z_b));
            scale = std(ret_w(:,m));
            mu_hat(m) = z_b * scale;
        end
    end
end

T       = size(ret_w,1);
Sigma_f = cov(F_w);
eps2    = chi2inv(p.chi2_alpha, K) * trace(Sigma_f) / T;
epsilon = sqrt(eps2);
Phi     = Sigma_f + 1e-6*eye(K);

%tikhonov matrix: path-graph laplacian for 4 maturities (row sums = 0)
L_tikh = diag([1 2 2 1]) - diag([1 1 1],1) - diag([1 1 1],-1);

%cvx robust
try
    cvx_begin quiet
        variable w_r(M)
        minimize( -w_r'*mu_hat + ...
                  epsilon*norm(Phi^(0.5)*B'*w_r) + ...
                  p.gamma*quad_form(w_r,Sigma) + ...
                  p.lambda_T*quad_form(w_r,L_tikh) + ...
                  p.kappa*norm(w_r-w_prev,1) )
        subject to
            sum(w_r) == 1
            w_r >= p.w_min
            w_r <= p.w_max
    cvx_end
    if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
        bm.w_rob = double(w_r);
    else
        bm.w_rob = ones(M,1)/M;
    end
catch
    bm.w_rob = ones(M,1)/M;
end

%mean-variance uses rolling mean
try
    cvx_begin quiet
        variable w_mv(M)
        minimize( -w_mv'*mu_mv + p.gamma*quad_form(w_mv,Sigma) )
        subject to
            sum(w_mv) == 1
            w_mv >= p.w_min
            w_mv <= p.w_max
    cvx_end
    if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
        bm.w_mv = double(w_mv);
    else
        bm.w_mv = ones(M,1)/M;
    end
catch
    bm.w_mv = ones(M,1)/M;
end

bm.w_eq = ones(M,1)/M;
end