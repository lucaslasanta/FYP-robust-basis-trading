function out = layer2_portfolio(ret_w, F_w, basis_win, w_prev, state_prev, p)
%robust portfolio allocation using lpnn (layer 2)
%expected return signal: basis z-score per maturity (mean reversion)
%TU uses factor fallback when basis not available (pre July 2021)

M = 4;
K = p.K;

Sigma = cov(ret_w);
Sigma = lw_shrink(ret_w, Sigma);

B = ((F_w'*F_w) + p.ridge_lam*eye(K)) \ (F_w'*ret_w);
B = B';

%factor fallback signal
n_smooth = min(4, size(F_w,1));
f_smooth = mean(F_w(end-n_smooth+1:end,:))';
mu_factor = B * f_smooth;

%basis z-score signal per maturity
mu_hat = mu_factor;   %initialise with factor signal
n_basis_rows = size(basis_win,1);
basis_now    = basis_win(end,:);

for m = 1:M
    if ~isnan(basis_now(m)) && n_basis_rows >= 8
        %sufficient history: compute z-score for this maturity
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
    %else keep factor signal for this maturity
end

%ellipsoidal uncertainty
T       = size(ret_w,1);
Sigma_f = cov(F_w);
eps2    = chi2inv(p.chi2_alpha, K) * trace(Sigma_f) / T;
epsilon = sqrt(eps2);
Phi     = Sigma_f + 1e-6*eye(K);

%tikhonov matrix: path-graph laplacian for 4 maturities (row sums = 0)
L_tikh = diag([1 2 2 1]) - diag([1 1 1],1) - diag([1 1 1],-1);

w_k   = w_prev;
lam_k = state_prev.lambda;
mu_k  = state_prev.mu;
tau   = p.tau;
c     = p.c;

for iter = 1:p.n_iter
    grad_ret  = -mu_hat;

    BPBt = B*Phi*B';
    nv   = sqrt(max(w_k'*BPBt*w_k,0));
    if nv > 1e-10
        grad_rob = epsilon*BPBt*w_k/nv;
    else
        grad_rob = zeros(M,1);
    end

    grad_risk  = 2*p.gamma*Sigma*w_k;
    grad_tikh  = 2*p.lambda_T*L_tikh*w_k;
    grad_tc    = p.kappa*sign(w_k - w_prev);
    grad_f     = grad_ret + grad_rob + grad_risk + grad_tikh + grad_tc;

    h_bud       = sum(w_k) - 1;
    grad_bud    = (lam_k + c*h_bud)*ones(M,1);
    g_lo        = p.w_min - w_k;
    g_hi        = w_k - p.w_max;
    act_lo      = g_lo > 0;
    act_hi      = g_hi > 0;
    grad_bnds   = -(mu_k(1:4)+c*g_lo).*act_lo + (mu_k(5:8)+c*g_hi).*act_hi;

    w_k       = w_k - tau*(grad_f + grad_bud + grad_bnds);
    w_k       = max(p.w_min, min(p.w_max, w_k));
    lam_k     = lam_k + tau*h_bud;
    mu_k(1:4) = mu_k(1:4) + tau*(g_lo.*act_lo);
    mu_k(5:8) = mu_k(5:8) + tau*(g_hi.*act_hi);
end

h_f  = sum(w_k)-1;
g_lf = p.w_min-w_k; alf = g_lf>0;
g_hf = w_k-p.w_max; ahf = g_hf>0;
gt   = grad_f + grad_bud + grad_bnds;
E    = 0.5*norm(gt)^2 + 0.5*h_f^2 + ...
       0.5*sum(g_lf(alf).^2) + 0.5*sum(g_hf(ahf).^2);

%count how many maturities are using basis signal
n_basis_active = sum(~isnan(basis_now) & ...
    arrayfun(@(m) length(basis_win(~isnan(basis_win(:,m)),m))>=8, 1:M));

out.w              = w_k;
out.lambda         = lam_k;
out.mu             = mu_k;
out.lyap_energy    = E;
out.epsilon        = epsilon;
out.mu_hat         = mu_hat;
out.signal_source  = sprintf('basis_%d/4', n_basis_active);
end