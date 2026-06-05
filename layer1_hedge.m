function h = layer1_hedge(r_basis, f_level, h_prev, lam_prev, p)
%robust hedge of basis return against level factor
%finds h minimising worst-case var(r_basis - h*f_level) over uncertainty in cov estimate

n      = length(f_level);
sig2_f = max(var(f_level), 1e-10);
cv     = cov(r_basis, f_level);
sig_bf = cv(1,2);

%95% confidence interval on covariance estimate
se_bf = sqrt(max(var(r_basis)*sig2_f + sig_bf^2, 0)) / sqrt(n);
delta = 1.96 * se_bf;

%allow both signs: level factor can help or hurt the basis
h_min = -2.0;
h_max =  2.0;

h_k    = h_prev;
lam_lo = lam_prev(1);
lam_hi = lam_prev(2);

%adaptive step: hessian of objective is 2*sig2_f; need tau < 1/(2*sig2_f) for stability
tau1 = min(p.tau, 0.9 / (2*sig2_f));

for iter = 1:p.n_iter
    %gradient of worst-case objective: d/dh [sig2_f*h^2 - 2h*(sig_bf - delta)]
    grad_obj = 2*sig2_f*h_k - 2*(sig_bf - delta);

    %separate multipliers for lower and upper bound constraints
    g_lo = h_min - h_k;
    g_hi = h_k   - h_max;

    grad_lo = -(lam_lo + p.c*max(g_lo, 0));
    grad_hi =  (lam_hi + p.c*max(g_hi, 0));

    h_k = h_k - tau1*(grad_obj + grad_lo + grad_hi);
    h_k = max(h_min, min(h_max, h_k));

    lam_lo = max(0, lam_lo + tau1*max(g_lo, 0));
    lam_hi = max(0, lam_hi + tau1*max(g_hi, 0));
end

h.ratio  = h_k;
h.lambda = [lam_lo; lam_hi];
end