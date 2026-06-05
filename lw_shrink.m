function S = lw_shrink(X, S_sample)
%ledoit-wolf analytical shrinkage towards scaled identity
%x: (T x N) data matrix, S_sample: sample covariance
[T, N] = size(X);

%target: mu * I where mu = trace(S)/N
mu   = trace(S_sample) / N;
F_tgt = mu * eye(N);

%optimal shrinkage intensity (ledoit-wolf 2004 analytical formula)
X_c  = X - mean(X);
sum_asym = 0;
for t = 1:T
    x_t   = X_c(t,:)';
    sum_asym = sum_asym + norm(x_t*x_t' - S_sample, 'fro')^2;
end
delta2 = sum_asym / (T^2);

gamma2 = norm(S_sample - F_tgt, 'fro')^2;
if gamma2 < 1e-14
    alpha = 0;
else
    alpha = min(1, delta2 / gamma2);
end

S = (1 - alpha) * S_sample + alpha * F_tgt;
end
