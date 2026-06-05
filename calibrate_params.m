function [gamma_opt, kappa_opt] = calibrate_params(returns_w, F, dates_w, basis_levels, p)
%grid search on validation period 2017-2018

val_start = datetime(2017,1,1);
val_end   = datetime(2018,12,31);

gamma_grid = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5];
kappa_grid = [0.0001, 0.0005, 0.001, 0.005, 0.01];

best_sr   = -inf;
gamma_opt = 0.01;
kappa_opt = 0.001;

fprintf('Calibrating on validation period...\n');

for gi = 1:length(gamma_grid)
    for ki = 1:length(kappa_grid)
        p_tmp       = p;
        p_tmp.gamma = gamma_grid(gi);
        p_tmp.kappa = kappa_grid(ki);

        port_ret = run_validation(returns_w, F, dates_w, basis_levels, ...
            val_start, val_end, p_tmp);

        if length(port_ret) > 4
            sr = mean(port_ret)/std(port_ret)*sqrt(52);
            if sr > best_sr
                best_sr   = sr;
                gamma_opt = gamma_grid(gi);
                kappa_opt = kappa_grid(ki);
            end
        end
    end
end
fprintf('Best validation Sharpe = %.3f\n', best_sr);
end


function port_ret = run_validation(returns_w, F, dates_w, basis_levels, ...
    val_start, val_end, p)
window = p.window;
bwin   = p.basis_window;
M      = 4;

idx_s = find(dates_w >= val_start, 1, 'first');
idx_e = find(dates_w <= val_end,   1, 'last');

if isempty(idx_s)||isempty(idx_e)||idx_s<=window
    port_ret = []; return;
end

w_k      = ones(M,1)/M;
state    = struct('lambda',0,'mu',zeros(8,1));
n_val    = idx_e - idx_s + 1;
port_ret = zeros(n_val,1);

for k = 1:n_val
    t         = idx_s + k - 1;
    win_idx   = (t-window):(t-1);
    b_start   = max(1, t-bwin);
    basis_win = basis_levels(b_start:t-1,:);

    out          = layer2_portfolio(returns_w(win_idx,:), F(win_idx,:), ...
                                    basis_win, w_k, state, p);
    w_k          = out.w;
    state.lambda = out.lambda;
    state.mu     = out.mu;
    port_ret(k)  = returns_w(t,:)*w_k;
end
end
