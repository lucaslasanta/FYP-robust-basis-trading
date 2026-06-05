function results = backtest(returns_w, F, dates_w, dv01, basis_levels, p)
%rolling backtest over test period

test_start = p.test_start;
window     = p.window;
bwin       = p.basis_window;
M          = 4;
T_test     = length(dates_w) - test_start + 1;

fprintf('Running backtest: %s to %s (%d weeks)\n', ...
    datestr(dates_w(test_start)), datestr(dates_w(end)), T_test);

dates_test  = dates_w(test_start:end);
w_lpnn          = zeros(T_test,M);
w_rob           = zeros(T_test,M);
w_mv            = zeros(T_test,M);
w_eq            = repmat(ones(1,M)/M, T_test,1);
ret_lpnn        = zeros(T_test,1);
ret_rob         = zeros(T_test,1);
ret_mv          = zeros(T_test,1);
ret_eq          = zeros(T_test,1);
lyap_E          = zeros(T_test,1);
h_ratios        = zeros(T_test,M);
eps_track       = zeros(T_test,1);
mu_hat_log      = zeros(T_test,M);
sig_source      = cell(T_test,1);
realized_r_mat  = zeros(T_test,M);

w_k    = ones(M,1)/M;
state  = struct('lambda',0,'mu',zeros(8,1));
for m = 1:M
    h_state(m).ratio  = 0.0;   %start at zero: no level hedge initially
    h_state(m).lambda = [0;0]; %separate lo/hi multipliers
end
w_prev_rob = ones(M,1)/M;
w_prev_mv  = ones(M,1)/M;

for i = 1:T_test
    t       = test_start + i - 1;
    win_idx = (t-window):(t-1);
    ret_win = returns_w(win_idx,:);
    F_win   = F(win_idx,:);

    %basis window: history up to t-1
    b_start   = max(1, t-bwin);
    basis_win = basis_levels(b_start:t-1,:);

    %layer 1: robust hedge of each basis return against the level factor
    %residual removes the common yield-curve level component, making
    %layer-2 inputs more idiosyncratic and better diversified
    f_level      = F_win(:,1);
    r_hedged_win = ret_win;
    for m = 1:M
        h_out = layer1_hedge(ret_win(:,m), f_level, ...
                             h_state(m).ratio, h_state(m).lambda, p);
        h_state(m).ratio  = h_out.ratio;
        h_state(m).lambda = h_out.lambda;
        r_hedged_win(:,m) = ret_win(:,m) - h_out.ratio * f_level;
    end

    %layer 2: lpnn on hedged residual returns
    out  = layer2_portfolio(r_hedged_win, F_win, basis_win, w_k, state, p);
    w_k  = out.w;
    state.lambda = out.lambda;
    state.mu     = out.mu;

    %benchmarks also use hedged returns for fair comparison
    bm = benchmarks(r_hedged_win, F_win, basis_win, w_prev_rob, p);

    w_lpnn(i,:)     = w_k';
    w_rob(i,:)      = bm.w_rob';
    w_mv(i,:)       = bm.w_mv';
    h_ratios(i,:)   = arrayfun(@(s) s.ratio, h_state);
    lyap_E(i)       = out.lyap_energy;
    eps_track(i)    = out.epsilon;
    mu_hat_log(i,:) = out.mu_hat';
    sig_source{i}   = out.signal_source;

    %realized hedged return: basis return minus level-factor hedge leg
    h_vec           = arrayfun(@(s) s.ratio, h_state)';
    r_t_raw         = returns_w(t,:)';
    r_t             = r_t_raw - h_vec * F(t,1);
    realized_r_mat(i,:) = r_t';
    ret_lpnn(i)     = w_k'*r_t;
    ret_rob(i)   = bm.w_rob'*r_t;
    ret_mv(i)    = bm.w_mv'*r_t;
    ret_eq(i)    = mean(r_t);

    w_prev_rob = bm.w_rob;
    w_prev_mv  = bm.w_mv;

    if mod(i,52)==0
        n_basis = sum(cellfun(@(s) contains(s,'basis'), sig_source(1:i)));
        fprintf('  Year %d complete (basis signal active: %d/%d weeks)\n', ...
            round(i/52), n_basis, i);
    end
end

results.dates       = dates_test;
results.ret_lpnn    = ret_lpnn;
results.ret_rob     = ret_rob;
results.ret_mv      = ret_mv;
results.ret_eq      = ret_eq;
results.w_lpnn      = w_lpnn;
results.w_rob       = w_rob;
results.w_mv        = w_mv;
results.w_eq        = w_eq;
results.h_ratios    = h_ratios;
results.lyap_E      = lyap_E;
results.epsilon     = eps_track;
results.mu_hat_log  = mu_hat_log;
results.sig_source      = sig_source;
results.realized_r_mat  = realized_r_mat;
end