%main.m - master script
clear; clc; close all;
addpath(genpath(pwd));

[returns_w, factors, dates_w, dv01, fut_px, basis_levels] = load_data();

%train/val/test split
train_end  = find(dates_w <= datetime(2016,12,31), 1, 'last');
val_end    = find(dates_w <= datetime(2018,12,31), 1, 'last');
test_start = find(dates_w >= datetime(2019,1,1),   1, 'first');

fprintf('Training ends:   %s (idx %d)\n', datestr(dates_w(train_end)),  train_end);
fprintf('Validation ends: %s (idx %d)\n', datestr(dates_w(val_end)),    val_end);
fprintf('Test starts:     %s (idx %d)\n', datestr(dates_w(test_start)), test_start);

%pre-test pca on training period only
[F, pca_model] = run_pca(factors, 4, train_end);

%parameters
p.window       = 104;
p.K            = 4;
p.tau          = 0.1;
p.c            = 1;
p.n_iter       = 2500;
p.lambda_T     = 0.01;
p.w_min        = -1;
p.w_max        =  1;
p.chi2_alpha   = 0.95;
p.ridge_lam    = 1e-4;
p.basis_window = 52;
p.train_end    = train_end;
p.val_end      = val_end;
p.test_start   = test_start;

[p.gamma, p.kappa] = calibrate_params(returns_w, F, dates_w, basis_levels, p);
fprintf('Calibrated: gamma=%.4f  kappa=%.6f\n', p.gamma, p.kappa);

results = backtest(returns_w, F, dates_w, dv01, basis_levels, p);

evaluate(results, dates_w, p);

%discretisation_analysis(returns_w, F, dates_w, basis_levels, p);    

%adaptive-step verification backtest: the discretisation analysis recommends
%tau_1 = 0.9*tau_c = 0.45 for redeployment. production results above retain
%tau = 0.10 throughout; this run is reported only as a robustness check in
%the discretisation section, confirming the residual LPNN-CVX gap at
%tau = 0.10 is finite-iteration solver error
p_adapt     = p;
p_adapt.tau = 0.45;
res_a = backtest(returns_w, F, dates_w, dv01, basis_levels, p_adapt);
r_a   = res_a.ret_lpnn;
r_p   = results.ret_lpnn;
fprintf('\n========== ADAPTIVE-STEP VERIFICATION (tau = 0.45 = 0.9*tau_c) ==========\n');
fprintf('%-32s %12s %12s\n', '', 'tau = 0.45', 'tau = 0.10');
fprintf('%-32s %12.4f %12.4f\n', 'LPNN Sharpe', ...
    mean(r_a)*52/(std(r_a)*sqrt(52)), mean(r_p)*52/(std(r_p)*sqrt(52)));
fprintf('%-32s %12.6f %12.6f\n', 'Corr(LPNN, CVX)', ...
    corr(r_a, res_a.ret_rob), corr(r_p, results.ret_rob));
fprintf('%-32s %12.4f %12.4f\n', 'Avg weekly turnover', ...
    mean(sum(abs(diff(res_a.w_lpnn)),2)), mean(sum(abs(diff(results.w_lpnn)),2)));
fprintf('%-32s %12.3e %12.3e\n', 'Max Lyapunov energy', ...
    max(res_a.lyap_E), max(results.lyap_E));
fprintf('==========================================================================\n');

%iteration-count convergence table: LPNN approach to the CVX solution as
%n_iter increases (production value 2500 already computed above)
fprintf('\n========== CONVERGENCE vs ITERATION COUNT ==========\n');
fprintf('%-10s %12s %12s %16s\n', 'n_iter', 'LPNN Sharpe', 'CVX Sharpe', 'Corr(LPNN,CVX)');
for n_it = [400 800 1500]
    p_tmp        = p;
    p_tmp.n_iter = n_it;
    res_n = backtest(returns_w, F, dates_w, dv01, basis_levels, p_tmp);
    fprintf('%-10d %12.4f %12.4f %16.4f\n', n_it, ...
        mean(res_n.ret_lpnn)*52/(std(res_n.ret_lpnn)*sqrt(52)), ...
        mean(res_n.ret_rob)*52/(std(res_n.ret_rob)*sqrt(52)), ...
        corr(res_n.ret_lpnn, res_n.ret_rob));
end
fprintf('%-10d %12.4f %12.4f %16.4f   (production run above)\n', 2500, ...
    mean(r_p)*52/(std(r_p)*sqrt(52)), ...
    mean(results.ret_rob)*52/(std(results.ret_rob)*sqrt(52)), ...
    corr(r_p, results.ret_rob));
fprintf('%-10d %12.4f %12.4f %16.4f   (adaptive tau=0.45 run above)\n', 2500, ...
    mean(r_a)*52/(std(r_a)*sqrt(52)), ...
    mean(res_a.ret_rob)*52/(std(res_a.ret_rob)*sqrt(52)), ...
    corr(r_a, res_a.ret_rob));
fprintf('====================================================\n');


%transaction_cost_analysis(results);
%sensitivity_analysis(returns_w, F, dates_w, dv01, basis_levels, p);
%tikhonov_frontier(returns_w, F, dates_w, dv01, basis_levels, p);
%plot_signal_flow_graph();