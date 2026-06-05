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

%plot_signal_flow_graph();

%sensitivity_analysis(returns_w, F, dates_w, dv01, basis_levels, p);
%tikhonov_frontier(returns_w, F, dates_w, dv01, basis_levels, p);
%transaction_cost_analysis(results);
discretisation_analysis(returns_w, F, dates_w, basis_levels, p);