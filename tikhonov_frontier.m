function tikhonov_frontier(returns_w, F, dates_w, dv01, basis_levels, p_base)
%sweeps lambda_T and reports performance/risk metrics for each value
%shows how tikhonov maturity-smoothness regularisation affects the portfolio

lambda_grid = [0, 0.0001, 0.0003, 0.001, 0.003, 0.005, 0.01, 0.02, 0.05, 0.1];
n           = length(lambda_grid);

net_cost_bps = 0.1; %institutional treasury futures cost (large bank desks)

fprintf('\n========== TIKHONOV FRONTIER ==========\n');
fprintf('(NetSharpe at %.1f bps one-way cost)\n', net_cost_bps);
fprintf('%-10s %8s %8s %8s %10s %8s %8s %8s %8s %8s %8s\n', ...
    'lambda_T','AnnRet','AnnVol','Sharpe','NetSharpe','MaxDD','VaR95','Turnover','GrossExp','DistEq','CorrEq');

metrics = zeros(n, 10);

for i = 1:n
    p_tmp          = p_base;
    p_tmp.lambda_T = lambda_grid(i);
    res = backtest(returns_w, F, dates_w, dv01, basis_levels, p_tmp);
    r   = res.ret_lpnn;

    ann_ret = mean(r)*52;
    ann_vol = std(r)*sqrt(52);
    sharpe  = ann_ret / ann_vol;
    cum     = cumprod(1+r);
    maxdd   = max((cummax(cum) - cum) ./ cummax(cum));
    var95   = -quantile(r, 0.05)*sqrt(52);
    to_vec  = [0; sum(abs(diff(res.w_lpnn)), 2)];
    to      = mean(to_vec(2:end));
    gross   = mean(sum(abs(res.w_lpnn), 2));
    dist_eq = mean(sum(abs(res.w_lpnn - 0.25), 2));
    corr_eq = corr(r, res.ret_eq);
    r_net   = r - (net_cost_bps/10000)*to_vec;
    net_sr  = mean(r_net)*52 / (std(r_net)*sqrt(52));

    metrics(i,:) = [ann_ret, ann_vol, sharpe, maxdd, var95, to, gross, dist_eq, corr_eq, net_sr];

    fprintf('%-10.4f %8.4f %8.4f %8.4f %10.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', ...
        lambda_grid(i), ann_ret, ann_vol, sharpe, net_sr, maxdd, var95, to, gross, dist_eq, corr_eq);
end
fprintf('========================================\n');

%for log-scale plot replace lambda=0 with small sentinel
lam_plot = lambda_grid;
lam_plot(lam_plot == 0) = 1e-5;
baseline = 0.01;

figure('Name','Tikhonov Frontier','Position',[100 100 900 600]);

subplot(2,2,1);
semilogx(lam_plot, metrics(:,3), 'b-o', 'LineWidth', 1.2); grid on;
xline(baseline, '--r', 'Baseline');
xlabel('lambda_T'); ylabel('Sharpe ratio');
title('Sharpe vs Tikhonov \lambda_T');

subplot(2,2,2);
semilogx(lam_plot, metrics(:,7), 'b-o', 'LineWidth', 1.2); grid on;
xline(baseline, '--r', 'Baseline');
xlabel('lambda_T'); ylabel('Gross exposure');
title('Gross Exposure vs \lambda_T');

subplot(2,2,3);
semilogx(lam_plot, metrics(:,6), 'b-o', 'LineWidth', 1.2); grid on;
xline(baseline, '--r', 'Baseline');
xlabel('lambda_T'); ylabel('Turnover');
title('Turnover vs \lambda_T');

subplot(2,2,4);
semilogx(lam_plot, metrics(:,4), 'b-o', 'LineWidth', 1.2); grid on;
xline(baseline, '--r', 'Baseline');
xlabel('lambda_T'); ylabel('Max drawdown');
title('MaxDD vs \lambda_T');

saveas(gcf, 'fig_tikhonov_frontier.pdf');
end
