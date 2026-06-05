function transaction_cost_analysis(results)
%ex-post transaction cost stress test using realistic costs for treasury futures
%
%cost basis (one-way, bps per unit of portfolio turnover):
%  0.10 bps: large institutional desks (major banks) with direct market access
%  0.20 bps: typical best-execution for institutional treasury futures
%  0.30 bps: slightly larger trades or less liquid sessions
%  0.50 bps: conservative institutional estimate
%  0.75 bps: mid-tier institutional / intraday liquidity constraints
%  1.00 bps: upper bound for institutional; entry-level for smaller firms
%
%source: rates trading desks at major banks;
%confirmed by cme group treasury futures market data and federal reserve (2024).
%note: retail brokers and smaller trading firms typically face 5-20 bps,
%which is outside the scope of this institutional strategy.

cost_bps = [0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0];
T        = length(results.ret_lpnn);

%per-week turnover vectors (T x 1); week 1 assumed zero
to_lpnn = [0; sum(abs(diff(results.w_lpnn)), 2)];
to_rob  = [0; sum(abs(diff(results.w_rob)),  2)];

strat_names = {'LPNN Robust', 'CVX Robust'};
strat_rets  = {results.ret_lpnn, results.ret_rob};
strat_to    = {to_lpnn, to_rob};

net_sharpe = zeros(length(cost_bps), 2);
net_annret = zeros(length(cost_bps), 2);
net_maxdd  = zeros(length(cost_bps), 2);

fprintf('\n========== TRANSACTION COST STRESS TEST ==========\n');
fprintf('(Realistic range for institutional treasury futures: 0.1-1.0 bps one-way)\n');

for s = 1:2
    fprintf('\n%s:\n', strat_names{s});
    fprintf('%-14s %12s %12s %12s\n', 'Cost (bps)', 'Net AnnRet', 'Net Sharpe', 'Net MaxDD');
    r_gross = strat_rets{s};
    to_vec  = strat_to{s};

    for j = 1:length(cost_bps)
        cost_t  = (cost_bps(j) / 10000) * to_vec;
        r_net   = r_gross - cost_t;
        ann_ret = mean(r_net)*52;
        ann_vol = std(r_net)*sqrt(52);
        sharpe  = ann_ret / ann_vol;
        cum     = cumprod(1+r_net);
        maxdd   = max((cummax(cum) - cum) ./ cummax(cum));

        net_sharpe(j,s) = sharpe;
        net_annret(j,s) = ann_ret;
        net_maxdd(j,s)  = maxdd;

        fprintf('%-14s %12.4f %12.4f %12.4f\n', ...
            sprintf('%.1f bps', cost_bps(j)), ann_ret, sharpe, maxdd);
    end
end
fprintf('===================================================\n');

figure('Name','Transaction Cost Analysis','Position',[100 100 900 380]);

subplot(1,2,1);
plot(cost_bps, net_sharpe(:,1), 'b-o', 'LineWidth', 1.2); hold on;
plot(cost_bps, net_sharpe(:,2), 'r--s', 'LineWidth', 1.2);
xline(0.1, ':k', '0.1 bp'); xline(0.5, ':k', '0.5 bps');
legend({'LPNN Robust','CVX Robust'}, 'Location', 'northeast');
xlabel('Transaction cost (bps)'); ylabel('Net Sharpe ratio');
title('Net Sharpe vs Transaction Cost'); grid on; hold off;

subplot(1,2,2);
plot(cost_bps, net_annret(:,1)*100, 'b-o', 'LineWidth', 1.2); hold on;
plot(cost_bps, net_annret(:,2)*100, 'r--s', 'LineWidth', 1.2);
xline(0.1, ':k', '0.1 bp'); xline(0.5, ':k', '0.5 bps');
legend({'LPNN Robust','CVX Robust'}, 'Location', 'northeast');
xlabel('Transaction cost (bps)'); ylabel('Net Annual Return (%)');
title('Net Annual Return vs Transaction Cost'); grid on; hold off;

saveas(gcf, 'fig_transaction_costs.pdf');
end
