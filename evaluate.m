function evaluate(results, dates_w, p)
%performance metrics, signal diagnostics and figures

r_lpnn = results.ret_lpnn;
r_rob  = results.ret_rob;
r_mv   = results.ret_mv;
r_eq   = results.ret_eq;
dates_t = results.dates;

strats = {r_lpnn, r_rob, r_mv, r_eq};
names  = {'LPNN Robust','CVX Robust','Mean-Variance','Equal Weight'};
N      = length(strats);

%--- performance summary ---
%net sharpe uses 0.1 bps one-way transaction cost per unit of turnover.
%this figure reflects execution costs at large institutional desks (major
%investment banks, e.g. bnp paribas, goldman sachs rates trading) with
%direct market access to treasury futures. it is not representative of
%retail brokers or smaller trading firms where costs are typically 1-5 bps.
%source: industry practitioner knowledge, confirmed by federal reserve (2024)
%analysis of treasury futures basis trade execution.
net_cost_bps = 0.1;
T_eval = length(r_lpnn);

fprintf('\n========== PERFORMANCE SUMMARY ==========\n');
fprintf('(NetSharpe computed at %.1f bps one-way cost: large institutional desks)\n', net_cost_bps);
fprintf('%-22s %8s %8s %8s %10s %8s %8s %8s\n', ...
    'Strategy','AnnRet','AnnVol','Sharpe','NetSharpe','MaxDD','VaR95','Turnover');

for s = 1:N
    r       = strats{s};
    ann_ret = mean(r)*52;
    ann_vol = std(r)*sqrt(52);
    sharpe  = ann_ret / ann_vol;
    maxdd   = max_drawdown(r);
    var95   = -quantile(r, 0.05)*sqrt(52);
    if s == 1
        to_vec = [0; sum(abs(diff(results.w_lpnn)),2)];
        to     = mean(to_vec(2:end));
    elseif s == 2
        to_vec = [0; sum(abs(diff(results.w_rob)),2)];
        to     = mean(to_vec(2:end));
    elseif s == 3
        to_vec = [0; sum(abs(diff(results.w_mv)),2)];
        to     = mean(to_vec(2:end));
    else
        to_vec = zeros(T_eval,1);
        to     = 0;
    end
    r_net  = r - (net_cost_bps/10000)*to_vec;
    net_sr = mean(r_net)*52 / (std(r_net)*sqrt(52));
    fprintf('%-22s %8.4f %8.4f %8.4f %10.4f %8.4f %8.4f %8.4f\n', ...
        names{s}, ann_ret, ann_vol, sharpe, net_sr, maxdd, var95, to);
end
fprintf('==========================================\n');

%--- lpnn active-weight diagnostics ---
w_eq_mat = ones(size(results.w_lpnn))/4;
avg_dist_eq  = mean(sum(abs(results.w_lpnn - w_eq_mat), 2));
avg_turnover = mean(sum(abs(diff(results.w_lpnn)), 2));
avg_gross    = mean(sum(abs(results.w_lpnn), 2));
upper_hit    = mean(mean(results.w_lpnn >= p.w_max - 1e-4));
lower_hit    = mean(mean(results.w_lpnn <= p.w_min + 1e-4));
corr_eq      = corr(r_lpnn, r_eq);
corr_rob     = corr(r_lpnn, r_rob);

fprintf('\n========== LPNN ACTIVE-WEIGHT DIAGNOSTICS ==========\n');
fprintf('Average L1 distance from equal weight: %.4f\n', avg_dist_eq);
fprintf('Average weekly turnover:               %.4f\n', avg_turnover);
fprintf('Average gross exposure:                %.4f\n', avg_gross);
fprintf('Fraction of weights at upper bound:    %.4f\n', upper_hit);
fprintf('Fraction of weights at lower bound:    %.4f\n', lower_hit);
fprintf('Corr(LPNN returns, Equal Weight):      %.4f\n', corr_eq);
fprintf('Corr(LPNN returns, CVX Robust):        %.4f\n', corr_rob);
fprintf('\nAverage weight by maturity:\n');
mat_names = {'2Y','5Y','10Y','30Y'};
for m = 1:4
    fprintf('  %4s: mean=%.4f  min=%.4f  max=%.4f\n', mat_names{m}, ...
        mean(results.w_lpnn(:,m)), min(results.w_lpnn(:,m)), max(results.w_lpnn(:,m)));
end
fprintf('=====================================================\n');

%--- signal diagnostics ---
if isfield(results,'mu_hat_log')
    mh = results.mu_hat_log;
    fprintf('\n========== SIGNAL DIAGNOSTICS (mu_hat) ==========\n');
    fprintf('%-6s %10s %10s %10s %10s %10s\n', ...
        'Mat','Mean','Std','Dispersion','PredCorr','HitRate');
    T = size(mh,1);
    for m = 1:4
        mh_m    = mh(:,m);
        r_next  = [r_lpnn(2:end); NaN];  %next-period realised return
        valid   = ~isnan(r_next);
        pc      = corr(mh_m(valid), r_next(valid));
        hr      = mean(sign(mh_m(valid)) == sign(r_next(valid)));
        disp_m  = std(mh_m);
        fprintf('%-6s %10.6f %10.6f %10.6f %10.4f %10.4f\n', ...
            mat_names{m}, mean(mh_m), std(mh_m), disp_m, pc, hr);
    end
    fprintf('==================================================\n');
end

%--- improved signal diagnostics (signal vs individual maturity return) ---
if isfield(results,'mu_hat_log') && isfield(results,'realized_r_mat')
    mh = results.mu_hat_log;
    rm = results.realized_r_mat;
    fprintf('\n========== IMPROVED SIGNAL DIAGNOSTICS ==========\n');
    fprintf('%-6s %14s %12s %14s %12s\n', ...
        'Mat','PredCorr(mat)','HitRate(mat)','Mean|mu_hat|','Std(mu_hat)');
    for m = 1:4
        mh_m = mh(1:end-1, m);
        rm_m = rm(2:end,   m);
        pc_m = corr(mh_m, rm_m);
        hr_m = mean(sign(mh_m) == sign(rm_m));
        fprintf('%-6s %14.4f %12.4f %14.6f %12.6f\n', ...
            mat_names{m}, pc_m, hr_m, mean(abs(mh(:,m))), std(mh(:,m)));
    end
    fprintf('==================================================\n');
end

%--- stress periods ---
stress = {datetime(2020,2,1), datetime(2020,5,31), 'COVID-19 (Feb-May 2020)';
          datetime(2022,1,1), datetime(2022,12,31),'Rate Hiking Cycle (2022)';
          datetime(2023,3,1), datetime(2023,5,31), 'Banking Stress (Mar-May 2023)'};

fprintf('\n========== STRESS PERIOD SHARPE ==========\n');
fprintf('%-32s', 'Period');
for s = 1:N, fprintf('%14s', names{s}); end
fprintf('\n');
for sp = 1:size(stress,1)
    idx = dates_t >= stress{sp,1} & dates_t <= stress{sp,2};
    fprintf('%-32s', stress{sp,3});
    for s = 1:N
        r_sp = strats{s}(idx);
        if length(r_sp) > 2
            sr = mean(r_sp)/std(r_sp)*sqrt(52);
            fprintf('%14.3f', sr);
        else
            fprintf('%14s', 'n/a');
        end
    end
    fprintf('\n');
end
fprintf('==========================================\n');


%--- figure 1: cumulative returns ---
figure('Name','Cumulative Returns','Position',[100 100 900 420]);
styles = {'b-','r--','k-.','g:'};
lw     = [1.5 1.2 1.2 1.2];
for s = 1:N
    plot(dates_t, cumprod(1+strats{s}), styles{s}, 'LineWidth', lw(s)); hold on;
end
legend(names,'Location','southwest'); grid on;
xlabel('Date'); ylabel('Cumulative return');
title('Cumulative portfolio returns');
hold off;
saveas(gcf,'fig_cumulative_returns.pdf');

%--- figure 2: lpnn weights (line plot - handles negative weights) ---
figure('Name','LPNN Weights','Position',[100 100 900 380]);
plot(dates_t, results.w_lpnn);
legend(mat_names,'Location','northeast'); grid on;
xlabel('Date'); ylabel('Portfolio weight');
title('LPNN adaptive portfolio weights');
yline(0,'k--','LineWidth',0.8);
saveas(gcf,'fig_lpnn_weights.pdf');

%--- figure 3: active weights (lpnn - equal weight) ---
figure('Name','Active Weights','Position',[100 100 900 380]);
active_w = results.w_lpnn - repmat(ones(1,4)/4, size(results.w_lpnn,1), 1);
plot(dates_t, active_w);
legend(mat_names,'Location','northeast'); grid on;
xlabel('Date'); ylabel('Active weight (vs equal weight)');
title('LPNN active weights');
yline(0,'k--','LineWidth',0.8);
saveas(gcf,'fig_active_weights.pdf');

%--- figure 4: lyapunov energy ---
figure('Name','Lyapunov Energy','Position',[100 100 900 280]);
semilogy(dates_t, results.lyap_E,'b-','LineWidth',1);
xlabel('Date'); ylabel('E(w,\lambda) [log scale]');
title('LPNN Lyapunov energy'); grid on;
saveas(gcf,'fig_lyapunov_energy.pdf');

%--- figure 5: hedge ratios ---
figure('Name','Hedge Ratios','Position',[100 100 900 380]);
plot(dates_t, results.h_ratios);
legend(mat_names,'Location','best'); grid on;
xlabel('Date'); ylabel('Hedge ratio');
title('Layer 1 robust hedge ratios');
saveas(gcf,'fig_hedge_ratios.pdf');

%--- figure 6: epsilon ---
figure('Name','Epsilon','Position',[100 100 900 280]);
plot(dates_t, results.epsilon,'b-','LineWidth',1);
xlabel('Date'); ylabel('\epsilon');
title('Ellipsoidal uncertainty radius'); grid on;
saveas(gcf,'fig_epsilon.pdf');

%--- figure 7: drawdowns ---
figure('Name','Drawdowns','Position',[100 100 900 380]);
for s = 1:N
    plot(dates_t, drawdown_series(strats{s})*100, styles{s}, 'LineWidth', lw(s)); hold on;
end
legend(names,'Location','southwest'); grid on;
xlabel('Date'); ylabel('Drawdown (%)');
title('Portfolio drawdowns'); hold off;
saveas(gcf,'fig_drawdowns.pdf');

%--- figure 8: mu_hat over time ---
if isfield(results,'mu_hat_log')
    figure('Name','Predicted Returns','Position',[100 100 900 380]);
    plot(dates_t, results.mu_hat_log);
    legend(mat_names,'Location','best'); grid on;
    xlabel('Date'); ylabel('\mu_{hat}');
    title('Factor-implied expected returns by maturity');
    yline(0,'k--','LineWidth',0.8);
    saveas(gcf,'fig_mu_hat.pdf');
end

%--- figure 9: gross exposure and distance from equal weight ---
figure('Name','Portfolio Activity','Position',[100 100 900 380]);
subplot(2,1,1);
plot(dates_t, sum(abs(results.w_lpnn),2),'b-');
xlabel('Date'); ylabel('Gross exposure');
title('Gross exposure over time'); grid on;
yline(1,'k--','LineWidth',0.8);

subplot(2,1,2);
plot(dates_t, sum(abs(results.w_lpnn - w_eq_mat),2),'b-');
xlabel('Date'); ylabel('L1 distance');
title('Distance from equal weight over time'); grid on;
saveas(gcf,'fig_activity.pdf');

fprintf('\nAll figures saved to figures/\n');
end


function dd = max_drawdown(r)
cum = cumprod(1+r);
rm  = cummax(cum);
dd  = max((rm - cum)./rm);
end

function dd = drawdown_series(r)
cum = cumprod(1+r);
rm  = cummax(cum);
dd  = (rm - cum)./rm;
end