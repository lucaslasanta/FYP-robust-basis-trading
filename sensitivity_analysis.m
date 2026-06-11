function sensitivity_analysis(returns_w, F, dates_w, dv01, basis_levels, p_base)
%runs named experiment configurations and prints comparison table
%purpose is diagnostic - not parameter tuning on the test set

fprintf('\n========== EXPERIMENT CONFIGURATIONS ==========\n');

configs = struct();

%baseline
configs(1).name    = 'Baseline';
configs(1).p       = p_base;

%low robustness
p_tmp = p_base; p_tmp.chi2_alpha = 0.80;
configs(2).name = 'Low Robustness (chi2=0.80)';
configs(2).p    = p_tmp;

%no tikhonov
p_tmp = p_base; p_tmp.lambda_T = 0;
configs(3).name = 'No Tikhonov';
configs(3).p    = p_tmp;

%no transaction cost
p_tmp = p_base; p_tmp.kappa = 0;
configs(4).name = 'No TC Penalty';
configs(4).p    = p_tmp;

%reduced friction: low robustness + no tikhonov + no TC
p_tmp = p_base;
p_tmp.chi2_alpha = 0.80;
p_tmp.lambda_T   = 0;
p_tmp.kappa      = 0;
configs(5).name  = 'Reduced Friction';
configs(5).p     = p_tmp;

%looser bounds
p_tmp = p_base; p_tmp.w_min = -1.5; p_tmp.w_max = 1.5;
configs(6).name = 'Looser Bounds (-1.5,1.5)';
configs(6).p    = p_tmp;

%higher gamma (baseline calibrated to 0.001, so test 10x risk aversion)
p_tmp = p_base; p_tmp.gamma = 0.01;
configs(7).name = 'High Gamma (0.01)';
configs(7).p    = p_tmp;

net_cost_bps = 0.1; %institutional treasury futures cost (large bank desks)

%print header
fprintf('(NetSharpe at %.1f bps one-way cost)\n', net_cost_bps);
fprintf('%-32s %8s %8s %10s %8s %8s %8s %8s %8s\n', ...
    'Config','AnnRet','Sharpe','NetSharpe','MaxDD','Turnover','GrossExp','DistEq','CorrEq');

for ci = 1:length(configs)
    res = backtest(returns_w, F, dates_w, dv01, basis_levels, configs(ci).p);
    r   = res.ret_lpnn;

    ann_ret = mean(r)*52;
    ann_vol = std(r)*sqrt(52);
    sharpe  = ann_ret/ann_vol;
    maxdd   = max_dd(r);
    to_vec  = [0; sum(abs(diff(res.w_lpnn)), 2)];
    to      = mean(to_vec(2:end));
    gross   = mean(sum(abs(res.w_lpnn),2));
    dist_eq = mean(sum(abs(res.w_lpnn - repmat(ones(1,4)/4, size(res.w_lpnn,1),1)),2));
    corr_eq = corr(r, res.ret_eq);
    r_net   = r - (net_cost_bps/10000)*to_vec;
    net_sr  = mean(r_net)*52 / (std(r_net)*sqrt(52));

    fprintf('%-32s %8.4f %8.4f %10.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', ...
        configs(ci).name, ann_ret, sharpe, net_sr, maxdd, to, gross, dist_eq, corr_eq);
end
fprintf('================================================\n');
end

function dd = max_dd(r)
cum = cumprod(1+r);
rm  = cummax(cum);
dd  = max((rm-cum)./rm);
end