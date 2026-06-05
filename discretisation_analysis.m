function discretisation_analysis(returns_w, F, dates_w, basis_levels, p)
%DISCRETISATION_ANALYSIS
%Compares four discretisation schemes for the Layer-2 LPNN primal-dual ODE:
%
%   dw/dt   = -nabla_w L_aug(w, lambda, mu)   [primal]
%   dlam/dt =  h_bud(w) = sum(w) - 1          [budget multiplier]
%   dmu/dt  =  g(w) = [g_lo.*act_lo;          [bound multipliers]
%                       g_hi.*act_hi]
%
%Four methods with DISTINCT step sizes so all curves are visible:
%  1. Fixed Euler       tau = p.tau = 0.1       (production baseline)
%  2. RK4               tau = 0.25              (larger step; wider RK4 stability region)
%  3. Adaptive Euler    tau_in = 0.6, capped to 0.9*tau_c (demonstrates adaptation)
%  4. Euler above bound tau = 0.6 > tau_c       (smooth-Hessian bound violated; converges
%                                                due to non-smooth implicit stabilisation)
%
%Layer-1 hedge ratios computed fresh for each date (deterministic given window).
%CVX robust solution used as reference w* for ||w_k - w*|| tracking.
%
%Outputs: printed convergence table + fig_discretisation.pdf

fprintf('\n========== DISCRETISATION ANALYSIS ==========\n');

M = 4;  K = p.K;

labels     = {'Normal (26 Sep 2019)', 'Stress (27 Mar 2020)'};
dates_str  = {'26-Sep-2019',          '27-Mar-2020'};

fig = figure('Name','Discretisation Analysis',...
    'Position',[50 50 1200 530],'Color','white');

for di = 1:2

    %--- locate date and extract rolling window ---
    tgt   = datetime(dates_str{di}, 'InputFormat','dd-MMM-yyyy');
    t     = find(dates_w >= tgt, 1, 'first');
    assert(~isempty(t) && t > p.window, 'Date not found or insufficient history.');

    win_idx   = (t - p.window):(t-1);
    ret_win   = returns_w(win_idx, :);
    F_win     = F(win_idx, :);
    b_start   = max(1, t - p.basis_window);
    basis_win = basis_levels(b_start:t-1, :);

    fprintf('\n--- %s  (index %d, %s) ---\n', labels{di}, t, datestr(dates_w(t)));

    %--- Layer 1: run hedge fresh from zero (Option A) ---
    f_level  = F_win(:,1);
    r_hedged = ret_win;
    for m = 1:M
        h_out = layer1_hedge(ret_win(:,m), f_level, 0.0, [0;0], p);
        r_hedged(:,m) = ret_win(:,m) - h_out.ratio * f_level;
    end

    %--- compute optimisation inputs ONCE (mirrors layer2_portfolio.m exactly) ---
    Sigma = cov(r_hedged);
    Sigma = lw_shrink(r_hedged, Sigma);

    B = ((F_win'*F_win) + p.ridge_lam*eye(K)) \ (F_win'*r_hedged);
    B = B';

    n_smooth = min(4, size(F_win,1));
    f_smooth = mean(F_win(end-n_smooth+1:end,:))';
    mu_hat   = B * f_smooth;                     %initialise with factor signal

    basis_now = basis_win(end,:);
    for m = 1:M
        if ~isnan(basis_now(m)) && size(basis_win,1) >= 8
            bvals = basis_win(:,m);
            bvals = bvals(~isnan(bvals));
            if length(bvals) >= 8
                mu_b  = mean(bvals);
                sig_b = max(std(bvals), 1e-8);
                z_b   = max(-2, min(2, (basis_now(m)-mu_b)/sig_b));
                mu_hat(m) = z_b * std(r_hedged(:,m));
            end
        end
    end

    T_w     = size(r_hedged,1);
    Sigma_f = cov(F_win);
    eps2    = chi2inv(p.chi2_alpha, K) * trace(Sigma_f) / T_w;
    epsilon = sqrt(eps2);
    Phi     = Sigma_f + 1e-6*eye(K);

    L_tikh      = 2*diag([1 2 2 1]) - diag([1 1 1],1) - diag([1 1 1],-1);
    L_tikh(1,1) = 1;  L_tikh(4,4) = 1;

    %--- smooth Hessian and stability analysis ---
    %  nabla^2_w L_aug (smooth part) = 2*gamma*Sigma + 2*lambda_T*L_T + c*ones(M,M)
    %  The augmented Lagrangian budget penalty (c/2)*(sum(w)-1)^2 adds c*ones(M,M)
    %  to the Hessian. lambda_max(c*ones(M,M)) = c*M (dominant for c=1, M=4).
    %  Euler stability condition: tau < tau_c = 2 / lambda_max(H_smooth).
    H_smooth = 2*p.gamma*Sigma + 2*p.lambda_T*L_tikh + p.c*ones(M,M);
    lam_max  = max(real(eig(H_smooth)));
    lam_max  = max(lam_max, 1e-8);
    tau_c    = 2 / lam_max;

    %step sizes (four distinct values so all curves are visually separate)
    tau_fe    = p.tau;                         % 0.10 -- production baseline
    tau_rk4   = 0.25;                          % 0.25 -- larger step (within Euler margin)
    tau_in    = 0.60;                          % 0.60 -- input tau for adaptive demo
    tau_ae    = min(tau_in, 0.9 * tau_c);      %       capped at 90% of tau_c
    tau_above = tau_in;                        % 0.60 -- same input, no cap applied

    fprintf('  lambda_max(H_smooth)  = %.6f\n', lam_max);
    fprintf('  tau_c = 2/lambda_max  = %.4f\n', tau_c);
    fprintf('  [1] Fixed Euler        tau = %.2f  (%.1f%% of tau_c)\n', tau_fe,    100*tau_fe/tau_c);
    fprintf('  [2] RK4                tau = %.2f  (%.1f%% of tau_c)\n', tau_rk4,   100*tau_rk4/tau_c);
    fprintf('  [3] Adaptive Euler     tau_in = %.1f  ->  tau_1 = min(%.1f, 0.9x%.4f) = %.4f  (%.1f%% of tau_c)\n', ...
        tau_in, tau_in, tau_c, tau_ae, 100*tau_ae/tau_c);
    fprintf('  [4] Euler above bound  tau = %.2f  (%.1f%% of tau_c > 100%%)\n', tau_above, 100*tau_above/tau_c);

    %--- CVX global optimum as reference w* ---
    bm     = benchmarks(r_hedged, F_win, basis_win, ones(M,1)/M, p);
    w_star = bm.w_rob;
    fprintf('  CVX w* = [%.4f  %.4f  %.4f  %.4f]\n', w_star);

    %--- pack optimization context (precompute BPBt to avoid repeating M^2*K ops) ---
    ctx.mu_hat  = mu_hat;
    ctx.Sigma   = Sigma;
    ctx.B       = B;
    ctx.BPBt    = B * Phi * B';
    ctx.epsilon = epsilon;
    ctx.L_tikh  = L_tikh;
    ctx.p       = p;
    ctx.w_prev  = ones(M,1)/M;   %TC reference = starting point
    ctx.M       = M;

    %--- initial conditions ---
    w0   = ones(M,1)/M;
    lam0 = 0;
    mu0  = zeros(2*M,1);

    %--- solver specifications ---
    method_names = {'Fixed Euler',          'RK4 (4 evals/step)', ...
                    'Adaptive Euler',        'Euler above smooth bound'};
    method_taus  = [tau_fe, tau_rk4, tau_ae, tau_above];
    method_rk4   = [false,  true,    false,  false];
    n_iter       = p.n_iter;

    all_E    = zeros(n_iter, 4);
    all_dist = zeros(n_iter, 4);
    all_wfin = zeros(M, 4);
    all_tms  = zeros(1, 4);
    all_it3  = zeros(1, 4);
    all_it5  = zeros(1, 4);

    for mi = 1:4
        tau_m = method_taus(mi);
        rk4_m = method_rk4(mi);
        w_k   = w0;  lam_k = lam0;  mu_k = mu0;
        it3   = Inf; it5   = Inf;
        t0    = tic;

        for k = 1:n_iter
            if rk4_m
                %=== RK4: four gradient evaluations per step ===
                %Stage 1
                [dw1, dl1, dm1] = ode_rhs(w_k, lam_k, mu_k, ctx);

                %Stage 2: half step using stage-1 slope
                w2 = max(p.w_min, min(p.w_max, w_k   + (tau_m/2)*dw1));
                l2 = lam_k + (tau_m/2)*dl1;
                m2 = mu_k  + (tau_m/2)*dm1;
                [dw2, dl2, dm2] = ode_rhs(w2, l2, m2, ctx);

                %Stage 3: half step using stage-2 slope
                w3 = max(p.w_min, min(p.w_max, w_k   + (tau_m/2)*dw2));
                l3 = lam_k + (tau_m/2)*dl2;
                m3 = mu_k  + (tau_m/2)*dm2;
                [dw3, dl3, dm3] = ode_rhs(w3, l3, m3, ctx);

                %Stage 4: full step using stage-3 slope
                w4 = max(p.w_min, min(p.w_max, w_k   + tau_m*dw3));
                l4 = lam_k + tau_m*dl3;
                m4 = mu_k  + tau_m*dm3;
                [dw4, dl4, dm4] = ode_rhs(w4, l4, m4, ctx);

                %Weighted RK4 update
                w_k   = w_k   + (tau_m/6)*(dw1 + 2*dw2 + 2*dw3 + dw4);
                w_k   = max(p.w_min, min(p.w_max, w_k));   %project at final step only
                lam_k = lam_k + (tau_m/6)*(dl1 + 2*dl2 + 2*dl3 + dl4);
                mu_k  = mu_k  + (tau_m/6)*(dm1 + 2*dm2 + 2*dm3 + dm4);
                mu_k  = max(0, mu_k);   %dual variables non-negative

            else
                %=== Euler: one gradient evaluation per step ===
                %Compute full gradient at CURRENT (w_k, lam_k, mu_k)
                [grad_w, h_bud, g_lo, g_hi, act_lo, act_hi] = ...
                    compute_grad(w_k, lam_k, mu_k, ctx);

                %Primal update + projection
                w_k = w_k - tau_m * grad_w;
                w_k = max(p.w_min, min(p.w_max, w_k));

                %Dual updates use h_bud, g_lo, g_hi from BEFORE primal update
                %(matches layer2_portfolio.m lines 84-86 exactly)
                lam_k         = lam_k + tau_m * h_bud;
                mu_k(1:M)     = mu_k(1:M)     + tau_m*(g_lo.*act_lo);
                mu_k(M+1:2*M) = mu_k(M+1:2*M) + tau_m*(g_hi.*act_hi);
            end

            %Lyapunov energy at updated state (w_{k+1}, lam_{k+1}, mu_{k+1})
            Ek = lyap_energy(w_k, lam_k, mu_k, ctx);
            all_E(k, mi)    = Ek;
            all_dist(k, mi) = norm(w_k - w_star);

            if isinf(it3) && Ek < 1e-3, it3 = k; end
            if isinf(it5) && Ek < 1e-5, it5 = k; end
        end

        all_wfin(:,mi) = w_k;
        all_tms(mi)    = toc(t0)*1000;
        all_it3(mi)    = it3;
        all_it5(mi)    = it5;
    end

    %--- convergence table ---
    fprintf('\n  Convergence table  (N_iter=%d):\n', n_iter);
    fprintf('  %-30s %6s %14s %14s %12s %11s %8s\n', ...
        'Method','tau','Iters(E<1e-3)','Iters(E<1e-5)','Final E','||w-w*||','ms');
    fprintf('  %s\n', repmat('-',1,99));
    for mi = 1:4
        fprintf('  %-30s %6.4f %14s %14s %12.3e %11.6f %8.1f\n', ...
            method_names{mi}, method_taus(mi), ...
            fmt_n(all_it3(mi),n_iter), fmt_n(all_it5(mi),n_iter), ...
            all_E(end,mi), all_dist(end,mi), all_tms(mi));
    end
    fprintf('  Note: Adaptive Euler tau_in=%.1f capped to tau_1=%.4f (90%% of tau_c=%.4f)\n',...
        tau_in, tau_ae, tau_c);
    fprintf('  Note: Euler above bound tau=%.2f lies %.1f%% above tau_c; converges due to non-smooth stabilisation\n',...
        tau_above, 100*(tau_above/tau_c - 1));
    fprintf('  RK4 efficiency:  %s total grad evals to E<1e-5  vs  Euler (fixed): %s\n',...
        fmt_rk4(all_it5(2),n_iter), fmt_n(all_it5(1),n_iter));

    %--- subplot: Lyapunov energy vs iteration ---
    col_fe    = [0.85 0.15 0.10];   %red
    col_rk4   = [0.10 0.60 0.20];   %green
    col_ae    = [0.10 0.40 0.80];   %blue
    col_above = [0.50 0.50 0.50];   %grey

    iters = (1:n_iter)';

    subplot(1, 2, di);
    %draw in reverse convergence order so faster methods appear on top
    h4 = semilogy(iters, all_E(:,4), ':',  'Color',col_above, 'LineWidth',2.2); hold on;
    h3 = semilogy(iters, all_E(:,3), '--', 'Color',col_ae,    'LineWidth',2.2);
    h2 = semilogy(iters, all_E(:,2), '-.', 'Color',col_rk4,   'LineWidth',2.2);
    h1 = semilogy(iters, all_E(:,1), '-',  'Color',col_fe,    'LineWidth',2.2);

    yline(1e-3,'--k','LineWidth',0.9);
    yline(1e-5,'-.k','LineWidth',0.9);
    text(30, 1e-3*2.5,'$10^{-3}$','Interpreter','latex','FontSize',8.5);
    text(30, 1e-5*2.5,'$10^{-5}$','Interpreter','latex','FontSize',8.5);

    xlabel('Iteration $k$','Interpreter','latex','FontSize',11);
    ylabel('Lyapunov Energy $E_k$','Interpreter','latex','FontSize',11);
    title(labels{di},'Interpreter','latex','FontSize',11,'FontWeight','bold');

    lstr = {sprintf('Fixed Euler ($\\tau=%.2f$, %.0f\\%% of $\\tau_c$)',      tau_fe,    100*tau_fe/tau_c),...
            sprintf('RK4 ($\\tau=%.2f$, $4{\\times}$ grad evals/step)',        tau_rk4),...
            sprintf('Adaptive Euler ($\\tau_{\\rm in}=%.1f \\to \\tau_1=%.3f$, 90\\%% of $\\tau_c$)', tau_in, tau_ae),...
            sprintf('Euler above smooth bound ($\\tau=%.1f > \\tau_c=%.3f$)',  tau_above, tau_c)};
    legend([h1 h2 h3 h4], lstr, 'Interpreter','latex','FontSize',8,'Location','southwest');
    grid on;  hold off;

end  %di

sgtitle({'Discretisation Analysis: Layer-2 LPNN Inner Optimisation',...
    'ODE: $\dot{w}=-\nabla_w\mathcal{L}_{\mathrm{aug}}$, $\;\dot{\lambda}=h_{\mathrm{bud}}(w)$, $\;\dot{\mu}=g(w)$'},...
    'Interpreter','latex','FontSize',11);

saveas(fig,'fig_discretisation.pdf');
fprintf('\nSaved: fig_discretisation.pdf\n');
fprintf('==============================================\n');
end


%% ================================================================
%  GRADIENT OF L_aug w.r.t. w  --  matches layer2_portfolio.m exactly
%  Returns:  grad_w   = full gradient (used for Euler primal update)
%            h_bud    = budget constraint value  sum(w)-1
%            g_lo     = lower bound slack  w_min - w  (positive = violated)
%            g_hi     = upper bound slack  w - w_max  (positive = violated)
%            act_lo   = logical: lower bound active
%            act_hi   = logical: upper bound active
%% ================================================================
function [grad_w, h_bud, g_lo, g_hi, act_lo, act_hi] = ...
        compute_grad(w_k, lam_k, mu_k, ctx)
    p  = ctx.p;
    M  = ctx.M;
    c  = p.c;

    %return signal (linear)
    grad_ret  = -ctx.mu_hat;

    %robust gradient (nonlinear): epsilon * BPBt * w / ||Phi^{1/2} B' w||_2
    nv = sqrt(max(w_k' * ctx.BPBt * w_k, 0));
    if nv > 1e-10
        grad_rob = ctx.epsilon * ctx.BPBt * w_k / nv;
    else
        grad_rob = zeros(M,1);
    end

    %risk and Tikhonov (both linear in w_k)
    grad_risk  = 2*p.gamma    * ctx.Sigma  * w_k;
    grad_tikh  = 2*p.lambda_T * ctx.L_tikh * w_k;

    %transaction cost (non-smooth: sub-gradient of kappa*||w - w_prev||_1)
    grad_tc    = p.kappa * sign(w_k - ctx.w_prev);

    %combined objective sub-gradient
    grad_f     = grad_ret + grad_rob + grad_risk + grad_tikh + grad_tc;

    %augmented Lagrangian: budget equality constraint  sum(w) = 1
    h_bud      = sum(w_k) - 1;
    grad_bud   = (lam_k + c*h_bud) * ones(M,1);

    %augmented Lagrangian: box inequality constraints
    g_lo       = p.w_min - w_k;          %positive when lower bound violated
    g_hi       = w_k - p.w_max;          %positive when upper bound violated
    act_lo     = g_lo > 0;
    act_hi     = g_hi > 0;
    grad_bnds  = -(mu_k(1:M)     + c*g_lo) .* act_lo + ...
                  (mu_k(M+1:2*M) + c*g_hi) .* act_hi;

    grad_w     = grad_f + grad_bud + grad_bnds;
end


%% ================================================================
%  LYAPUNOV ENERGY  --  matches layer2_portfolio.m lines 89-94
%  E_k = (1/2)||nabla_w L_aug||^2 + (1/2)||h_bud||^2
%         + (1/2)||g_lo(active)||^2 + (1/2)||g_hi(active)||^2
%% ================================================================
function E = lyap_energy(w_k, lam_k, mu_k, ctx)
    [gt, h_f, g_lf, g_hf, alf, ahf] = compute_grad(w_k, lam_k, mu_k, ctx);
    E = 0.5*norm(gt)^2 + 0.5*h_f^2 + ...
        0.5*sum(g_lf(alf).^2) + 0.5*sum(g_hf(ahf).^2);
end


%% ================================================================
%  ODE RHS FOR RK4  --  returns (dw, dlam, dmu) for the coupled ODE:
%    dw/dt   = -nabla_w L_aug
%    dlam/dt =  h_bud(w)
%    dmu/dt  = [g_lo.*act_lo; g_hi.*act_hi]
%% ================================================================
function [dw, dlam, dmu] = ode_rhs(w_k, lam_k, mu_k, ctx)
    [grad_w, h_bud, g_lo, g_hi, act_lo, act_hi] = compute_grad(w_k, lam_k, mu_k, ctx);
    M    = ctx.M;
    dw   = -grad_w;
    dlam =  h_bud;
    dmu  =  [g_lo.*act_lo; g_hi.*act_hi];
end


%% ================================================================
%  FORMATTING HELPERS
%% ================================================================
function s = fmt_n(n, n_max)
    if isinf(n), s = sprintf('>%d', n_max);
    else,         s = sprintf('%d', n); end
end

function s = fmt_rk4(n, n_max)
    if isinf(n), s = sprintf('>%d', 4*n_max);
    else,         s = sprintf('%d (= %d x 4)', 4*n, n); end
end
