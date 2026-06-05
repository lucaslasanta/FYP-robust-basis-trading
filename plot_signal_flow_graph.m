function plot_signal_flow_graph()
%LPNN signal flow graph -- clean vertical flow, stability-focused
%Fixes: no LaTeX commands outside math mode, no \dfrac/\tfrac/\searrow,
%       feedback arc within figure bounds, stability clearly labelled.
%Saves fig_signal_flow_graph.pdf

fig = figure('Name','LPNN Signal Flow Graph',...
    'Position',[30 30 1180 860],'Color','white');
axes('Position',[0 0 1 1],'Visible','off');
hold on;
xlim([0 10]);  ylim([0 10]);

xn = @(x) x/10;
yn = @(y) y/10;

%colour palette
c_input = [0.84 0.92 1.00];   %light blue   - inputs
c_lin   = [0.86 1.00 0.86];   %light green  - linear gradients
c_nonl  = [1.00 0.93 0.78];   %light amber  - nonlinear gradients
c_lagr  = [0.93 0.88 1.00];   %light purple - Lagrangian gradient sum
c_prim  = [0.80 0.92 1.00];   %blue-grey    - primal update
c_dual  = [1.00 0.97 0.80];   %light yellow - dual updates
c_lyap  = [0.78 1.00 0.86];   %light mint   - Lyapunov / output
col_arr = [0.20 0.20 0.20];   %dark grey    - arrows
col_fb  = [0.12 0.38 0.75];   %blue         - feedback arc
lw_box  = 1.3;
lw_arr  = 1.1;

% ================================================================
%  TITLE
% ================================================================
text(5, 9.72,'LPNN Coupled Primal--Dual Signal Flow Graph',...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','bold',...
    'Interpreter','latex');
text(5, 9.40,...
    '\textit{Linear path (green): $\nabla_{\mathrm{ret}},\nabla_{\mathrm{risk}},\nabla_{\mathrm{tikh}}$. \quad Nonlinear path (amber): $\nabla_{\mathrm{rob}}^{*},\nabla_{\mathrm{tc}}^{*}$. \quad Lyapunov certificate guarantees KKT convergence.}',...
    'HorizontalAlignment','center','FontSize',7.2,'Interpreter','latex',...
    'Color',[0.35 0.35 0.35]);

% ================================================================
%  ROW 1 -- INPUT PARAMETER BOXES   (y=8.52, h=0.58)
% ================================================================
r1y=8.52;  r1h=0.58;  fs1=8.5;

draw_box([0.10 r1y 1.50 r1h],'$\hat{\mu}$ (signal)',           c_input,fs1,lw_box);
draw_box([1.75 r1y 1.50 r1h],'$\Sigma,\;\gamma$ (risk)',        c_input,fs1,lw_box);
draw_box([3.40 r1y 1.60 r1h],'$L_T,\;\lambda_T$ (Tikhonov)',   c_input,fs1,lw_box);

%dashed separator between linear and nonlinear columns
sep_x = 5.18;
ann_line(fig,xn,yn, sep_x,r1y-0.10, sep_x,r1y+r1h+0.10, [0.50 0.50 0.50],0.9,'--');

draw_box([5.28 r1y 1.68 r1h],'$B,\Phi,\varepsilon$ (robust)',  c_nonl,fs1,lw_box);
draw_box([7.10 r1y 1.68 r1h],'$\kappa,\;w_{\mathrm{prev}}$ (TC)',c_nonl,fs1,lw_box);

%column header labels (placed at vertical midpoint of row 1)
text(2.55, r1y+r1h+0.14,'Linear inputs',...
    'FontSize',7,'Interpreter','latex','HorizontalAlignment','center',...
    'Color',[0.00 0.42 0.10],'FontWeight','bold');
text(6.85, r1y+r1h+0.14,'Nonlinear inputs',...
    'FontSize',7,'Interpreter','latex','HorizontalAlignment','center',...
    'Color',[0.60 0.32 0.00],'FontWeight','bold');

% ================================================================
%  ROW 2 -- GRADIENT NODES   (y=6.82, h=1.00)
% ================================================================
r2y=6.82;  r2h=1.00;  fs2=8.0;

draw_box([0.10 r2y 1.42 r2h],...
    {'$\nabla_{\mathrm{ret}}$','$= -\hat{\mu}$'},                         c_lin,fs2,lw_box);
draw_box([1.68 r2y 1.48 r2h],...
    {'$\nabla_{\mathrm{risk}}$','$= 2\gamma\Sigma w_k$'},                  c_lin,fs2,lw_box);
draw_box([3.30 r2y 1.76 r2h],...
    {'$\nabla_{\mathrm{tikh}}$','$= 2\lambda_T L_T w_k$'},                 c_lin,fs2,lw_box);

ann_line(fig,xn,yn, sep_x,r2y-0.10, sep_x,r2y+r2h+0.10, [0.50 0.50 0.50],0.9,'--');

draw_box([5.28 r2y 1.85 r2h],...
    {'$\nabla_{\mathrm{rob}}^{*}$','$=\varepsilon\frac{B\Phi B^{\top}w_k}{\|B^{\top}w_k\|}$'}, c_nonl,fs2,lw_box);
draw_box([7.25 r2y 1.75 r2h],...
    {'$\nabla_{\mathrm{tc}}^{*}$','$=\kappa\,\mathrm{sgn}(w_k - w_{\mathrm{prev}})$'},         c_nonl,fs2,lw_box);

%row label
text(2.55, r2y+r2h+0.13,'Linear gradient terms',...
    'FontSize',7,'Interpreter','latex','HorizontalAlignment','center',...
    'Color',[0.00 0.42 0.10],'FontWeight','bold');
text(6.85, r2y+r2h+0.13,'Nonlinear gradient terms ($*$)',...
    'FontSize',7,'Interpreter','latex','HorizontalAlignment','center',...
    'Color',[0.60 0.32 0.00],'FontWeight','bold');

%arrows row1 -> row2 (vertical, centred on each box)
arr(fig,xn,yn,  0.85, r1y,  0.81, r2y+r2h, col_arr,lw_arr);
arr(fig,xn,yn,  2.50, r1y,  2.42, r2y+r2h, col_arr,lw_arr);
arr(fig,xn,yn,  4.20, r1y,  4.18, r2y+r2h, col_arr,lw_arr);
arr(fig,xn,yn,  6.12, r1y,  6.20, r2y+r2h, col_arr,lw_arr);
arr(fig,xn,yn,  7.94, r1y,  8.12, r2y+r2h, col_arr,lw_arr);

% ================================================================
%  ROW 3 -- COMBINED LAGRANGIAN GRADIENT   (y=5.12, h=1.00)
% ================================================================
r3y=5.12;  r3h=1.00;  fs3=7.8;

draw_box([0.10 r3y 9.50 r3h],...
    {'Combined Augmented Lagrangian Gradient  $\nabla_w\mathcal{L}_{\mathrm{aug}}$',...
     '$\nabla_w\mathcal{L}_{\mathrm{aug}} = \nabla_{\mathrm{ret}} + \nabla_{\mathrm{risk}} + \nabla_{\mathrm{tikh}} + \nabla_{\mathrm{rob}}^{*} + \nabla_{\mathrm{tc}}^{*} + (\lambda_k + c\,h_{\mathrm{bud}})\,\mathbf{1} + \mu_k^{\mathrm{lo}}\,\mathbf{1}_{[w<w_{\min}]} - \mu_k^{\mathrm{hi}}\,\mathbf{1}_{[w>w_{\max}]}$'},...
    c_lagr,fs3,lw_box);

%arrows row2 -> row3
arr(fig,xn,yn,  0.81, r2y,  0.81, r3y+r3h, col_arr,lw_arr);
arr(fig,xn,yn,  2.42, r2y,  2.42, r3y+r3h, col_arr,lw_arr);
arr(fig,xn,yn,  4.18, r2y,  4.18, r3y+r3h, col_arr,lw_arr);
arr(fig,xn,yn,  6.20, r2y,  6.20, r3y+r3h, col_arr,lw_arr);
arr(fig,xn,yn,  8.12, r2y,  8.12, r3y+r3h, col_arr,lw_arr);

% ================================================================
%  ROW 4 -- PRIMAL UPDATE   (y=3.22, h=1.18)
%  NOTE: no LaTeX spacing commands (\quad) outside $...$
% ================================================================
r4y=3.22;  r4h=1.18;  fs4=8.0;

draw_box([0.10 r4y 9.50 r4h],...
    {'\bf Primal Update -- Euler discretisation of $\dot{w}=-\nabla_w\mathcal{L}_{\mathrm{aug}}$ (continuous-time ODE)',...
     'Euler step: $w_{k+1} = w_k - \tau\,\nabla_w\mathcal{L}_{\mathrm{aug}}(w_k,\lambda_k,\mu_k)$;  stability requires $\tau < 2/\lambda_{\max}(\nabla^2\mathcal{L})$',...
     'Adaptive step (nonlinear Layer 1): $\tau_1 = \min\!\left(\tau,\;0.9/(2\sigma_f^2)\right)$ -- prevents Euler instability under regime shifts',...
     'Projection: $w_{k+1} \leftarrow \mathrm{clip}(w_{k+1},\,w_{\min},\,w_{\max})$'},...
    c_prim,fs4,lw_box);

%single centred arrow row3 -> row4
arr(fig,xn,yn, 4.85, r3y, 4.85, r4y+r4h, col_arr,lw_arr+0.2);

% ================================================================
%  ROW 5 -- DUAL UPDATES + LYAPUNOV   (y=1.72, h=1.10)
% ================================================================
r5y=1.72;  r5h=1.10;  fs5=8.0;

draw_box([0.10 r5y 2.90 r5h],...
    {'Budget Multiplier Update',...
     '$\lambda_{k+1} = \lambda_k + \tau\,h_{\mathrm{bud}}(w_{k+1})$',...
     '$h_{\mathrm{bud}} = \sum_i w_i - 1$'},...
    c_dual,fs5,lw_box);

draw_box([3.20 r5y 3.15 r5h],...
    {'Bound Multiplier Updates',...
     '$\mu^{\mathrm{lo}}_{k+1} = \mu^{\mathrm{lo}}_k + \tau\max(w_{\min}-w_{k+1},\,0)$',...
     '$\mu^{\mathrm{hi}}_{k+1} = \mu^{\mathrm{hi}}_k + \tau\max(w_{k+1}-w_{\max},\,0)$'},...
    c_dual,fs5,lw_box);

%Lyapunov box -- stability certificate (prominent)
draw_box([6.55 r5y 3.05 r5h],...
    {'\bf Lyapunov Stability Certificate  $(\dot{E}\leq 0)$',...
     '$E_k = \frac{1}{2}\|\nabla_w\mathcal{L}\|^2 + \frac{1}{2}\|h_{\mathrm{bud}}\|^2 + \frac{1}{2}\|g_{\mathrm{lo}}\|^2 + \frac{1}{2}\|g_{\mathrm{hi}}\|^2$',...
     'Monotone decrease $\Rightarrow$ asymptotic convergence to KKT point'},...
    c_lyap,fs5,lw_box);

%arrows row4 -> row5
arr(fig,xn,yn, 1.55, r4y, 1.55, r5y+r5h, col_arr,lw_arr);
arr(fig,xn,yn, 4.78, r4y, 4.78, r5y+r5h, col_arr,lw_arr);
arr(fig,xn,yn, 8.08, r4y, 8.08, r5y+r5h, col_arr,lw_arr);

% ================================================================
%  ROW 6 -- OUTPUT   (y=0.22, h=0.90)
% ================================================================
r6y=0.22;  r6h=0.90;  fs6=8.5;

draw_box([2.50 r6y 5.00 r6h],...
    {'\bf Output: $w^* = w_k$ at $k=n_{\mathrm{iter}}$',...
     'Warm-start: $w_{\mathrm{prev}} \leftarrow w^*$ for next rebalancing date'},...
    c_lyap,fs6,lw_box);

arr(fig,xn,yn, 5.00, r5y, 5.00, r6y+r6h, col_arr,lw_arr+0.2);

% ================================================================
%  STATE FEEDBACK ARC  (right side, x=9.60)
%  w_{k+1},lambda_{k+1},mu_{k+1}  fed back as w_k,lambda_k,mu_k
%  Connects dual-update row (bottom) back to Lagrangian row (top)
% ================================================================
fx      = 9.60;
fb_bot  = r5y + r5h/2;
fb_top  = r3y + r3h/2;
fb_mid  = r4y + r4h/2;

%three horizontal stubs: from each active row to the vertical line
ann_line(fig,xn,yn, 9.60,fb_bot,  fx,fb_bot,  col_fb,1.5,'-');
ann_line(fig,xn,yn, 9.60,fb_mid,  fx,fb_mid,  col_fb,1.5,'-');
ann_line(fig,xn,yn, 9.60,fb_top,  fx,fb_top,  col_fb,1.5,'-');
%vertical backbone
ann_line(fig,xn,yn, fx,fb_bot, fx,fb_top, col_fb,1.5,'-');
%arrowhead into Lagrangian row (pointing left)
annotation(fig,'arrow',[xn(fx) xn(9.62)],[yn(fb_top) yn(fb_top)],...
    'Color',col_fb,'HeadLength',8,'HeadWidth',7,...
    'HeadStyle','vback2','LineWidth',1.5);

%label -- white background box so it reads clearly over any row colour
text(fx-0.12, (fb_top+fb_bot)/2,...
    {'State feedback','$w_k \leftarrow w_{k+1}$','$\lambda_k \leftarrow \lambda_{k+1}$','$\mu_k \leftarrow \mu_{k+1}$'},...
    'FontSize',7.0,'Interpreter','latex','Color',col_fb,...
    'HorizontalAlignment','right','VerticalAlignment','middle',...
    'BackgroundColor','white','EdgeColor',col_fb,'Margin',3,'LineWidth',0.7);

% ================================================================
%  LEGEND (bottom-left, clear of content)
% ================================================================
lx=0.12; ly=0.25; lw2=2.30; lh2=1.35;
rectangle('Position',[lx ly lw2 lh2],...
    'FaceColor',[0.98 0.98 0.98],'EdgeColor',[0.60 0.60 0.60],'LineWidth',0.8);
text(lx+lw2/2, ly+lh2-0.13,'Legend','FontSize',7,'FontWeight','bold',...
    'Interpreter','latex','HorizontalAlignment','center');
rectangle('Position',[lx+0.10 ly+0.85 0.20 0.20],'FaceColor',c_lin,...
    'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.7);
text(lx+0.38, ly+0.95,'Linear gradient term',...
    'FontSize',6.5,'Interpreter','latex','VerticalAlignment','middle');
rectangle('Position',[lx+0.10 ly+0.55 0.20 0.20],'FaceColor',c_nonl,...
    'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.7);
text(lx+0.38, ly+0.65,'Nonlinear gradient term ($*$)',...
    'FontSize',6.5,'Interpreter','latex','VerticalAlignment','middle');
rectangle('Position',[lx+0.10 ly+0.25 0.20 0.14],'FaceColor',col_fb+0.6*(1-col_fb),...
    'EdgeColor',col_fb,'LineWidth',0.7);
text(lx+0.38, ly+0.32,'State feedback arc',...
    'FontSize',6.5,'Interpreter','latex','VerticalAlignment','middle');

saveas(fig,'fig_signal_flow_graph.pdf');
fprintf('Saved: fig_signal_flow_graph.pdf\n');
end

% ================================================================
%  LOCAL HELPERS
% ================================================================
function draw_box(pos, txt, fc, fsize, lw)
    rectangle('Position',pos,'Curvature',0.08,...
        'FaceColor',fc,'EdgeColor',[0.28 0.28 0.28],'LineWidth',lw);
    text(pos(1)+pos(3)/2, pos(2)+pos(4)/2, txt,...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',fsize,'Interpreter','latex','BackgroundColor','none');
end

function arr(fig, xn, yn, x1, y1, x2, y2, col, lw)
    annotation(fig,'arrow',[xn(x1) xn(x2)],[yn(y1) yn(y2)],...
        'Color',col,'HeadLength',7,'HeadWidth',6,...
        'HeadStyle','vback2','LineWidth',lw);
end

function ann_line(fig, xn, yn, x1, y1, x2, y2, col, lw, ls)
    annotation(fig,'line',[xn(x1) xn(x2)],[yn(y1) yn(y2)],...
        'Color',col,'LineWidth',lw,'LineStyle',ls);
end
