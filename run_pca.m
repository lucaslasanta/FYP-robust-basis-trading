function [F, pca_model] = run_pca(factors, K, train_end_idx)
%applies pre-test pca: loadings estimated only on training period
%then applied to full sample using training-period standardisation
%this avoids look-ahead bias in the backtest
%
%inputs:
%  factors       - (T x 11) matrix: 10 yc changes + repo (all observations)
%  K             - number of factors (K-1 pca components + repo)
%  train_end_idx - last index of training period (before test period)
%
%outputs:
%  F         - (T x K) factor matrix for full sample, no look-ahead
%  pca_model - struct with loadings and standardisation params for reporting

yc_full  = factors(:, 1:10);
repo_full = factors(:, 11);

%estimate mean and std from training period only
yc_train  = yc_full(1:train_end_idx, :);
mu_train  = mean(yc_train);
std_train = std(yc_train);
std_train(std_train < 1e-10) = 1;

%standardise full sample using training-period parameters
yc_std_full  = (yc_full - mu_train) ./ std_train;
yc_std_train = yc_std_full(1:train_end_idx, :);

%pca on training period only
[coeff, ~, ~, ~, explained] = pca(yc_std_train);

n_pc     = K - 1;
B_pca    = coeff(:, 1:n_pc);
expl_var = sum(explained(1:n_pc));

fprintf('Pre-test PCA: first %d components explain %.1f%% of training-period variance\n', ...
    n_pc, expl_var);
fprintf('  PC1 (level):     %5.1f%%  -- dominant yield-curve factor used in Layer 1 hedge\n', explained(1));
fprintf('  PC2 (slope):     %5.1f%%\n', explained(2));
fprintf('  PC3 (curvature): %5.1f%%\n', explained(3));

%project full sample onto training-period loadings
F_yc = yc_std_full * B_pca;   %(T x n_pc)

%standardise repo using training period only
mu_repo  = mean(repo_full(1:train_end_idx));
std_repo = std(repo_full(1:train_end_idx));
if std_repo < 1e-10, std_repo = 1; end
repo_std = (repo_full - mu_repo) / std_repo;

%full factor matrix
F = [F_yc, repo_std];

%store model for reference
pca_model.coeff      = B_pca;
pca_model.mu_yc      = mu_train;
pca_model.std_yc     = std_train;
pca_model.mu_repo    = mu_repo;
pca_model.std_repo   = std_repo;
pca_model.expl_var      = expl_var;
pca_model.expl_by_pc    = explained(1:n_pc);
pca_model.n_pc          = n_pc;
end