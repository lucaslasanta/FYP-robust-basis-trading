function [returns_w, factors, dates_w, dv01, fut_px, basis_levels] = load_data(~)
%reads from csv files - returns normalised weekly returns and basis levels
%basis levels: FV TY US available from 2015, TU from July 2021 only

fprintf('Reading model_data2.csv...\n');
md2   = readmatrix('model_data2.csv','OutputType','string','NumHeaderLines',1);
nrows = size(md2,1);

dates_all = parse_dates(md2(:,1));
ret_w_all = str2dbl_mat(md2(:,[3 5 7 9]));
repo_all  = str2dbl_mat(md2(:,10));
yc_all    = str2dbl_mat(md2(:,[11 12 13 14 15 16 17 18 19 21]));

fprintf('Reading futures.csv...\n');
fut_raw    = readmatrix('futures.csv','OutputType','string','NumHeaderLines',4);
fut_dates  = parse_dates(fut_raw(:,1));
fut_px_all = str2dbl_mat(fut_raw(:,[2 3 4 5]));
for m = 1:4
    for t = 2:size(fut_px_all,1)
        if isnan(fut_px_all(t,m)), fut_px_all(t,m) = fut_px_all(t-1,m); end
    end
end

fprintf('Reading CTD_selected.csv...\n');
ctd_raw   = readmatrix('CTD_selected.csv','OutputType','string','NumHeaderLines',1);
ctd_dates = parse_dates(ctd_raw(:,1));
dv01_raw  = str2dbl_mat(ctd_raw(:,[9 18 27 36]));
for m = 1:4
    fv = find(~isnan(dv01_raw(:,m)),1,'first');
    if ~isempty(fv) && fv>1, dv01_raw(1:fv-1,m) = dv01_raw(fv,m); end
    for t = 2:size(dv01_raw,1)
        if isnan(dv01_raw(t,m)), dv01_raw(t,m) = dv01_raw(t-1,m); end
    end
end

fprintf('Reading basis.csv...\n');
%cols: 1=TU_basis(#N/A pre-Jul2021), 4=FV_basis, 7=TY_basis, 10=US_basis
basis_raw   = readmatrix('basis.csv','OutputType','string','NumHeaderLines',1);
basis_dates = parse_dates(basis_raw(:,1));
basis_px    = str2dbl_mat(basis_raw(:,[2 5 8 11]));

%TU basis is missing pre July 2021 - forward fill from first valid value
%for consistency keep as NaN before first observation (handled in signal)
tu_first = find(~isnan(basis_px(:,1)),1,'first');
fprintf('  TU_basis available from %s\n', datestr(basis_dates(tu_first)));
fprintf('  FV/TY/US basis available from %s\n', datestr(basis_dates(find(~isnan(basis_px(:,2)),1,'first'))));

fprintf('Filtering to Fridays and aligning...\n');
is_friday = weekday(dates_all) == 6;
has_data  = all(~isnan(ret_w_all),2) & all(~isnan(yc_all),2) & ~isnan(repo_all);

fut_px_aligned = nan(nrows,4);
fri_idx = find(is_friday & has_data);
for k = 1:length(fri_idx)
    t   = fri_idx(k);
    idx = find(fut_dates <= dates_all(t),1,'last');
    if ~isempty(idx), fut_px_aligned(t,:) = fut_px_all(idx,:); end
end
has_px = all(~isnan(fut_px_aligned),2);
valid  = is_friday & has_data & has_px;

px_fri    = fut_px_aligned(valid,:);
returns_w = -(ret_w_all(valid,:) ./ px_fri);
dates_w   = dates_all(valid);
factors   = [yc_all(valid,:), repo_all(valid)];
fut_px    = px_fri;

%align basis levels to weekly friday dates
n_w          = sum(valid);
basis_levels = nan(n_w,4);
for t = 1:n_w
    idx = find(basis_dates <= dates_w(t),1,'last');
    if ~isempty(idx)
        basis_levels(t,:) = basis_px(idx,:);
    end
end
%TU remains NaN before July 2021 - signal code handles this per maturity

%align dv01
dv01 = nan(n_w,4);
for t = 1:n_w
    idx = find(ctd_dates <= dates_w(t),1,'last');
    if ~isempty(idx), dv01(t,:) = dv01_raw(idx,:); end
end

fprintf('Done: %d weekly obs, %s to %s\n', n_w, ...
    datestr(dates_w(1)), datestr(dates_w(end)));
fprintf('Return std: TU=%.5f FV=%.5f TY=%.5f US=%.5f\n', ...
    std(returns_w(:,1)), std(returns_w(:,2)), ...
    std(returns_w(:,3)), std(returns_w(:,4)));
for m = 1:4
    n_b = sum(~isnan(basis_levels(:,m)));
    fprintf('  Basis obs mat %d: %d/%d weeks\n', m, n_b, n_w);
end
end

function d = parse_dates(raw)
n = length(raw); d = NaT(n,1);
base = datetime(1899,12,30);
for i = 1:n
    try
        v = strtrim(char(raw(i)));
        if isempty(v), continue; end
        num = str2double(v);
        if ~isnan(num) && num > 10000
            d(i) = base + days(round(num)); continue;
        end
        for fmt = {'dd/MM/yyyy','yyyy-MM-dd','MM/dd/yyyy'}
            try, d(i) = datetime(v,'InputFormat',fmt{1}); break; catch, end
        end
    catch, end
end
end

function out = str2dbl_mat(data)
out = nan(size(data));
for i = 1:numel(data)
    try
        v = strtrim(char(data(i)));
        n = str2double(v);
        if ~isnan(n), out(i) = n; end
    catch, end
end
end