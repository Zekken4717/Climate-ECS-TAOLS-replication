% =========================================================================
% Plot_Data_Unit_Root_Normality_Tests.m
%
% Pre-estimation diagnostics for the climate multicointegration study.
% Produces three figures and two printed tables.
%
% FIGURES
%   Figure 1 -- Time series of all four data series (1850-2014)
%               Left panel:  both radiative forcing series
%               Right panel: both surface temperature series
%
%   Figure 2 -- Histograms with fitted normal density (2 x 2 panel)
%               Data: AR(3) residuals of first differences.
%               JB statistic and p-value annotated on each panel.
%
%   Figure 3 -- Normal Q-Q plots of AR(3) residuals of first differences
%               (2 x 2 panel).
%               S-shaped deviation indicates heavy-tailed distribution.
%
% TABLES
%   Table I  -- Unit Root Tests on Level Series
%               Left:  ADF test (H0: unit root; ARD model, BIC lag).
%               Right: KPSS test (H0: trend stationarity; Andrews bandwidth).
%               Fail to reject ADF AND reject KPSS => strong evidence of I(1).
%
%   Table II -- Normality Tests on AR(3) Residuals of First Differences
%               Model: Dy_t = c + phi_1*Dy_{t-1} + phi_2*Dy_{t-2} +
%                             phi_3*Dy_{t-3} + e_t
%               Jarque-Bera test (H0: normality) on OLS residuals e_t.
%               Filtering out autocorrelation before JB gives a more
%               stringent test of the distributional assumption on innovations.
%               Reject H0 => non-normal innovations.
%
% Series tested:
%   TotalRF   -- full-efficacy radiative forcing
%   MarvelRF  -- partial-efficacy RF (Marvel et al. efficacy correction)
%   Berkeley  -- Berkeley Earth global surface temperature
%   HadCRUT   -- HadCRUT global surface temperature
%
% Data: OHC Data.xls, sheet 'Data', annual 1850-2014 (T = 165 obs).
%
% Outputs (saved to Results_Diagnostics\):
%   unit_root_normality_tests.out        -- printed tables (diary)
%   rf_temperature_combined.fig/.eps/.png -- time series plot (Figure 1)
%   fd_histograms.fig / .eps / .png      -- AR(3) resid histograms (Figure 2)
%   fd_qqplots.fig / .eps / .png         -- AR(3) resid Q-Q plots (Figure 3)
%
% Also saves rf_temperature_combined.png/.eps to ..\Paper_June_2026\ for LaTeX.
%
% Toolbox requirements:
%   Econometrics Toolbox            -- adftest, kpsstest
%   Statistics & Machine Learning   -- jbtest, skewness, kurtosis, qqplot
%
% Reference:
%   Bruns, Csereklyei & Stern (2020). A multicointegration model of global
%   climate change. Journal of Econometrics, 214(1), 175-197.
% =========================================================================

clear;
close all;

cd(fileparts(mfilename('fullpath')));

% -------------------------------------------------------------------------
% OUTPUT DIRECTORY
% -------------------------------------------------------------------------
results_dir = 'Results_Diagnostics';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

diary_file = fullfile(results_dir, 'unit_root_normality_tests.out');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);
diary on;

fprintf('=================================================================\n');
fprintf(' Pre-Estimation Diagnostic Tests\n');
fprintf(' Climate Multicointegration Analysis -- Justin Sun (2026)\n');
fprintf(' Data: OHC Data.xls, annual 1850-2014\n');
fprintf(' Run date: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf('=================================================================\n\n');

% -------------------------------------------------------------------------
% LOAD DATA
% -------------------------------------------------------------------------
data  = readtable('OHC Data.xls', 'Sheet', 'Data');
year  = (1850:1:2014)';
T_lv  = length(year);          % 165 annual observations in levels

f   = data.TotalRF;             % full-efficacy radiative forcing
f_m = data.MarvelRF;            % partial-efficacy RF (Marvel et al.)
s   = data.Berkeley;            % Berkeley Earth surface temperature
s_h = data.HADCRUT;             % HadCRUT surface temperature

% First differences (164 obs, covering 1851-2014)
df   = diff(f);
df_m = diff(f_m);
ds   = diff(s);
ds_h = diff(s_h);

year_fd = year(2:end);
T_fd    = length(year_fd);      % 164 observations

% Bundle for looping
series_lv = {f,   f_m,    s,         s_h};
series_fd = {df,  df_m,   ds,        ds_h};
lbl_lv    = {'TotalRF', 'MarvelRF', 'Berkeley', 'HadCRUT'};
lbl_fd    = {'D.TotalRF', 'D.MarvelRF', 'D.Berkeley', 'D.HadCRUT'};
n_ser     = 4;
p_max     = 3;      % max ADF lag (appropriate for annual data, T = 165)


% =========================================================================
% TABLE I.  UNIT ROOT TESTS ON LEVEL SERIES
% =========================================================================
%
%  ADF Test:  H0 = unit root.  Fail to reject H0 => evidence of unit root.
%             Specification: constant only (ARD model).
%             ADF lag chosen by BIC over 0 to p_max; lag shown in ().
%
%  KPSS Test: H0 = trend stationarity.  Reject H0 => evidence of unit root.
%             Bandwidth: Andrews (1991) automatic selection.
%
% =========================================================================
fprintf('=================================================================\n');
fprintf(' Table I.  Unit Root Tests -- Level Series\n');
fprintf('           T = %d annual observations, 1850-2014\n', T_lv);
fprintf('=================================================================\n\n');
fprintf(' ADF  H0: unit root.  Fail to reject => evidence of unit root.\n');
fprintf(' KPSS H0: trend stationarity.  Reject => evidence of unit root.\n');
fprintf(' If ADF fails to reject AND KPSS rejects: strong evidence of I(1).\n');
fprintf(' ADF lag (in parentheses) chosen by BIC; max lag = %d.\n', p_max);
fprintf(' KPSS bandwidth: Andrews (1991) automatic method.\n\n');

% ---------- Table header ----------
SEP1 = repmat('-', 1, 54);
fprintf('  %-12s | %-20s | %-16s\n', ...
    '', '  ADF: Const (ARD)', ' KPSS: Trend');
fprintf('  %-12s | %-10s  %-6s  | %-7s  %-6s\n', ...
    'Series', 'Stat(Lag)', 'p-val', 'Stat', 'p-val');
fprintf('  %s\n', SEP1);

% ---------- Table rows ----------
for i = 1:n_ser
    y   = series_lv{i};
    lbl = lbl_lv{i};

    % ADF -- constant only (ARD)
    p_A = adf_bic_select(y, p_max);
    [~, pvA, stA] = adftest(y, 'lags', p_A, 'model', 'ARD', 'alpha', 0.05);
    strA = sprintf('%7.3f(%d)', stA, p_A);    % e.g. " -1.234(1)"

    % KPSS -- trend stationarity
    [~, pvK, stK] = kpsstest(y, 'trend', true);

    fprintf('  %-12s | %-10s  %6.4f  | %7.4f  %6.4f\n', ...
        lbl, strA, pvA, stK, pvK);
end

fprintf('  %s\n\n', SEP1);

fprintf('  Notes:\n');
fprintf('  (a) ADF critical values (Dickey-Fuller distribution, ARD model):\n');
fprintf('      10%%: -2.57;  5%%: -2.87;  1%%: -3.43.\n');
fprintf('  (b) KPSS critical values (trend stationarity, KPSS 1992, Table 1):\n');
fprintf('      10%%: 0.119;  5%%: 0.146;  1%%: 0.216.\n');
fprintf('  (c) KPSS p-values are bounded in [0.01, 0.10] by the tabulated\n');
fprintf('      critical values.  p <= 0.01 means the stat exceeds the 1%% CV.\n');
fprintf('  (d) ADF lag in parentheses; selected to minimise BIC over 0-%d.\n\n', p_max);


% -------------------------------------------------------------------------
% Compute JB statistics on raw first differences (silent -- no table
% printed).  Results are stored in jb_store for histogram annotations.
% -------------------------------------------------------------------------
jb_store = zeros(n_ser, 4);   % [jbstat, pval, skewness, ex.kurtosis]

for i = 1:n_ser
    dy  = series_fd{i};
    sk  = skewness(dy);
    ek  = kurtosis(dy) - 3;
    [~, pval, jbstat] = jbtest(dy, 0.05);
    jb_store(i, :) = [jbstat, pval, sk, ek];
end

% =========================================================================
% TABLE II.  NORMALITY TESTS ON AR(3) RESIDUALS FROM FIRST DIFFERENCES
% =========================================================================
%
%  Model:  Dy_t = c + phi_1*Dy_{t-1} + phi_2*Dy_{t-2} + phi_3*Dy_{t-3} + e_t
%
%  OLS residuals e_t are tested for normality (Jarque-Bera).
%  Fitting AR(3) to each first-difference series removes autocorrelation
%  structure before testing the distributional assumption, providing a
%  more stringent and cleaner check of Gaussianity.
%
%  T_ar = T_fd - 3 = 161 (three lagged regressors consume three obs).
%
% =========================================================================
fprintf('=================================================================\n');
fprintf(' Table II.  Normality Tests -- AR(3) Residuals of First Differences\n');
fprintf('            Model: Dy_t = c + phi_1*Dy_{t-1} + phi_2*Dy_{t-2} + phi_3*Dy_{t-3} + e_t\n');
fprintf('            Residuals: T = %d obs (1854-2014, after 3-lag alignment)\n', T_fd - 3);
fprintf('=================================================================\n\n');
fprintf(' H0: AR(3) residual is normally distributed.\n');
fprintf(' JB = (n/6) * [Skewness^2 + (ExcessKurtosis^2)/4]  ~  chi^2(2).\n');
fprintf(' Critical values: 10%% = 4.61,  5%% = 5.99,  1%% = 9.21.\n\n');

p_ar = 3;
T_ar = T_fd - p_ar;    % 164 - 3 = 161

% ---------- Table header ----------
SEP3 = repmat('-', 1, 76);
fprintf('  %-12s  %9s  %9s  %10s  %8s  %6s  %6s  %6s\n', ...
    'Series', 'Skewness', 'Ex.Kurt', 'JB stat', 'p-value', 'Rej1%', 'Rej5%', 'Rej10%');
fprintf('  %s\n', SEP3);

% ---------- Table rows ----------
jb_ar3_store = zeros(n_ser, 4);   % [jbstat, pval, skewness, ex.kurtosis]
resid_ar3    = cell(n_ser, 1);    % AR(3) residuals for each series

for i = 1:n_ser
    dy  = series_fd{i};     % first-difference series (T_fd x 1)
    lbl = lbl_fd{i};

    % AR(3) OLS: design matrix columns = [1, Dy_{t-1}, Dy_{t-2}, Dy_{t-3}]
    %   Dependent  : dy(4 : end)     -- t = 4,...,164  (161 obs)
    %   Lag 1 col  : dy(3 : end-1)   -- Dy_{t-1}
    %   Lag 2 col  : dy(2 : end-2)   -- Dy_{t-2}
    %   Lag 3 col  : dy(1 : end-3)   -- Dy_{t-3}
    Y_dep = dy(p_ar+1 : end);
    X_ar  = [ones(T_ar, 1), ...
             dy(p_ar   : end-1), ...
             dy(p_ar-1 : end-2), ...
             dy(p_ar-2 : end-3)];

    b_ar  = X_ar \ Y_dep;
    resid = Y_dep - X_ar * b_ar;

    sk = skewness(resid);
    ek = kurtosis(resid) - 3;

    [~, pval, jbstat] = jbtest(resid, 0.05);
    h1  = jbtest(resid, 0.01);
    h5  = jbtest(resid, 0.05);
    h10 = jbtest(resid, 0.10);

    jb_ar3_store(i, :) = [jbstat, pval, sk, ek];
    resid_ar3{i}       = resid;

    r1  = 'No '; if h1,  r1  = 'Yes'; end
    r5  = 'No '; if h5,  r5  = 'Yes'; end
    r10 = 'No '; if h10, r10 = 'Yes'; end

    fprintf('  %-12s  %9.4f  %9.4f  %10.3f  %8.4f  %6s  %6s  %6s\n', ...
        lbl, sk, ek, jbstat, pval, r1, r5, r10);
end

fprintf('  %s\n\n', SEP3);

fprintf('  Notes:\n');
fprintf('  (a) AR(3) model fitted by OLS; T_ar = %d residuals (1854-2014).\n', T_ar);
fprintf('  (b) Excess kurtosis = sample kurtosis - 3 (normal distribution => 0).\n');
fprintf('  (c) Removing AR dynamics before JB provides a more stringent test:\n');
fprintf('      if non-normality persists in e_t, it cannot be attributed to\n');
fprintf('      residual autocorrelation.\n\n');

diary off;


% =========================================================================
% FIGURE 1:  Time series of all four data series (1 x 2 panel)
%            Left:  both radiative forcing series (full- and partial-efficacy)
%            Right: both surface temperature series (Berkeley Earth and HadCRUT)
% =========================================================================
figure(1)

subplot(1,2,1)
plot(year, f,   'r--', 'LineWidth', 1.5);
hold on
plot(year, f_m, 'k-',  'LineWidth', 1.5);
hold off
xlim([1850 2014]);
xticks(1850:25:2000);
legend({'Full-efficacy RF', 'Partial-efficacy RF'}, ...
       'Location', 'northwest', 'FontSize', 8.5);
xlabel('Year', 'FontSize', 10);
ylabel('Watts Per Square Meter', 'FontSize', 10);
title('(a) Radiative Forcing', 'FontSize', 10);
grid on

subplot(1,2,2)
plot(year, s_h, 'r--', 'LineWidth', 1.5);
hold on
plot(year, s,   'k-',  'LineWidth', 1.5);
hold off
xlim([1850 2014]);
xticks(1850:25:2000);
legend({'HadCRUT', 'Berkeley Earth'}, ...
       'Location', 'northwest', 'FontSize', 8.5);
xlabel('Year', 'FontSize', 10);
ylabel('Celsius', 'FontSize', 10);
title('(b) Surface Temperature', 'FontSize', 10);
grid on

%% Save Figure 1
figure(1)
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
width  = 7.0;
height = 3;
left   = (papersize(1) - width)  / 2;
bottom = (papersize(2)/2 - height) / 2;
set(gcf, 'PaperPosition', [left, bottom, width, height]);

saveas(gcf, fullfile(results_dir, 'rf_temperature_combined.fig'));
saveas(gcf, fullfile(results_dir, 'rf_temperature_combined.eps'), 'epsc');
print(gcf,  fullfile(results_dir, 'rf_temperature_combined.png'), '-dpng', '-r300');

% Also save to Paper_Justin3 for LaTeX inclusion
paper_dir = fullfile('..', 'Paper_June_2026');
if exist(paper_dir, 'dir')
    print(gcf, fullfile(paper_dir, 'rf_temperature_combined.png'), '-dpng', '-r300');
    print(gcf, fullfile(paper_dir, 'rf_temperature_combined.eps'), '-depsc', '-vector', '-r300');
end


% =========================================================================
% FIGURE 2:  Histograms of AR(3) residuals with fitted normal density
%            (2 x 2 panel).  JB statistic and p-value annotated on each
%            panel.  Residuals from: Dy_t = c + phi_1*Dy_{t-1} +
%            phi_2*Dy_{t-2} + phi_3*Dy_{t-3} + e_t.
% =========================================================================
figure(2);

ar3_hist_titles  = {'(a) AR(3) Resid: D.TotalRF',  '(b) AR(3) Resid: D.MarvelRF', ...
                    '(c) AR(3) Resid: D.Berkeley',  '(d) AR(3) Resid: D.HadCRUT'};
fd_ylabels = {'W m^{-2}', 'W m^{-2}', '^{\circ}C', '^{\circ}C'};

for i = 1:n_ser
    subplot(2, 2, i);

    dy = resid_ar3{i};
    mu = mean(dy);
    sg = std(dy);

    % Density-normalised histogram
    histogram(dy, 'Normalization', 'pdf', ...
              'FaceColor', [0.75, 0.75, 0.75], 'EdgeColor', 'w');
    hold on;
    % Fitted normal density
    xr = linspace(min(dy) - 0.5*sg, max(dy) + 0.5*sg, 300);
    plot(xr, normpdf(xr, mu, sg), 'r-', 'LineWidth', 1.5);
    hold off;

    title(ar3_hist_titles{i}, 'FontSize', 9);
    xlabel(fd_ylabels{i}, 'FontSize', 9);
    ylabel('Density', 'FontSize', 9);

    % JB annotation (upper-right corner) -- from AR(3) residuals
    ann = sprintf('JB = %.2f\np = %.4f', jb_ar3_store(i,1), jb_ar3_store(i,2));
    text(0.97, 0.97, ann, 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', [0.7 0.7 0.7]);

    legend({'Observed', 'Normal fit'}, 'Location', 'northwest', 'FontSize', 8);
    grid on; box on;
    set(gca, 'FontSize', 9);
end

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5, 0.75, 7.5, 5.5]);
saveas(gcf, fullfile(results_dir, 'fd_histograms.fig'));
print(gcf, fullfile(results_dir, 'fd_histograms.eps'), '-depsc', '-vector', '-r300');
print(gcf, fullfile(results_dir, 'fd_histograms.png'), '-dpng', '-r300');


% =========================================================================
% FIGURE 3:  Normal Q-Q plots of AR(3) residuals of first differences
%            (2 x 2 panel).  Systematic S-shaped deviation from the
%            reference line indicates heavy-tailed (leptokurtic)
%            distribution.
% =========================================================================
figure(3);

qq_titles = {'(a) Q-Q: AR(3) Resid: D.TotalRF',  '(b) Q-Q: AR(3) Resid: D.MarvelRF', ...
             '(c) Q-Q: AR(3) Resid: D.Berkeley',  '(d) Q-Q: AR(3) Resid: D.HadCRUT'};

for i = 1:n_ser
    subplot(2, 2, i);
    qqplot(resid_ar3{i});
    title(qq_titles{i}, 'FontSize', 9);
    xlabel('Standard Normal Quantiles', 'FontSize', 8);
    ylabel('Sample Quantiles', 'FontSize', 8);
    grid on;
    set(gca, 'FontSize', 9);
end

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5, 0.75, 7.5, 5.5]);
saveas(gcf, fullfile(results_dir, 'fd_qqplots.fig'));
print(gcf, fullfile(results_dir, 'fd_qqplots.eps'), '-depsc', '-vector', '-r300');
print(gcf, fullfile(results_dir, 'fd_qqplots.png'), '-dpng', '-r300');

fprintf('\nDone. Results saved to: %s\\\n', results_dir);
fprintf('  unit_root_normality_tests.out\n');
fprintf('  rf_temperature_combined.fig / .eps / .png\n');
fprintf('  fd_histograms.fig / .eps / .png\n');
fprintf('  fd_qqplots.fig / .eps / .png\n\n');


% =========================================================================
% LOCAL FUNCTION: BIC-based ADF lag selection
% =========================================================================
function opt_p = adf_bic_select(y, p_max)
%ADF_BIC_SELECT  Select ADF lag order by BIC (Bayesian Information Criterion).
%   ARD model: Δy_t = α + ρ*y_{t-1} + γ_1*Δy_{t-1} + ... + γ_p*Δy_{t-p} + ε_t
%
%   y     -- column vector, the level series (length T)
%   p_max -- maximum lag to evaluate (integer >= 0)
%
%   Returns opt_p in {0, 1, ..., p_max} that minimises BIC.
%   BIC = T_eff * log(RSS / T_eff) + k * log(T_eff),  k = number of regressors.

    T  = length(y);
    dy = diff(y);       % T-1 first differences

    best_bic = Inf;
    opt_p    = 0;

    for p = 0:p_max
        T_eff = T - p - 1;              % effective sample after lag alignment

        % Dependent variable: Δy_{p+2}, ..., Δy_T
        Y_dep = dy(p+1 : end);          % length T_eff

        % Lagged level regressor: y_{p+1}, ..., y_{T-1}
        X_lag = y(p+1 : end-1);         % length T_eff

        % Lagged-difference regressors: Δy_{t-j}, j = 1, ..., p
        if p > 0
            X_dlag = zeros(T_eff, p);
            for j = 1:p
                X_dlag(:, j) = dy(p+1-j : end-j);
            end
        else
            X_dlag = zeros(T_eff, 0);
        end

        % Build full design matrix (ARD: intercept + lagged level + lagged diffs)
        X = [ones(T_eff, 1), X_lag, X_dlag];

        % OLS residuals and BIC
        b   = X \ Y_dep;
        res = Y_dep - X * b;
        k   = size(X, 2);
        RSS = sum(res .^ 2);
        bic_val = T_eff * log(RSS / T_eff) + k * log(T_eff);

        if bic_val < best_bic
            best_bic = bic_val;
            opt_p    = p;
        end
    end
end
