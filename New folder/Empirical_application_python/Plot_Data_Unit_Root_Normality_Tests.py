# =========================================================================
# Plot_Data_Unit_Root_Normality_Tests.py
#
# Pre-estimation diagnostics for the climate multicointegration study.
# Produces three figures and two printed tables.
#
# FIGURES
#   Figure 1 -- Time series of all four data series (1850-2014)
#   Figure 2 -- Histograms of AR(3) residuals with fitted normal density (2x2 panel)
#   Figure 3 -- Normal Q-Q plots of AR(3) residuals (2x2 panel)
#
# TABLES
#   Table I  -- ADF (H0: unit root) and KPSS (H0: trend stationary) tests
#   Table II -- Jarque-Bera normality on AR(3) residuals of first differences
#
# Data: OHC Data.xls, sheet 'Data', annual 1850-2014 (T = 165 obs).
#
# Python dependencies: numpy, pandas, scipy, statsmodels, matplotlib, xlrd
#
# Reference:
#   Bruns, Csereklyei & Stern (2020). Journal of Econometrics, 214(1), 175-197.
# =========================================================================

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import jarque_bera
from statsmodels.tsa.stattools import adfuller, kpss
from datetime import datetime

script_dir  = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, 'Results_Diagnostics')
paper_dir   = os.path.join(os.path.dirname(script_dir), 'Paper_June_2026')
os.makedirs(results_dir, exist_ok=True)

log_path = os.path.join(results_dir, 'unit_root_normality_tests.out')
open(log_path, 'w').close()   # overwrite

def log(msg=''):
    print(msg)
    with open(log_path, 'a') as fh:
        fh.write(msg + '\n')

# -------------------------------------------------------------------------
# LOCAL FUNCTION: BIC-based ADF lag selection
# ARD model: Δy_t = α + ρ*y_{t-1} + γ_1*Δy_{t-1} + ... + γ_p*Δy_{t-p} + ε_t
# ARD model: Δy_t = α + ρ*y_{t-1} + γ_1*Δy_{t-1} + ... + γ_p*Δy_{t-p} + ε_t
# BIC = T_eff * log(RSS / T_eff) + k * log(T_eff)
# -------------------------------------------------------------------------
def adf_bic_select(y, p_max):
    T  = len(y)
    dy = np.diff(y)    # length T-1; dy[i] = y[i+1] - y[i]

    best_bic = np.inf
    opt_p    = 0

    for p in range(p_max + 1):
        T_eff = T - p - 1

        Y_dep = dy[p:]             # Δy at effective obs, length T_eff
        X_lag = y[p:T - 1]        # lagged level y_{t-1}, length T_eff

        if p > 0:
            X_dlag = np.column_stack([dy[p - j: T - 1 - j] for j in range(1, p + 1)])
            X = np.column_stack([np.ones(T_eff), X_lag, X_dlag])
        else:
            X = np.column_stack([np.ones(T_eff), X_lag])

        b   = np.linalg.lstsq(X, Y_dep, rcond=None)[0]
        res = Y_dep - X @ b
        k   = X.shape[1]
        RSS = float(res @ res)
        bic = T_eff * np.log(RSS / T_eff) + k * np.log(T_eff)

        if bic < best_bic:
            best_bic = bic
            opt_p    = p

    return opt_p

# -------------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------------
data = pd.read_excel(os.path.join(script_dir, 'OHC Data.xls'), sheet_name='Data')
year  = np.arange(1850, 2015)
T_lv  = len(year)    # 165

f   = data['TotalRF'].values
f_m = data['MarvelRF'].values
s   = data['Berkeley'].values
s_h = data['HADCRUT'].values

df   = np.diff(f)
df_m = np.diff(f_m)
ds   = np.diff(s)
ds_h = np.diff(s_h)

year_fd = year[1:]
T_fd    = len(year_fd)    # 164

series_lv  = [f,   f_m,    s,         s_h]
series_fd  = [df,  df_m,   ds,        ds_h]
lbl_lv     = ['TotalRF', 'MarvelRF', 'Berkeley', 'HadCRUT']
lbl_fd     = ['D.TotalRF', 'D.MarvelRF', 'D.Berkeley', 'D.HadCRUT']
n_ser      = 4
p_max      = 3

# =========================================================================
# TABLE I.  UNIT ROOT TESTS ON LEVEL SERIES
# =========================================================================
log('=================================================================')
log(' Pre-Estimation Diagnostic Tests')
log(' Climate Multicointegration Analysis -- Justin Sun (2026)')
log(' Data: OHC Data.xls, annual 1850-2014')
log(f' Run date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
log('=================================================================\n')

log('=================================================================')
log(f' Table I.  Unit Root Tests -- Level Series')
log(f'           T = {T_lv} annual observations, 1850-2014')
log('=================================================================\n')
log(' ADF  H0: unit root.  Fail to reject => evidence of unit root.')
log(' KPSS H0: trend stationarity.  Reject => evidence of unit root.')
log(' If ADF fails to reject AND KPSS rejects: strong evidence of I(1).')
log(f' ADF lag (in parentheses) chosen by BIC; max lag = {p_max}.')
log(' KPSS bandwidth: Newey-West (1994) automatic selection (statsmodels default).\n')

SEP1 = '-' * 54
log(f'  {"":12s} | {"  ADF: Const (ARD)":20s} | {" KPSS: Trend":16s}')
log(f'  {"Series":12s} | {"Stat(Lag)":10s}  {"p-val":6s}  | {"Stat":7s}  {"p-val":6s}')
log(f'  {SEP1}')

for i in range(n_ser):
    y   = series_lv[i]
    lbl = lbl_lv[i]

    # ADF: BIC lag selection then test with constant only
    p_A      = adf_bic_select(y, p_max)
    adf_stat, adf_pval, *_ = adfuller(y, maxlag=p_A, regression='c', autolag=None)
    strA     = f'{adf_stat:7.3f}({p_A})'

    # KPSS: trend stationarity null, auto bandwidth
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        kpss_stat, kpss_pval, *_ = kpss(y, regression='ct', nlags='auto')

    log(f'  {lbl:12s} | {strA:10s}  {adf_pval:6.4f}  | {kpss_stat:7.4f}  {kpss_pval:6.4f}')

log(f'  {SEP1}\n')
log('  Notes:')
log('  (a) ADF critical values (Dickey-Fuller distribution, ARD model):')
log('      10%: -2.57;  5%: -2.86;  1%: -3.43.')
log('  (b) KPSS critical values (trend stationarity, KPSS 1992, Table 1):')
log('      10%: 0.119;  5%: 0.146;  1%: 0.216.')
log('  (c) KPSS p-values are bounded in [0.01, 0.10] by the tabulated')
log('      critical values.  p <= 0.01 means the stat exceeds the 1% CV.')
log(f'  (d) ADF lag in parentheses; selected to minimise BIC over 0-{p_max}.\n')

# =========================================================================
# TABLE II.  NORMALITY TESTS ON AR(3) RESIDUALS OF FIRST DIFFERENCES
# =========================================================================
log('=================================================================')
log(' Table II.  Normality Tests -- AR(3) Residuals of First Differences')
log('            Model: Dy_t = c + phi_1*Dy_{t-1} + phi_2*Dy_{t-2} + phi_3*Dy_{t-3} + e_t')
log(f'            Residuals: T = {T_fd - 3} obs (1854-2014, after 3-lag alignment)')
log('=================================================================\n')
log(' H0: AR(3) residual is normally distributed.')
log(' JB = (n/6) * [Skewness^2 + (ExcessKurtosis^2)/4]  ~  chi^2(2).')
log(' Critical values: 10% = 4.61,  5% = 5.99,  1% = 9.21.\n')

p_ar = 3
T_ar = T_fd - p_ar   # 161

SEP3 = '-' * 76
log(f'  {"Series":12s}  {"Skewness":>9s}  {"Ex.Kurt":>9s}  {"JB stat":>10s}  '
    f'{"p-value":>8s}  {"Rej1%":>6s}  {"Rej5%":>6s}  {"Rej10%":>6s}')
log(f'  {SEP3}')

jb_ar3_store = np.zeros((n_ser, 4))
ar3_residuals = []   # AR(3) residuals stored for Figures 2 and 3

for i in range(n_ser):
    dy  = series_fd[i]
    lbl = lbl_fd[i]

    Y_dep = dy[p_ar:]
    X_ar  = np.column_stack([
        np.ones(T_ar),
        dy[p_ar - 1: -1],
        dy[p_ar - 2: -2],
        dy[p_ar - 3: -3],
    ])

    b     = np.linalg.lstsq(X_ar, Y_dep, rcond=None)[0]
    resid = Y_dep - X_ar @ b
    ar3_residuals.append(resid)

    sk      = stats.skew(resid)
    ek      = stats.kurtosis(resid, fisher=True)
    jb_stat, jb_pval = jarque_bera(resid)

    chi2_crit = {0.01: 9.21, 0.05: 5.99, 0.10: 4.61}
    r1  = 'Yes' if jb_stat > chi2_crit[0.01] else 'No '
    r5  = 'Yes' if jb_stat > chi2_crit[0.05] else 'No '
    r10 = 'Yes' if jb_stat > chi2_crit[0.10] else 'No '

    jb_ar3_store[i] = [jb_stat, jb_pval, sk, ek]

    log(f'  {lbl:12s}  {sk:9.4f}  {ek:9.4f}  {jb_stat:10.3f}  '
        f'{jb_pval:8.4f}  {r1:6s}  {r5:6s}  {r10:6s}')

log(f'  {SEP3}\n')
log('  Notes:')
log(f'  (a) AR(3) model fitted by OLS; T_ar = {T_ar} residuals (1854-2014).')
log('  (b) Excess kurtosis = sample kurtosis - 3 (normal distribution => 0).')
log('  (c) Removing AR dynamics before JB provides a more stringent test:')
log('      if non-normality persists in e_t, it cannot be attributed to')
log('      residual autocorrelation.\n')

# =========================================================================
# FIGURE 1:  Time series (1 x 2 panel)
# =========================================================================
fig1, axes = plt.subplots(1, 2, figsize=(7.0, 3))

ax = axes[0]
ax.plot(year, f,   'r--', linewidth=1.5, label='Full-efficacy RF')
ax.plot(year, f_m, 'k-',  linewidth=1.5, label='Partial-efficacy RF')
ax.set_xlim(1850, 2014)
ax.set_xticks(range(1850, 2001, 25))
ax.legend(loc='upper left', fontsize=8.5)
ax.set_xlabel('Year', fontsize=10)
ax.set_ylabel('Watts Per Square Meter', fontsize=10)
ax.set_title('(a) Radiative Forcing', fontsize=10)
ax.grid(True)

ax = axes[1]
ax.plot(year, s_h, 'r--', linewidth=1.5, label='HadCRUT')
ax.plot(year, s,   'k-',  linewidth=1.5, label='Berkeley Earth')
ax.set_xlim(1850, 2014)
ax.set_xticks(range(1850, 2001, 25))
ax.legend(loc='upper left', fontsize=8.5)
ax.set_xlabel('Year', fontsize=10)
ax.set_ylabel('Celsius', fontsize=10)
ax.set_title('(b) Surface Temperature', fontsize=10)
ax.grid(True)

fig1.tight_layout()

for ext in ['eps', 'png']:
    fig1.savefig(os.path.join(results_dir, f'rf_temperature_combined.{ext}'),
                 format=ext, dpi=300)

if os.path.isdir(paper_dir):
    for ext in ['eps', 'png']:
        shutil.copy(os.path.join(results_dir, f'rf_temperature_combined.{ext}'),
                    os.path.join(paper_dir,   f'rf_temperature_combined.{ext}'))

plt.close(fig1)

# =========================================================================
# FIGURE 2:  Histograms of AR(3) residuals with fitted normal density (2 x 2 panel)
# =========================================================================
ar3_titles  = ['(a) AR(3) Res.: D.TotalRF',  '(b) AR(3) Res.: D.MarvelRF',
               '(c) AR(3) Res.: D.Berkeley', '(d) AR(3) Res.: D.HadCRUT']
ar3_ylabels = [r'W m$^{-2}$', r'W m$^{-2}$', r'$^\circ$C', r'$^\circ$C']

fig2, axes = plt.subplots(2, 2, figsize=(7.5, 5.5))

for i, ax in enumerate(axes.ravel()):
    resid = ar3_residuals[i]
    mu = resid.mean()
    sg = resid.std(ddof=1)

    ax.hist(resid, bins='auto', density=True,
            color=(0.75, 0.75, 0.75), edgecolor='white')
    xr = np.linspace(resid.min() - 0.5 * sg, resid.max() + 0.5 * sg, 300)
    ax.plot(xr, stats.norm.pdf(xr, mu, sg), 'r-', linewidth=1.5)

    ax.set_title(ar3_titles[i], fontsize=9)
    ax.set_xlabel(ar3_ylabels[i], fontsize=9)
    ax.set_ylabel('Density', fontsize=9)

    ann = f'JB = {jb_ar3_store[i, 0]:.2f}\np = {jb_ar3_store[i, 1]:.4f}'
    ax.text(0.97, 0.97, ann, transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='square,pad=0.3', facecolor='white',
                      edgecolor=(0.7, 0.7, 0.7)))

    ax.legend(['Observed', 'Normal fit'], loc='upper left', fontsize=8)
    ax.grid(True)
    ax.tick_params(labelsize=9)

fig2.tight_layout()
for ext in ['eps', 'png']:
    fig2.savefig(os.path.join(results_dir, f'fd_histograms.{ext}'),
                 format=ext, dpi=300)
plt.close(fig2)

# =========================================================================
# FIGURE 3:  Normal Q-Q plots of AR(3) residuals (2 x 2 panel)
# =========================================================================
qq_titles = ['(a) Q-Q: AR(3) Res. D.TotalRF',  '(b) Q-Q: AR(3) Res. D.MarvelRF',
             '(c) Q-Q: AR(3) Res. D.Berkeley', '(d) Q-Q: AR(3) Res. D.HadCRUT']

fig3, axes = plt.subplots(2, 2, figsize=(7.5, 5.5))

for i, ax in enumerate(axes.ravel()):
    resid = ar3_residuals[i]
    (osm, osr), (slope, intercept, _) = stats.probplot(resid, dist='norm')
    ax.plot(osm, osr,           'o', markersize=3, color='steelblue')
    ax.plot(osm, slope * np.array(osm) + intercept, 'r-', linewidth=1.2)
    ax.set_title(qq_titles[i], fontsize=9)
    ax.set_xlabel('Standard Normal Quantiles', fontsize=8)
    ax.set_ylabel('Sample Quantiles', fontsize=8)
    ax.grid(True)
    ax.tick_params(labelsize=9)

fig3.tight_layout()
for ext in ['eps', 'png']:
    fig3.savefig(os.path.join(results_dir, f'fd_qqplots.{ext}'),
                 format=ext, dpi=300)
plt.close(fig3)

print(f'\nDone. Results saved to: {results_dir}')
print('  unit_root_normality_tests.out')
print('  rf_temperature_combined.eps / .png')
print('  fd_histograms.eps / .png')
print('  fd_qqplots.eps / .png\n')
