# =========================================================================
# Empirical_Analysis_Climate_trendType_1_in_mc.py
#
# Purpose:
#   Estimates the Equilibrium Climate Sensitivity (ECS) and the
#   multicointegration parameter phi using the Transformed and Augmented OLS
#   (TAOLS) applied to the multicointegration framework of Bruns, Csereklyei,
#   and Stern (JOE 2020).
#
#   Two models are estimated sequentially:
#     (1) Multicointegration model  -- regresses cumulated RF (F) on
#         cumulated temperature (S), temperature (s), and its first
#         difference (ds), and a trend, recovering lambda and phi.
#     (2) Regular cointegration model -- regresses RF (f) on temperature
#         (s), its first difference (ds), and a (lower order) trend,
#         recovering lambda.
#
#   ECS = (5.35 * log(2)) / lambda_hat
#   % surface warming = 100 * 0.31 / phi_hat
#
# Data:
#   'OHC Data.xls', sheet 'Data' -- annual data 1850-2014 (165 obs).
#
# Reference:
#   Bruns, S. B., Csereklyei, Z., & Stern, D. I. (2020). A multicointegration
#   model of global climate change. Journal of Econometrics, 214(1), 175-197.
# =========================================================================

import os
try:
    os.chdir(os.path.dirname(os.path.abspath(__file__)))  # mirror MATLAB's cd(fileparts(mfilename('fullpath')))
except NameError:
    pass  # running interactively; assume CWD is already correct

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')          # non-interactive backend; swap to 'TkAgg' if you want pop-up windows
import matplotlib.pyplot as plt
from scipy.stats import t as t_dist
from scipy.linalg import lstsq, inv

# -------------------------------------------------------------------------
# USER SETTINGS
# -------------------------------------------------------------------------

# 0 = constant only, 1 = constant + linear trend
trend_type_coint  = 0
trend_type_mcoint = trend_type_coint + 1   # automatically one order higher

# 'I', 'II', or 'III'  (Table 3, Bruns et al. 2020)
model = 'III'

# -------------------------------------------------------------------------
# PHYSICAL CONSTANTS
# -------------------------------------------------------------------------
DELTA_F2X = 5.35 * np.log(2)   # radiative forcing from CO2 doubling [W/m^2]
HC_ATMOS  = 0.31                # atmospheric heat capacity per unit surface [W-yr/m^2/degC]

# -------------------------------------------------------------------------
# OUTPUT DIRECTORY & LOG FILE
# -------------------------------------------------------------------------
results_dir = f'Results_Model_{model}_cointTrendType_{trend_type_coint}'
os.makedirs(results_dir, exist_ok=True)

log_path = os.path.join(results_dir, f'Model_{model}_cointTrendType_{trend_type_coint}.out')

def log(msg=''):
    print(msg)
    with open(log_path, 'a') as fh:
        fh.write(msg + '\n')

# Overwrite any previous log
open(log_path, 'w').close()

log(f'coint Trend Type: {trend_type_coint}')
log(f'Mcoint Trend Type: {trend_type_mcoint}')
log(f'Model: {model}')

# -------------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------------
data = pd.read_excel('OHC Data.xls', sheet_name='Data')

f   = data['TotalRF'].values
F   = data['CTotalRF'].values
f_m = data['MarvelRF'].values
F_m = data['CMarvelRF'].values
s   = data['Berkeley'].values
S   = data['CBerkeley'].values
s_h = data['HADCRUT'].values
S_h = data['CHadcrut'].values

# Ocean heat content series (1940-2014)
OHC700    = data['OHCCheng700'].values   # ocean heat content (Cheng et al. 2017)
OHC2000   = data['OHCCheng2000'].values  # ocean heat content (Cheng et al. 2017)
OHCMarvel = data['OHCMarvel'].values     # ocean heat content (Marvel et al. 2016)

start_ohc = 1940 - 1850                 # index 90 (0-based)
OHC700    = OHC700[start_ohc:]
OHC2000   = OHC2000[start_ohc:]
OHCMarvel = OHCMarvel[start_ohc:]
OHCMarvel = OHCMarvel - OHCMarvel[0]    # de-mean from 1940 value

year = np.arange(1940, 2015)            # for OHC plot (75 observations)

# -------------------------------------------------------------------------
# INITIAL OHC FIGURE (TAOLS predictions added inside the loop below)
# -------------------------------------------------------------------------
fig_ohc, ax_ohc = plt.subplots(figsize=(5, 3.5))
ax_ohc.plot(year, OHC2000,   'k-',  linewidth=1.5, markersize=3,
            label='Cheng et al. (2017) OHC 0-2000m')
ax_ohc.plot(year, OHC700,    'ro-', linewidth=1.5, markersize=3,
            label='Cheng et al. (2017) OHC 0-700m')
ax_ohc.plot(year, OHCMarvel, 'b+-', linewidth=1.5, markersize=3,
            label='Marvel et al. (2016)')
ax_ohc.set_ylabel('$10^{22}$ Joules')
ax_ohc.set_xlabel('Year')
ax_ohc.set_xticks(range(1940, 2011, 10))
ax_ohc.grid(True)
ax_ohc.set_axisbelow(True)

# -------------------------------------------------------------------------
# PRE-PROCESS THE MAIN DATA
# -------------------------------------------------------------------------
T_full = len(f)

ds   = s[1:]   - s[:-1]
ds_h = s_h[1:] - s_h[:-1]

# Drop first observation to align with first-differenced series
F   = F[1:];   f   = f[1:]
S   = S[1:];   s   = s[1:]
F_m = F_m[1:]; f_m = f_m[1:]
S_h = S_h[1:]; s_h = s_h[1:]

T = T_full - 1                    # effective sample (164); index 0 = year 1851
t = np.arange(1, T + 1).reshape(-1, 1)   # column vector 1..T

# -------------------------------------------------------------------------
# COSINE BASIS MATRIX
# -------------------------------------------------------------------------
K_end   = T - 10
K_start = 8 + 2 * trend_type_mcoint
K_step  = 2
nK      = (K_end - K_start) // K_step + 1

Phi = np.zeros((T, K_end))
for k in range(1, K_end + 1):
    Phi[:, k - 1] = np.sqrt(2) * np.sin(np.pi * t.ravel() / T * (k - 0.5)) / np.sqrt(T)

# -------------------------------------------------------------------------
# PROJECT SERIES ONTO COSINE BASIS
# -------------------------------------------------------------------------
if model == 'I':
    V_F  = Phi.T @ F.reshape(-1, 1)
    V_S  = Phi.T @ S.reshape(-1, 1)
    V_f  = Phi.T @ f.reshape(-1, 1)
    V_s  = Phi.T @ s.reshape(-1, 1)
    V_ds = Phi.T @ ds.reshape(-1, 1)
    lambda_MLE = 1.709
    phi_MLE    = 41.482
    filename   = f'Full_RF_BerkeleyT_trendType_{trend_type_mcoint}_mc'
    F_use = F;  S_use = S;  s_use = s

elif model == 'II':
    V_F  = Phi.T @ F_m.reshape(-1, 1)
    V_S  = Phi.T @ S.reshape(-1, 1)
    V_f  = Phi.T @ f_m.reshape(-1, 1)
    V_s  = Phi.T @ s.reshape(-1, 1)
    V_ds = Phi.T @ ds.reshape(-1, 1)
    lambda_MLE = 1.326
    phi_MLE    = 33.355
    filename   = f'Partial_RF_BerkeleyT_trendType_{trend_type_mcoint}_mc'
    F_use = F_m;  S_use = S;  s_use = s

else:  # model == 'III'
    V_F  = Phi.T @ F_m.reshape(-1, 1)
    V_S  = Phi.T @ S_h.reshape(-1, 1)
    V_f  = Phi.T @ f_m.reshape(-1, 1)
    V_s  = Phi.T @ s_h.reshape(-1, 1)
    V_ds = Phi.T @ ds_h.reshape(-1, 1)
    lambda_MLE = 1.567
    phi_MLE    = 30.695
    filename   = f'Partial_RF_HADCRUT_trendType_{trend_type_mcoint}_mc'
    F_use = F_m;  S_use = S_h;  s_use = s_h

# Trend regressors for the multicointegration model
if trend_type_mcoint == 1:
    trend = np.hstack([np.ones((T, 1)), t])
else:
    trend = np.hstack([np.ones((T, 1)), t, t ** 2])

V_ell = Phi.T @ trend   # shape: (K_end, 1 + trend_type_mcoint)

# =========================================================================
# ESTIMATION: MULTICOINTEGRATION MODEL
# =========================================================================
d_x = 1
n_params_mc = 3 * d_x + trend_type_mcoint + 1

TA_OLS_mc    = np.zeros((nK, n_params_mc))
TA_OLS_se_mc = np.zeros((nK, n_params_mc))
TA_OLS_cv_mc = np.zeros(nK)

# OHC index: effective sample starts at year 1851 (index 0); 1940 is index 89
ohc_start = 1940 - 1851  # = 89

for idx, k in enumerate(range(K_start, K_end + 1, K_step)):
    V_F_k    = V_F[:k]
    V_S_k    = V_S[:k]
    V_s_k    = V_s[:k]
    V_ds_k   = V_ds[:k]
    V_ell_k  = V_ell[:k]

    reg = np.hstack([V_S_k, V_s_k, V_ds_k, V_ell_k])

    theta_hat, _, _, _ = lstsq(reg, V_F_k)

    residual = V_F_k - reg @ theta_hat
    dof      = k - 3 * d_x - 1 - trend_type_mcoint
    omega_ee = float(residual.T @ residual) / dof

    se = np.sqrt(np.diag(omega_ee * inv(reg.T @ reg)))

    TA_OLS_mc[idx, :]    = theta_hat.ravel()
    TA_OLS_se_mc[idx, :] = se
    TA_OLS_cv_mc[idx]    = t_dist.ppf(0.975, dof)

    if idx == 0:
        assert theta_hat.ravel().shape[0] == n_params_mc, \
            f'theta_hat length {theta_hat.ravel().shape[0]} != n_params_mc {n_params_mc}'

    # --- Predicted system heat content (OHC comparison) ---
    th = theta_hat.ravel()
    Q1 = F_use - th[0] * S_use
    Q2 = th[1] * s_use
    Q  = (Q1 + Q2) / 2
    Q  = Q[ohc_start:]              # align to 1940-2014
    Q  = Q - Q[0]                   # set 1940 value to zero
    Q  = Q * 1.609 * 0.81           # W-yr/m^2 -> 10^22 J, then scale to top-2000m OHC

    if k == K_start:
        ax_ohc.plot(year, Q, ':', color=[1, 0.5, 0], linewidth=1.5,
                    label=f'TAOLS Predictions for K={K_start}:{K_step}:{K_end}')
    else:
        ax_ohc.plot(year, Q, ':', color=[1, 0.5, 0], linewidth=1.5)

# -------------------------------------------------------------------------
# SAVE FIGURE: OHC comparison
# -------------------------------------------------------------------------
ax_ohc.legend(loc='upper left', fontsize=9)
fig_ohc.tight_layout()
fig_ohc.savefig(os.path.join(results_dir, f'{filename}_OHC_comparison.eps'), format='eps', dpi=300)
fig_ohc.savefig(os.path.join(results_dir, f'{filename}_OHC_comparison.pdf'), dpi=300)
plt.close(fig_ohc)

# -------------------------------------------------------------------------
# SUMMARY STATISTICS: multicointegration
# -------------------------------------------------------------------------
log('\n==============================================')
log(' Results for the multicointegration model ')

log('Summary statistics for lambda')
lam = TA_OLS_mc[:, 0]
log(f' Min: {lam.min():.4f},\n Max: {lam.max():.4f},\n Mean: {lam.mean():.4f},\n Median: {np.median(lam):.4f}')
log()

log('Summary statistics for standard errors of lambda')
_se_lam_log = TA_OLS_se_mc[:, 0]
log(f' Min: {_se_lam_log.min():.4f},\n Max: {_se_lam_log.max():.4f},\n Mean: {_se_lam_log.mean():.4f},\n Median: {np.median(_se_lam_log):.4f}')
log()

ECS_mc = DELTA_F2X / TA_OLS_mc[:, 0]
log('Summary statistics for ECS')
log(f' Min: {ECS_mc.min():.4f},\n Max: {ECS_mc.max():.4f},\n Mean: {ECS_mc.mean():.4f},\n Median: {np.median(ECS_mc):.4f}')
log()

log('Summary statistics for phi')
phi = TA_OLS_mc[:, 1]
log(f' Min: {phi.min():.4f},\n Max: {phi.max():.4f},\n Mean: {phi.mean():.4f},\n Median: {np.median(phi):.4f}')
log()

pct = 100 * HC_ATMOS / TA_OLS_mc[:, 1]
log('Summary statistics for percent (100 * HC_ATMOS / phi)')
log(f' Min: {pct.min():.4f},\n Max: {pct.max():.4f},\n Mean: {pct.mean():.4f},\n Median: {np.median(pct):.4f}')

# -------------------------------------------------------------------------
# BUILD PLOT ARRAYS: lambda and ECS (multicointegration)
# -------------------------------------------------------------------------
Ks = np.arange(K_start, K_end + 1, K_step)

cv = TA_OLS_cv_mc
se_lam = TA_OLS_se_mc[:, 0]

TA_OLS_mc_lambda = np.column_stack([
    TA_OLS_mc[:, 0],
    TA_OLS_mc[:, 0] + cv * se_lam,
    TA_OLS_mc[:, 0] - cv * se_lam,
])

se_adj_ECS = DELTA_F2X / (TA_OLS_mc[:, 0] ** 2)
TA_OLS_mc_ECS = np.column_stack([
    DELTA_F2X / TA_OLS_mc[:, 0],
    DELTA_F2X / TA_OLS_mc[:, 0] + cv * se_lam * se_adj_ECS,
    DELTA_F2X / TA_OLS_mc[:, 0] - cv * se_lam * se_adj_ECS,
])

# -------------------------------------------------------------------------
# FIGURE 1: lambda and ECS (multicointegration)
# -------------------------------------------------------------------------
fig1, axes = plt.subplots(1, 2, figsize=(7.5, 3))

ax = axes[0]
ax.fill_between(Ks, TA_OLS_mc_lambda[:, 1], TA_OLS_mc_lambda[:, 2],
                alpha=0.2, color='b', linewidth=0)
ax.plot(Ks, TA_OLS_mc_lambda[:, 0], 'k-o', markevery=5, linewidth=1.5,
        label=r'$\hat{\lambda}_{\rm{TAOLS}}$')
ax.plot(Ks, TA_OLS_mc_lambda[:, 1], 'r:o', markevery=5, linewidth=1.5,
        label='95% CI')
ax.plot(Ks, TA_OLS_mc_lambda[:, 2], 'r:o', markevery=5, linewidth=1.5)
ax.axhline(lambda_MLE, color='b', linestyle=':', linewidth=1.5,
           label=r'$\hat{\lambda}_{\rm{MLE}}$')
ax.legend(fontsize=9, loc='lower right')
ax.set_xlabel('$K$', fontsize=9)
ax.set_title(r'(a) TAOLS of $\lambda$ and 95% CI', fontsize=9)
ax.set_xticks(range(10, 151, 20))
ax.set_yticks(np.arange(1, 2.41, 0.1))
ax.axis([10, 154, 1.0, 2.4])
ax.grid(True); ax.set_box_aspect(None)
ax.tick_params(labelsize=9)

ax = axes[1]
ax.fill_between(Ks, TA_OLS_mc_ECS[:, 1], TA_OLS_mc_ECS[:, 2],
                alpha=0.2, color='b', linewidth=0)
ax.plot(Ks, TA_OLS_mc_ECS[:, 0], 'k-o', markevery=5, linewidth=1.5,
        label='TAOLS of ECS')
ax.plot(Ks, TA_OLS_mc_ECS[:, 1], 'r:o', markevery=5, linewidth=1.5,
        label='95% CI')
ax.plot(Ks, TA_OLS_mc_ECS[:, 2], 'r:o', markevery=5, linewidth=1.5)
ax.axhline(DELTA_F2X / lambda_MLE, color='b', linestyle=':', linewidth=1.5,
           label='MLE of ECS')
ax.legend(fontsize=9, loc='lower right')
ax.set_xlabel('$K$', fontsize=9)
ax.set_title(r'(b) TAOLS of ECS and 95% CI', fontsize=9)
ax.set_xticks(range(10, 151, 20))
ax.set_yticks(np.arange(0.5, 3.51, 0.25))
ax.axis([10, 154, 0.5, 3.5])
ax.grid(True)
ax.tick_params(labelsize=9)

fig1.tight_layout()
fig1.savefig(os.path.join(results_dir, f'{filename}_lambda_ECS.eps'), format='eps', dpi=300)
fig1.savefig(os.path.join(results_dir, f'{filename}_lambda_ECS.pdf'), dpi=300)

np.savez(os.path.join(results_dir, f'{filename}_lambda_ECS.npz'),
         Ks=Ks, TA_OLS_mc_lambda=TA_OLS_mc_lambda, TA_OLS_mc_ECS=TA_OLS_mc_ECS)

plt.close(fig1)

# -------------------------------------------------------------------------
# BUILD PLOT ARRAYS: phi and percent (multicointegration)
# -------------------------------------------------------------------------
se_phi = TA_OLS_se_mc[:, 1]

TA_OLS_mc_phi = np.column_stack([
    TA_OLS_mc[:, 1],
    TA_OLS_mc[:, 1] + cv * se_phi,
    TA_OLS_mc[:, 1] - cv * se_phi,
])

se_adj_pct = HC_ATMOS / (TA_OLS_mc[:, 1] ** 2)
TA_OLS_mc_percent = 100 * np.column_stack([
    HC_ATMOS / TA_OLS_mc[:, 1],
    HC_ATMOS / TA_OLS_mc[:, 1] + cv * se_phi * se_adj_pct,
    HC_ATMOS / TA_OLS_mc[:, 1] - cv * se_phi * se_adj_pct,
])

# -------------------------------------------------------------------------
# FIGURE 2: phi and percent (multicointegration)
# -------------------------------------------------------------------------
fig2, axes = plt.subplots(1, 2, figsize=(7.5, 3))

ax = axes[0]
ax.fill_between(Ks, TA_OLS_mc_phi[:, 1], TA_OLS_mc_phi[:, 2],
                alpha=0.2, color='b', linewidth=0)
ax.plot(Ks, TA_OLS_mc_phi[:, 0], 'k-o', markevery=5, linewidth=1.5,
        label=r'$\hat{\phi}_{\rm{TAOLS}}$')
ax.plot(Ks, TA_OLS_mc_phi[:, 1], 'r:o', markevery=5, linewidth=1.5,
        label='95% CI')
ax.plot(Ks, TA_OLS_mc_phi[:, 2], 'r:o', markevery=5, linewidth=1.5)
ax.axhline(phi_MLE, color='b', linestyle=':', linewidth=1.5,
           label=r'$\hat{\phi}_{\rm{MLE}}$')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlabel('$K$', fontsize=9)
ax.set_title(r'(a) TAOLS of $\phi$ and 95% CI', fontsize=9)
ax.set_xticks(range(10, 151, 20))
ax.axis([10, 154, -5, 52])
ax.grid(True)
ax.tick_params(labelsize=9)

ax = axes[1]
ax.fill_between(Ks, TA_OLS_mc_percent[:, 1], TA_OLS_mc_percent[:, 2],
                alpha=0.2, color='b', linewidth=0)
ax.plot(Ks, TA_OLS_mc_percent[:, 0], 'k-o', markevery=5, linewidth=1.5,
        label='% Heat Content Directed to Warming the Atmosphere')
ax.plot(Ks, TA_OLS_mc_percent[:, 1], 'r:o', markevery=5, linewidth=1.5,
        label='95% CI')
ax.plot(Ks, TA_OLS_mc_percent[:, 2], 'r:o', markevery=5, linewidth=1.5)
ax.axhline(100 * HC_ATMOS / phi_MLE, color='b', linestyle=':', linewidth=1.5,
           label='MLE of %')
ax.legend(fontsize=9, loc='lower right')
ax.set_xlabel('$K$', fontsize=9)
ax.set_title('(b) % of Total Heat Content towards Surface Warming', fontsize=9, x=0.43)
ax.set_xticks(range(10, 151, 20))
ax.set_xlim(10, 154)
ax.grid(True)
ax.tick_params(labelsize=9)

fig2.tight_layout()
fig2.savefig(os.path.join(results_dir, f'{filename}_phi_percent.eps'), format='eps', dpi=300)
fig2.savefig(os.path.join(results_dir, f'{filename}_phi_percent.pdf'), dpi=300)

np.savez(os.path.join(results_dir, f'{filename}_phi_percent.npz'),
         Ks=Ks, TA_OLS_mc_phi=TA_OLS_mc_phi, TA_OLS_mc_percent=TA_OLS_mc_percent)

plt.close(fig2)

# -------------------------------------------------------------------------
# FIGURE 3: trend coefficients (multicointegration)
# -------------------------------------------------------------------------
if trend_type_mcoint == 1:
    ind = 3 * d_x + 1   # 0-based index of linear trend coefficient
    pt    = TA_OLS_mc[:, ind]
    se    = TA_OLS_se_mc[:, ind]
    upper = pt + cv * se
    lower = pt - cv * se

    fig3, ax = plt.subplots(figsize=(5, 3.5))
    ax.fill_between(Ks, upper, lower, alpha=0.2, color='b', linewidth=0)
    ax.plot(Ks, pt,    'k-o', markevery=5, linewidth=1.5,
            label='TAOLS of linear trend coeff')
    ax.plot(Ks, upper, 'r:o', markevery=5, linewidth=1.5, label='95% CI')
    ax.plot(Ks, lower, 'r:o', markevery=5, linewidth=1.5)
    ax.legend(fontsize=9)
    ax.set_xlabel('$K$', fontsize=9)
    ax.set_ylabel('TAOLS', fontsize=9)
    ax.set_title(r'TAOLS of linear trend coefficient and 95% CI', fontsize=9)
    ax.set_xticks(range(10, 151, 20))
    ax.set_xlim(10, 154)
    ax.grid(True)
    ax.tick_params(labelsize=9)
    fig3.tight_layout()
    fig3.savefig(os.path.join(results_dir, f'{filename}_trend.pdf'), dpi=300)
    plt.close(fig3)

elif trend_type_mcoint == 2:
    fig3, axes = plt.subplots(1, 2, figsize=(7.5, 3))

    for i, (label, title) in enumerate([
        ('TAOLS of linear trend coeff',    '(a) Linear trend coefficient'),
        ('TAOLS of quadratic trend coeff', '(b) Quadratic trend coefficient'),
    ]):
        ind   = 3 * d_x + 1 + i
        pt    = TA_OLS_mc[:, ind]
        se    = TA_OLS_se_mc[:, ind]
        upper = pt + cv * se
        lower = pt - cv * se

        ax = axes[i]
        ax.fill_between(Ks, upper, lower, alpha=0.2, color='b', linewidth=0)
        ax.plot(Ks, pt,    'k-o', markevery=5, linewidth=1.5, label=label)
        ax.plot(Ks, upper, 'r:o', markevery=5, linewidth=1.5, label='95% CI')
        ax.plot(Ks, lower, 'r:o', markevery=5, linewidth=1.5)
        ax.legend(fontsize=9)
        ax.set_xlabel('$K$', fontsize=9)
        ax.set_ylabel('TAOLS', fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.set_xticks(range(10, 151, 20))
        ax.set_xlim(10, 154)
        ax.grid(True)
        ax.tick_params(labelsize=9)

    fig3.tight_layout()
    fig3.savefig(os.path.join(results_dir, f'{filename}_trend.pdf'), dpi=300)
    plt.close(fig3)

# =========================================================================
# ESTIMATION: REGULAR COINTEGRATION MODEL
# =========================================================================
if trend_type_coint == 0:
    trend = np.ones((T, 1))
else:
    trend = np.hstack([np.ones((T, 1)), t])

V_ell = Phi.T @ trend

if model == 'I':
    filename_coint = f'Full_RF_BerkeleyT_trendType_{trend_type_coint}_coint'
elif model == 'II':
    filename_coint = f'Partial_RF_BerkeleyT_trendType_{trend_type_coint}_coint'
else:
    filename_coint = f'Partial_RF_HADCRUT_trendType_{trend_type_coint}_coint'

n_params_c = 2 * d_x + trend_type_coint + 1

TA_OLS    = np.zeros((nK, n_params_c))
TA_OLS_se = np.zeros((nK, n_params_c))
TA_OLS_cv = np.zeros(nK)

for idx, k in enumerate(range(K_start, K_end + 1, K_step)):
    V_f_k   = V_f[:k]
    V_s_k   = V_s[:k]
    V_ds_k  = V_ds[:k]
    V_ell_k = V_ell[:k]

    reg = np.hstack([V_s_k, V_ds_k, V_ell_k])

    theta_hat, _, _, _ = lstsq(reg, V_f_k)

    residual = V_f_k - reg @ theta_hat
    dof      = k - 2 * d_x - 1 - trend_type_coint
    omega_ee = float(residual.T @ residual) / dof

    se = np.sqrt(np.diag(omega_ee * inv(reg.T @ reg)))

    TA_OLS[idx, :]    = theta_hat.ravel()
    TA_OLS_se[idx, :] = se
    TA_OLS_cv[idx]    = t_dist.ppf(0.975, dof)

# -------------------------------------------------------------------------
# SUMMARY STATISTICS: cointegration
# -------------------------------------------------------------------------
log('\n==============================================')
log('Results for the regular cointegration model ')

lam_c = TA_OLS[:, 0]
log('Summary statistics for lambda')
log(f' Min: {lam_c.min():.4f},\n Max: {lam_c.max():.4f},\n Mean: {lam_c.mean():.4f},\n Median: {np.median(lam_c):.4f}')
log()

se_lam_c = TA_OLS_se[:, 0]
log('Summary statistics for standard errors of lambda')
log(f' Min: {se_lam_c.min():.4f},\n Max: {se_lam_c.max():.4f},\n Mean: {se_lam_c.mean():.4f},\n Median: {np.median(se_lam_c):.4f}')
log()

ECS_c = DELTA_F2X / TA_OLS[:, 0]
log('Summary statistics for ECS')
log(f' Min: {ECS_c.min():.4f},\n Max: {ECS_c.max():.4f},\n Mean: {ECS_c.mean():.4f},\n Median: {np.median(ECS_c):.4f}')

# -------------------------------------------------------------------------
# BUILD PLOT ARRAYS: lambda and ECS (cointegration)
# -------------------------------------------------------------------------
cv_c = TA_OLS_cv

TA_OLS_lambda = np.column_stack([
    TA_OLS[:, 0],
    TA_OLS[:, 0] + cv_c * se_lam_c,
    TA_OLS[:, 0] - cv_c * se_lam_c,
])

se_adj_ECS_c = DELTA_F2X / (TA_OLS_lambda[:, 0] ** 2)
TA_OLS_ECS = np.column_stack([
    DELTA_F2X / TA_OLS[:, 0],
    DELTA_F2X / TA_OLS[:, 0] + cv_c * se_lam_c * se_adj_ECS_c,
    DELTA_F2X / TA_OLS[:, 0] - cv_c * se_lam_c * se_adj_ECS_c,
])

# -------------------------------------------------------------------------
# FIGURE 4: lambda and ECS (cointegration)
# -------------------------------------------------------------------------
fig4, axes = plt.subplots(1, 2, figsize=(7.5, 3))

ax = axes[0]
ax.fill_between(Ks, TA_OLS_lambda[:, 1], TA_OLS_lambda[:, 2],
                alpha=0.2, color='b', linewidth=0)
ax.plot(Ks, TA_OLS_lambda[:, 0], 'k-o', markevery=5, linewidth=1.5,
        label=r'$\hat{\lambda}_{\rm{TAOLS}}$')
ax.plot(Ks, TA_OLS_lambda[:, 1], 'r:o', markevery=5, linewidth=1.5,
        label='95% CI')
ax.plot(Ks, TA_OLS_lambda[:, 2], 'r:o', markevery=5, linewidth=1.5)
ax.legend(fontsize=9, loc='lower right')
ax.set_xlabel('$K$', fontsize=9)
ax.set_title(r'(a) TAOLS of $\lambda$ and 95% CI', fontsize=9)
ax.set_xticks(range(10, 151, 20))
ax.set_yticks(np.arange(1, 3.31, 0.2))
ax.axis([10, 154, 1.5, 3.3])
ax.grid(True)
ax.tick_params(labelsize=9)

ax = axes[1]
ax.fill_between(Ks, TA_OLS_ECS[:, 1], TA_OLS_ECS[:, 2],
                alpha=0.2, color='b', linewidth=0)
ax.plot(Ks, TA_OLS_ECS[:, 0], 'k-o', markevery=5, linewidth=1.5,
        label='TAOLS of ECS')
ax.plot(Ks, TA_OLS_ECS[:, 1], 'r:o', markevery=5, linewidth=1.5,
        label='95% CI')
ax.plot(Ks, TA_OLS_ECS[:, 2], 'r:o', markevery=5, linewidth=1.5)
ax.legend(fontsize=9, loc='lower right')
ax.set_xlabel('$K$', fontsize=9)
ax.set_title(r'(b) TAOLS of ECS and 95% CI', fontsize=9)
ax.set_xticks(range(10, 151, 20))
ax.set_yticks(np.arange(0.25, 3.01, 0.25))
ax.axis([10, 154, 0.25, 3])
ax.grid(True)
ax.tick_params(labelsize=9)

fig4.tight_layout()
fig4.savefig(os.path.join(results_dir, f'{filename_coint}_lambda_ECS.eps'), format='eps', dpi=300)
fig4.savefig(os.path.join(results_dir, f'{filename_coint}_lambda_ECS.pdf'), dpi=300)

np.savez(os.path.join(results_dir, f'{filename_coint}_lambda_ECS.npz'),
         Ks=Ks, TA_OLS_lambda=TA_OLS_lambda, TA_OLS_ECS=TA_OLS_ECS)

plt.close(fig4)

# -------------------------------------------------------------------------
# FIGURE 5: linear trend coefficient (cointegration, only when trend=1)
# -------------------------------------------------------------------------
if trend_type_coint == 1:
    ind   = 2 * d_x + 1   # 0-based index of linear trend coeff in TA_OLS
    pt    = TA_OLS[:, ind]
    se    = TA_OLS_se[:, ind]
    upper = pt + cv_c * se
    lower = pt - cv_c * se

    fig5, ax = plt.subplots(figsize=(5, 3.5))
    ax.fill_between(Ks, upper, lower, alpha=0.2, color='b', linewidth=0)
    ax.plot(Ks, pt,    'k-o', markevery=5, linewidth=1.5,
            label='TAOLS of linear trend coeff')
    ax.plot(Ks, upper, 'r:o', markevery=5, linewidth=1.5, label='95% CI')
    ax.plot(Ks, lower, 'r:o', markevery=5, linewidth=1.5)
    ax.legend(fontsize=9)
    ax.set_xlabel('$K$', fontsize=9)
    ax.set_ylabel('TAOLS', fontsize=9)
    ax.set_title(r'TAOLS of linear trend coefficient and 95% CI', fontsize=9)
    ax.set_xticks(range(10, 151, 20))
    ax.set_xlim(10, 154)
    ax.grid(True)
    ax.tick_params(labelsize=9)
    fig5.tight_layout()
    fig5.savefig(os.path.join(results_dir, f'{filename_coint}_trend.pdf'), dpi=300)
    plt.close(fig5)

log('\nDone.')
