# extract_ranges_and_combine.py
#
# Run from Empirical_application_python\ to:
#   1. Extract TAOLS numerical ranges for Models I, II, III
#   2. Generate combined lambda/ECS and phi/percent figures
#   3. Export combined figures as PNG to ..\Paper_June_2026\
#
# Prerequisite: Empirical_Analysis_Climate_trendType_1_in_mc.py must have
#   been run with model = 'I', 'II', and 'III' (producing all six .npz files).

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root  = os.path.dirname(script_dir)

dir_I   = os.path.join(script_dir, 'Results_Model_I_cointTrendType_0')
dir_II  = os.path.join(script_dir, 'Results_Model_II_cointTrendType_0')
dir_III = os.path.join(script_dir, 'Results_Model_III_cointTrendType_0')
out_dir = os.path.join(repo_root, 'Paper_June_2026')
os.makedirs(out_dir, exist_ok=True)

# -------------------------------------------------------------------------
# 1. Load lambda/ECS npz files
# -------------------------------------------------------------------------
d = np.load(os.path.join(dir_I,   'Full_RF_BerkeleyT_trendType_1_mc_lambda_ECS.npz'))
Ks      = d['Ks']
lambda1 = d['TA_OLS_mc_lambda']   # N x 3: [estimate, upper CI, lower CI]
ECS1    = d['TA_OLS_mc_ECS']

d = np.load(os.path.join(dir_II,  'Partial_RF_BerkeleyT_trendType_1_mc_lambda_ECS.npz'))
lambda2 = d['TA_OLS_mc_lambda']
ECS2    = d['TA_OLS_mc_ECS']

d = np.load(os.path.join(dir_III, 'Partial_RF_HADCRUT_trendType_1_mc_lambda_ECS.npz'))
lambda3 = d['TA_OLS_mc_lambda']
ECS3    = d['TA_OLS_mc_ECS']

# -------------------------------------------------------------------------
# 2. Load phi/percent npz files
# -------------------------------------------------------------------------
d = np.load(os.path.join(dir_I,   'Full_RF_BerkeleyT_trendType_1_mc_phi_percent.npz'))
phi1 = d['TA_OLS_mc_phi']
pct1 = d['TA_OLS_mc_percent']

d = np.load(os.path.join(dir_II,  'Partial_RF_BerkeleyT_trendType_1_mc_phi_percent.npz'))
phi2 = d['TA_OLS_mc_phi']
pct2 = d['TA_OLS_mc_percent']

d = np.load(os.path.join(dir_III, 'Partial_RF_HADCRUT_trendType_1_mc_phi_percent.npz'))
phi3 = d['TA_OLS_mc_phi']
pct3 = d['TA_OLS_mc_percent']

# -------------------------------------------------------------------------
# 3. Print numerical ranges
# -------------------------------------------------------------------------
print('\n=== TAOLS Numerical Ranges ===')

print('\n--- Model I (Full RF, Berkeley Earth) ---')
print(f'lambda: [{lambda1[:,0].min():.3f}, {lambda1[:,0].max():.3f}],  mean = {lambda1[:,0].mean():.3f}')
print(f'ECS:    [{ECS1[:,0].min():.2f}, {ECS1[:,0].max():.2f}] C,  mean = {ECS1[:,0].mean():.2f} C')
print(f'phi:    [{phi1[:,0].min():.2f}, {phi1[:,0].max():.2f}] W-yr/m2')
print(f'pct:    [{pct1[:,0].min():.2f}, {pct1[:,0].max():.2f}] %')

print('\n--- Model II (Partial RF, Berkeley Earth) ---')
print(f'lambda: [{lambda2[:,0].min():.3f}, {lambda2[:,0].max():.3f}],  mean = {lambda2[:,0].mean():.3f}')
print(f'ECS:    [{ECS2[:,0].min():.2f}, {ECS2[:,0].max():.2f}] C,  mean = {ECS2[:,0].mean():.2f} C')
print(f'phi:    [{phi2[:,0].min():.2f}, {phi2[:,0].max():.2f}] W-yr/m2')
print(f'pct:    [{pct2[:,0].min():.2f}, {pct2[:,0].max():.2f}] %')

print('\n--- Model III (Partial RF, HadCRUT) ---')
print(f'lambda: [{lambda3[:,0].min():.3f}, {lambda3[:,0].max():.3f}],  mean = {lambda3[:,0].mean():.3f}')
print(f'ECS:    [{ECS3[:,0].min():.2f}, {ECS3[:,0].max():.2f}] C,  mean = {ECS3[:,0].mean():.2f} C')
print(f'phi:    [{phi3[:,0].min():.2f}, {phi3[:,0].max():.2f}] W-yr/m2')
print(f'pct:    [{pct3[:,0].min():.2f}, {pct3[:,0].max():.2f}] %')

print('\n--- Cross-model summary ---')
all_lambda = np.concatenate([lambda1[:,0], lambda2[:,0], lambda3[:,0]])
all_ECS    = np.concatenate([ECS1[:,0],    ECS2[:,0],    ECS3[:,0]])
all_phi    = np.concatenate([phi1[:,0],    phi2[:,0],    phi3[:,0]])
all_pct    = np.concatenate([pct1[:,0],    pct2[:,0],    pct3[:,0]])
print(f'lambda (all models): [{all_lambda.min():.3f}, {all_lambda.max():.3f}]')
print(f'ECS    (all models): [{all_ECS.min():.2f}, {all_ECS.max():.2f}] C')
print(f'phi    (all models): [{all_phi.min():.2f}, {all_phi.max():.2f}] W-yr/m2')
print(f'pct    (all models): [{all_pct.min():.2f}, {all_pct.max():.2f}] %')

# -------------------------------------------------------------------------
# Colour palette
# -------------------------------------------------------------------------
C1 = (0.8, 0.2, 0.2)   # Model I  -- red
C2 = (0.2, 0.4, 0.8)   # Model II -- blue
C3 = (0.2, 0.7, 0.2)   # Model III -- green

def fill_ci(ax, Ks, arr, color):
    ax.fill_between(Ks, arr[:, 1], arr[:, 2], color=color, alpha=0.18, linewidth=0)

# -------------------------------------------------------------------------
# 4. Combined lambda/ECS figure
# -------------------------------------------------------------------------
fig1, axes = plt.subplots(1, 2, figsize=(9, 3.5))

# --- Panel (a): lambda ---
ax = axes[0]
fill_ci(ax, Ks, lambda1, C1)
fill_ci(ax, Ks, lambda2, C2)
fill_ci(ax, Ks, lambda3, C3)
p1, = ax.plot(Ks, lambda1[:, 0], '-',  color=C1, linewidth=1.5, label='Model I')
p2, = ax.plot(Ks, lambda2[:, 0], '-s', color=C2, linewidth=1.5, markevery=5,
              markersize=4, label='Model II')
p3, = ax.plot(Ks, lambda3[:, 0], '-o', color=C3, linewidth=1.5, markevery=5,
              markersize=4, label='Model III')
ax.legend(fontsize=8, loc='lower right')
ax.set_xlabel('$K$', fontsize=10)
ax.set_ylabel(r'$\hat{\lambda}$', fontsize=10)
ax.set_title('(a) TAOLS of $\\lambda$ and 95% CI', fontsize=10)
ax.set_xticks(range(10, 151, 20))
ax.set_xlim(10, 154)
ax.grid(True)

# --- Panel (b): ECS ---
ax = axes[1]
fill_ci(ax, Ks, ECS1, C1)
fill_ci(ax, Ks, ECS2, C2)
fill_ci(ax, Ks, ECS3, C3)
q1, = ax.plot(Ks, ECS1[:, 0], '-',  color=C1, linewidth=1.5, label='Model I')
q2, = ax.plot(Ks, ECS2[:, 0], '-s', color=C2, linewidth=1.5, markevery=5,
              markersize=4, label='Model II')
q3, = ax.plot(Ks, ECS3[:, 0], '-o', color=C3, linewidth=1.5, markevery=5,
              markersize=4, label='Model III')
ax.axhline(2.80, color='k', linestyle=':', linewidth=1.2)  # MLE ECS for Model II: DELTA_F2X / 1.326 ≈ 2.80 (Bruns et al. 2020, Table 3)
ax.legend(fontsize=8, loc='upper right')
ax.set_xlabel('$K$', fontsize=10)
ax.set_ylabel(r'ECS ($^\circ$C)', fontsize=10)
ax.set_title('(b) TAOLS of ECS and 95% CI', fontsize=10)
ax.set_xticks(range(10, 151, 20))
ax.set_xlim(10, 154)
ax.grid(True)

fig1.tight_layout()
outfile1 = os.path.join(out_dir, 'combined_mc_lambda_ECS.png')
fig1.savefig(outfile1, dpi=300)
print(f'\nSaved: {outfile1}')
plt.close(fig1)

# -------------------------------------------------------------------------
# 5. Combined phi/percent figure
# -------------------------------------------------------------------------
fig2, axes = plt.subplots(1, 2, figsize=(9, 3.5))

# --- Panel (a): phi ---
ax = axes[0]
fill_ci(ax, Ks, phi1, C1)
fill_ci(ax, Ks, phi2, C2)
fill_ci(ax, Ks, phi3, C3)
r1, = ax.plot(Ks, phi1[:, 0], '-',  color=C1, linewidth=1.5, label='Model I')
r2, = ax.plot(Ks, phi2[:, 0], '-s', color=C2, linewidth=1.5, markevery=5,
              markersize=4, label='Model II')
r3, = ax.plot(Ks, phi3[:, 0], '-o', color=C3, linewidth=1.5, markevery=5,
              markersize=4, label='Model III')
ax.legend(fontsize=8, loc='upper right')
ax.set_xlabel('$K$', fontsize=10)
ax.set_ylabel(r'$\hat{\phi}$ (W-yr/m$^2$)', fontsize=10)
ax.set_title('(a) TAOLS of $\\phi$ and 95% CI', fontsize=10)
ax.set_xticks(range(10, 151, 20))
ax.set_xlim(10, 154)
ax.grid(True)

# --- Panel (b): percent ---
ax = axes[1]
fill_ci(ax, Ks, pct1, C1)
fill_ci(ax, Ks, pct2, C2)
fill_ci(ax, Ks, pct3, C3)
s1, = ax.plot(Ks, pct1[:, 0], '-',  color=C1, linewidth=1.5, label='Model I')
s2, = ax.plot(Ks, pct2[:, 0], '-s', color=C2, linewidth=1.5, markevery=5,
              markersize=4, label='Model II')
s3, = ax.plot(Ks, pct3[:, 0], '-o', color=C3, linewidth=1.5, markevery=5,
              markersize=4, label='Model III')
ax.legend(fontsize=8, loc='upper right')
ax.set_xlabel('$K$', fontsize=10)
ax.set_ylabel('% of Total Heat Content', fontsize=10)
ax.set_title('(b) % of Total Heat Content toward Surface Warming', fontsize=10)
ax.set_xticks(range(10, 151, 20))
ax.set_xlim(10, 154)
ax.grid(True)

fig2.tight_layout()
outfile2 = os.path.join(out_dir, 'combined_mc_phi_percent.png')
fig2.savefig(outfile2, dpi=300)
print(f'Saved: {outfile2}')
plt.close(fig2)

print('\nDone.')
