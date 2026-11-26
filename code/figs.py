import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

a = np.load('arrays.npz', allow_pickle=True)
res = json.load(open('robust.json'))
BLUE, ORANGE, INK, INK2, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.edgecolor': INK2,
                     'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2})

# Figure 1: log noise vs log exposure index
I, noise = a['I'], a['noise']
fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=300)
jit = np.exp(np.random.default_rng(1).normal(0, 0.012, len(I)))
ax.scatter(I * jit, noise, s=22, color=BLUE, alpha=0.55, edgecolor='white', linewidth=0.6, zorder=3)
xs = np.linspace(0.72, 1.8, 50)
b = np.polyfit(np.log(I), np.log(noise), 1)
ax.plot(xs, np.exp(b[1]) * xs ** b[0], color=BLUE, lw=2, zorder=4)
ref = np.exp(b[1]) * xs ** -0.5
ax.plot(xs, ref, color=INK2, lw=1.4, ls='--', zorder=2)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xticks([0.75, 1.0, 1.25, 1.5, 1.75]); ax.set_xticklabels(['0.75', '1.00', '1.25', '1.50', '1.75'])
ax.set_yticks([0.7, 0.8, 0.9, 1.0, 1.2, 1.4]); ax.set_yticklabels(['0.7', '0.8', '0.9', '1.0', '1.2', '1.4'])
ax.minorticks_off()
ax.set_xlabel('Relative exposure index I (standard protocol = 1)')
ax.set_ylabel('Noise σ (grey levels, log scale)')
ph = res['physics']
ax.text(1.78, 1.36, f"measured slope {ph['slope']:.2f}\n(95% CI {ph['lo']:.2f} to {ph['hi']:.2f})",
        color=BLUE, ha='right', va='top', fontsize=9)
ax.text(1.78, np.exp(b[1]) * 1.78 ** -0.5 * 0.96, 'quantum-noise limit, slope −0.50', color=INK2, ha='right', va='top', fontsize=9)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.grid(True, color=GRID, lw=0.6, zorder=0)
fig.tight_layout(); fig.savefig('../figures/Fig2_noise_vs_exposure.png'); plt.close(fig)

# Figure 2: estimated DAP per protocol, with published adult reference values
prot = [('64 kV / 4 mA / 13.9 s', 77.7, 12), ('66 kV / 4 mA / 13.9 s', 82.7, 5), ('72 kV / 4 mA / 14.7 s (standard)', 104.0, 105),
        ('68 kV / 5 mA / 13.9 s', 109.7, 5), ('72 kV / 5 mA / 14.7 s', 130.1, 1), ('68 kV / 6 mA / 14.7 s', 139.2, 11),
        ('72 kV / 6 mA / 14.7 s', 156.1, 1), ('74 kV / 6 mA / 14.7 s', 164.9, 13), ('72 kV / 7 mA / 14.7 s', 182.1, 4)]
fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=300)
y = np.arange(len(prot))[::-1]
vals = np.array([p[1] for p in prot])
ax.errorbar(vals, y, xerr=[vals - vals * 0.6, vals * 1.6 - vals], fmt='none', ecolor=GRID, elinewidth=5, capsize=0, zorder=1)
ax.scatter(vals, y, s=46, color=BLUE, edgecolor='white', linewidth=1.5, zorder=3)
for yi, (lab, v, n) in zip(y, prot):
    ax.text(v, yi + 0.28, f'{v:.0f}', ha='center', va='bottom', fontsize=8, color=INK)
ax.set_yticks(y); ax.set_yticklabels([f'{p[0]}  (n={p[2]})' for p in prot], fontsize=8.5)
ax.set_xlim(30, 310); ax.set_ylim(-0.7, len(prot) + 0.2)
ax.set_xlabel('Estimated DAP (mGy\u00b7cm\u00b2)')
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.grid(True, axis='x', color=GRID, lw=0.6, zorder=0)
fig.tight_layout(); fig.savefig('../figures/Fig3_DAP_by_protocol.png'); plt.close(fig)
print('ok')
