import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p = json.load(open('phase2.json'))
BLUE, INK, INK2, GRID = '#2a78d6', '#0b0b0b', '#52514e', '#e4e3df'
rows = [r for r in p['fom'] if r['n'] >= 3]
fig, ax = plt.subplots(figsize=(7.2, 3.8), dpi=300)
y = np.arange(len(rows))[::-1]
for yi, r in zip(y, rows):
    ax.plot([r['fom_rel_lo'], r['fom_rel_hi']], [yi, yi], color=GRID, lw=5, solid_capstyle='round', zorder=1)
    ax.scatter(r['fom_rel'], yi, s=46, color=BLUE, edgecolor='white', linewidth=1.5, zorder=3)
    ax.text(r['fom_rel'], yi + 0.25, f"{r['fom_rel']:.2f}", ha='center', va='bottom', fontsize=8, color=INK)
labels = []
for r in rows:
    lab = r['protocol'].replace('72 kV / 4 mA / 14.7 s', '72 kV / 4 mA / 14.7 s (standard)')
    if r['protocol'] == '72 kV / 7 mA / 14.7 s':
        lab += ' *'
    labels.append(f"{lab}  — DAP {r['dap']:.0f}  (n={r['n']})")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
ax.axvline(1.0, color=INK2, lw=1, ls='--', zorder=0)
ax.set_xlabel('Dose efficiency CNR\u00b2/DAP (standard protocol = 1)')
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.grid(True, axis='x', color=GRID, lw=0.6, zorder=0)
ax.set_xlim(0.5, 2.0)
fig.tight_layout(); fig.savefig('../figures/Fig4_dose_efficiency.png'); plt.close(fig)
print('ok')
