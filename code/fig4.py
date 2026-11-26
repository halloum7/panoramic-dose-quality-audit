import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p = json.load(open('phase3.json'))['nps']
f = np.array(p['freq'])
cols = {'Bas': '#2a78d6', 'Standard': '#52514e', 'Eleve': '#eb6834'}
names = {'Bas': 'Low exposure (n=17)', 'Standard': 'Standard (n=110)', 'Eleve': 'High exposure (n=30)'}
fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=300)
m = f > 0.2
for k in ['Bas', 'Standard', 'Eleve']:
    y = np.array(p[f'profile_{k}'])
    ax.plot(f[m], y[m], color=cols[k], lw=2, label=names[k])
ax.set_yscale('log')
ax.set_xlabel('Spatial frequency (cycles/mm)')
ax.set_ylabel('Residual noise power (relative units, log scale)')
ax.axvspan(3.0, 5.0, color='#e4e3df', alpha=0.5, lw=0)
ax.text(4.0, ax.get_ylim()[1] * 0.6, 'high-frequency band', ha='center', va='top', fontsize=8, color='#52514e')
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.grid(True, color='#e4e3df', lw=0.6)
ax.legend(frameon=False, fontsize=8.5, loc='lower left')
fig.tight_layout(); fig.savefig('../figures/Fig5_noise_power_spectrum.png'); plt.close(fig)
print('ok')
