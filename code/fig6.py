import json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r = json.load(open('phase3b.json'))
rows = [('Noise σ, whole region', -0.29, -0.39, -0.19, -0.5),
        ('Noise σ, ramus and condyle', r['noise_ramus_condyle']['slope'], r['noise_ramus_condyle']['lo'], r['noise_ramus_condyle']['hi'], -0.5),
        ('Noise σ, molars', r['noise_molars']['slope'], r['noise_molars']['lo'], r['noise_molars']['hi'], -0.5),
        ('Noise σ, premolars and canines', r['noise_premolars_canines']['slope'], r['noise_premolars_canines']['lo'], r['noise_premolars_canines']['hi'], -0.5),
        ('Noise σ, incisors', r['noise_incisors']['slope'], r['noise_incisors']['lo'], r['noise_incisors']['hi'], -0.5),
        ('Noise power, low band', r['nps_low']['slope'], r['nps_low']['lo'], r['nps_low']['hi'], -1.0),
        ('Noise power, high band', r['nps_high']['slope'], r['nps_high']['lo'], r['nps_high']['hi'], -1.0)]
BLUE, ORANGE, INK, GRID = '#2a78d6', '#eb6834', '#28272a', '#e4e3df'
fig, ax = plt.subplots(figsize=(6.6, 3.9), dpi=300)
for i, (lab, m, lo, hi, q) in enumerate(rows[::-1]):
    ax.plot([lo, hi], [i, i], color=BLUE, lw=2.2, solid_capstyle='round')
    ax.scatter([m], [i], s=44, color=BLUE, edgecolor='white', linewidth=1.4, zorder=3)
    ax.scatter([q], [i], s=46, marker='D', color='white', edgecolor=ORANGE, linewidth=1.6, zorder=3)
    ax.text(0.02, i, f'{m:.2f}', va='center', fontsize=8, color=INK)
ax.axvline(0, color='#9c9a92', lw=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r_[0] for r_ in rows[::-1]], fontsize=8.5)
ax.axhline(1.5, color=GRID, lw=0.8)
ax.set_xlim(-1.1, 0.12); ax.set_xlabel('Log–log slope against the exposure index (95% CI)')
ax.scatter([], [], s=44, color=BLUE, label='Measured slope'); ax.scatter([], [], s=46, marker='D', color='white', edgecolor=ORANGE, linewidth=1.6, label='Quantum-noise expectation')
ax.legend(frameon=False, fontsize=8, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2)
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.tick_params(axis='y', length=0); ax.grid(True, axis='x', color=GRID, lw=0.6, zorder=0)
fig.tight_layout(); fig.savefig('../figures/Fig6_slopes_vs_quantum.png'); plt.close(fig); print('ok')
