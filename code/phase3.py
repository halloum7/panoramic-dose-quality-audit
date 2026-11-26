import os
"""Phase 3: noise power spectrum, regional analysis, effective dose, local typical DAP, sample size."""
import json, sys, os
import numpy as np
from scipy import stats, ndimage
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metrics import load_pixels, immerkaer_sigma

D = os.environ.get('DICOM_DIR', 'dicom')
a = np.load('arrays.npz', allow_pickle=True)
files = json.load(open('files157.json'))
I, DAP, age, sex, g = a['I'], a['DAP'], a['age'], a['sex'], a['g']
out = {}

def ols(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b; n, p = X.shape
    se = np.sqrt(np.diag((e @ e / (n - p)) * np.linalg.inv(X.T @ X)))
    return b, se, n - p
def adj(y):
    X = np.column_stack([np.ones_like(I), I, age, sex]); b, se, df = ols(y, X)
    return dict(beta=b[1], lo=b[1] - 1.96 * se[1], hi=b[1] + 1.96 * se[1], p=2 * stats.t.sf(abs(b[1] / se[1]), df),
                rho=stats.spearmanr(I, y)[0])

# ---------- 1. NPS of the high-pass residual in central region
NB = 32
fr = None
nps_all = []
regions = {'ramus_condyle': [(0.00, 0.15), (0.85, 1.00)], 'molars': [(0.15, 0.33), (0.67, 0.85)],
           'premolars_canines': [(0.33, 0.45), (0.55, 0.67)], 'incisors': [(0.45, 0.55)]}
reg = {k: {'noise': [], 'cnr': []} for k in regions}
for f in files:
    px = load_pixels(os.path.join(D, f)).astype(float)
    H, W = px.shape
    roi = px[int(0.15 * H):int(0.85 * H), int(0.10 * W):int(0.90 * W)]
    # residual after removing anatomy (gaussian sigma 4 px)
    res = roi - ndimage.gaussian_filter(roi, 4)
    # tile 128x128 patches, keep flattest 25% by local gradient
    ps = 128; spectra = []; flat = []
    grad = ndimage.gaussian_gradient_magnitude(ndimage.gaussian_filter(roi, 2), 1)
    for y0 in range(0, res.shape[0] - ps + 1, ps):
        for x0 in range(0, res.shape[1] - ps + 1, ps):
            flat.append((grad[y0:y0 + ps, x0:x0 + ps].mean(), y0, x0))
    flat.sort(); keep = flat[:max(4, len(flat) // 4)]
    win = np.outer(np.hanning(ps), np.hanning(ps))
    for _, y0, x0 in keep:
        p = res[y0:y0 + ps, x0:x0 + ps]; p = (p - p.mean()) * win
        spectra.append(np.abs(np.fft.fftshift(np.fft.fft2(p))) ** 2)
    S = np.mean(spectra, 0) * (0.096 ** 2) / (ps * ps)  # mm^2 scaling (relative units)
    yy, xx = np.indices(S.shape) - ps // 2
    r = np.sqrt(xx ** 2 + yy ** 2) / ps / 0.096  # cycles/mm
    bins = np.linspace(0, 1 / (2 * 0.096), NB + 1)
    idx = np.digitize(r.ravel(), bins) - 1
    prof = np.array([S.ravel()[idx == i].mean() for i in range(NB)])
    fr = 0.5 * (bins[1:] + bins[:-1])
    nps_all.append(prof)
    # regional metrics on rows 15-85 %
    band = px[int(0.15 * H):int(0.85 * H)]
    for k, spans in regions.items():
        parts = [band[:, int(s0 * W):int(s1 * W)] for s0, s1 in spans]
        sig = np.mean([immerkaer_sigma(p) for p in parts])
        con = np.mean([np.percentile(p, 95) - np.percentile(p, 5) for p in parts])
        reg[k]['noise'].append(sig); reg[k]['cnr'].append(con / sig)
nps_all = np.array(nps_all)
lo_f = (fr > 0.3) & (fr < 1.5); hi_f = (fr > 3.0) & (fr < 5.0)
hf_frac = nps_all[:, hi_f].sum(1) / nps_all[:, (fr > 0.3)].sum(1)
total = nps_all[:, lo_f].mean(1); high = nps_all[:, hi_f].mean(1)
out['nps'] = {
    'freq_bands_cyc_per_mm': {'low': [0.3, 1.5], 'high': [3.0, 5.0]},
    'low_band_vs_I': adj(np.log(total)), 'high_band_vs_I': adj(np.log(high)), 'hf_fraction_vs_I': adj(hf_frac),
    'group_median_hf_fraction': {gg: float(np.median(hf_frac[g == gg])) for gg in ['Bas', 'Standard', 'Eleve']},
}
for gg in ['Bas', 'Standard', 'Eleve']:
    out['nps'][f'profile_{gg}'] = np.median(nps_all[g == gg], 0).tolist()
out['nps']['freq'] = fr.tolist()

# ---------- 2. regional analysis
out['regions'] = {}
for k in regions:
    n_ = np.array(reg[k]['noise']); c_ = np.array(reg[k]['cnr'])
    out['regions'][k] = {'noise_median': float(np.median(n_)), 'noise_vs_I': adj(n_), 'cnr_vs_I': adj(c_)}

# ---------- 3. effective dose (E/DAP = 0.08 uSv per mGy cm2; Helmrot & Alm Carlsson 2005 via Lee et al.)
CC = 0.080
E = DAP * CC
out['effective_dose'] = {gg: dict(median=float(np.median(E[g == gg])), min=float(E[g == gg].min()), max=float(E[g == gg].max())) for gg in ['Bas', 'Standard', 'Eleve']}
out['effective_dose']['all_median'] = float(np.median(E))

# ---------- 4. local typical values (median) and 75th percentile of DAP by age class
adult = age >= 18
out['local_dap'] = {
    'adults': dict(n=int(adult.sum()), median=float(np.median(DAP[adult])), p75=float(np.percentile(DAP[adult], 75))),
    'children_lt18': dict(n=int((~adult).sum()), median=float(np.median(DAP[~adult])), p75=float(np.percentile(DAP[~adult], 75))),
}

# ---------- 5. sample size for future observer study (weighted kappa precision & group comparison)
# (a) precision of kappa: n for 95% CI half-width 0.10 at kappa 0.6 (approx. SE ~ sqrt((1-k^2)/n) heuristic)
def n_for_halfwidth(k, hw):
    return int(np.ceil((1.96 ** 2) * (1 - k ** 2) / hw ** 2))
out['sample_size'] = {'kappa_0.6_hw_0.10': n_for_halfwidth(0.6, 0.10), 'kappa_0.6_hw_0.15': n_for_halfwidth(0.6, 0.15)}
# (b) detect difference in proportion acceptable (score<=2) between standard (90%) and low (75%) exposure, alpha .05 power .8, 2-sided
from math import ceil
p1, p2 = 0.90, 0.75; za, zb = 1.96, 0.8416; pbar = (p1 + p2) / 2
n_pg = ceil(((za * np.sqrt(2 * pbar * (1 - pbar)) + zb * np.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2) / (p1 - p2) ** 2)
out['sample_size']['per_group_90_vs_75pct'] = n_pg

json.dump(out, open('phase3.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in out.items() if k != 'nps'}, indent=1, default=lambda x: round(float(x), 4)))
print(json.dumps({k: v for k, v in out['nps'].items() if not k.startswith('profile') and k != 'freq'}, indent=1, default=lambda x: round(float(x), 5)))
