import os
"""Robustness, sensitivity and physics analyses for the panoramic exposure-quality study."""
import json, sys, os
import numpy as np
from scipy import stats, ndimage
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metrics import load_pixels, immerkaer_sigma

D = os.environ.get('DICOM_DIR', 'dicom')
rng = np.random.default_rng(20261121)

img = {x[0]: x for x in json.load(open('img.json'))}
img['99.dcm'] = ['99.dcm', 'P155', '2026', 44, 'M', 72, 4, 14.7]
img['74.dcm'][4] = 'M'
files = json.load(open('files157.json'))
met = json.load(open('met.json'))

kv = np.array([img[f][5] for f in files], float)
ma = np.array([img[f][6] for f in files], float)
t = np.array([img[f][7] for f in files], float)
age = np.array([img[f][3] for f in files], float)
sex = np.array([1.0 if img[f][4] == 'M' else 0.0 for f in files])
pid = np.array([img[f][1] for f in files])
I = (kv / 72) ** 2 * (ma / 4) * (t / 14.7)
DAP = 0.05 * I * 4 * 14.7 * 35.38944  # mGy cm2, central assumption

keys = ['mean', 'noise_sigma', 'contrast_p95_p5', 'snr', 'cnr', 'sharpness_tenengrad', 'entropy', 'sat_white_pct', 'sat_black_pct']
V = {k: np.array([met[f][k] for f in files]) for k in keys}

def ols(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    n, p = X.shape
    s2 = e @ e / (n - p)
    cov = s2 * np.linalg.inv(X.T @ X)
    return b, np.sqrt(np.diag(cov)), n - p

def adj(y, x, mask=None):
    m = np.ones_like(y, bool) if mask is None else mask
    X = np.column_stack([np.ones(m.sum()), x[m], age[m], sex[m]])
    b, se, df = ols(y[m], X)
    tval = b[1] / se[1]
    return b[1], b[1] - 1.96 * se[1], b[1] + 1.96 * se[1], 2 * stats.t.sf(abs(tval), df)

def boot_rho(x, y, n=2000):
    idx = rng.integers(0, len(x), (n, len(x)))
    r = [stats.spearmanr(x[i], y[i])[0] for i in idx]
    return np.percentile(r, [2.5, 97.5])

def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p)
    out = np.empty(m); run = 0
    for r, i in enumerate(o):
        run = max(run, min(1, (m - r) * p[i])); out[i] = run
    return out

res = {'main': {}}
pvals = []
for k in keys:
    rho, p = stats.spearmanr(I, V[k])
    lo, hi = boot_rho(I, V[k])
    b, bl, bh, pa = adj(V[k], I)
    res['main'][k] = dict(rho=rho, rho_lo=lo, rho_hi=hi, beta=b, beta_lo=bl, beta_hi=bh, p_adj=pa)
    pvals.append(pa)
for k, ph in zip(keys, holm(pvals)):
    res['main'][k]['p_holm'] = ph

# one image per patient (keep first occurrence)
seen = set(); keep = []
for i, p in enumerate(pid):
    keep.append(p not in seen); seen.add(p)
keep = np.array(keep)

# sensitivity subsets
subsets = {
    'one_per_patient': keep,
    'adults_only': age >= 18,
    'program_14_7s_only': t > 14.5,
}
res['sens'] = {}
for name, m in subsets.items():
    res['sens'][name] = {'n': int(m.sum())}
    for k in ['noise_sigma', 'snr', 'cnr', 'mean']:
        rho, p = stats.spearmanr(I[m], V[k][m])
        b, bl, bh, pa = adj(V[k], I, m)
        res['sens'][name][k] = dict(rho=rho, beta=b, beta_lo=bl, beta_hi=bh, p_adj=pa)

# alternative noise estimators and ROI on raw pixels
def wavelet_mad_sigma(img):
    # Donoho: sigma = median(|HH|)/0.6745 using Haar diagonal detail
    a = img[0:img.shape[0] // 2 * 2, 0:img.shape[1] // 2 * 2].astype(float)
    hh = (a[0::2, 0::2] - a[0::2, 1::2] - a[1::2, 0::2] + a[1::2, 1::2]) / 2.0
    return np.median(np.abs(hh)) / 0.6745

def flat_local_std(img, w=7, q=10):
    f = img.astype(float)
    m = ndimage.uniform_filter(f, w); m2 = ndimage.uniform_filter(f * f, w)
    sd = np.sqrt(np.clip(m2 - m * m, 0, None))
    g = ndimage.gaussian_gradient_magnitude(ndimage.gaussian_filter(f, 2), 1)
    mask = (g <= np.percentile(g, q)) & (f > 20) & (f < 235)
    return np.median(sd[mask])

alt = {'wavelet_mad': [], 'flat_local_std': [], 'immerkaer_roi2': []}
for f in files:
    px = load_pixels(os.path.join(D, f))
    H, W = px.shape
    roi = px[int(0.15 * H):int(0.85 * H), int(0.10 * W):int(0.90 * W)]
    roi2 = px[int(0.25 * H):int(0.75 * H), int(0.20 * W):int(0.80 * W)]
    alt['wavelet_mad'].append(wavelet_mad_sigma(roi))
    alt['flat_local_std'].append(flat_local_std(roi))
    alt['immerkaer_roi2'].append(immerkaer_sigma(roi2.astype(float)))
res['alt_noise'] = {}
for k, v in alt.items():
    v = np.array(v)
    rho, p = stats.spearmanr(I, v)
    b, bl, bh, pa = adj(v, I)
    res['alt_noise'][k] = dict(rho=rho, p=p, beta=b, beta_lo=bl, beta_hi=bh, p_adj=pa,
                               corr_with_immerkaer=stats.spearmanr(v, V['noise_sigma'])[0])

# physics: log(noise) vs log(I); quantum-limited expectation slope -0.5
X = np.column_stack([np.ones(len(I)), np.log(I), age, sex])
b, se, df = ols(np.log(V['noise_sigma']), X)
res['physics'] = dict(slope=b[1], lo=b[1] - 1.96 * se[1], hi=b[1] + 1.96 * se[1],
                      p_vs_0=2 * stats.t.sf(abs(b[1] / se[1]), df),
                      p_vs_minus05=2 * stats.t.sf(abs((b[1] + 0.5) / se[1]), df))
b2, se2, _ = ols(np.log(V['noise_sigma']), np.column_stack([np.ones(len(I)), np.log(I)]))
res['physics']['slope_unadjusted'] = b2[1]

# DAP per group
g = np.where(I < 0.9, 'Bas', np.where(I < 1.2, 'Standard', 'Eleve'))
res['dap_groups'] = {gg: dict(n=int((g == gg).sum()), median=float(np.median(DAP[g == gg])),
                              min=float(DAP[g == gg].min()), max=float(DAP[g == gg].max())) for gg in ['Bas', 'Standard', 'Eleve']}
res['age_by_group'] = {gg: float(np.median(age[g == gg])) for gg in ['Bas', 'Standard', 'Eleve']}

np.savez('arrays.npz', I=I, DAP=DAP, age=age, sex=sex, noise=V['noise_sigma'], cnr=V['cnr'], snr=V['snr'], mean=V['mean'], g=g)
json.dump(res, open('robust.json', 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=lambda x: round(float(x), 4)))
