"""Phase 2: dose-efficiency FOM, kV vs mAs model, Monte-Carlo DAP uncertainty."""
import json, numpy as np
from scipy import stats

a = np.load('arrays.npz', allow_pickle=True)
img = {x[0]: x for x in json.load(open('img.json'))}
img['99.dcm'] = ['99.dcm', 'P155', '2026', 44, 'M', 72, 4, 14.7]
files = json.load(open('files157.json'))
met = json.load(open('met.json'))
kv = np.array([img[f][5] for f in files], float)
ma = np.array([img[f][6] for f in files], float)
t = np.array([img[f][7] for f in files], float)
age, sex = a['age'], a['sex']
noise, cnr, snr = a['noise'], a['cnr'], a['snr']
contrast = np.array([met[f]['contrast_p95_p5'] for f in files])
mas = ma * t
rng = np.random.default_rng(7)

S1M = 0.6 * 14.7456 * (100 / 50) ** 2
def dap(kv_, mas_, Y=0.05, n=2.0, S=S1M):
    return Y * (kv_ / 72) ** n * mas_ * S
DAP = dap(kv, mas)
out = {}

# ---- 1. FOM = CNR^2 / DAP per protocol
prot = [f"{int(k)} kV / {int(m)} mA / {tt:.1f} s" for k, m, tt in zip(kv, ma, t)]
prot = np.array(prot)
fom = cnr ** 2 / DAP
ref = np.median(fom[prot == '72 kV / 4 mA / 14.7 s'])
rows = []
for p in sorted(set(prot), key=lambda s: dap(float(s.split()[0]), float(s.split('/')[1].split()[0]) * float(s.split('/')[2].split()[0]))):
    m = prot == p
    f = fom[m]
    if m.sum() >= 3:
        bs = [np.median(rng.choice(f, m.sum())) for _ in range(2000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
    else:
        lo = hi = np.nan
    rows.append(dict(protocol=p, n=int(m.sum()), dap=float(DAP[m][0]), cnr_med=float(np.median(cnr[m])),
                     fom_med=float(np.median(f)), fom_rel=float(np.median(f) / ref),
                     fom_rel_lo=float(lo / ref) if m.sum() >= 3 else None, fom_rel_hi=float(hi / ref) if m.sum() >= 3 else None))
out['fom'] = rows
# FOM vs exposure index: adjusted regression of log FOM on log I
I = a['I']
X = np.column_stack([np.ones_like(I), np.log(I), age, sex])
def ols(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b; n, p = X.shape
    cov = (e @ e / (n - p)) * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
    return b, se, n - p
b, se, df = ols(np.log(fom), X)
out['fom_slope'] = dict(slope=b[1], lo=b[1] - 1.96 * se[1], hi=b[1] + 1.96 * se[1], p=2 * stats.t.sf(abs(b[1] / se[1]), df))
rho, p = stats.spearmanr(I, fom); out['fom_slope'].update(rho=rho, p_rho=p)

# ---- 2. kV vs mAs model on log scale
Xk = np.column_stack([np.ones_like(kv), np.log(kv), np.log(mas), age, sex])
out['collinearity_r_logkv_logmas'] = float(np.corrcoef(np.log(kv), np.log(mas))[0, 1])
vif = 1 / (1 - out['collinearity_r_logkv_logmas'] ** 2)
out['vif'] = vif
out['kv_mas'] = {}
for name, y in [('noise', noise), ('cnr', cnr), ('contrast', contrast), ('snr', snr)]:
    b, se, df = ols(np.log(y), Xk)
    out['kv_mas'][name] = {
        'b_logkV': b[1], 'lo_kV': b[1] - 1.96 * se[1], 'hi_kV': b[1] + 1.96 * se[1], 'p_kV': 2 * stats.t.sf(abs(b[1] / se[1]), df),
        'b_logmAs': b[2], 'lo_mAs': b[2] - 1.96 * se[2], 'hi_mAs': b[2] + 1.96 * se[2], 'p_mAs': 2 * stats.t.sf(abs(b[2] / se[2]), df)}
out['n_unique_kv'] = sorted(set(kv.tolist())); out['n_unique_mas'] = sorted(set(np.round(mas, 1).tolist()))

# ---- 3. Monte-Carlo DAP uncertainty
N = 100000
Y = np.exp(rng.uniform(np.log(0.03), np.log(0.08), N))       # tube output mGy/mAs @1m, 72 kV (log-uniform)
width = rng.uniform(0.5, 0.7, N)                                # beam width at detector, cm
sid = rng.uniform(45, 55, N)                                    # focus-detector distance, cm
n_exp = rng.uniform(1.8, 2.5, N)                                # kV exponent
height = 14.7456
mc = []
for r in rows:
    k = float(r['protocol'].split()[0]); m_ = float(r['protocol'].split('/')[1].split()[0]); tt = float(r['protocol'].split('/')[2].split()[0])
    d = Y * (k / 72) ** n_exp * m_ * tt * width * height * (100 / sid) ** 2
    mc.append(dict(protocol=r['protocol'], median=float(np.median(d)), lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5))))
out['mc_dap'] = mc
# share of variance from each input (standard protocol, log-scale)
d0 = Y * 4 * 14.7 * width * height * (100 / sid) ** 2
ld = np.log(d0)
out['mc_variance_share'] = {nm: float(np.corrcoef(np.log(v), ld)[0, 1] ** 2) for nm, v in [('tube_output', Y), ('beam_width', width), ('distance', sid)]}

json.dump(out, open('phase2.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=lambda x: round(float(x), 4)))
