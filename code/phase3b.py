import re, json
src = open('phase3.py').read()
# keep everything up to the regional analysis output, then add log-log fits
cut = src.index('# ---------- 2. regional analysis')
exec(src[:cut])
def loglog(y):
    X = np.column_stack([np.ones_like(I), np.log(I), age, sex]); b, se, df = ols(np.log(y), X)
    rng = np.random.default_rng(20261121); bs = []
    for _ in range(2000):
        i = rng.integers(0, len(I), len(I)); bb, *_ = np.linalg.lstsq(X[i], np.log(y)[i], rcond=None); bs.append(bb[1])
    return dict(slope=b[1], lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), p=2 * stats.t.sf(abs(b[1] / se[1]), df))
r = {'nps_low': loglog(total), 'nps_high': loglog(high)}
for k in regions:
    r['noise_' + k] = loglog(np.array(reg[k]['noise'])); r['cnr_' + k] = loglog(np.array(reg[k]['cnr']))
json.dump(r, open('phase3b.json', 'w'), indent=1, default=float)
for k, v in r.items(): print(k, {a: round(float(b), 4) if a != 'p' else float('%.2g' % b) for a, b in v.items()})
