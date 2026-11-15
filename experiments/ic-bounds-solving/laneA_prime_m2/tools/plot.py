"""data/sd_vs_k.svg: sd (modal over unit instances, large-prime rows), d_lf (measured / recorded-law), k+2 and 2k+2 vs k.
[pcc:unverified]"""
import sys, os, json, collections
here = os.path.dirname(__file__)
sys.path.append(os.path.join(here, '_pylib'))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
D = os.path.join(here, '..', 'data')
rows = json.load(open(D + '/all_instances.json'))
ks = sorted(set(r['k'] for r in rows))
sd_mode, sd_min = {}, {}
for k in ks:
    c = collections.Counter(r['sd_inverse'] for r in rows if r['k'] == k)
    sd_mode[k] = c.most_common(1)[0][0]; sd_min[k] = min(c)
meas = {k: collections.Counter(r['d_lf'] for r in rows if r['k'] == k and r['d_lf'] is not None) for k in ks}
fig, ax = plt.subplots(figsize=(7, 4.5))
kk = list(range(min(ks), max(ks) + 1))
ax.plot(kk, [2 * k + 2 for k in kk], '--', color='#999', label='2k+2 (dimension-type bound)')
ax.plot(kk, [k + 2 for k in kk], ':', color='#555', label='k+2 (U1 bound on d_lf)')
ax.plot(ks, [sd_mode[k] for k in ks], 'o-', color='#c0392b', label='sd = 2+maxwt(g^-1) (modal)')
ax.plot(ks, [sd_min[k] for k in ks], 'x', color='#c0392b', label='sd (smallest seen, exceptions)')
mk = [k for k in ks if meas[k]]
ax.plot(mk, [meas[k].most_common(1)[0][0] for k in mk], 's', color='#2471a3', label='d_lf measured here')
ax.plot(kk, [(k + 14) // 4 for k in kk], '-', color='#2471a3', alpha=.4, label='d_lf law floor((k+14)/4) (recorded, not measured here)')
ax.set_xlabel('depth k'); ax.set_ylabel('degree'); ax.legend(fontsize=7); ax.grid(alpha=.3)
ax.set_title('Two-summand prime chains: static solving degree vs last fall degree [pcc:unverified]', fontsize=9)
fig.tight_layout(); fig.savefig(D + '/sd_vs_k.svg')
print('wrote', D + '/sd_vs_k.svg')
