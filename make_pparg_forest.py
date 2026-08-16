from pathlib import Path
import csv
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
src = ROOT / '6.ML_reanalysis_20260809' / 'rank_normalized_candidate_expression.csv'
out = ROOT / 'submission_package' / 'figures' / 'standalone'
out.mkdir(parents=True, exist_ok=True)

rows = []
with src.open(newline='', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        if r['set'] == 'discovery':
            r['label'] = int(r['label'])
            r['PPARG'] = float(r['PPARG'])
            rows.append(r)

cohorts = ['GSE73661', 'GSE75214', 'GSE87466', 'GSE107499']
results = []
for cohort in cohorts:
    vals = [r for r in rows if r['cohort'] == cohort]
    uc = [r['PPARG'] for r in vals if r['label'] == 1]
    hc = [r['PPARG'] for r in vals if r['label'] == 0]
    diff = sum(uc) / len(uc) - sum(hc) / len(hc)
    var_uc = sum((x - sum(uc)/len(uc))**2 for x in uc) / (len(uc)-1)
    var_hc = sum((x - sum(hc)/len(hc))**2 for x in hc) / (len(hc)-1)
    se = math.sqrt(var_uc / len(uc) + var_hc / len(hc))
    results.append((cohort, diff, diff - 1.96 * se, diff + 1.96 * se, len(hc), len(uc)))

plt.rcParams.update({'font.family': 'Arial', 'font.size': 10, 'axes.titlesize': 14,
                     'axes.labelsize': 11, 'xtick.labelsize': 10, 'ytick.labelsize': 10})
fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=300)
y = list(range(len(results)-1, -1, -1))
color = '#147D7E'
for yi, (cohort, diff, low, high, nhc, nuc) in zip(y, results):
    ax.plot([low, high], [yi, yi], color=color, linewidth=2.0, solid_capstyle='round')
    ax.plot([low, low], [yi-0.09, yi+0.09], color=color, linewidth=1.2)
    ax.plot([high, high], [yi-0.09, yi+0.09], color=color, linewidth=1.2)
    ax.scatter(diff, yi, s=55, color=color, edgecolor='white', linewidth=0.8, zorder=3)

ax.axvline(0, color='#8A98A6', linestyle=(0, (4, 3)), linewidth=1)
ax.set_yticks(y, [r[0] for r in results])
ax.set_xlabel('PPARG rank-normalized mean difference (UC − HC)')
ax.set_title('Cross-cohort PPARG suppression')
ax.text(0.0, 1.02, 'Points: effect estimates; bars: 95% CI', transform=ax.transAxes,
        ha='left', va='bottom', color='#566A78', fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.grid(axis='x', color='#E8EEF3', linewidth=0.7)
ax.set_axisbelow(True)
fig.text(0.99, 0.015, 'Discovery cohorts; HC and UC sample sizes shown in source table',
         ha='right', va='bottom', fontsize=8, color='#566A78')
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ('png', 'pdf', 'tif'):
    fig.savefig(out / f'PPARG_discovery_effect_forest.{ext}', dpi=600 if ext == 'tif' else 300,
                bbox_inches='tight')
plt.close(fig)

with (out / 'PPARG_discovery_effect_forest_values.csv').open('w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['cohort', 'effect_UC_minus_HC', 'CI_low', 'CI_high', 'n_HC', 'n_UC'])
    w.writerows(results)
