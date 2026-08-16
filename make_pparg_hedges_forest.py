from pathlib import Path
import csv
import math

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
manifest_path = ROOT / '6.ML_reanalysis_20260809' / 'sample_manifest.csv'
data_dir = next(ROOT.glob('2*/'))
out_dir = ROOT / 'submission_package' / 'figures' / 'standalone'
out_dir.mkdir(parents=True, exist_ok=True)

cohorts = ['GSE73661', 'GSE75214', 'GSE87466', 'GSE107499']
with manifest_path.open(encoding='utf-8-sig', newline='') as f:
    manifest = list(csv.DictReader(f))

results = []
for cohort in cohorts:
    source = data_dir / f'{cohort}.csv'
    with source.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        pparg_row = next(row for row in reader if row[0].upper() == 'PPARG')
    values = {sample: float(pparg_row[i]) for i, sample in enumerate(header[1:], start=1)}
    cohort_meta = [r for r in manifest if r['cohort'] == cohort and r['set'] == 'discovery']
    hc = [values[r['sample']] for r in cohort_meta if r['label'] == '0']
    uc = [values[r['sample']] for r in cohort_meta if r['label'] == '1']
    n_hc, n_uc = len(hc), len(uc)
    mean_hc, mean_uc = sum(hc) / n_hc, sum(uc) / n_uc
    var_hc = sum((x - mean_hc) ** 2 for x in hc) / (n_hc - 1)
    var_uc = sum((x - mean_uc) ** 2 for x in uc) / (n_uc - 1)
    df = n_hc + n_uc - 2
    pooled_sd = math.sqrt(((n_hc - 1) * var_hc + (n_uc - 1) * var_uc) / df)
    cohen_d = (mean_uc - mean_hc) / pooled_sd
    correction = 1 - 3 / (4 * df - 1)
    hedges_g = correction * cohen_d
    var_g = correction**2 * ((n_hc + n_uc) / (n_hc * n_uc) + cohen_d**2 / (2 * df))
    se_g = math.sqrt(var_g)
    results.append({
        'cohort': cohort, 'g': hedges_g, 'low': hedges_g - 1.96 * se_g,
        'high': hedges_g + 1.96 * se_g, 'n_HC': n_hc, 'n_UC': n_uc,
        'mean_HC': mean_hc, 'mean_UC': mean_uc,
    })

# Publication-friendly raster figure. PDF/TIFF are exported alongside PNG.
# A narrower, taller canvas matches panel F in the final two-row Figure 2.
W, H = 2050, 2200
img = Image.new('RGB', (W, H), 'white')
draw = ImageDraw.Draw(img)
font_path = Path('C:/Windows/Fonts/arial.ttf')
bold_path = Path('C:/Windows/Fonts/arialbd.ttf')
font = lambda size, bold=False: ImageFont.truetype(str(bold_path if bold else font_path), size)
title_font, subtitle_font = font(58, True), font(34)
label_font, tick_font = font(35), font(31)
dark, teal, grey, grid = '#233746', '#147D7E', '#8A98A6', '#E8EEF3'

draw.text((W / 2, 105), 'Cross-cohort PPARG suppression', fill=dark,
          font=title_font, anchor='ma')
draw.text((W / 2, 182), "Hedges' g (UC - HC) with 95% confidence intervals",
          fill='#566A78', font=subtitle_font, anchor='ma')

left, right, top, bottom = 650, 1500, 430, 1860
xmin, xmax = -3.5, 0.5
def xpix(x):
    return left + (x - xmin) / (xmax - xmin) * (right - left)
for tick in [-3, -2, -1, 0]:
    xp = xpix(tick)
    draw.line((xp, top - 20, xp, bottom), fill=grid, width=2)
    draw.text((xp, bottom + 28), f'{tick}', fill=dark, font=tick_font, anchor='ma')
draw.line((xpix(0), top - 15, xpix(0), bottom), fill=grey, width=3)
draw.line((left, bottom, right, bottom), fill=dark, width=3)
draw.text(((left + right) / 2, bottom + 98), "Hedges' g (UC - HC)", fill=dark,
          font=label_font, anchor='ma')
draw.text((left - 42, top - 74), 'Discovery cohort', fill=dark, font=label_font, anchor='ra')

ys = [top + i * ((bottom - top - 70) / (len(results) - 1)) for i in range(len(results))]
for result, y in zip(results, ys):
    label = f"{result['cohort']}  (HC={result['n_HC']}, UC={result['n_UC']})"
    draw.text((left - 35, y), label, fill=dark, font=label_font, anchor='ra')
    lo, hi, g = xpix(result['low']), xpix(result['high']), xpix(result['g'])
    draw.line((lo, y, hi, y), fill=teal, width=8)
    draw.line((lo, y - 17, lo, y + 17), fill=teal, width=5)
    draw.line((hi, y - 17, hi, y + 17), fill=teal, width=5)
    draw.ellipse((g - 16, y - 16, g + 16, y + 16), fill=teal, outline='white', width=4)
    draw.text((right + 25, y), f"{result['g']:.2f} [{result['low']:.2f}, {result['high']:.2f}]",
              fill=dark, font=tick_font, anchor='lm')

png = out_dir / 'Fig2F_PPARG_discovery_Hedges_g_forest.png'
tif = out_dir / 'Fig2F_PPARG_discovery_Hedges_g_forest.tif'
pdf = out_dir / 'Fig2F_PPARG_discovery_Hedges_g_forest.pdf'
img.save(png, dpi=(600, 600))
img.save(tif, dpi=(600, 600), compression='tiff_lzw')
img.save(pdf, 'PDF', resolution=300.0)

mirror_dir = ROOT / '图片1200dpi'
mirror_dir.mkdir(parents=True, exist_ok=True)
img.save(mirror_dir / png.name, dpi=(1200, 1200))

with (out_dir / 'Fig2F_PPARG_discovery_Hedges_g_values.csv').open('w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['cohort', 'n_HC', 'n_UC', 'mean_HC', 'mean_UC', 'hedges_g', 'CI_low', 'CI_high'])
    writer.writeheader()
    for r in results:
        writer.writerow({'cohort': r['cohort'], 'n_HC': r['n_HC'], 'n_UC': r['n_UC'],
                         'mean_HC': f"{r['mean_HC']:.6f}", 'mean_UC': f"{r['mean_UC']:.6f}",
                         'hedges_g': f"{r['g']:.6f}", 'CI_low': f"{r['low']:.6f}", 'CI_high': f"{r['high']:.6f}"})
print('Wrote:', png, tif, pdf)
for r in results:
    print(r['cohort'], f"g={r['g']:.4f}", f"CI=[{r['low']:.4f}, {r['high']:.4f}]")
