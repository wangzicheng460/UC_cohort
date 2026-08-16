from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"
OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)

DARK = "#233746"
TEAL = "#168487"
ORANGE = "#E56B4A"
BLUE = "#397FB5"
PURPLE = "#8D63B7"
GREEN = "#38B86D"
GREY = "#AAB7C1"
LIGHT = "#E7EDF2"
RED = "#D95A4E"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.titlesize": 15,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.edgecolor": DARK,
    "axes.labelcolor": DARK,
    "xtick.color": DARK,
    "ytick.color": DARK,
    "text.color": DARK,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save_all(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tif", dpi=600, bbox_inches="tight",
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(MIRROR / f"{stem}.png", dpi=1200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1)
    return out


def hedges_g(x_uc: np.ndarray, x_hc: np.ndarray) -> tuple[float, float]:
    x_uc, x_hc = np.asarray(x_uc, float), np.asarray(x_hc, float)
    n1, n0 = len(x_uc), len(x_hc)
    df = n1 + n0 - 2
    v1, v0 = np.var(x_uc, ddof=1), np.var(x_hc, ddof=1)
    sp2 = ((n1 - 1) * v1 + (n0 - 1) * v0) / df
    if not np.isfinite(sp2) or sp2 <= 0:
        return np.nan, np.nan
    d = (np.mean(x_uc) - np.mean(x_hc)) / math.sqrt(sp2)
    J = 1 - 3 / (4 * df - 1)
    g = J * d
    var_g = J**2 * ((n1 + n0) / (n1 * n0) + d**2 / (2 * df))
    return g, var_g


def random_effects(effects: np.ndarray, variances: np.ndarray) -> tuple[float, float, float, float]:
    ok = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    y, v = effects[ok], variances[ok]
    w = 1 / v
    fixed = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - fixed) ** 2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (Q - (len(y) - 1)) / c) if c > 0 else 0.0
    wr = 1 / (v + tau2)
    pooled = np.sum(wr * y) / np.sum(wr)
    se = math.sqrt(1 / np.sum(wr))
    return pooled, pooled - 1.96 * se, pooled + 1.96 * se, tau2


def prepare_metadata(ids: pd.Series) -> pd.DataFrame:
    manifest = pd.read_csv(ROOT / "6.ML_reanalysis_20260809" / "sample_manifest.csv")
    manifest = manifest.set_index("sample")
    base = ids.str.replace(r"_(con|treat)$", "", regex=True)
    meta = manifest.loc[base].reset_index(drop=False)
    meta["id"] = ids.to_numpy()
    return meta


def make_panel_a() -> pd.DataFrame:
    cib = pd.read_csv(ROOT / "7.CIBERSORT" / "CIBERSORT-Results.txt", sep="\t")
    meta = prepare_metadata(cib["id"])
    dat = cib.merge(meta[["id", "cohort", "label"]], on="id", how="inner")
    cell_cols = [c for c in cib.columns if c != "id"]
    rows = []
    for cell in cell_cols:
        cohort_effects = []
        for cohort, sub in dat.groupby("cohort"):
            # CIBERSORT fractions are bounded and frequently zero-inflated.
            # The arcsine-square-root transform stabilizes their variance before SMD estimation.
            uc = np.arcsin(np.sqrt(np.clip(sub.loc[sub.label == 1, cell].to_numpy(), 0, 1)))
            hc = np.arcsin(np.sqrt(np.clip(sub.loc[sub.label == 0, cell].to_numpy(), 0, 1)))
            if len(uc) < 2 or len(hc) < 2:
                continue
            g, var = hedges_g(uc, hc)
            cohort_effects.append((cohort, g, var, len(hc), len(uc)))
        if len(cohort_effects) < 2:
            continue
        pooled, lo, hi, tau2 = random_effects(
            np.array([x[1] for x in cohort_effects]),
            np.array([x[2] for x in cohort_effects]),
        )
        se = (hi - lo) / (2 * 1.96)
        p = 2 * stats.norm.sf(abs(pooled / se))
        rows.append({"cell_type": cell, "g": pooled, "low": lo, "high": hi,
                     "p": p, "tau2": tau2, "k": len(cohort_effects)})
    res = pd.DataFrame(rows)
    res["FDR"] = bh_fdr(res.p.to_numpy())
    res = res.sort_values("g")
    res.to_csv(OUT / "Fig6A_CIBERSORT_random_effects_values.csv", index=False)

    significant = res.loc[res.FDR < 0.05].copy()
    if len(significant) > 12:
        significant = significant.iloc[np.argsort(np.abs(significant.g))[-12:]].sort_values("g")
    elif len(significant) < 8:
        significant = res.iloc[np.argsort(np.abs(res.g))[-10:]].sort_values("g")

    fig, ax = plt.subplots(figsize=(7.7, 6.4))
    y = np.arange(len(significant))
    colors = np.where(significant.g >= 0, ORANGE, BLUE)
    ax.hlines(y, significant.low, significant.high, color=colors, lw=2.2)
    ax.scatter(significant.g, y, c=colors, s=58, zorder=3, edgecolor="white", linewidth=0.8)
    ax.axvline(0, color="#94A4B0", ls="--", lw=1.2)
    ax.set_yticks(y, significant.cell_type)
    ax.set_xlabel("Pooled Hedges' g (UC - HC), 95% CI")
    ax.set_title("Cross-cohort immune-cell alterations", loc="center", pad=14)
    ax.grid(axis="x", color=LIGHT, lw=1)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    span = max(abs(significant.low.min()), abs(significant.high.max()))
    ax.set_xlim(-span * 1.08, span * 1.08)
    for yi, (_, r) in enumerate(significant.iterrows()):
        label = f"FDR {r.FDR:.3f}" if r.FDR >= 0.001 else "FDR < 0.001"
        ax.text(1.025, yi, label, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=8.5, clip_on=False)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE,
                              label="Lower in UC", markersize=8),
                       Line2D([0], [0], marker="o", color="none", markerfacecolor=ORANGE,
                              label="Higher in UC", markersize=8)],
              loc="lower center", bbox_to_anchor=(0.5, -0.19), ncol=2, frameon=False)
    fig.subplots_adjust(left=0.36, right=0.83, top=0.90, bottom=0.19)
    save_all(fig, "Fig6A_cross_cohort_immune_forest")
    return res


def partial_spearman(x: np.ndarray, y: np.ndarray, label: np.ndarray) -> tuple[float, float, float, float]:
    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    X = np.column_stack([np.ones(len(label)), label])
    rx = xr - X @ np.linalg.lstsq(X, xr, rcond=None)[0]
    ry = yr - X @ np.linalg.lstsq(X, yr, rcond=None)[0]
    rho = np.corrcoef(rx, ry)[0, 1]
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1 / math.sqrt(max(len(x) - 4, 1))
    return rho, math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se), se**2


def make_panel_b() -> pd.DataFrame:
    xcell = pd.read_csv(ROOT / "8.IOBR" / "xCell.txt", sep="\t", index_col=0)
    expr = pd.read_csv(ROOT / "7.CIBERSORT" / "normalize.txt", sep="\t", index_col=0)
    samples = pd.Series(xcell.columns, name="id")
    meta = prepare_metadata(samples)
    meta = meta.set_index("id")
    pparg = expr.loc["PPARG", samples].astype(float)
    outcomes = ["Epithelial cells", "Plasma cells"]
    rows = []
    for outcome in outcomes:
        score = xcell.loc[outcome, samples].astype(float)
        zs, vs = [], []
        for cohort in ["GSE73661", "GSE75214", "GSE87466", "GSE107499"]:
            mask = meta.loc[samples, "cohort"].to_numpy() == cohort
            lab = meta.loc[samples[mask], "label"].to_numpy(float)
            rho, lo, hi, var_z = partial_spearman(
                pparg.loc[samples[mask]].to_numpy(), score.loc[samples[mask]].to_numpy(), lab)
            z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
            zs.append(z); vs.append(var_z)
            rows.append({"cell_type": outcome, "cohort": cohort, "n": int(mask.sum()),
                         "rho": rho, "low": lo, "high": hi, "pooled": False})
        pooled_z, lo_z, hi_z, _ = random_effects(np.array(zs), np.array(vs))
        pooled = math.tanh(pooled_z)
        lo, hi = math.tanh(lo_z), math.tanh(hi_z)
        rows.append({"cell_type": outcome, "cohort": "Pooled", "n": int(len(samples)),
                     "rho": pooled, "low": lo, "high": hi, "pooled": True})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "Fig6B_xCell_partial_Spearman_values.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    cohort_order = ["GSE73661", "GSE75214", "GSE87466", "GSE107499"]
    cohort_colors = {
        "GSE73661": "#4E79A7",
        "GSE75214": "#76B7B2",
        "GSE87466": "#B07AA1",
        "GSE107499": "#F28E2B",
    }
    x_positions = [0, 1]
    for cohort in cohort_order:
        epithelial = res[(res.cell_type == "Epithelial cells") & (res.cohort == cohort)].iloc[0]
        plasma = res[(res.cell_type == "Plasma cells") & (res.cohort == cohort)].iloc[0]
        vals = [epithelial.rho, plasma.rho]
        ax.plot(x_positions, vals, color=cohort_colors[cohort], lw=2.0, alpha=0.82,
                marker="o", markersize=7, markeredgecolor="white", markeredgewidth=0.8,
                label=cohort, zorder=3)

    pooled_ep = res[(res.cell_type == "Epithelial cells") & res.pooled].iloc[0]
    pooled_pl = res[(res.cell_type == "Plasma cells") & res.pooled].iloc[0]
    pooled_vals = [pooled_ep.rho, pooled_pl.rho]
    ax.plot(x_positions, pooled_vals, color=DARK, lw=4.2, marker="D", markersize=10,
            markerfacecolor=DARK, markeredgecolor="white", markeredgewidth=1.2,
            label="Pooled", zorder=5)
    ax.errorbar(0, pooled_ep.rho,
                yerr=[[pooled_ep.rho - pooled_ep.low], [pooled_ep.high - pooled_ep.rho]],
                fmt="none", ecolor=TEAL, elinewidth=2.4, capsize=6, capthick=2.0, zorder=4)
    ax.errorbar(1, pooled_pl.rho,
                yerr=[[pooled_pl.rho - pooled_pl.low], [pooled_pl.high - pooled_pl.rho]],
                fmt="none", ecolor=PURPLE, elinewidth=2.4, capsize=6, capthick=2.0, zorder=4)

    ax.axhline(0, color="#94A4B0", ls="--", lw=1.2)
    ax.set_xlim(-0.38, 1.38)
    ax.set_ylim(-0.82, 0.82)
    ax.set_xticks(x_positions, ["Epithelial cells", "Plasma cells"])
    ax.set_ylabel("Disease-adjusted partial Spearman ρ")
    ax.set_title("Opposing PPARG–xCell associations", loc="center", pad=14)
    ax.grid(axis="y", color=LIGHT, lw=1)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0, labelsize=12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=5,
              frameon=False, fontsize=8.2, handlelength=1.6, columnspacing=0.9)
    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.24, top=0.87)
    save_all(fig, "Fig6B_PPARG_xCell_adjusted_associations")
    return res


def rounded_box(ax, xy, width, height, text, fc, ec=DARK, lw=1.5,
                fontsize=11, weight="normal"):
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.018,rounding_size=0.018",
                           facecolor=fc, edgecolor=ec, linewidth=lw)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=weight, color=DARK)
    return patch


def make_panel_c() -> None:
    rel = pd.read_csv(ROOT / "submission_package" / "supplementary" / "prepared_csv" / "S10_GutMGene.csv")
    fig, ax = plt.subplots(figsize=(9.0, 5.7))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    microbes = rel.microbe.drop_duplicates().tolist()
    micro_y = {microbes[0]: 0.66, microbes[1]: 0.27}
    micro_color = {microbes[0]: "#CDEEDB", microbes[1]: "#E4D8F0"}
    mets = rel.metabolite.tolist()
    ys = np.linspace(0.87, 0.13, len(mets))
    met_y = dict(zip(mets, ys))
    for m in microbes:
        rounded_box(ax, (0.04, micro_y[m] - 0.055), 0.27, 0.11, m,
                    micro_color[m], ec=GREEN if m == "Enterococcus faecalis" else PURPLE,
                    lw=2, fontsize=11, weight="bold")
    for met in mets:
        selected = met == "Lariciresinol"
        rounded_box(ax, (0.71, met_y[met] - 0.043), 0.24, 0.086, met,
                    "#F5C7BF" if selected else "#EAF0F4",
                    ec=ORANGE if selected else "#738897", lw=3 if selected else 1.3,
                    fontsize=11.5 if selected else 10.5, weight="bold" if selected else "normal")
    for _, r in rel.iterrows():
        selected = r.metabolite == "Lariciresinol"
        color = ORANGE if selected else (GREEN if r.microbe == "Enterococcus faecalis" else PURPLE)
        arrow = FancyArrowPatch((0.31, micro_y[r.microbe]), (0.71, met_y[r.metabolite]),
                                arrowstyle="-", connectionstyle="arc3,rad=0.0",
                                linewidth=4.2 if selected else 1.6,
                                color=color, alpha=1 if selected else 0.68)
        ax.add_patch(arrow)
    ax.text(0.175, 0.965, "Gut microbe", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.83, 0.965, "Database-linked metabolite", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.5, 1.055, "GutMGene-derived microbe–metabolite network", ha="center",
            va="center", fontsize=15, fontweight="bold", transform=ax.transAxes)
    ax.annotate("Prioritized for structural testing", xy=(0.71, met_y["Lariciresinol"]),
                xytext=(0.50, met_y["Lariciresinol"] - 0.11), fontsize=10.5, fontweight="bold",
                color=ORANGE, arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.5))
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.05)
    save_all(fig, "Fig6C_GutMGene_node_link_network")


def make_panel_d() -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.7))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 1.035, "Independent evidence streams converge on structural testing",
            ha="center", va="center", fontsize=15, fontweight="bold", transform=ax.transAxes)

    rounded_box(ax, (0.055, 0.70), 0.38, 0.18,
                "Reproducible PPARG suppression\nBulk transcriptomics",
                "#E1F0F1", ec=TEAL, lw=2, fontsize=11.5, weight="bold")
    rounded_box(ax, (0.055, 0.43), 0.38, 0.18,
                "Epithelial localization\nSingle-cell and spatial data",
                "#E1F0F1", ec=TEAL, lw=2, fontsize=11.5, weight="bold")
    rounded_box(ax, (0.565, 0.70), 0.38, 0.18,
                "E. faecalis–lariciresinol\nGutMGene-supported hypothesis",
                "#FBE9E4", ec=ORANGE, lw=2, fontsize=11.5, weight="bold")
    rounded_box(ax, (0.565, 0.43), 0.38, 0.18,
                "Lariciresinol\nPrioritized microbiota-derived compound",
                "#FBE9E4", ec=ORANGE, lw=2, fontsize=11.5, weight="bold")
    rounded_box(ax, (0.25, 0.10), 0.50, 0.17,
                "PPARG–lariciresinol\nDocking and molecular dynamics",
                "#FFF3CF", ec=DARK, lw=2, fontsize=13, weight="bold")

    for start, end, color in [
        ((0.245, 0.70), (0.245, 0.61), TEAL),
        ((0.755, 0.70), (0.755, 0.61), ORANGE),
        ((0.245, 0.43), (0.43, 0.27), TEAL),
        ((0.755, 0.43), (0.57, 0.27), ORANGE),
    ]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                                     lw=2, color=color,
                                     linestyle="--" if start[1] == 0.43 else "-"))
    ax.text(0.205, 0.645, "PPARG target", ha="right", va="center",
            fontsize=9.5, color=TEAL)
    ax.text(0.795, 0.645, "Candidate ligand", ha="left", va="center",
            fontsize=9.5, color=ORANGE)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.04)
    save_all(fig, "Fig6D_evidence_convergence_schematic")


if __name__ == "__main__":
    a = make_panel_a()
    b = make_panel_b()
    make_panel_c()
    make_panel_d()
    print("Fig6A selected signals:")
    print(a.sort_values("FDR").head(12).to_string(index=False))
    print("\nFig6B associations:")
    print(b.to_string(index=False))
    print(f"\nWrote standalone panels to: {OUT}")
