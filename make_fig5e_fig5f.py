from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
SPATIAL = ROOT / "空转" / "GSE189184_spatial_validation_20260810"
TABLES = SPATIAL / "tables"
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"
OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)

DARK = "#243746"
HC = "#3C78A8"
UC = "#D64B45"
GRID = "#E4E9ED"

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def style_ax(ax, grid_axis="y"):
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(DARK)
    ax.spines["bottom"].set_color(DARK)
    ax.tick_params(colors=DARK)


def save_multi(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(
        OUT / f"{stem}.tif", dpi=1200, bbox_inches="tight", facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(MIRROR / f"{stem}.png", dpi=1200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {stem}.pdf/.png/.tif")


def fig5e():
    d = pd.read_csv(TABLES / "sample_PPARG_activity.csv")
    groups = ["HC", "UC"]
    colors = [HC, UC]
    values = [
        d.loc[d["group"] == group, "epithelial_weighted_mean"].to_numpy()
        for group in groups
    ]

    fig, ax = plt.subplots(figsize=(6.2, 5.25))
    bp = ax.boxplot(
        values, positions=[0, 1], widths=0.46, patch_artist=True,
        showfliers=False, medianprops={"linewidth": 2.3},
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set(facecolor="white", edgecolor=color, linewidth=2.0)
    for whisker, color in zip(bp["whiskers"], [HC, HC, UC, UC]):
        whisker.set(color=color, linewidth=1.7)
    for cap, color in zip(bp["caps"], [HC, HC, UC, UC]):
        cap.set(color=color, linewidth=1.7)
    for median, color in zip(bp["medians"], colors):
        median.set(color=color)

    offsets = {
        "B10": -0.055, "C5": 0.055,
        "B12": -0.080, "B13": -0.040, "B4": 0.000,
        "B5": 0.040, "C2": 0.080,
    }
    for xi, (group, color) in enumerate(zip(groups, colors)):
        sub = d.loc[d["group"] == group].copy()
        for _, row in sub.iterrows():
            x = xi + offsets[row["sample"]]
            y = row["epithelial_weighted_mean"]
            ax.scatter(x, y, s=66, color=color, edgecolor="white", linewidth=0.8, zorder=3)
            ax.annotate(
                row["sample"], (x, y), xytext=(5, 5), textcoords="offset points",
                fontsize=10.2, color=DARK, ha="left", va="bottom",
            )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(groups)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Epithelial-weighted PPARG regulon activity")
    ax.set_title("Epithelial-weighted PPARG regulon activity", loc="center", pad=12)
    style_ax(ax, "y")
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.13, top=0.88)
    save_multi(fig, "Fig5E_epithelial_weighted_PPARG_activity")


def fig5f():
    d = pd.read_csv(TABLES / "spatial_PPARG_glycolysis_correlations_by_section.csv")
    d = d[["sample", "group", "rho_epithelial_top_quartile"]].copy()
    d = d.sort_values("rho_epithelial_top_quartile", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    colors = d["group"].map({"HC": HC, "UC": UC}).to_numpy()

    fig, ax = plt.subplots(figsize=(6.5, 5.25))
    ax.axvline(0, color="#92999E", linestyle="--", linewidth=1.2, zorder=0)
    for yi, (_, row), color in zip(y, d.iterrows(), colors):
        rho = row["rho_epithelial_top_quartile"]
        ax.hlines(yi, 0, rho, color=color, linewidth=2.5, alpha=0.78)
        ax.scatter(rho, yi, s=82, color=color, edgecolor="white", linewidth=0.9, zorder=3)
        ax.text(rho + 0.018, yi, f"{rho:.2f}", va="center", ha="left",
                fontsize=11, color=DARK)

    ax.set_yticks(y)
    ax.set_yticklabels(d["sample"])
    ax.set_xlim(-0.03, 0.62)
    ax.set_xlabel("Spearman ρ in epithelial-enriched spots")
    ax.set_ylabel("Visium section")
    ax.set_title("Section-level PPARG–glycolysis spatial association",
                 loc="center", pad=12)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=HC, label="HC", markersize=7),
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=UC, label="UC", markersize=7),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", ncol=2,
              columnspacing=1.0, handletextpad=0.4)
    style_ax(ax, "x")
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.14, top=0.88)
    save_multi(fig, "Fig5F_PPARG_glycolysis_lollipop")


if __name__ == "__main__":
    fig5e()
    fig5f()
