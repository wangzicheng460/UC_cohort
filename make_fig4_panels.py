from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parent
SC = ROOT / "11.单细胞"
GLY = SC / "GSE214695_PPARG_glycolysis_20260809"
REG = SC / "GSE214695_composition_dorothea_20260810"
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"
OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)

DARK = "#243746"
TEAL = "#168587"
HC = "#4C78A8"
UC = "#D64B4B"
ORANGE = "#D55E00"
GREY = "#4D4D4D"
GRID = "#E7EBEE"

CELL_ORDER = [
    "Epithelial cells", "Stromal cells", "Mast cells", "Enteric glia",
    "Endothelial cells", "Neutrophils", "Myeloid cells", "B cells",
    "T/NK/ILC", "Plasma cells",
]

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 12,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def style_ax(ax):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(DARK)
    ax.spines["bottom"].set_color(DARK)
    ax.tick_params(colors=DARK)


def save_multi(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tif", dpi=1200, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(MIRROR / f"{stem}.png", dpi=1200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {stem}.pdf/.png/.tif")


def panel_a():
    """Add standard UMAP axes to the original broad-cell-type raster."""
    source = GLY / "figures_png" / "Fig1_cell_atlas_by_celltype.png"
    im = Image.open(source).convert("RGB")
    # Remove the original raster title before adding a centered title and axes.
    im = im.crop((0, 88, im.width, im.height))

    left_margin, top_margin, right_margin, bottom_margin = 180, 110, 35, 150
    canvas = Image.new(
        "RGB",
        (im.width + left_margin + right_margin, im.height + top_margin + bottom_margin),
        "white",
    )
    canvas.paste(im, (left_margin, top_margin))
    draw = ImageDraw.Draw(canvas)
    bold = Path("C:/Windows/Fonts/arialbd.ttf")
    regular = Path("C:/Windows/Fonts/arial.ttf")
    title_font = ImageFont.truetype(str(bold), 68)
    tick_font = ImageFont.truetype(str(regular), 34)
    axis_font = ImageFont.truetype(str(bold), 42)

    # The source atlas uses the same UMAP embedding and limits as the companion
    # PPARG UMAPs (ticks at -15, -10, ..., 15 and -10, 0, 10).  The axes are
    # placed along the data panel, leaving the right-hand cell-type legend intact.
    plot_left = left_margin + 60
    plot_right = left_margin + 1480
    plot_top = top_margin + 8
    plot_bottom = top_margin + im.height - 28
    axis_color = "#1F2933"
    draw.text(((plot_left + plot_right) / 2, 24), "Broad cell-type atlas",
              font=title_font, fill=DARK, anchor="ma")
    draw.line((plot_left, plot_top, plot_left, plot_bottom), fill=axis_color, width=4)
    draw.line((plot_left, plot_bottom, plot_right, plot_bottom), fill=axis_color, width=4)

    x_ticks = [-15, -10, -5, 0, 5, 10, 15]
    for value, x in zip(x_ticks, np.linspace(plot_left, plot_right, len(x_ticks))):
        x = int(round(x))
        draw.line((x, plot_bottom, x, plot_bottom + 12), fill=axis_color, width=3)
        draw.text((x, plot_bottom + 20), str(value), font=tick_font,
                  fill=axis_color, anchor="ma")
    y_ticks = [-10, 0, 10]
    y_positions = np.linspace(plot_bottom - 0.83 * (plot_bottom - plot_top),
                              plot_bottom - 0.17 * (plot_bottom - plot_top), len(y_ticks))
    for value, y in zip(y_ticks, y_positions[::-1]):
        y = int(round(y))
        draw.line((plot_left - 12, y, plot_left, y), fill=axis_color, width=3)
        draw.text((plot_left - 22, y), str(value), font=tick_font,
                  fill=axis_color, anchor="rm")
    draw.text(((plot_left + plot_right) / 2, canvas.height - 38), "UMAP 1",
              font=axis_font, fill=axis_color, anchor="ms")
    y_label = Image.new("RGBA", (260, 60), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_label)
    y_draw.text((130, 30), "UMAP 2", font=axis_font, fill=axis_color, anchor="mm")
    y_label = y_label.rotate(90, expand=True)
    canvas.paste(y_label, (20, int((plot_top + plot_bottom - y_label.height) / 2)), y_label)
    # Slightly compress the horizontal footprint for multi-panel assembly while
    # retaining the original white background and all plotted points/labels.
    compressed_width = int(round(canvas.width * 0.78))
    canvas = canvas.resize((compressed_width, canvas.height), Image.Resampling.LANCZOS)
    stem = "Fig4A_broad_celltype_UMAP"
    canvas.save(OUT / f"{stem}.png", dpi=(600, 600))
    canvas.save(OUT / f"{stem}.tif", dpi=(1200, 1200), compression="tiff_lzw")
    canvas.save(OUT / f"{stem}.pdf", "PDF", resolution=300.0)
    canvas.save(MIRROR / f"{stem}.png", dpi=(1200, 1200))
    print(f"Wrote {stem}.pdf/.png/.tif")


def panel_b():
    d = pd.read_csv(REG / "tables" / "broad_celltype_composition_tests.csv")
    d = d.set_index("cell_type_broad").reindex(CELL_ORDER).reset_index()
    y = np.arange(len(d))[::-1]
    sig = d["FDR"] < 0.05
    colors = np.where(sig, ORANGE, GREY)

    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    ax.axvline(0, color="#93999D", linestyle="--", linewidth=1.2, zorder=0)
    for yi, row, color in zip(y, d.to_dict("records"), colors):
        ax.errorbar(row["proportion_difference"] * 100, yi,
                    xerr=[[((row["proportion_difference"] - row["bootstrap_CI_low"]) * 100)],
                          [((row["bootstrap_CI_high"] - row["proportion_difference"]) * 100)]],
                    fmt="o", color=color, capsize=4, linewidth=2.0, markersize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(d["cell_type_broad"])
    ax.set_xlabel("UC − HC mean proportion (bootstrap 95% CI)")
    ax.set_title("Donor-level cell composition differences", loc="center", pad=12)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=GREY, label="FDR ≥ 0.05"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=ORANGE, label="FDR < 0.05"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left", ncol=1,
              bbox_to_anchor=(1.015, 0.98), fontsize=10.5,
              borderaxespad=0, labelspacing=0.9)
    ax.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=100, decimals=0))
    style_ax(ax)
    fig.subplots_adjust(left=0.23, right=0.79, bottom=0.14, top=0.88)
    save_multi(fig, "Fig4B_donor_cell_composition_forest")


def panel_c():
    d = pd.read_csv(GLY / "tables" / "broad_celltype_pparg_glycolysis_summary.csv")
    order = [
        "Enteric glia", "Endothelial cells", "Stromal cells", "Mast cells",
        "Neutrophils", "Myeloid cells", "Plasma cells", "B cells",
        "T/NK/ILC", "Epithelial cells",
    ]
    d["cell_type_broad"] = pd.Categorical(d["cell_type_broad"], categories=order, ordered=True)
    d = d.sort_values(["cell_type_broad", "group"])
    ymap = {name: len(order) - 1 - i for i, name in enumerate(order)}
    xmap = {"HC": 0, "UC": 1}
    x = d["group"].map(xmap).to_numpy()
    y = d["cell_type_broad"].map(ymap).astype(float).to_numpy()
    size = 18 + d["PPARG_pct"].to_numpy() / max(d["PPARG_pct"].max(), 1) * 550
    norm = mpl.colors.Normalize(vmin=0, vmax=d["PPARG_mean"].max())
    cmap = mpl.colors.LinearSegmentedColormap.from_list("pparg", ["#FEE5D9", "#A50F15"])

    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    ax.scatter(x, y, s=size, c=d["PPARG_mean"], cmap=cmap, norm=norm,
               edgecolor="white", linewidth=0.7)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["HC", "UC"])
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order[::-1])
    ax.set_xlim(-0.6, 1.65); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_title("PPARG across broad cell types", loc="center", pad=12)
    style_ax(ax)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.036, pad=0.04)
    cbar.set_label("Average expression")
    pct_values = [5, 10, 15]
    handles = [ax.scatter([], [], s=18 + v / max(d["PPARG_pct"].max(), 1) * 550,
                          facecolor="white", edgecolor=DARK) for v in pct_values]
    ax.legend(handles, [f"{v}%" for v in pct_values], title="PPARG-positive",
              frameon=False, loc="center left", bbox_to_anchor=(1.18, 0.72),
              labelspacing=1.2)
    fig.subplots_adjust(left=0.27, right=0.78, bottom=0.10, top=0.91)
    save_multi(fig, "Fig4C_PPARG_celltype_dotplot")


def panel_d():
    d = pd.read_csv(GLY / "tables" / "UC_broad_celltype_cocolocalization_ranking.csv")
    palette = {
        "Epithelial cells": "#1787C0", "T/NK/ILC": "#E69F00", "B cells": "#56B4E9",
        "Plasma cells": "#CC79A7", "Myeloid cells": "#D55E00", "Neutrophils": "#F0C808",
        "Mast cells": "#009E73", "Stromal cells": "#7A5195",
        "Endothelial cells": "#2F4B7C", "Enteric glia": "#8C564B",
    }
    fig, ax = plt.subplots(figsize=(8.7, 6.5))
    ax.axvline(d["PPARG_pct"].median(), color="#94999C", linestyle="--", linewidth=1.1)
    ax.axhline(d["Hallmark_Glycolysis"].median(), color="#94999C", linestyle="--", linewidth=1.1)
    sizes = 55 + d["n_cells"] / d["n_cells"].max() * 420
    # Manual label offsets avoid overlaps among low-PPARG cell types while
    # preserving the plotted coordinates and quantitative interpretation.
    label_pos = {
        # Labels are staggered around the dense upper-left cluster.  All are
        # connected to their point so their positions remain unambiguous.
        "Plasma cells": (-1.22, 0.06980, "left"),
        "Stromal cells": (2.28, 0.06965, "left"),
        "Myeloid cells": (1.30, 0.06465, "left"),
        "B cells": (-0.82, 0.05305, "left"),
        "Mast cells": (-0.82, 0.04865, "left"),
        "T/NK/ILC": (1.28, 0.04665, "left"),
        "Endothelial cells": (5.85, 0.05135, "left"),
        "Enteric glia": (0.18, 0.03535, "left"),
        "Neutrophils": (0.18, 0.02030, "left"),
        "Epithelial cells": (13.95, 0.06715, "right"),
    }
    for (_, row), s in zip(d.iterrows(), sizes):
        ax.scatter(row["PPARG_pct"], row["Hallmark_Glycolysis"], s=s,
                   color=palette[row["cell_type_broad"]], alpha=0.93,
                   edgecolor="white", linewidth=0.8)
        tx, ty, ha = label_pos[row["cell_type_broad"]]
        ax.annotate(row["cell_type_broad"],
                    (row["PPARG_pct"], row["Hallmark_Glycolysis"]),
                    xytext=(tx, ty), textcoords="data", fontsize=12.5,
                    ha=ha, va="center",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.6),
                    arrowprops=dict(arrowstyle="->", color=DARK, lw=1.0,
                                    shrinkA=3, shrinkB=4, mutation_scale=11,
                                    connectionstyle="arc3,rad=0.04"))
    ax.set_xlabel("Median donor PPARG-positive cells (%)")
    ax.set_ylabel("Median donor Hallmark glycolysis UCell score")
    ax.set_title("PPARG–glycolysis cellular co-localization in UC",
                 loc="center", pad=12, fontsize=16, fontweight="bold")
    style_ax(ax)
    ax.set_xlim(-1.55, d["PPARG_pct"].max() + 0.7)
    ax.set_ylim(d["Hallmark_Glycolysis"].min() - 0.002, d["Hallmark_Glycolysis"].max() + 0.002)
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.90)
    save_multi(fig, "Fig4D_PPARG_glycolysis_cocolocalization")


def panel_e():
    d = pd.read_csv(REG / "tables" / "dorothea_regulon_activity_by_donor.csv")
    d = d[(d["source"] == "PPARG") & (d["lineage"] == "Epithelial cells")].copy()
    d["group"] = pd.Categorical(d["group"], categories=["HC", "UC"], ordered=True)

    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    groups = [d.loc[d["group"] == g, "score"].to_numpy() for g in ["HC", "UC"]]
    bp = ax.boxplot(groups, positions=[0, 1], widths=0.48, patch_artist=True,
                    showfliers=False, medianprops=dict(linewidth=2.0))
    for patch, color in zip(bp["boxes"], [HC, UC]):
        patch.set(facecolor="white", edgecolor=color, linewidth=1.8)
    for whisker, color in zip(bp["whiskers"], [HC, HC, UC, UC]):
        whisker.set(color=color, linewidth=1.5)
    for cap, color in zip(bp["caps"], [HC, HC, UC, UC]):
        cap.set(color=color, linewidth=1.5)
    for med, color in zip(bp["medians"], [HC, UC]):
        med.set(color=color)
    rng = np.random.default_rng(20260813)
    for xi, (vals, color) in enumerate(zip(groups, [HC, UC])):
        jitter = rng.uniform(-0.07, 0.07, len(vals))
        ax.scatter(xi + jitter, vals, color=color, s=48, edgecolor="white",
                   linewidth=0.6, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["HC", "UC"])
    ax.set_ylabel("Predicted PPARG regulon activity\n(normalized weighted mean)")
    ax.set_title("Epithelial PPARG regulon activity",
                 loc="center", pad=12, fontsize=16, fontweight="bold")
    style_ax(ax)
    fig.subplots_adjust(left=0.20, right=0.98, bottom=0.14, top=0.88)
    save_multi(fig, "Fig4E_epithelial_PPARG_regulon_activity")


if __name__ == "__main__":
    panel_a()
    panel_b()
    panel_c()
    panel_d()
    panel_e()
