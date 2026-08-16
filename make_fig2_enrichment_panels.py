from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "4.富集分析"
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"
OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)

DARK = "#243746"
GRID = "#E6EBEF"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.labelsize": 13,
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 12,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def ratio_value(value):
    a, b = str(value).split("/")
    return float(a) / float(b)


def save_multi(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tif", dpi=1200, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(MIRROR / f"{stem}.png", dpi=1200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {stem}.pdf/.png/.tif")


def bubble_legend(ax, counts):
    values = np.linspace(min(counts), max(counts), 4).round().astype(int)
    handles = [ax.scatter([], [], s=24 + (v - min(counts)) / max(max(counts)-min(counts), 1) * 130,
                          facecolor="white", edgecolor="#667580", linewidth=1.0)
               for v in values]
    return handles, [str(v) for v in values]


def make_go():
    go = pd.read_csv(DATA / "GO.txt", sep="\t")
    selected = {
        "BP": [
            "response to molecule of bacterial origin",
            "response to lipopolysaccharide",
            "leukocyte migration",
            "regulation of inflammatory response",
        ],
        "CC": [
            "apical plasma membrane",
            "external side of plasma membrane",
            "collagen-containing extracellular matrix",
            "microvillus",
        ],
        "MF": [
            "cytokine activity",
            "cytokine receptor binding",
            "immune receptor activity",
            "chemokine activity",
        ],
    }
    pieces = []
    for ontology, terms in selected.items():
        sub = go[(go["ONTOLOGY"] == ontology) & go["Description"].isin(terms)].copy()
        sub["order"] = sub["Description"].map({v: i for i, v in enumerate(terms)})
        pieces.append(sub.sort_values("order"))
    d = pd.concat(pieces, ignore_index=True)
    d["GeneRatioValue"] = d["GeneRatio"].map(ratio_value)

    heights = [len(selected[k]) for k in ("BP", "CC", "MF")]
    fig, axes = plt.subplots(3, 1, figsize=(9.6, 6.8),
                             gridspec_kw={"height_ratios": heights, "hspace": 0.06})
    norm = mpl.colors.Normalize(vmin=d["qvalue"].min(), vmax=d["qvalue"].max())
    cmap = mpl.colors.LinearSegmentedColormap.from_list("q", ["#E66561", "#4A86B8"])
    cmin, cmax = d["Count"].min(), d["Count"].max()
    for ax, ontology in zip(axes, ("BP", "CC", "MF")):
        sub = d[d["ONTOLOGY"] == ontology].iloc[::-1]
        y = np.arange(len(sub))
        sizes = 42 + (sub["Count"] - cmin) / max(cmax - cmin, 1) * 150
        ax.scatter(sub["GeneRatioValue"], y, s=sizes, c=sub["qvalue"], cmap=cmap,
                   norm=norm, edgecolor="#4B5860", linewidth=0.8, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([textwrap.fill(x, 42) for x in sub["Description"]], fontsize=12.5)
        # Add vertical padding so the first and last bubbles are not clipped
        # by the facet borders after high-resolution export.
        ax.set_ylim(-0.35, len(sub) - 0.65)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.set_xlim(0.015, 0.11)
        ax.tick_params(colors=DARK)
        for spine in ax.spines.values():
            spine.set_color("#8B979E")
        ax.text(1.012, 0.5, ontology, transform=ax.transAxes, rotation=-90,
                ha="left", va="center", fontsize=11.5, color=DARK,
                bbox=dict(boxstyle="square,pad=0.35", fc="#D7D7D7", ec="#8B8B8B", lw=0.8))
    axes[-1].set_xlabel("GeneRatio", color=DARK)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = fig.add_axes([0.82, 0.61, 0.022, 0.27])
    cbar = fig.colorbar(sm, cax=cax)
    go_ticks = np.linspace(d["qvalue"].min(), d["qvalue"].max(), 5)
    cbar.set_ticks(go_ticks)
    cbar.set_ticklabels([f"{x / 1e-7:.1f}" for x in go_ticks])
    cbar.ax.yaxis.get_offset_text().set_visible(False)
    cbar.set_label(r"q-value ($\times 10^{-7}$)", fontsize=12)
    cbar.ax.tick_params(labelsize=10.5)
    handles, labels = bubble_legend(axes[0], d["Count"].to_numpy())
    fig.legend(handles, labels, title="Count", frameon=False, fontsize=10.5,
               title_fontsize=11.5, loc="upper left", bbox_to_anchor=(0.80, 0.47),
               labelspacing=0.9, handletextpad=0.8)
    fig.subplots_adjust(left=0.39, right=0.76, bottom=0.11, top=0.98)
    save_multi(fig, "Fig2B_refined_UC_relevant_GO")


def make_kegg():
    kegg = pd.read_csv(DATA / "KEGG.txt", sep="\t")
    # Preserve the ten pathways shown in the existing main panel, while enlarging type.
    terms = [
        "Cytokine-cytokine receptor interaction",
        "Viral protein interaction with cytokine and cytokine receptor",
        "Rheumatoid arthritis",
        "Staphylococcus aureus infection",
        "Hematopoietic cell lineage",
        "Complement and coagulation cascades",
        "IL-17 signaling pathway",
        "Amoebiasis",
        "Pertussis",
        "Leishmaniasis",
    ]
    d = kegg[kegg["Description"].isin(terms)].copy()
    d["order"] = d["Description"].map({v: i for i, v in enumerate(terms)})
    d = d.sort_values("order").iloc[::-1]
    d["GeneRatioValue"] = d["GeneRatio"].map(ratio_value)
    norm = mpl.colors.Normalize(vmin=d["qvalue"].min(), vmax=d["qvalue"].max())
    cmap = mpl.colors.LinearSegmentedColormap.from_list("q", ["#E66561", "#4A86B8"])
    cmin, cmax = d["Count"].min(), d["Count"].max()

    fig, ax = plt.subplots(figsize=(7.0, 6.8))
    y = np.arange(len(d))
    sizes = 42 + (d["Count"] - cmin) / max(cmax - cmin, 1) * 150
    ax.scatter(d["GeneRatioValue"], y, s=sizes, c=d["qvalue"], cmap=cmap,
               norm=norm, edgecolor="#4B5860", linewidth=0.8, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([textwrap.fill(x, 48) for x in d["Description"]], fontsize=12.5)
    ax.set_xlabel("GeneRatio", color=DARK)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(0.035, 0.105)
    for spine in ax.spines.values():
        spine.set_color("#8B979E")
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.06, shrink=0.46, anchor=(0, 0.95))
    # Keep the scientific multiplier in the colorbar label rather than as
    # offset text above the bar, where it would overlap the main panel.
    kegg_ticks = np.linspace(d["qvalue"].min(), d["qvalue"].max(), 4)
    cbar.set_ticks(kegg_ticks)
    cbar.set_ticklabels([f"{x / 1e-6:.1f}" for x in kegg_ticks])
    cbar.ax.yaxis.get_offset_text().set_visible(False)
    cbar.set_label(r"q-value ($\times 10^{-6}$)", fontsize=12)
    cbar.ax.tick_params(labelsize=10.5)
    handles, labels = bubble_legend(ax, d["Count"].to_numpy())
    ax.legend(handles, labels, title="Count", frameon=False, fontsize=10.5,
              title_fontsize=11.5, loc="center left", bbox_to_anchor=(1.20, 0.35),
              labelspacing=0.9, handletextpad=0.8)
    fig.subplots_adjust(left=0.47, right=0.78, bottom=0.11, top=0.98)
    save_multi(fig, "Fig2C_refined_KEGG")


if __name__ == "__main__":
    make_go()
    make_kegg()
