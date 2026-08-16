from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"
OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)

GREEN = "#7FB685"
BLUE = "#9EC5E5"
PINK = "#F29AA0"
DARK = "#243746"
ARROW = "#6B7780"


def setup_axis(ax, xlim, ylim):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_three_set(ax):
    setup_axis(ax, (0.0, 7.0), (0.0, 5.45))
    radius = 1.52
    centers = {
        "glycolysis": (2.25, 3.18),
        "wgcna": (4.75, 3.18),
        "degs": (3.50, 1.78),
    }
    for key, color in (("glycolysis", GREEN), ("wgcna", BLUE), ("degs", PINK)):
        ax.add_patch(Circle(centers[key], radius, facecolor=color,
                            edgecolor="none", alpha=0.66, zorder=1))

    ax.text(1.35, 5.16, "Glycolysis-related genes\n(n = 911)", ha="center",
            va="center", fontsize=14, fontweight="bold", color=GREEN, zorder=5)
    ax.text(5.65, 5.16, "WGCNA genes\n(n = 294)", ha="center",
            va="center", fontsize=14, fontweight="bold", color=BLUE, zorder=5)
    ax.text(3.50, 0.04, "DEGs (n = 681)", ha="center", va="bottom",
            fontsize=14, fontweight="bold", color=PINK, zorder=5)

    # Region counts reproduce the source three-set intersection exactly.
    labels = [
        (1.37, 3.43, "857"),
        (5.63, 3.43, "23"),
        (3.50, 1.05, "380"),
        (3.50, 3.64, "0"),
        (2.43, 2.30, "30"),
        (4.57, 2.30, "247"),
    ]
    for x, y, value in labels:
        ax.text(x, y, value, ha="center", va="center", fontsize=15,
                fontweight="bold", color=DARK, zorder=5)

    ax.text(3.50, 2.57, "24", ha="center", va="center", fontsize=18,
            fontweight="bold", color=DARK,
            bbox=dict(boxstyle="circle,pad=0.34", facecolor="white",
                      edgecolor=DARK, linewidth=1.4, alpha=0.96))


def draw_transition(ax):
    setup_axis(ax, (0.0, 7.0), (0.0, 1.0))
    arrow = FancyArrowPatch((3.50, 0.92), (3.50, 0.08), arrowstyle="-|>",
                            mutation_scale=20, linewidth=1.7, color=ARROW)
    ax.add_patch(arrow)
    ax.text(3.72, 0.50, "24 candidate genes", ha="left", va="center",
            fontsize=13, fontweight="bold", color=DARK)


def draw_two_set(ax):
    setup_axis(ax, (0.0, 7.0), (0.0, 3.65))
    radius = 1.28
    left = (2.60, 1.62)
    right = (4.40, 1.62)
    ax.add_patch(Circle(left, radius, facecolor=GREEN, edgecolor="none", alpha=0.66,
                        zorder=1))
    ax.add_patch(Circle(right, radius, facecolor=BLUE, edgecolor="none", alpha=0.66,
                        zorder=1))

    ax.text(1.70, 3.48, "Candidate genes\n(n = 24)", ha="center", va="top",
            fontsize=16, fontweight="bold", color=GREEN, zorder=5)
    ax.text(5.30, 3.48, "Gut microbiota-related genes\n(n = 117)",
            ha="center", va="top", fontsize=16, fontweight="bold", color=BLUE, zorder=5)
    ax.text(1.82, 1.62, "23", ha="center", va="center", fontsize=17,
            fontweight="bold", color=DARK)
    ax.text(5.18, 1.62, "116", ha="center", va="center", fontsize=17,
            fontweight="bold", color=DARK)
    ax.text(3.50, 1.72, "PPARG", ha="center", va="center", fontsize=21,
            fontweight="bold", color=DARK, zorder=5)
    ax.text(3.50, 1.35, "n = 1", ha="center", va="center", fontsize=13.5,
            color=DARK, zorder=5)


def build_figure():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 12,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.2, 9.5), facecolor="white")
    grid = fig.add_gridspec(3, 1, height_ratios=[5.45, 0.85, 3.65], hspace=0.00)
    draw_three_set(fig.add_subplot(grid[0]))
    draw_transition(fig.add_subplot(grid[1]))
    draw_two_set(fig.add_subplot(grid[2]))
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.99)
    return fig


def save():
    fig = build_figure()
    stem = "Fig2E_combined_candidate_convergence_PPARG_venn"
    outputs = {
        "pdf": OUT / f"{stem}.pdf",
        "png": OUT / f"{stem}.png",
        "tif": OUT / f"{stem}.tif",
    }
    fig.savefig(outputs["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(outputs["png"], dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs["tif"], dpi=1200, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(MIRROR / f"{stem}.png", dpi=1200, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("Wrote " + ", ".join(path.name for path in outputs.values()))


if __name__ == "__main__":
    save()
