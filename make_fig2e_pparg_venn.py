from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"
OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "candidate": "#7FB685",
    "microbiota": "#9EC5E5",
    "dark": "#243746",
    "overlap": "#4E8F96",
}


def build_figure():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 13,
        "axes.titleweight": "bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(8.4, 4.8), facecolor="white")
    ax.set_aspect("equal")
    ax.axis("off")

    radius = 1.48
    left_center = (2.42, 2.02)
    right_center = (4.18, 2.02)

    ax.add_patch(Circle(left_center, radius, facecolor=COLORS["candidate"],
                        edgecolor="none", alpha=0.88))
    ax.add_patch(Circle(right_center, radius, facecolor=COLORS["microbiota"],
                        edgecolor="none", alpha=0.88))

    # Set names include their total sizes so the diagram remains compact.
    ax.text(1.87, 3.66, "Candidate genes\n(n = 24)", ha="center", va="bottom",
            fontsize=15, fontweight="bold", color=COLORS["candidate"])
    ax.text(4.73, 3.66, "Gut microbiota-related genes\n(n = 117)",
            ha="center", va="bottom", fontsize=15, fontweight="bold",
            color=COLORS["microbiota"])

    ax.text(1.60, 2.02, "23", ha="center", va="center", fontsize=18,
            fontweight="bold", color=COLORS["dark"])
    ax.text(5.00, 2.02, "116", ha="center", va="center", fontsize=18,
            fontweight="bold", color=COLORS["dark"])
    ax.text(3.30, 2.10, "PPARG", ha="center", va="center", fontsize=18,
            fontweight="bold", color=COLORS["dark"])
    ax.text(3.30, 1.78, "n = 1", ha="center", va="center", fontsize=12,
            color=COLORS["dark"])

    ax.set_xlim(0.45, 6.15)
    ax.set_ylim(0.30, 4.05)
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.04, top=0.98)
    return fig


def save():
    fig = build_figure()
    stem = "Fig2E_PPARG_candidate_GM_venn"
    pdf = OUT / f"{stem}.pdf"
    png = OUT / f"{stem}.png"
    tif = OUT / f"{stem}.tif"

    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(tif, dpi=1200, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(MIRROR / f"{stem}.png", dpi=1200, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"Wrote {pdf.name}, {png.name}, and {tif.name}")


if __name__ == "__main__":
    save()
