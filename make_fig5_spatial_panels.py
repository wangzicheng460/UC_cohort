from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin  # noqa: F401
from matplotlib import colormaps


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "空转" / "GSE189184_spatial_validation_20260810" / "figures_individual_png"
OUT = ROOT / "submission_package" / "figures" / "standalone"
MIRROR = ROOT / "图片1200dpi"

OUT.mkdir(parents=True, exist_ok=True)
MIRROR.mkdir(parents=True, exist_ok=True)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(str(Path("C:/Windows/Fonts") / name), size=size)


def centered_text(draw, box, text, text_font, fill="#23384A"):
    x0, y0, x1, y1 = box
    bounds = draw.textbbox((0, 0), text, font=text_font)
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    x = x0 + ((x1 - x0) - width) / 2
    y = y0 + ((y1 - y0) - height) / 2 - bounds[1]
    draw.text((x, y), text, font=text_font, fill=fill)


PANELS = [
    {
        "letter": "A",
        "stem": "SpatialAnnotations",
        "title": "Spatial domain annotation",
        "top": 160,
        "left_right": 1495,
        "plot_right": 1605,
        "legend": "domains",
    },
    {
        "letter": "B",
        "stem": "EpithelialFraction",
        "title": "Epithelial fraction",
        "top": 100,
        "left_right": 1550,
        "plot_right": 1605,
        "legend": "continuous",
        "ticks": [0.00, 0.25, 0.50, 0.75, 1.00],
        "tick_decimals": 2,
    },
    {
        "letter": "C",
        "stem": "PPARGActivity",
        "title": "PPARG regulon activity",
        "top": 100,
        "left_right": 1660,
        "plot_right": 1605,
        "legend": "continuous",
        "ticks": [-2, -1, 0, 1, 2, 3],
        "tick_decimals": 0,
    },
    {
        "letter": "D",
        "stem": "GlycolysisActivity",
        "title": "Hallmark glycolysis activity",
        "top": 125,
        "left_right": 1660,
        "plot_right": 1605,
        "legend": "continuous",
        "ticks": [0.025, 0.050, 0.075, 0.100, 0.125],
        "tick_decimals": 3,
    },
]


def draw_domain_legend(draw, x, center_y):
    entries = [
        ("Epithelial-rich", "#F28E2B"),
        ("Plasma-rich", "#EDC948"),
        ("Stromal-rich", "#76B7B2"),
        ("Endothelial-rich", "#59A14F"),
    ]
    text_font = font(32)
    row_height = 62
    y = center_y - row_height * len(entries) / 2
    for label, color in entries:
        cy = round(y + row_height / 2)
        draw.ellipse((x, cy - 10, x + 20, cy + 10), fill=color)
        draw.text((x + 42, cy), label, font=text_font, fill="#111111", anchor="lm")
        y += row_height


def format_tick(value, decimals):
    return f"{value:.{decimals}f}"


def draw_continuous_legend(draw, x, center_y, ticks, decimals):
    bar_width = 42
    bar_height = 360
    top = round(center_y - bar_height / 2)
    cmap = colormaps["plasma"]
    for i in range(bar_height):
        t = 1 - i / (bar_height - 1)
        rgba = cmap(t)
        color = tuple(round(channel * 255) for channel in rgba[:3])
        draw.line((x, top + i, x + bar_width, top + i), fill=color, width=1)
    draw.rectangle((x, top, x + bar_width, top + bar_height), outline="#666666", width=1)
    tick_font = font(31)
    low, high = min(ticks), max(ticks)
    for value in ticks:
        frac = (value - low) / (high - low)
        y = top + bar_height - frac * bar_height
        draw.line((x + bar_width, y, x + bar_width + 12, y), fill="white", width=2)
        draw.text(
            (x + bar_width + 24, y), format_tick(value, decimals),
            font=tick_font, fill="#111111", anchor="lm"
        )


def make_panel(spec):
    hc = Image.open(SOURCE / f"{spec['stem']}_C5.png").convert("RGB")
    uc = Image.open(SOURCE / f"{spec['stem']}_B4.png").convert("RGB")

    hc = hc.crop((0, spec["top"], spec["left_right"], hc.height))
    uc = uc.crop((0, spec["top"], spec["plot_right"], uc.height))

    target_height = 1500
    hc = hc.resize((round(hc.width * target_height / hc.height), target_height), Image.Resampling.LANCZOS)
    uc = uc.resize((round(uc.width * target_height / uc.height), target_height), Image.Resampling.LANCZOS)

    margin = 44
    gap = 150
    title_height = 92
    group_height = 66
    bottom = 30
    legend_width = 430 if spec["legend"] == "domains" else 250
    width = margin * 2 + hc.width + gap + uc.width + legend_width
    height = title_height + group_height + target_height + bottom
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    title_font = font(55, bold=True)
    group_font = font(43, bold=True)
    centered_text(draw, (margin, 0, width - margin, title_height), spec["title"], title_font)

    hc_x = margin
    uc_x = margin + hc.width + gap
    centered_text(draw, (hc_x, title_height, hc_x + hc.width, title_height + group_height), "C5 (HC)", group_font)
    centered_text(draw, (uc_x, title_height, uc_x + uc.width, title_height + group_height), "B4 (UC)", group_font)

    y = title_height + group_height
    canvas.paste(hc, (hc_x, y))
    canvas.paste(uc, (uc_x, y))
    legend_x = uc_x + uc.width + 34
    legend_center_y = y + target_height / 2
    if spec["legend"] == "domains":
        draw_domain_legend(draw, legend_x, legend_center_y)
    else:
        draw_continuous_legend(
            draw, legend_x, legend_center_y, spec["ticks"], spec["tick_decimals"]
        )

    base = f"Fig5{spec['letter']}_{spec['stem']}"
    png = OUT / f"{base}.png"
    pdf = OUT / f"{base}.pdf"
    tif = OUT / f"{base}.tif"
    canvas.save(png, dpi=(600, 600), optimize=True)
    canvas.save(pdf, resolution=600.0)
    canvas.save(tif, dpi=(600, 600), compression="tiff_lzw")
    canvas.save(MIRROR / f"{base}.png", dpi=(600, 600), optimize=True)
    print(png)


for panel in PANELS:
    make_panel(panel)
