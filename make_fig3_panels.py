from __future__ import annotations

from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parent
ML = ROOT / "6.ML_reanalysis_20260809"
OUT = ROOT / "submission_package" / "figures" / "standalone"
OUT.mkdir(parents=True, exist_ok=True)
COPY_PNG = ROOT / "图片1200dpi"
COPY_PNG.mkdir(parents=True, exist_ok=True)

COLORS = {
    "dark": "#233746",
    "muted": "#566A78",
    "grid": "#E8EEF3",
    "teal": "#147D7E",
    "coral": "#E76F51",
    "blue": "#3B82B8",
    "grey": "#A9B5BE",
    "light": "#F4F7FA",
}

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.labelcolor": COLORS["dark"],
    "axes.edgecolor": COLORS["dark"],
    "xtick.color": COLORS["dark"],
    "ytick.color": COLORS["dark"],
    "text.color": COLORS["dark"],
    "axes.titleweight": "bold",
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color=COLORS["grid"], linewidth=0.8)
    ax.set_axisbelow(True)


def save_multi(fig, stem: str):
    png = OUT / f"{stem}.png"
    tif = OUT / f"{stem}.tif"
    pdf = OUT / f"{stem}.pdf"
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(tif, dpi=600, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    # Keep a high-resolution PNG mirror with the other 1200-dpi figure assets.
    fig.savefig(COPY_PNG / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {stem}: {pdf.name}, {tif.name}, {png.name}")


def fig3a():
    d = pd.read_csv(ML / "workflow_summary.csv")
    d = d[np.isfinite(d["mean_auc"]) & np.isfinite(d["min_auc"])].copy()
    selected = d["workflow"].eq("all + glmnet_0")
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    norm = Normalize(vmin=float(d["sd_auc"].min()), vmax=float(d["sd_auc"].max()))
    sc = ax.scatter(
        d["mean_auc"], d["min_auc"], c=d["sd_auc"], cmap="cividis", norm=norm,
        s=54, alpha=0.82, edgecolors="white", linewidths=0.35, zorder=2,
    )
    ax.scatter(
        d.loc[selected, "mean_auc"], d.loc[selected, "min_auc"],
        s=170, marker="*", color=COLORS["coral"], edgecolors=COLORS["dark"],
        linewidths=0.8, zorder=4,
    )
    row = d.loc[selected].iloc[0]
    ax.annotate(
        "Selected Ridge\nmean AUC = %.3f\nminimum AUC = %.3f" % (row["mean_auc"], row["min_auc"]),
        xy=(row["mean_auc"], row["min_auc"]), xytext=(0.947, 0.992), textcoords="data",
        fontsize=11, color=COLORS["dark"], ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=COLORS["coral"], linewidth=0.8),
        arrowprops=dict(arrowstyle="->", color=COLORS["coral"], linewidth=1.0,
                        connectionstyle="arc3,rad=-0.12"),
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.025, fraction=0.045)
    cbar.set_label("Across-cohort SD\n(lower = more stable)", color=COLORS["dark"])
    cbar.ax.tick_params(colors=COLORS["dark"])
    ax.set_xlim(0.943, 0.999)
    ax.set_ylim(0.862, 1.003)
    ax.set_xlabel("Mean LOCO AUC")
    ax.set_ylabel("Minimum cohort AUC")
    ax.set_title("117-workflow performance landscape", loc="center", pad=12, fontsize=16)
    style_ax(ax)
    fig.tight_layout()
    save_multi(fig, "Fig3A_117_workflow_performance_landscape")


def fig3b():
    d = pd.read_csv(ML / "workflow_auc_by_outer_cohort.csv")
    workflows = [
        ("all + glmnet_0", "24-gene Ridge"),
        ("all + glmnet_1", "LASSO"),
        ("all + glmnet_0.5", "Elastic Net"),
        ("all + svm_linear", "Linear SVM"),
        ("all + rf", "Random forest"),
    ]
    keep = [x[0] for x in workflows]
    labels = dict(workflows)
    d = d[d["workflow"].isin(keep)].copy()
    d["model"] = d["workflow"].map(labels)
    model_order = [x[1] for x in workflows]
    cohorts = ["GSE73661", "GSE75214", "GSE87466", "GSE107499"]
    palette = {
        "24-gene Ridge": COLORS["teal"], "LASSO": "#6E8797", "Elastic Net": "#8AA1AF",
        "Linear SVM": "#A5B4BE", "Random forest": "#C0CBD1",
    }
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 4.9), sharex=True, sharey=True)
    for ax, cohort in zip(axes, cohorts):
        sub = d[d["heldout_cohort"].eq(cohort)].set_index("model").reindex(model_order).reset_index()
        y = np.arange(len(sub))[::-1]
        ax.hlines(y, 0.90, sub["auc"], color=COLORS["grid"], linewidth=2.3, zorder=1)
        ax.scatter(sub["auc"], y, s=52, c=[palette[m] for m in sub["model"]],
                   edgecolors="white", linewidths=0.6, zorder=3)
        for yi, value, model in zip(y, sub["auc"], sub["model"]):
            ax.text(value + 0.0013, yi, f"{value:.3f}", va="center", ha="left", fontsize=10,
                    color=COLORS["dark"])
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.set_xlim(0.90, 1.012)
        ax.set_ylim(-0.7, len(model_order) - 0.3)
        ax.set_title(cohort, fontsize=13, pad=10)
        style_ax(ax)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticklabels(model_order)
    axes[0].tick_params(axis="y", labelleft=True)
    axes[0].set_ylabel("Model")
    fig.supxlabel("AUC in the held-out discovery cohort", y=0.03)
    fig.suptitle("24-gene Ridge and transparent baselines across held-out cohorts", x=0.5, ha="center",
                 fontsize=16, fontweight="bold", y=1.04)
    fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[m],
                               markeredgecolor="white", markersize=8, label=m) for m in model_order],
               loc="upper center", bbox_to_anchor=(0.72, 1.00), ncol=3, frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.20, right=0.995, bottom=0.17, top=0.79, wspace=0.10)
    save_multi(fig, "Fig3B_heldout_cohort_model_auc")


def bootstrap_auc_ci(y, score, n_boot=5000, seed=20260813):
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    observed = roc_auc_score(y, score)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        ii = np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]
        try:
            values.append(roc_auc_score(y[ii], score[ii]))
        except ValueError:
            pass
    low, high = np.quantile(values, [0.025, 0.975])
    return float(observed), float(low), float(high)


def fit_ridge_cv(x, y, seed=20260813):
    scaler = StandardScaler().fit(x)
    z = scaler.transform(x)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    Cs = np.logspace(-3, 3, 61)
    best_c, best_score = None, -np.inf
    for C in Cs:
        fold_scores = []
        for train, test in cv.split(z, y):
            model = LogisticRegression(C=float(C), solver="liblinear", max_iter=5000)
            model.fit(z[train], y[train])
            fold_scores.append(roc_auc_score(y[test], model.predict_proba(z[test])[:, 1]))
        score = float(np.mean(fold_scores))
        if score > best_score:
            best_score, best_c = score, float(C)
    model = LogisticRegression(C=best_c, solver="liblinear", max_iter=5000)
    model.fit(z, y)
    return model, scaler, best_c, best_score


def external_model_data():
    expr = pd.read_csv(ML / "rank_normalized_candidate_expression.csv")
    discovery = expr[expr["set"].eq("discovery")].copy()
    external = expr[expr["set"].eq("external")].copy()
    genes = [c for c in expr.columns if c not in {"sample", "cohort", "set", "label", "subgroup"}]
    stable = pd.read_csv(ML / "gene_stability_summary.csv")
    stable_genes = stable.loc[stable["stable_core"], "gene"].tolist()

    # Use the locked predictions produced by the original leakage-controlled run.
    locked = pd.read_csv(ML / "external_predictions.csv")
    locked = locked[["sample", "cohort", "label", "score"]].rename(columns={"score": "8-gene Ridge"})

    ridge24, scaler24, best_c, cv_auc = fit_ridge_cv(discovery[genes].to_numpy(), discovery["label"].to_numpy())
    ext24 = ridge24.predict_proba(scaler24.transform(external[genes]))[:, 1]

    # Orient PPARG in the discovery direction, then fit a univariate logistic calibration model.
    direction = stable.loc[stable["gene"].eq("PPARG"), "discovery_direction"].iloc[0]
    sign = 1.0 if direction == "UC_up" else -1.0
    p_train = (sign * discovery["PPARG"].to_numpy()).reshape(-1, 1)
    p_ext = (sign * external["PPARG"].to_numpy()).reshape(-1, 1)
    p_scaler = StandardScaler().fit(p_train)
    p_model = LogisticRegression(C=1e6, solver="liblinear", max_iter=5000)
    p_model.fit(p_scaler.transform(p_train), discovery["label"].to_numpy())
    ext_pparg = p_model.predict_proba(p_scaler.transform(p_ext))[:, 1]

    scores = locked.merge(external[["sample", "cohort", "label"]], on=["sample", "cohort", "label"], how="left")
    # Attach the model scores by sample without relying on row order.
    score_map24 = dict(zip(external["sample"], ext24))
    score_map_pparg = dict(zip(external["sample"], ext_pparg))
    scores["24-gene Ridge"] = scores["sample"].map(score_map24)
    scores["PPARG single gene"] = scores["sample"].map(score_map_pparg)
    return scores, stable_genes, {"best_C": best_c, "cv_auc": cv_auc}


def fig3e():
    scores, stable_genes, fit_info = external_model_data()
    models = ["8-gene Ridge", "24-gene Ridge", "PPARG single gene"]
    cohorts = ["GSE47908", "GSE13367"]
    rows = []
    for cohort in cohorts:
        sub = scores[scores["cohort"].eq(cohort)]
        for model in models:
            auc, low, high = bootstrap_auc_ci(sub["label"], sub[model], seed=20260813 + len(rows))
            rows.append({"cohort": cohort, "model": model, "n": len(sub), "auc": auc,
                         "ci_low": low, "ci_high": high,
                         "brier": float(np.mean((sub[model] - sub["label"]) ** 2)),
                         "ci_method": "bootstrap percentile, 5000 resamples"})
    metrics = pd.DataFrame(rows)
    metrics["fit_note"] = f"24-gene ridge selected C={fit_info['best_C']:.6g}; 5-fold discovery CV AUC={fit_info['cv_auc']:.3f}"
    metrics.to_csv(OUT / "Fig3E_external_model_metrics.csv", index=False)

    palette = {"8-gene Ridge": COLORS["teal"], "24-gene Ridge": COLORS["coral"], "PPARG single gene": COLORS["blue"]}
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.0), sharey=True)
    for ax, cohort in zip(axes, cohorts):
        sub = metrics[metrics["cohort"].eq(cohort)].set_index("model").reindex(models).reset_index()
        y = np.arange(len(models))[::-1]
        ax.axvline(0.5, color=COLORS["grid"], linewidth=0.8, zorder=0)
        for yi, row in zip(y, sub.to_dict("records")):
            ax.errorbar(row["auc"], yi, xerr=[[row["auc"] - row["ci_low"]], [row["ci_high"] - row["auc"]]],
                        fmt="o", color=palette[row["model"]], markeredgecolor="white", markeredgewidth=0.7,
                        capsize=3.5, linewidth=1.7, markersize=6.5, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([])
        # Keep the lower CI endpoint of GSE13367/PPARG visible instead of
        # clipping it at the 0.50 axis boundary.
        ax.set_xlim(0.45, 1.08)
        ax.set_xticks(np.arange(0.5, 1.01, 0.1))
        ax.set_ylim(-0.7, len(models) - 0.3)
        ax.set_title(cohort, fontsize=13, pad=10)
        style_ax(ax)
        for yi, row in zip(y, sub.to_dict("records")):
            ax.text(1.02, yi, f"AUC {row['auc']:.3f}\nBrier {row['brier']:.3f}",
                    transform=ax.get_yaxis_transform(), va="center", ha="left",
                    fontsize=10, color=COLORS["dark"], clip_on=False)
    axes[0].set_yticklabels(models)
    axes[0].tick_params(axis="y", labelleft=True)
    axes[0].set_ylabel("External model")
    fig.supxlabel("AUC with 95% CI", y=0.045)
    fig.suptitle("External model comparison", x=0.5, ha="center", fontsize=16, fontweight="bold", y=1.04)
    fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[m],
                               markeredgecolor="white", markersize=8, label=m) for m in models],
               loc="upper center", bbox_to_anchor=(0.75, 1.0), ncol=3, frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.22, right=0.995, bottom=0.18, top=0.78, wspace=0.24)
    save_multi(fig, "Fig3E_external_model_comparison")
    return scores


def calibration_metrics(scores):
    rows = []
    for cohort, sub in scores.groupby("cohort"):
        y = sub["label"].to_numpy(dtype=float)
        p = np.clip(sub["8-gene Ridge"].to_numpy(dtype=float), 1e-6, 1 - 1e-6)
        lp = np.log(p / (1 - p))
        intercept_fit = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit()
        slope_fit = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit()
        i_ci = intercept_fit.conf_int()[0]
        s_ci = slope_fit.conf_int()[1]
        rows.extend([
            {"cohort": cohort, "metric": "Calibration-in-the-large", "estimate": float(intercept_fit.params[0]),
             "ci_low": float(i_ci[0]), "ci_high": float(i_ci[1]), "reference": 0.0,
             "method": "Binomial GLM with logit prediction offset; Wald 95% CI"},
            {"cohort": cohort, "metric": "Calibration slope", "estimate": float(slope_fit.params[1]),
             "ci_low": float(s_ci[0]), "ci_high": float(s_ci[1]), "reference": 1.0,
             "method": "Binomial GLM of outcome on logit prediction; Wald 95% CI"},
        ])
    return pd.DataFrame(rows)


def fig3f(scores):
    metrics = calibration_metrics(scores)
    metrics.to_csv(OUT / "Fig3F_external_calibration_metrics.csv", index=False)
    cohorts = ["GSE47908", "GSE13367"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.5), sharey=True)
    specs = [("Calibration-in-the-large", 0.0, (-1.6, 2.4)), ("Calibration slope", 1.0, (0.0, 2.8))]
    for ax, (metric, ref, xlim) in zip(axes, specs):
        sub = metrics[metrics["metric"].eq(metric)].set_index("cohort").reindex(cohorts).reset_index()
        y = np.arange(len(sub))[::-1]
        label_x = float(sub["ci_high"].max()) + 0.12
        ax.axvline(ref, color=COLORS["grey"], linewidth=1.0, linestyle="--", zorder=0)
        for yi, row in zip(y, sub.to_dict("records")):
            ax.errorbar(row["estimate"], yi,
                        xerr=[[row["estimate"] - row["ci_low"]], [row["ci_high"] - row["estimate"]]],
                        fmt="o", color=COLORS["teal"], markeredgecolor="white", markeredgewidth=0.7,
                        capsize=4, linewidth=1.8, markersize=7, zorder=3)
            ax.text(label_x, yi,
                    f"{row['estimate']:.2f} [{row['ci_low']:.2f}, {row['ci_high']:.2f}]",
                    va="center", ha="left", fontsize=10, color=COLORS["dark"])
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.6, len(cohorts) - 0.4)
        ax.set_title(metric, fontsize=13, pad=10)
        ax.set_xlabel("Estimate (95% CI)")
        style_ax(ax)
    axes[0].set_yticklabels(cohorts)
    axes[0].tick_params(axis="y", labelleft=True)
    axes[0].set_ylabel("Locked external cohort")
    fig.suptitle("Simple calibration of the locked 8-gene Ridge", x=0.5, ha="center", fontsize=16,
                 fontweight="bold", y=1.04)
    # Reserve a clear central gutter for the enlarged numeric CI labels.
    fig.subplots_adjust(left=0.17, right=0.99, bottom=0.23, top=0.77, wspace=0.34)
    save_multi(fig, "Fig3F_external_calibration")


def fig3c():
    d = pd.read_csv(ML / "gene_stability_summary.csv")
    d = d.sort_values("selection_frequency", ascending=True).copy()
    d["stable_core"] = d["stable_core"].astype(bool)
    fig, ax = plt.subplots(figsize=(7.1, 7.0))
    y = np.arange(len(d))
    colors = np.where(d["stable_core"], COLORS["teal"], COLORS["grey"])
    ax.barh(y, d["selection_frequency"], color=colors, height=0.68, edgecolor="none")
    ax.axvline(0.5, color=COLORS["coral"], linestyle="--", linewidth=1.2)
    for yi, value in zip(y, d["selection_frequency"]):
        ax.text(min(value + 0.018, 1.02), yi, f"{value * 100:.1f}%", va="center", ha="left",
                fontsize=10, color=COLORS["dark"])
    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"], fontsize=10.5)
    ax.set_xlim(0, 1.10)
    ax.set_xlabel("Selection frequency across sparse selectors")
    ax.set_ylabel("Candidate gene")
    ax.set_title("Stable gene selection frequency", loc="center", pad=12, fontsize=16)
    ax.legend(handles=[
        Patch(facecolor=COLORS["teal"], label="Stable core (8 genes)"),
        Patch(facecolor=COLORS["grey"], label="Other candidates"),
        Line2D([0], [0], color=COLORS["coral"], linestyle="--", label="Selection threshold = 50%"),
    ], loc="lower right", frameon=False, fontsize=10)
    style_ax(ax)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    save_multi(fig, "Fig3C_gene_selection_frequency")


def fig3d():
    from sklearn.metrics import roc_curve

    pred = pd.read_csv(ML / "external_predictions.csv")
    auc_table = pd.read_csv(ML / "external_validation_auc.csv").set_index("cohort")
    colors = {"GSE47908": COLORS["teal"], "GSE13367": COLORS["coral"]}
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    handles = []
    for cohort in ["GSE47908", "GSE13367"]:
        sub = pred[pred["cohort"].eq(cohort)]
        fpr, tpr, _ = roc_curve(sub["label"], sub["score"])
        row = auc_table.loc[cohort]
        line, = ax.plot(fpr, tpr, color=colors[cohort], linewidth=2.4,
                        label=(f"{cohort} (n={int(row['n'])}); "
                               f"AUC {row['auc']:.3f} [{row['ci_low']:.3f}, {row['ci_high']:.3f}]"))
        handles.append(line)
    ax.plot([0, 1], [0, 1], linestyle="--", color=COLORS["grey"], linewidth=1.0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("1 - specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_title("Locked external validation", loc="center", pad=12, fontsize=16)
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=10,
              handlelength=2.0, borderaxespad=0.3)
    style_ax(ax)
    fig.tight_layout()
    save_multi(fig, "Fig3D_locked_external_roc")


def fig3d_summary_table():
    pred = pd.read_csv(ML / "external_predictions.csv")
    auc_table = pd.read_csv(ML / "external_validation_auc.csv")
    rows = []
    for _, row in auc_table.iterrows():
        sub = pred[pred["cohort"].eq(row["cohort"])]
        rows.append([
            row["cohort"], int(row["n"]), f"{row['auc']:.3f}",
            f"{row['ci_low']:.3f}–{row['ci_high']:.3f}",
            f"{np.mean((sub['score'] - sub['label']) ** 2):.3f}",
        ])
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    ax.axis("off")
    ax.set_title("External performance summary", loc="center", pad=16, fontsize=14, fontweight="bold")
    table = ax.table(
        cellText=rows,
        colLabels=["Cohort", "N", "AUC", "95% CI", "Brier score"],
        cellLoc="center", colLoc="center", loc="center",
        colWidths=[0.24, 0.12, 0.16, 0.28, 0.20],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    table.scale(1.0, 1.85)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_linewidth(0.6)
        if r == 0:
            cell.set_facecolor(COLORS["teal"])
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("white" if r % 2 else COLORS["light"])
            cell.get_text().set_color(COLORS["dark"])
    fig.tight_layout()
    save_multi(fig, "Fig3D_external_performance_summary_table")


if __name__ == "__main__":
    fig3a()
    fig3b()
    fig3c()
    fig3d()
    fig3d_summary_table()
    scores = fig3e()
    fig3f(scores)
    print("Stable genes in locked model:", ", ".join(pd.read_csv(ML / "gene_stability_summary.csv").query("stable_core")["gene"].tolist()))
