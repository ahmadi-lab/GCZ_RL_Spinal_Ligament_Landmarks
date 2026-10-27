"""LaTeX table bodies of the revision analyses (rows only; the table environments and captions are
in the manuscript). Reads paper/outputs/revision and paper/outputs, writes manuscript/tables/*.tex.

Usage: python paper/revision_tables.py [out_dir]
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REV = ROOT / "paper" / "outputs" / "revision"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "manuscript" / "tables"
LIGS = ["ALL", "PLL", "LF", "ISL", "SSL", "ITL", "CL"]


def f3(x):
    return f"{x:.3f}"


def bold_min(vals, fmt=f3, high=False):
    best = max(vals) if high else min(vals)
    return [f"\\textbf{{{fmt(v)}}}" if abs(v - best) < 1e-12 else fmt(v) for v in vals]


def pfmt(p):
    return "$<$0.001" if p < 0.001 else (f"{p:.3f}" if p < 0.01 else f"{p:.2f}")


def write(name, lines):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text("\n".join(lines) + "\n", encoding="utf-8")


def overall():
    """Primary comparison: matched-pair and matching-independent metrics side by side."""
    n = json.loads((ROOT / "paper" / "outputs" / "paper_numbers.json").read_text())
    rows = []
    for m in ("GCZ", "SLD", "GCZ-RL"):
        d = n[m]
        diff = "--" if m == "GCZ-RL" else (  # GCZ-RL minus this method, with its 95% bootstrap interval
            f"${n[f'GCZRL_vs_{m}']['macro_mae_difference_mm']:.3f}$ "
            f"[${n[f'GCZRL_vs_{m}']['bootstrap_ci95_mm'][0]:.3f}$, ${n[f'GCZRL_vs_{m}']['bootstrap_ci95_mm'][1]:.3f}$]")
        p = "--" if m == "GCZ-RL" else ("$<$0.001" if n[f"GCZRL_vs_{m}"]["wilcoxon_p"] < 0.001
                                        else f"{n[f'GCZRL_vs_{m}']['wilcoxon_p']:.2f}")
        rows.append([m, d["macro_mae"], diff, p, d["pass5_case_avg"], d["closest_point_pred_to_ref"],
                     d["closest_point_ref_to_pred"], d["symmetric_nn"], d["ref_within_5mm_pct"], d["pooled_coverage_pct"]])
    cols = list(zip(*rows))
    mae = bold_min(cols[1])
    p5 = bold_min(cols[4], lambda v: f"{v:.1f}", high=True)
    p2r, r2p, sym = bold_min(cols[5]), bold_min(cols[6]), bold_min(cols[7])
    r5 = bold_min(cols[8], lambda v: f"{v:.1f}", high=True)
    lines = [f"{r[0]} & {mae[i]} & {r[2]} & {r[3]} & {p5[i]} & {p2r[i]} & {r2p[i]} & {sym[i]} & {r5[i]} & {r[9]:.1f} \\\\"
             for i, r in enumerate(rows)]
    write("tab_overall.tex", lines)


def ablation():
    a = pd.read_csv(REV / "table_ablation_components.csv")
    names = {"GCZ": "GCZ (geometry only)",
             "GCZ + projection only (no residual)": "GCZ + projection, no residual",
             "GCZ + raw residual": "GCZ + raw residual",
             "GCZ + bounded residual": "GCZ + bounded residual",
             "GCZ + bounded residual + projection": "GCZ + bounded residual + projection ($\\alpha_L=1$)",
             "GCZ-RL (complete)": "Complete GCZ-RL (nested $\\alpha_L$)"}
    lines = []
    for _, r in a.iterrows():
        p = "--" if pd.isna(r.get("wilcoxon_p_vs_GCZ")) else pfmt(r.wilcoxon_p_vs_GCZ)
        lines.append(f"{names[r.configuration]} & {r.macro_mae_mm:.3f} & {p} & {r.pass5_case_avg_pct:.1f} & "
                     f"{r.ref_to_pred_mm:.3f} & {r.ref_within_5mm_pct:.1f} & {r.MAE_ITL:.3f} & "
                     f"{r.admissible_vertex_pct:.1f} & {r.mean_distance_to_surface_mm:.2f} \\\\")
    write("tab_ablation.tex", lines)
    ns = pd.read_csv(REV / "table_ablation_no_side_constraint.csv").iloc[0]
    lines = [f"{names[r.configuration]} & " + " & ".join(f"{r[f'MAE_{l}']:.3f}" for l in LIGS) + " \\\\"
             for _, r in a.iterrows()]
    lines.append("\\midrule")
    lines.append("GCZ-RL without side constraint (nested re-run) & " + " & ".join(f"{ns[f'MAE_{l}']:.3f}" for l in LIGS) + " \\\\")
    write("tab_ablation_class.tex", lines)


def compactness():
    c = pd.read_csv(REV / "table_region_compactness.csv")
    lines = []
    for _, r in c.iterrows():
        lines.append(f"{r.ligament} & {r.n_ref} & {r.inside_selection_pct:.1f} & {int(r.projection_vertices_median)} & "
                     f"{r.projection_vertex_fraction_pct} & {r.projection_area_fraction_pct} \\\\".replace("-", "--"))
    write("tab_compactness.tex", lines)
    per = pd.read_csv(REV / "region_compactness_per_vertebra.csv")
    p = per[per.zone == "projection"]
    lines = []
    for lig in LIGS:
        x = p[p.ligament == lig]
        cells = []
        for h in ("S01", "S02", "S03"):
            y = x[x.held_out == h]
            cells += [f"{100 * y.n_inside.sum() / y.n_ref.sum():.1f}", f"{y.area_fraction_pct.median():.1f}"]
        lines.append(f"{lig} & " + " & ".join(cells) + " \\\\")
    write("tab_compactness_by_subject.tex", lines)


def atlas():
    a = pd.read_csv(REV / "table_atlas_icp.csv")
    lines = []
    for _, r in a.iterrows():
        lines.append(f"{r.target.replace('_', ' ')} & {r.atlas.replace('_', ' ')} & {'same' if r.same_level else 'cross'} & "
                     f"{r.icp_residual_mm:.2f} & {r.GCZ_mm:.3f} & {r.GCZRL_mm:.3f} & "
                     f"{r.second_env_atlas.replace('_', ' ')}{'$^{*}$' if r.atlas_changed else ''} & "
                     f"{r.second_env_residual_mm:.2f} & {r.second_env_GCZ_mm:.3f} & {r.second_env_GCZRL_mm:.3f} \\\\")
    write("tab_atlas.tex", lines)


def residuals():
    t = pd.read_csv(REV / "table_residual_agreement.csv")
    lines = [f"{r.ligament} & {r.points} & {r.true_residual_mean_mm:.2f} & {r.pred_clipped_mean_mm:.2f} & "
             f"{r.median_cosine:.2f} & {r.pct_cosine_positive:.1f} & {r.pct_closer_after_clipped:.1f} \\\\"
             for _, r in t.iterrows()]
    write("tab_residual_agreement.tex", lines)
    perm = pd.read_csv(REV / "feature_importance_permutation.csv")
    mdi = pd.read_csv(REV / "feature_importance_impurity.csv")
    g = perm.groupby("group").increase_mm
    order = g.mean().sort_values(ascending=False).index
    lines = [f"{k} & {g.mean()[k]:.3f} & {g.min()[k]:.3f}--{g.max()[k]:.3f} & {mdi[k].mean():.3f} \\\\" for k in order]
    write("tab_importance.tex", lines)


if __name__ == "__main__":
    overall()
    ablation()
    compactness()
    atlas()
    residuals()
    print("tables written to", OUT)
