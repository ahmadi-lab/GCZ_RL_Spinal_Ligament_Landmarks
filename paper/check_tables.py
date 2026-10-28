"""Check every cell of the tables typed directly in manuscript/main.tex and supplementary.tex against
the code outputs.

Tables 6-8 and S10-S18 are written by the code (manuscript/tables/*.tex, copied and checked by
paper/build_manuscript_assets.py). This script covers the tables typed in the LaTeX source: Tables 1-5
and 9, and Supplementary Tables S2-S4 and S6-S9. Each cell is rebuilt from the output files with the
rounding of the manuscript and compared with the printed cell, and the bold cells are compared with
the best values of their row. Tables S1 and S5 list features and settings and hold no results.

Usage: python paper/check_tables.py        Writes paper/outputs/revision/table_cell_check.csv;
                                            exit code 1 if any cell differs.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
P1 = REPO / "part1_dissertation_57" / "outputs"
P2 = REPO / "part2_patient_loso_14" / "outputs"
NB = REPO / "part1_dissertation_57" / "original_notebook" / "outputs" / "tables"
OUT = REPO / "paper" / "outputs"
MAIN = (REPO / "manuscript" / "main.tex").read_text(encoding="utf-8")
SUPP = (REPO / "manuscript" / "supplementary.tex").read_text(encoding="utf-8")
LIGS = ["ALL", "PLL", "LF", "ISL", "SSL", "ITL", "CL"]


def f(x, d=3):
    return f"{x:.{d}f}"


def clean(c):
    while True:
        n = re.sub(r"\\(?:textbf|rev|textit)\{([^{}]*)\}", r"\1", c)
        if n == c:
            break
        c = n
    return " ".join(c.replace("\\%", "%").replace("~", " ").replace("$", "").split())


def printed(tex, label):
    """Data rows of the table with this label: {row label: (cells, bold flags)}."""
    i = tex.index(f"\\label{{{label}}}")
    end = min(j for j in (tex.find("\\end{tabular}", i), tex.find("\\end{longtable}", i)) if j > 0)
    body = tex[tex.index("\\midrule", i):end]
    rows = {}
    for seg in body.split("\\\\"):
        seg = re.sub(r"\\(?:mid|bottom|top)rule|\\cmidrule(?:\([a-z]*\))?\{[^}]*\}", "", seg).strip()
        if not seg or "\\multicolumn" in seg:
            continue
        cells = [c.strip() for c in seg.split("&")]
        rows[clean(cells[0])] = ([clean(c) for c in cells[1:]], ["\\textbf" in c for c in cells[1:]])
    return rows


def best(cells, groups):
    """Indices that should be bold: the lowest (min) or highest (max) printed value of each group."""
    out = set()
    for cols, how in groups:
        v = {k: float(cells[k].split()[0].rstrip("%")) for k in cols}
        t = (min if how == "min" else max)(v.values())
        out |= {k for k, x in v.items() if x == t}
    return out


def expected_tables():
    m = pd.read_csv(P2 / "per_vertebra_per_ligament_metrics.csv")
    nums = json.loads((OUT / "paper_numbers.json").read_text())
    folds = pd.read_csv(P2 / "table_P2_1_folds.csv").set_index("held_out_subject")
    ds = pd.read_csv(P1 / "table_4_1_dataset.csv").set_index("split")
    T = {}

    # Table 1: counts of the study design
    gcz = m[m.method == "GCZ"]
    T["tab:dataset"] = ("main", {
        "Annotated vertebrae": [str(ds.loc["Total", "vertebrae"]), None],
        "Patient-level cohort (primary, LOSO)": [str(folds.test_vertebrae.sum()), None],
        "Independent subjects (primary)": [str(len(folds)), None],
        "Vertebra--class cases (primary)": [str(len(gcz)), None],
        "Reference points (primary)": [str(gcz.n_ref.sum()), None],
        "Calibration vertebrae (development)": [str(ds.loc["Calibration", "vertebrae"]), None],
        "Held-out vertebrae (development)": [str(ds.loc["Held-out test", "vertebrae"]), None],
        "Ligament classes": [str(gcz.ligament.nunique()), None]}, [])

    # Table 2: alignment quality
    a, d = nums["atlas_residual_mm"], pd.read_csv(NB / "table_2_alignment_per_case.csv").mean_nnd_mm
    T["tab:alignment"] = ("main", {
        "Primary (LOSO), selected atlas": [f"{folds.test_vertebrae.sum()} lumbar", f"{f(a['mean'], 2)} ({f(a['min'], 2)}--{f(a['max'], 2)})"],
        "Development, after ICP": [f"{len(d)} lumbar", f"{f(d.mean(), 2)} ({f(d.min(), 2)}--{f(d.max(), 2)})"]}, [])

    # Table 3 and Supplementary Table S6: calibration regressor comparison
    fam = ["Ridge", "KNN", "Random forest", "Gradient boosting", "Extra trees"]
    t3 = {v: pd.read_csv(OUT / f"table3_{v}.csv").set_index("ligament") for v in ("as_published", "rows_fixed", "rows_fixed_ctx")}
    T["tab:residual_model_selection"] = ("main", {r: [f(t3["rows_fixed_ctx"].loc[r, c]) for c in fam] for r in LIGS + ["Mean"]},
                                         [(range(5), "min")])
    S6 = {"Mean, published protocol": [f(t3["as_published"].loc["Mean", c]) for c in fam],
          "Mean, validation rows corrected": [f(t3["rows_fixed"].loc["Mean", c]) for c in fam],
          "Mean, rows and fold quantities corrected (Table 3)": [f(t3["rows_fixed_ctx"].loc["Mean", c]) for c in fam]}
    for r in ["LF", "ISL", "ITL", "CL"]:
        S6[r] = [f"{f(t3['rows_fixed_ctx'].loc[r, c])} ({f(t3['as_published'].loc[r, c])})" for c in fam]
    T["tab:regcorr"] = ("supp", S6, [])

    # Table 4: primary evaluation by class
    lig = pd.read_csv(P2 / "table_P2_3_per_ligament.csv").set_index(["method", "ligament"])
    p5 = pd.read_csv(OUT / "table4_primary_per_ligament.csv").set_index(["method", "ligament"])
    meth = ["GCZ", "SLD", "GCZ-RL"]
    t4 = {r: [f(lig.loc[(k, r), "MAE_mm"]) for k in meth] + [f(p5.loc[(k, r), "Pass5_case_avg"], 1) for k in meth] for r in LIGS}
    t4["Macro mean"] = [f(nums[k]["macro_mae"]) for k in meth] + [f(nums[k]["pass5_case_avg"], 1) for k in meth]
    T["tab:loso_ligament"] = ("main", t4, [((0, 1, 2), "min"), ((3, 4, 5), "max")])

    # Table 5: primary evaluation by held-out subject
    sub = pd.read_csv(P2 / "table_P2_4_per_subject.csv").set_index(["method", "held_out_subject"])
    T["tab:loso_subject"] = ("main", {s: [str(folds.loc[s, "training_vertebrae"]), str(folds.loc[s, "test_vertebrae"])]
                                      + [f(sub.loc[(k, s), "MAE_mm"]) for k in meth] for s in folds.index},
                             [((2, 3, 4), "min")])

    # Table 9: development analysis
    dev = pd.read_csv(P1 / "table_4_4_heldout_results.csv").set_index(["method", "ligament"])
    T["tab:literature"] = ("main", {r: [c for k in ["GCZ", "GCZ-RL", "SLD"] for c in
                                        (f(dev.loc[(k, r), "MAE_mm"]), f(dev.loc[(k, r), "Pass5_percent"], 1) + "%")]
                                    for r in LIGS + ["Mean"]},
                           [((0, 2, 4), "min"), ((1, 3, 5), "max")])

    # Supplementary Table S2: point accounting
    order = ["GCZ", "GCZ-RL", "SLD"]

    def acc(g):
        ref = g[g.method == "GCZ"].n_ref.sum()
        out = [str(ref)]
        for k in order:
            h = g[g.method == k]
            out += [str(h.n_pred.sum()), f(100 * h.n_matched.sum() / ref, 1)]
        return out
    S2 = {r: acc(m[m.ligament == r]) for r in LIGS}
    S2.update({s: acc(m[m.subject == s]) for s in sorted(m.subject.unique())})
    S2["All"] = acc(m)
    T["tab:accounting"] = ("supp", S2, [])

    # Supplementary Table S3: matching-independent metrics
    ind = pd.read_csv(OUT / "table_S2_S3_point_accounting_and_matching_independent.csv").set_index(["method", "ligament"])
    cols = [("P_to_R", 3), ("R_to_P", 3), ("symmetric", 3), ("ref_within_5mm_pct", 1)]
    S3 = {r: [f(ind.loc[(k, r), c], d) for c, d in cols for k in order] for r in LIGS}
    keys = [("closest_point_pred_to_ref", 3), ("closest_point_ref_to_pred", 3), ("symmetric_nn", 3), ("ref_within_5mm_pct", 1)]
    S3["Mean over cases"] = [f(nums[k][c], d) for c, d in keys for k in order]
    T["tab:indep"] = ("supp", S3, [((0, 1, 2), "min"), ((3, 4, 5), "min"), ((6, 7, 8), "min"), ((9, 10, 11), "max")])

    # Supplementary Table S4: candidate-region feasibility
    z = pd.read_csv(OUT / "table_S4_candidate_region_feasibility.csv").set_index("ligament")
    size = lambda x: f"{x:g}"  # noqa: E731
    n = len(folds.index) and folds.test_vertebrae.sum()
    S4 = {r: [str(z.loc[r, "ref"]), f(z.loc[r, "inside_pct"], 1), size(z.loc[r, "selection_median"]), size(z.loc[r, "projection_median"]),
              f"{z.loc[r, 'fallback']}/{n}", f"{z.loc[r, 'fallback']}/{n}", str(z.loc[r, "side_empty"])] for r in LIGS}
    S4["All"] = [str(z.ref.sum()), f(nums["zones_inside_pct"], 1), "--", "--", f"{z.fallback.sum()}/{n * len(LIGS)}",
                 f"{z.fallback.sum()}/{n * len(LIGS)}", str(z.side_empty.sum())]
    T["tab:zones"] = ("supp", S4, [])

    # Supplementary Tables S7 and S8: correspondence diagnostic (original notebook)
    cs = pd.read_csv(NB / "table_3_correspondence_strategies__dissertation_table_4_2.csv")
    p1 = pd.read_csv(P1 / "per_vertebra_per_ligament_metrics.csv")
    nine = p1[(p1.method == p1.method.iloc[0]) & p1.case.isin(pd.read_csv(NB / "table_2_alignment_per_case.csv").case)]
    nref = nine.n_ref.sum()
    name = {"Final anatomy-constrained GCZ": "Final anatomy-constrained"}
    T["tab:matching_versions"] = ("supp", {name.get(r.strategy, r.strategy): [str(r.pairs), f(100 * r.pairs / nref, 1), f(r.MAE, 2), f(r.SD, 2),
                                                                             f(r.RMSE, 2), f(r.Pass5_percent, 1), str(r.N_gt_10mm)]
                                           for r in cs.itertuples()}, [])
    pl = pd.read_csv(NB / "table_4_per_ligament_mae.csv").set_index("ligament")
    T["tab:per_ligament"] = ("supp", {r: [f(pl.loc[r, "initial_MAE"], 2), f(pl.loc[r, "final_MAE"], 2), str(pl.loc[r, "N_gt_10mm_final"])]
                                      for r in LIGS}, [])

    # Supplementary Table S9: selected models
    sel = pd.read_csv(P2 / "table_P2_6_selected_models.csv").set_index(["held_out_subject", "ligament"])
    T["tab:selected"] = ("supp", {r: [c for s in folds.index for c in
                                      (f"{sel.loc[(s, r), 'extra_trees']}, {'0' if sel.loc[(s, r), 'alpha'] == 0 else f(sel.loc[(s, r), 'alpha'], 2)}",
                                       f(sel.loc[(s, r), "inner_loso_mae_mm"]))] for r in LIGS}, [])
    return T, {"Supplementary text: reference landmarks of the nine diagnostic vertebrae": (str(nref), "1145")}


def main():
    tables, extra = expected_tables()
    names = {"tab:dataset": "Table 1", "tab:alignment": "Table 2", "tab:residual_model_selection": "Table 3",
             "tab:loso_ligament": "Table 4", "tab:loso_subject": "Table 5", "tab:literature": "Table 9",
             "tab:accounting": "Table S2", "tab:indep": "Table S3", "tab:zones": "Table S4", "tab:regcorr": "Table S6",
             "tab:matching_versions": "Table S7", "tab:per_ligament": "Table S8", "tab:selected": "Table S9"}
    rows = []
    for label, (where, exp, groups) in tables.items():
        got = printed(MAIN if where == "main" else SUPP, label)
        for r, cells in exp.items():
            if r not in got:
                rows.append({"table": names[label], "row": r, "column": "", "printed": "(row missing)", "code": "", "ok": False})
                continue
            pc, pb = got[r]
            for k, c in enumerate(cells):
                if c is not None:
                    rows.append({"table": names[label], "row": r, "column": k + 1, "printed": pc[k] if k < len(pc) else "(missing)",
                                 "code": c, "ok": k < len(pc) and pc[k] == c})
            if groups:
                want, have = best(pc, groups), {k for k, b in enumerate(pb) if b}
                rows.append({"table": names[label], "row": r, "column": "bold", "printed": sorted(have), "code": sorted(want),
                             "ok": want == have})
        for r in set(got) - set(exp):
            if label not in ("tab:dataset", "tab:alignment"):   # descriptive rows (protocol, acceptance criterion)
                rows.append({"table": names[label], "row": r, "column": "", "printed": "(row not checked)", "code": "", "ok": False})
    for what, (code, text) in extra.items():
        rows.append({"table": "text", "row": what, "column": "", "printed": text, "code": code,
                     "ok": code == text and f"carry {text} reference landmarks" in SUPP})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "revision" / "table_cell_check.csv", index=False)
    for r in tab[~tab.ok].itertuples():
        print(f"DIFFERS  {r.table}, row {r.row}, column {r.column}: printed {r.printed}, code {r.code}")
    cells = tab[tab.column != "bold"]
    print(f"{int(cells.ok.sum())}/{len(cells)} table cells and {int(tab[tab.column == 'bold'].ok.sum())}/{(tab.column == 'bold').sum()} "
          f"bold markings of {tab.table.nunique() - ('text' in set(tab.table))} typed tables agree with the code outputs")
    sys.exit(0 if tab.ok.all() else 1)


if __name__ == "__main__":
    main()
