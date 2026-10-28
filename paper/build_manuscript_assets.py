"""Collect every figure and table of the manuscript from the outputs of the code, and check that
every \\includegraphics and every table input of manuscript/main.tex and supplementary.tex resolves.

GENERATED maps each file used by the manuscript to the output of the script that produces it; the
file is copied from that output, so the manuscript always shows what the code computes. STATIC lists
the drawn schematic and illustration figures of the original manuscript, which no script produces;
they are kept as files in manuscript/figs or manuscript/supplementary_figs. Any other figure or
table, or a missing source, is an error.

Run after the analysis scripts (run_all.sh does both).
Usage: python paper/build_manuscript_assets.py [--pdf]
  --pdf  also compiles the clean and highlighted manuscript and the supplement with latexmk
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MS = REPO / "manuscript"
F = "paper/outputs/figures/"
GENERATED = {  # manuscript file -> produced by
    "figs/fig5_patient_level.pdf": (F + "fig5_patient_level.pdf", "paper/figures_main.py fig5"),
    "figs/fig_qualitative_insets.pdf": (F + "fig_qualitative_insets.pdf", "paper/figures_anatomy.py"),
    "figs/fig_regions_margins.pdf": (F + "fig_regions_margins.pdf", "paper/figures_anatomy.py"),
    "figs/fig7_components.pdf": (F + "fig7_components.pdf", "paper/figures_main.py fig7"),
    "figs/fig8_icp_seeds.pdf": (F + "fig8_icp_seeds.pdf", "paper/figures_main.py fig8"),
    "figs/fig11_development_results.pdf": (F + "fig11_development_results.pdf", "paper/fig11_development.py"),
    "supplementary_figs/fig_matching_versions_clean.png": (
        "part1_dissertation_57/original_notebook/outputs/figs/fig_matching_versions_clean.png",
        "part1_dissertation_57/original_notebook/run_original_notebook.py"),
    "supplementary_figs/figS3_residual_vectors.png": (F + "figS3_residual_vectors.png", "paper/figures_main.py figS3"),
    "supplementary_figs/fig1_pipeline.pdf": (F + "fig1_pipeline.pdf", "paper/figures_main.py fig1"),
    "supplementary_figs/fig4_residual_correction.pdf": (F + "fig4_residual_correction.pdf", "paper/figures_main.py fig4"),
    "supplementary_figs/fig_margins_ALL_PLL.pdf": ("paper/outputs/margin_verification/fig_margins_ALL_PLL.pdf",
                                                  "paper/margin_verification.py"),
    "supplementary_figs/fig_qualitative_loso_S03_L4.png": ("paper/outputs/fig6_qualitative_loso_S03_L4.png",
                                                          "paper/fig6_qualitative_loso.py"),
    "supplementary_figs/fig_search_regions.pdf": (F + "fig_search_regions.pdf", "paper/figures_anatomy.py"),
    "tables/tab_region_verification.tex": ("paper/outputs/region_verification/tab_region_verification.tex",
                                           "paper/region_verification.py"),
    "tables/tab_icp_seeds.tex": ("paper/outputs/icp_seed_sensitivity/tab_icp_seeds.tex",
                                 "paper/icp_seed_sensitivity.py summary"),
    "tables/tab_itl_trace.tex": ("paper/outputs/revision/tab_itl_trace.tex", "paper/itl_trace.py"),
    "tables/reproducibility_table.tex": ("paper/outputs/revision/reproducibility_table.tex",
                                         "paper/reproducibility_check.py"),
    "graphical_abstract.png": ("paper/outputs/graphical_abstract.png", "paper/graphical_abstract.py"),
    "graphical_abstract.pdf": ("paper/outputs/graphical_abstract.pdf", "paper/graphical_abstract.py"),
}
# rows of the remaining tables are written into manuscript/tables by paper/revision_tables.py
TABLES_BY_REVISION_TABLES = "paper/revision_tables.py"
STATIC = {  # drawn figures of the original manuscript (no data are plotted in them)
    "fig_gcz_pipeline_benchmark_style4.png": "Fig. 1, workflow schematic",
    "fig_preprocessing_montage.png": "Fig. 2, preprocessing illustration",
    "fig_threshold_candidate_combined1.png": "Fig. 3, calibration schematic",
    "fig_prediction_pipeline.png": "Fig. 4, prediction workflow schematic",
    "gcz_residual_learning_framework.png": "Fig. 5, correction schematic",
    "figure4_cropped_labeled1.png": "Fig. S2, illustration of the matching diagnostic",
}


def main():
    errors = []
    for dest, (src, script) in GENERATED.items():
        s = REPO / src
        if not s.exists():
            errors.append(f"missing {src} (run {script})")
            continue
        (MS / dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(s, MS / dest)
    for tex, dirs in (("main.tex", ["figs"]), ("supplementary.tex", ["supplementary_figs", "figs"])):
        text = (MS / tex).read_text(encoding="utf-8")
        text = "\n".join(line.split("%")[0] if not line.lstrip().startswith("\\%") else line
                         for line in text.splitlines())
        for name in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text):
            hit = next((d for d in dirs if (MS / d / name).exists()), None)
            if hit is None:
                errors.append(f"{tex}: figure {name} not found")
            elif f"{hit}/{name}" not in GENERATED and name not in STATIC:
                errors.append(f"{tex}: figure {name} is neither generated by a script nor a listed drawn figure")
        for name in re.findall(r"\\inputrows\{([^}]+)\}", text):
            if not (MS / name).exists():
                errors.append(f"{tex}: table {name} not found (run {TABLES_BY_REVISION_TABLES})")
    if errors:
        print("\n".join("ERROR " + e for e in errors))
        sys.exit(1)
    print(f"manuscript assets: {len(GENERATED)} files copied from the code outputs; "
          f"every figure and table of main.tex and supplementary.tex resolves")
    if "--pdf" in sys.argv:
        compile_pdfs()


def compile_pdfs():
    build = MS / "_build"
    for name, highlight in (("main", False), ("main_highlighted", True)):
        d = build / name
        if d.exists():
            shutil.rmtree(d)
        # the compiled PDFs of an earlier build stay out of the build folder; figure PDFs are copied
        top_pdfs = {"main.pdf", "main_highlighted.pdf", "supplementary.pdf"}
        shutil.copytree(MS, d, ignore=lambda folder, names: [n for n in names if n == "_build" or
                                                            (Path(folder) == MS and n in top_pdfs)])
        if highlight:
            t = (d / "main.tex").read_text(encoding="utf-8")
            assert "\n\\highlightfalse\n" in t
            (d / "main.tex").write_text(t.replace("\n\\highlightfalse\n", "\n\\highlighttrue\n"), encoding="utf-8")
        latex(d, "main.tex")
        shutil.copy2(d / "main.pdf", MS / f"{name}.pdf")
        if not highlight:
            latex(d, "supplementary.tex")
            shutil.copy2(d / "supplementary.pdf", MS / "supplementary.pdf")
    print("compiled manuscript/main.pdf, main_highlighted.pdf, supplementary.pdf")


def latex(d, tex):
    r = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex], cwd=d,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        sys.exit(f"LaTeX failed for {tex}")
    log = (d / tex.replace(".tex", ".log")).read_text(encoding="utf-8", errors="replace")
    undefined = len(re.findall(r"undefined", log, re.IGNORECASE))
    if undefined:
        sys.exit(f"{tex}: {undefined} undefined references or citations")


if __name__ == "__main__":
    main()
