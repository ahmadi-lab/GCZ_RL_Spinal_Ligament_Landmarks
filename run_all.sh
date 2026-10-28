#!/usr/bin/env bash
# GCZ-RL - reproduce the paper with one command (Linux/macOS).
#   bash run_all.sh            Part 1 + Part 2 + checks + paper numbers
#   FRESH=1 bash run_all.sh    also recompute ICP geometry and SLD
#   SLOW=1  bash run_all.sh    also Table 3 (about 25 min), the original notebook, and the ICP seed analysis (about 100 min)
#   PDF=1   bash run_all.sh    also compile the manuscript PDFs (needs latexmk)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-.venv/bin/python}
export PYTHONIOENCODING=utf-8 MPLBACKEND=Agg PYTHONHASHSEED=0
if [ ! -x "$PY" ]; then python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt; fi
RAW=${GCZRL_RAW_ROOT:-dataset/raw}
[ -d "$RAW/VerSe_Ligaments_GT" ] || "$PY" scripts/fetch_data.py --out "$RAW"
F=""
if [ "${FRESH:-0}" = 1 ]; then F="--fresh"; [ -d code/third_party/SpineLib ] || "$PY" scripts/fetch_spinelib.py; fi
"$PY" part1_dissertation_57/run_part1.py $F
if [ "${SLOW:-0}" = 1 ]; then
  "$PY" part1_dissertation_57/original_notebook/run_original_notebook.py
  "$PY" part1_dissertation_57/original_notebook/extract_tables.py
  "$PY" paper/table3_regressor_comparison.py
  "$PY" paper/icp_seed_sensitivity.py run 0 1 2 3 4 5 6 7 8 9
  "$PY" paper/icp_seed_sensitivity.py summary
fi
"$PY" part2_patient_loso_14/run_part2.py $F
"$PY" part2_patient_loso_14/run_leakage_test.py
"$PY" code/verify_all.py
"$PY" paper/make_paper_numbers.py
"$PY" paper/fig6_qualitative_loso.py
"$PY" paper/check_fig6.py
"$PY" paper/graphical_abstract.py
"$PY" paper/revision_analyses.py
"$PY" paper/figures_main.py
"$PY" paper/region_verification.py
"$PY" paper/margin_verification.py
"$PY" paper/figures_anatomy.py
"$PY" paper/fig11_development.py
"$PY" paper/itl_trace.py
"$PY" paper/icp_seed_sensitivity.py summary
"$PY" paper/revision_tables.py
"$PY" paper/reproducibility_check.py
"$PY" paper/check_text_numbers.py
"$PY" paper/check_tables.py
if [ "${PDF:-0}" = 1 ]; then "$PY" paper/build_manuscript_assets.py --pdf; else "$PY" paper/build_manuscript_assets.py; fi
"$PY" -m unittest discover -s tests
echo DONE
