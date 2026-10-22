# Primary evaluation with the ICP geometry recomputed in a second environment

These files come from `part2_patient_loso_14/run_part2.py` after the three fold checkpoints
`checkpoints/fold_held_S0*_gcz.pkl` were deleted, so that the GCZ geometry (ICP atlas selection and
alignment) was recomputed from scratch. The run used the same code, the pinned library versions,
and the input files verified by `scripts/fetch_data.py`; the SLD output was read from its
checkpoint.

Changes against the paper:

- The atlas changed for 2 of the 14 held-out vertebrae, S01 L3 and S02 L3.
- The inner model selection changed in the S01 and S02 folds.
- Macro-averaged MAE was 3.208 mm for GCZ and 2.418 mm for GCZ-RL, against 2.879 and 2.150 mm in
  the paper. SLD was unchanged at 2.224 mm.

Supplementary Section S5 discusses the result, and `docs/CODE_PAPER_AUDIT.md` (finding 14)
records it.
