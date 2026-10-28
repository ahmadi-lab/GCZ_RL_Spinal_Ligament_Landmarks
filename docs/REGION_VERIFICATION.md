# Anatomical verification of the GCZ candidate regions

Script: `paper/region_verification.py`; outputs in `paper/outputs/region_verification/`. The
regions are rebuilt exactly as in the method (`gcz.build_thresholds`, `gcz.candidate_vertices`,
projection regions without the last fallback) from the training subjects of each fold, on the
14 held-out lumbar vertebrae of the primary evaluation. Nothing was shrunk, moved, or redrawn.
Values are medians over the 14 vertebrae (containment: pooled over all reference points).

## Where the broad coloured surfaces come from

1. **Not the visualization.** The renderer colours exactly the vertices of each region mask
   (automatic check: identical vertex sets for all seven classes).
2. **The masks themselves.** The regions are the vertices inside three calibrated bands
   (μ ± 2σ of `y_n`, `|x_n − 0.5|`, and `z_n`, Eq. 4). For a class whose annotations have two
   separate parts, such as superior and inferior rims (`ALL_S`/`ALL_I`, `PLL_S`/`PLL_I`, the four
   CL quadrants), the pooled band covers both parts and the surface between them. Because `z_n` is
   measured over the whole vertebra, the parts also sit at different normalized heights in
   different vertebrae (the inferior ALL margin lies between `z_n` 0.05 and 0.50), which widens the
   band further. In the fold holding out S03, the `z_n` band of ALL is 0.00–1.00, of PLL 0.04–1.00,
   and of CL 0.00–1.00, so these regions extend over the whole height of the vertebra.
3. **This is the calibration as designed and as described in the paper.** The regions are search
   regions that contain the footprints, not footprint segmentations (Sections 3.4 and 4.4).

## Class by class

Footprint: mesh surface within 2.5 mm of an expert point of the class.

| Class | Expert annotations (anatomy) | GCZ candidate region (what the masks contain) | Inside (%) | Region area (%) | Footprint area (%) | Region / footprint | Footprint in region (%) | Region on footprint (%) | Overlap with other regions (%) |
|---|---|---|---|---|---|---|---|---|---|
| ALL | Points along the anterior superior and inferior margins of the vertebral body | Anterior 21 % of the antero-posterior extent over the full height: the anterior body surface and the anterior parts of both endplates | 96.7 | 20.3 | 3.3 | 6.3 | 96.7 | 15.5 | 0.0 |
| PLL | Points along the posterior margins of the vertebral body, canal side | Posterior body slab (`y_n` 0.33–0.46) over almost the whole height, including the posterior body wall | 83.7 | 9.7 | 1.9 | 5.2 | 88.9 | 16.4 | 31.0 |
| LF | Points along the laminar margins bordering the canal | Laminar and pedicle zone (`y_n` 0.42–0.66), overlapping the CL and PLL regions | 87.4 | 17.7 | 3.0 | 6.0 | 94.2 | 15.6 | 84.3 |
| ISL | Points along the superior and inferior borders of the spinous process | Mid-sagittal slab (`|x_n − 0.5|` ≤ 0.05) covering the spinous process and its sides | 94.5 | 8.1 | 2.5 | 3.1 | 93.2 | 28.1 | 20.6 |
| SSL | One point at the tip of the spinous process | Posterior tip of the spinous process (`y_n` 0.98–1.00), shared with the ISL region | 85.7 | 0.4 | 0.1 | 2.5 | 86.2 | 32.5 | 78.0 |
| ITL | One point per side at the transverse process tips | The most lateral vertices only (`|x_n − 0.5|` = 0.50): a few vertices at each tip, smaller than the footprint | 92.9 | 0.1 | 0.3 | 0.3 | 29.9 | 100.0 | 0.0 |
| CL | Points along the margins of the superior and inferior articular processes (facet capsules) | Postero-lateral band over the full height (`z_n` 0.00–1.00) covering the articular processes and parts of the pedicles and laminae | 88.1 | 26.9 | 6.6 | 4.0 | 89.0 | 21.9 | 38.3 |

GCZ-RL landmarks: 100 % are vertices of the projection region for every class except ITL
(92.9 %; the two side-preserving fallback cases of S01 L5 and S02 L4).

## Discrepancies between the regions and the anatomical descriptions

- **ALL, PLL, CL.** The regions cover the surface between the superior and inferior attachment
  margins, which is not attachment anatomy (for ALL, the anterior halves of the endplates and the
  anterior body wall). The cause is the pooling of two annotation parts into one band.
- **LF.** The region overlaps the CL and PLL regions by 84 % of its area. Classes therefore share
  search space, and the class of a landmark is fixed by the transferred atlas landmark, not by
  the region.
- **ITL.** The lateral band collapses to the extreme lateral vertices (`|x_n − 0.5|` = 0.50), so the
  region is narrower than the annotated tip footprint and contains only 30 % of its area. Because
  containment is evaluated at the nearest vertex of each reference point, 92.9 % of the ITL points
  still count as inside.
- **Across all classes.** About 67–85 % of each region's area (ITL excepted) lies more than
  2.5 mm from any expert point. This is the price of 86–97 % footprint coverage with axis-aligned
  bands in normalized coordinates.

## ALL and PLL: superior and inferior body margins

Script: `paper/margin_verification.py`; outputs in `paper/outputs/margin_verification/`
(`margin_check.csv` per vertebra, `margin_zn_by_vertebra.csv`, `margin_summary.json`, and the
figure `fig_margins_ALL_PLL`, Supplementary Fig. S6). The question was whether the ALL and PLL
regions mark the middle of the vertebral body because of the camera, the normalization, the mask
computation, or the drawing code.

| Possible cause | Check | Result |
|---|---|---|
| Camera or axis direction | Superior (`_S`) expert files above inferior (`_I`) files along +z (LPS) in all 14 vertebrae; on S03 L4, the body margins measured on the mesh alone agree with the expert parts | Correct: superior above inferior in 28/28 vertebra–class cases (median 27.5 mm ALL, 27.4 mm PLL); mesh margins within 1.6 mm of the expert parts |
| Normalization | `gcz.aabb_norm` against (z − z_min)/(z_max − z_min) over the whole vertebra | Identical for every expert point |
| Drawing code | Rendered colours against the mask vertices (`region_verification.py`) | Identical vertex sets |
| Mask computation | Share of the region area in the middle half between the two margins of the same wall | **The masks do cover the middle of the body wall**: median 28 % (ALL) and 35 % (PLL) of the region area |

**Cause.** The method calibrates one `z_n` band per class from all its training points
(μ ± 2σ, Eq. 4), and `z_n` is measured over the whole vertebra, including the posterior
elements. The body margins therefore sit at different normalized heights in different vertebrae:
inferior ALL margin 0.05–0.50, superior 0.50–0.97. A band that contains the margins of all
training vertebrae spans 0.00–1.00 (ALL) and 0.04–1.00 (PLL). This is the method as designed and
as described in the paper, not an implementation error, and the figures draw it correctly.

**Effect on the existing results.** None of the reported numbers change; they were computed with
these masks. The superior–inferior position of each landmark comes from the atlas transfer, and
the region constrains only its antero-posterior and lateral position: none of the 140 ALL and 140
PLL landmarks of GCZ, and none of GCZ-RL, lies in the middle half between the margins (expert
points: none either). For these two classes the region is a weaker constraint than its appearance
suggests, but it does not move landmarks off the margins.

**Separate superior and inferior bands (counterfactual, not applied).** In the same normalized
coordinates, separate bands would leave the ALL region unchanged in the folds holding out S02 and
S03 (the two bands together span 0.00–1.00) and would remove about 45 % of the PLL region, but 8
ALL and 18 PLL landmarks of GCZ-RL, all at the margins of their vertebrae, would fall outside the
separate bands, because some held-out margins lie beyond the range of the training vertebrae.
Regions that follow the margins would need a body-relative height coordinate, for example `z`
normalized between the endplates of the vertebral body. That is a change of the method that
would require a full re-run of the nested evaluation; it has not been made.

## Consequences for the manuscript

- The original Fig. 3 (threshold calibration) is unchanged.
- A new figure shows the expert annotations, the candidate regions, and the GCZ-RL landmarks
  separately for every class on held-out S03 L4, with one camera per class. Its caption states
  that the regions are search regions containing the footprints.
- Results 4.4 and Supplementary Section S14 report the class-by-class numbers above.
- The reported containment, area fractions, and predictions are unchanged.
- A tighter region definition, such as separate bands for the superior and inferior parts, would
  be a method change. It is noted as future work, not applied.
