# GCZ-RL manuscript: Medical Image Analysis submission package

`main.tex` uses Elsevier's `elsarticle` class with author-year citations (`elsarticle-harv`), as
the MedIA Guide for Authors requires. Overleaf compiles it as is: pdfLaTeX → BibTeX →
pdfLaTeX ×2. Line numbers are on for review. For a single-column double-spaced copy, swap the two
`\documentclass` lines at the top.

## Files

| File | Purpose | MedIA requirement |
|---|---|---|
| `main.tex`, `references.bib`, `figs/` | Manuscript (15 pp., abstract 249 words, 7 keywords) | Editable source; abstract ≤ 250 words; 1–7 keywords |
| `supplementary.tex`, `supplementary_figs/` | Supplementary material (S1–S7) | — |
| `highlights.txt` | 5 bullets, each ≤ 85 characters | Required, separate file with "highlights" in its name |
| `graphical_abstract.pdf` / `.png` | Graphical abstract, 2600 px wide | Required, ≥ 1328 px wide, readable at 13 cm |
| `submission/Statement_of_Significance.docx` | Draft answers to all 18 questions | **Required. A missing or incomplete statement means desk rejection.** |
| `submission/cover_letter.tex` | Cover letter draft | — |
| `submission/author_vitae.docx` | Placeholders, ≤ 100 words each | Requested |

## Before submission: authors must decide or confirm

0. **Most important: robustness of the ICP atlas selection** (Supplementary S5; audit finding 14).
   The stored alignments reproduce every number exactly. Recomputing them in another environment
   changed the atlas of 2 of 14 vertebrae and gave GCZ-RL 2.418 mm instead of 2.150 mm (SLD
   2.224 mm). The current draft discloses this. The stronger option is to make atlas selection
   robust (pre-alignment and multi-start ICP), re-run all analyses, and report variability over
   seeds. Decide before submission.
1. **Authorship.** Dr. Lin's agreement and his approval of the manuscript. He remains the
   corresponding author, as in the previous draft; if that changes, move `\corref`/`\ead` together.
2. **Postal addresses and departments.** Check them in the `\affiliation` blocks.
3. **CRediT roles.** Every author should confirm their own roles. Dr. Vrionis's "Clinical
   guidance" is not a CRediT role.
4. **Journal declarations.** Answer question 16 of the Statement of Significance and the
   declarations of the submission system according to the journal's policies.
5. **Statement of Significance.** Fill the items marked [AUTHORS], in particular the 3–5
   Associate Editors (non-US, no conflict of interest).
6. **Code access during review.** Decide whether reviewers get access to the GitHub repository
   (it is private until publication).

## What changed relative to the previous draft

### Format

- **Template.** The manuscript uses the MedIA template, author-year citations, and a BibTeX file.
- **References.** All 27 DOIs were checked against Crossref, and author names now come from the
  Crossref records (several given names were wrong). SLD is cited as published at ISBI 2026,
  with its DOI.
- **Front matter.** New title, chosen by the authors: "GCZ-RL: Geometry-Calibrated Residual
  Learning for Spinal Ligament Landmark Localization on Vertebral Surface Meshes". Full postal
  addresses were added.
- **Required items added.** Highlights within the limits, graphical abstract, and Statement of
  Significance.

### Agreement with the code

Details are in `docs/CODE_PAPER_AUDIT.md` of the code repository.

- **Correction equations.** The equations now match the implementation: the correction is
  projected twice, and the blend weight is applied after the first projection.
- **Implementation details.** The text now names the mirrored classes and the mirroring plane,
  and states the snap radius, the side and quadrant snapping, the farthest-point sampling to ten
  landmarks, and the Hungarian pairing of residual targets.
- **Table 3 (regressor comparison).** It was recomputed after a validation-indexing bug was fixed;
  the old values are in Supplementary Table S6. Extra Trees is now best for all seven classes.
- **Development analysis.** The text now states that the GCZ baseline column and the GCZ input
  of GCZ-RL are different runs, and reports the input's own metrics.
- **Matching diagnostic.** The analysis (old Tables 3–4 and two figures) moved to Supplementary
  S6. A disclosure states that its prediction count was set from the reference count.

### New analyses and figures

- **Statistics.** Paired statistics were added: Wilcoxon tests and a vertebra-level cluster
  bootstrap.
- **Closest-point metrics.** Closest-point metrics in both directions were added, including
  SLD's own protocol, together with a Discussion paragraph on why SLD's errors here are higher
  than its published 0.7 mm.
- **Alignment check.** Alignment quality is now reported for the primary analysis.
- **Fig. 6 (new).** It shows all three methods on a held-out vertebra in three orthogonal views,
  as asked by question 9 of the Statement of Significance.
- **Settings table.** Supplementary Table S5 lists every setting, and Table S9 the models
  selected in each fold.
