# GCZ-RL manuscript: Medical Image Analysis submission package

`main.tex` uses Elsevier's `elsarticle` class with author-year citations (`elsarticle-harv`), as
the MedIA Guide for Authors requires. Overleaf compiles it as is: pdfLaTeX → BibTeX →
pdfLaTeX ×2. Line numbers are on for review. For a single-column double-spaced copy, swap the two
`\documentclass` lines at the top.

**Clean and highlighted versions.** Both are built from this one source. A switch near the top
of `main.tex` selects the version: `\highlightfalse` gives the clean version, and `\highlighttrue` marks the text,
captions, and tables that are new or revised in this revision in blue. Nothing else differs.

## Files

| File | Purpose | MedIA requirement |
|---|---|---|
| `main.tex`, `references.bib`, `figs/`, `tables/` | Manuscript (21 pp., abstract 249 words, 7 keywords). `tables/` holds the rows of the generated tables (Tables 6–8), written by `paper/revision_tables.py` of the code release | Editable source; abstract ≤ 250 words; 1–7 keywords |
| `supplementary.tex`, `supplementary_figs/`, `tables/` | Supplementary material, Sections S1–S19 (Tables S1–S18, Figs. S1–S8), in the separate supplementary source package | — |
| `highlights.txt` | 5 bullets, each ≤ 85 characters | Required, separate file with "highlights" in its name |
| `graphical_abstract.pdf` / `.png` | Graphical abstract, 2600 px wide | Required, ≥ 1328 px wide, readable at 13 cm |
| `submission/Statement_of_Significance.md` / `.docx` | Draft answers to all 18 questions; the Word file is generated from the Markdown with pandoc, so edit the Markdown and regenerate | **Required. A missing or incomplete statement means desk rejection.** |
| `submission/cover_letter.tex` / `.pdf` | Cover letter draft (rebuild the PDF after editing the `.tex`) | — |
| `submission/author_vitae.md` / `.docx` | Placeholders, ≤ 100 words each (Word file generated from the Markdown) | Requested |

## Before submission: authors must decide or confirm

1. **Authorship.** Dr. Lin's agreement and his approval of the manuscript. He remains the
   corresponding author, as in the previous draft; if that changes, move `\corref`/`\ead` together.
2. **Postal addresses and departments.** Check them in the `\affiliation` blocks.
3. **CRediT roles.** Every author should confirm their own roles. Dr. Vrionis's "Clinical
   guidance" is not a CRediT role.
4. **Generative-AI declaration.** The current `main.tex` has none. Check Elsevier's policy and add
   an accurate declaration if one applies.
5. **Statement of Significance.** Fill the items marked [AUTHORS], in particular the 3–5
   Associate Editors (non-US, no conflict of interest).
6. **Code access during review.** Decide whether reviewers get access to the GitHub repository
   (it is private until publication). The revision is on the branch `revision-lin-tang`; merging it
   into `main` awaits approval.
7. **ICP atlas selection.** The stored alignments reproduce every number exactly. Recomputing all
   alignments with ten different ICP seeds (Supplementary S13) changed 2–4 of the 14 atlases per
   run: GCZ-RL improved on GCZ in every run and for every subject (22.8–33.4 %), while its absolute
   error ranged from 2.10 to 2.92 mm (mean 2.36 mm; reported 2.150 mm). This is stated in Results
   4.4, the Discussion, and the Limitations. A more robust atlas initialization (pre-alignment,
   multi-start ICP) is named as future work.

## Changes in this revision (responses to Dr. Lin and Dr. Tang)

Every result was first reproduced from the public data with the released code
(Supplementary S12, Table S15: 68 claims verified; the S02 reduction is now 10.0 %, previously
computed from rounded means). The method and the reported primary results are unchanged.

- **Annotation provenance** from the source publications (Section 3.1, Limitations).
- **Small-data motivation**: generic landmark detection versus ligament attachment
  localization (Introduction).
- **Units of analysis**: subjects, vertebrae, vertebra–class cases, and points listed
  separately (Table 1); significance statements refer to paired comparisons within the
  evaluated set, with the per-subject consistency reported alongside (Sections 3.8, 4.4).
- **Roles of the components** (Sections 3.2–3.7) and the side-constraint exception, in which 2 of
  742 points keep their GCZ position (Sections 3.7, 3.10).
- **Component analyses** (Section 4.4, Tables 7–8, Fig. 9; Supplementary S8–S11): ablation,
  region coverage and compactness, residual interpretability, and atlas and ICP sensitivity.
- **Matching-independent metrics** in the main text (Table 6), with what each metric measures
  (Section 3.8).
- **Finite element modeling and external validation** in the Discussion and Conclusion.
- **Figures.** All original figures are kept. New: Fig. 6 (by subject and vertebra), Fig. 7
  (GCZ and GCZ-RL with magnified close-ups of the applied correction), Fig. 8 (expert attachment
  footprints, annotations, GCZ-RL landmarks, ALL and PLL margins), Fig. 9 (component analyses),
  Fig. 10 (ten ICP seeds); Supplementary Figs. S4–S5 (method on real data), Section S14 (region
  verification), Section S16 with Fig. S6 (ALL and PLL against the body margins). The previous
  qualitative figure (three methods, seven classes, three views) is unchanged as Supplementary
  Fig. S7; the calibrated search regions are Supplementary Fig. S8 (Section S18). Audit: `docs/FIGURE_AUDIT.md` and `docs/REGION_VERIFICATION.md` of the code release.
- **Highlights**: two reworded (constrained refinement; predictions on the bone surface).

## Changes made earlier, relative to the dissertation draft

- **Template and references.** MedIA template, author-year citations, BibTeX; all DOIs checked
  against Crossref; SLD cited as published at ISBI 2026.
- **Agreement with the code** (`docs/CODE_PAPER_AUDIT.md` of the code repository). The equations
  match the implementation, the implementation details are stated, and Table 3 was recomputed
  after a validation-indexing bug was fixed (old values in Supplementary Table S6). The matching
  diagnostic moved to Supplementary S6 with a disclosure.
- **New analyses.** Paired statistics, closest-point metrics in both directions, alignment
  quality of the primary analysis, Fig. 6 (three methods in three views), and Supplementary
  Tables S5 and S9 (settings and selected models).
