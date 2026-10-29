# Manuscript and submission files

LaTeX source of the *Medical Image Analysis* submission, together with the compiled PDFs and the
submission documents.

| File | Content |
|---|---|
| `main.tex`, `references.bib`, `figs/` | Manuscript (Elsevier `elsarticle`, author-year citations) |
| `main.pdf` | Compiled manuscript |
| `supplementary.tex`, `supplementary_figs/` | Supplementary material |
| `supplementary.pdf` | Compiled supplementary material |
| `highlights.txt`, `graphical_abstract.{pdf,png}` | Required highlights and graphical abstract |
| `submission/cover_letter.{tex,pdf}` | Cover letter |
| `submission/author_vitae.{md,docx}` | Author biographies (placeholders) |
| `README_SUBMISSION.md` | Submission checklist, open decisions, and changes from the previous draft |

**Build.** Run pdfLaTeX, then BibTeX, then pdfLaTeX twice on `main.tex`. Overleaf compiles the
folder as it is: zip this folder and use New Project → Upload Project.

**Before the repository is made public**, decide whether to keep this folder. The cover letter
and `README_SUBMISSION.md` are internal submission documents; the Statement of Significance is
kept with the submission files, not in this repository.
