---
title: "Statement of Significance"
subtitle: "GCZ-RL: Geometry-Calibrated Residual Learning for Spinal Ligament Landmark Localization on Vertebral Surface Meshes"
---

*Draft prepared for the authors. Items marked **[AUTHORS]** must be completed or confirmed by the authors before upload. The journal desk-rejects submissions with a missing or incomplete statement.*

## Question 1 — Peer review history

**[AUTHORS: confirm]** No. This manuscript has not been submitted to this or any other journal. The method and part of the development analysis appear in the first author's doctoral dissertation (Florida Atlantic University, 2026), which is not a peer-reviewed publication. Relative to the dissertation, the manuscript adds a patient-level leave-one-subject-out evaluation with nested model selection and a data-leakage test, a common-frame comparison with the released SLD implementation, paired statistics, matching-independent metrics, and a corrected calibration regressor comparison.

## Question 2 — Prior peer-reviewed version

**[AUTHORS: confirm]** No. No version of this work has been published as a peer-reviewed conference or journal paper.

## Question 3 — Suggested Associate Editors

**[AUTHORS]** Select 3–5 Associate Editors from the current editorial board (https://www.sciencedirect.com/journal/medical-image-analysis/about/editorial-board). Each must be based in a country other than the USA, must have no prior collaboration with any author, and must not be a Senior Editor or member of the Executive Committee. Suitable expertise: landmark detection, shape analysis and statistical shape models, spine image analysis, and point-cloud or mesh registration. Give a one-sentence rationale for each.

## Question 4 — Novelty

GCZ-RL localizes sparse spinal ligament attachment landmarks on vertebral surface meshes by learning annotation-calibrated admissible regions and a bounded residual correction of rigid atlas transfer. Every corrected landmark is projected onto the admissible region of its class, with a side-preserving fallback for bilateral landmarks, and all final predictions remain on the target vertebral surface; the learned component therefore refines a geometrically valid estimate (Introduction, contributions 1 and 2; Sections 3.4–3.10).

## Question 5 — Difference from the authors' prior publications

- Ahmadi, M., Chen, H., Lin, M., Biswas, D., Doulgeris, J., Tang, Y., Engeberg, E.D., Hashemi, J., Pires, G., Vrionis, F.D., 2025. Streamlined and efficient patient-specific modeling for lumbar spine segmentation and finite element analysis. Sci. Rep. 15, 35619. https://doi.org/10.1038/s41598-025-19664-6. That work builds patient-specific lumbar finite element models and places ligaments by fixed geometric rules. The present manuscript learns ligament attachment landmarks from expert annotations and evaluates them against those annotations.
- Ahmadi, M., Zhang, X., Lin, M., Tang, Y., Engeberg, E.D., Hashemi, J., Vrionis, F.D., 2025. Automated finite element modeling of the lumbar spine. World Neurosurg. 201, 124236. https://doi.org/10.1016/j.wneu.2025.124236. That work automates finite element modeling and simulation. Ligament attachment localization is not learned or evaluated there.
- Biswas, D., Zhang, X., Ahmadi, M., Lin, M., Tang, Y., Vrionis, F., Pires, G., Sharp, R., 2025. A unified computational framework for anatomical modeling of six lumbar ligaments. SSRN 5461109. https://doi.org/10.2139/ssrn.5461109. That preprint models ligament geometry from rule-based attachment definitions. The present work replaces fixed rules with calibrated regions and learned bounded correction, validated at the subject level.

## Question 6 — Closest published work

1. SLD: Segmentation-based landmark detection for spinal ligaments. DOI 10.1109/ISBI61048.2026.11515825. First author Lara Blomenkamp, last author Theresa Schöche. IEEE ISBI, 2026. It uses a rule-based detector with no training. We compare against its released implementation on the same vertebrae with the same evaluator.
2. Spinal ligaments detection on vertebrae meshes using registration and 3D edge detection. DOI 10.1109/EMBC53108.2024.10782838. First author Ivanna Kramer, last author Dietrich Paulus. IEEE EMBC, 2024. It combines registration with edge detection. GCZ-RL instead learns admissible regions and a bounded correction from annotations.
3. Automatic detection of cervical spine ligaments origin and insertion points. DOI 10.1109/ISBI.2019.8759223. First author Ibraheem Al-Dhamari, last author Dietrich Paulus. IEEE ISBI, 2019. It uses atlas registration only. GCZ-RL adds calibrated constraints and learned residual correction after atlas transfer.
4. ALPACA: A fast and accurate computer vision approach for automated landmarking of three-dimensional biological structures. DOI 10.1111/2041-210X.13689. First author Arthur Porto, last author A. Murat Maga. Methods in Ecology and Evolution, 2021. It is a generic point-cloud landmark transfer method. GCZ-RL encodes ligament-specific anatomical constraints.

## Question 7 — Comparison with the state of the art

GCZ-RL is better than the geometry-only GCZ baseline for each of the three held-out subjects, with paired tests over the evaluated vertebra–class cases supporting the difference within this cohort. It is not significantly different from SLD, the dedicated state-of-the-art detector.

- **Metrics.** Mean Euclidean error after one-to-one Hungarian matching (MAE) and Pass@5. Point coverage and matching-independent closest-point metrics in both directions are also reported.
- **Data and protocol.** 14 lumbar vertebrae of 3 subjects, evaluated by leave-one-subject-out with nested model selection.
- **Results.** Macro MAE was 2.150 mm for GCZ-RL, 2.879 mm for GCZ, and 2.224 mm for SLD.
  - Against GCZ: difference −0.729 mm, 95% vertebra-level cluster-bootstrap CI −1.078 to −0.381 mm, Wilcoxon signed-rank p < 0.001 over 98 vertebra–class cases. Improvement held for every held-out subject.
  - Against SLD: difference −0.074 mm, CI −0.386 to 0.228 mm, p = 0.35.
- **Closest-point protocol.** Under SLD's own protocol (each prediction to its nearest reference point), SLD had the lower error (1.787 vs 2.140 mm). We therefore describe the two methods as comparable, with complementary class-specific strengths (Section 4.3, Supplementary Table S3).

## Question 8 — Full resolution and field of view

The method operates on the complete vertebral surface mesh at its native resolution (median 33,445 vertices in the primary cohort), with no cropping or patch extraction.

- Candidate regions and projections are computed over all vertices.
- Random point subsets (1200 or 5000 points) are used only to estimate the rigid ICP transformation.
- The input is a surface mesh derived from CT segmentation, not an image. Images are not processed by the method.

## Question 9 — Findings across the full spatial extent

Yes. Fig. 6 shows the predictions of all three methods and the expert landmarks on a complete held-out vertebra in three orthogonal views (axial, sagittal, coronal). Figs. 2–4 show complete vertebral surfaces.

## Question 10 — Novelty information in the manuscript

Yes. Novelty and its quantification are stated in the following places:

- Introduction (contribution list)
- Section 2 (related work)
- Sections 4.3 and 4.6 (quantified comparison)
- Discussion

## Question 11 — Reproducibility

- **Multiple datasets, at least one public.** Yes. Both sources are public:
  - the ligament-landmark annotations released with SLD, which cover VerSe CT subjects and a cervical dataset;
  - the VerSe CT data.
- **Cross-validation and external testing.** The primary evaluation is patient-level leave-one-subject-out on every annotated lumbar vertebra. With only three lumbar subjects, stratification is not meaningful, so each fold holds out one complete subject. No external annotated test set with the same seven ligament classes was available.
  - Primary analysis: S01, S02 and S03 (VerSe).
  - Supplementary vertebra-wise development analysis: also includes the cervical subjects.
- **Data leakage.** All data of a subject are confined to one fold. Calibration, atlas, residual targets, correction bounds, model selection (inner leave-one-subject-out) and blend weights are computed from the training subjects only. A leakage test recomputed the S02 fold with every S02 landmark corrupted by 25 mm Gaussian noise; all 424 predictions were unchanged. Correlated mirrored samples are grouped with their vertebra, and only original samples are used for validation.
- **Reproducibility.** Seeds are fixed. From the stored ICP alignments, every reported number is regenerated exactly from the public data. Recomputing the alignments in a second environment changed the atlas of 2 of 14 vertebrae; this is disclosed with its effect in Supplementary Section S5. All code, trained models, predictions and table-generating scripts are released.

## Question 12 — Datasets and samples

- **VerSe-based ligament landmarks.** Public: VerSe CT and the SLD annotation release. Three subjects, 47 vertebrae (C3–L5), all annotated vertebrae used.
  - Primary evaluation: the 14 lumbar vertebrae, cross-validated by subject; training folds of 9, 10 and 9 vertebrae.
  - Development analysis: 20 calibration vertebrae and 27 scored held-out vertebrae.
- **Cervical ligament landmarks.** Public SLD release. Two subjects, 10 vertebrae (C3–C7). Used only in the development analysis, where they are predicted but not scored, because the original evaluation scored VerSe vertebrae only.
- **Not used.** The artificial Sawbone vertebrae of the SLD release, because they are synthetic models rather than patient data.

## Question 13 — Potential biases

- **Few subjects.** The cohort has only three CT subjects for the lumbar evaluation, so inter-subject variability is under-sampled. Results are reported per subject and with confidence intervals, and claims are limited to this cohort.
- **Single annotation protocol.** Annotations come from a single protocol and group, which may favor methods that learn from them (GCZ-RL) over rule-based ones (SLD). SLD is therefore run with its released implementation and assessed with several metrics, including its own closest-point protocol.
- **Unknown demographics and pathology.** Subject demographics and pathology are unknown.
- **Incomplete scan.** One scan ends at L4.
- **Sparse classes.** SSL and ITL have one point per vertebra or side, so their errors rest on few points.
- **Post hoc design decision.** One design decision (the lateral side constraint) was introduced after inspecting an earlier evaluation. It is disclosed and its effect on ITL is discussed.

## Question 14 — Data and code availability

Yes.

- **Code, trained models, predicted landmarks and per-vertebra metrics.** Released at https://github.com/ahmadi-lab/GCZ_RL_Spinal_Ligament_Landmarks and archived with a DOI on publication. **[AUTHORS: decide whether reviewers get access during review.]**
- **Source meshes and annotations.** Public at their original repositories. They are not redistributed; a script retrieves them from a pinned version and verifies their SHA-256 checksums.

## Question 15 — Reproducibility information in the manuscript

Yes. It is given in the following places:

- Sections 3.1, 3.7, 3.10 and 3.11
- Supplementary Section S5
- The Data and code availability section

## Question 16 — Use of large language models

**[AUTHORS: confirm and edit]** Yes. Claude (Anthropic) was used to:

- reformat the manuscript to the journal template;
- check that reported numbers agree with the released code and saved outputs;
- improve language.

Every number in the manuscript was regenerated from the released code. Every reference was checked against its DOI record (Crossref). The authors reviewed all text. The statement appears in the section "Declaration of generative AI and AI-assisted technologies in the manuscript preparation process".

## Question 17 — Patient-identifying information

Yes, we confirm that the submission includes no patient-identifying information. All data are public, de-identified vertebral surface meshes and landmark coordinates, and subjects are referred to by anonymized codes (S01–S05).

## Question 18 — Additional comments

**[AUTHORS]** Optional. For example, confirm the corresponding author and note the overlap with the dissertation (Question 1).
