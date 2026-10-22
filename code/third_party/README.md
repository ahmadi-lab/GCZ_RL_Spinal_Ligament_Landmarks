# Third-party code

`SpineLib/` is the library of the released SLD implementation, part of the SpineToolkit
3D Slicer extension: <https://github.com/VisSim-UniKO/SpineToolkit> (Blomenkamp et al., University
of Koblenz). The upstream repository has no licence file, so SpineLib is **not stored here**.
Fetch the 21 files used by `code/ligament/sld_official.py`, unmodified and verified by their git
blob SHA-1, from the pinned commit:

```bash
python scripts/fetch_spinelib.py
```

SpineLib is needed only to recompute SLD (`run_part2.py --fresh`). Without it, the stored SLD
output in `part2_patient_loso_14/checkpoints/sld_official_medial.pkl` is used.

`slicer_stub/` was written for this project. It replaces the 3D Slicer objects that SpineLib
expects, so SLD can run in plain Python. `slicer_stub/ExtractCenterline.py` replaces the VMTK
centerline with a cross-section medial approximation; this is the one methodological
substitution in the SLD comparison, and the paper documents it.
