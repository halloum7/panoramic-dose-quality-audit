# Diminishing returns of dose in digital panoramic radiography

Code and anonymised data for the article *Diminishing returns of dose in digital panoramic radiography: a DICOM-based physics audit of noise texture, dose efficiency and paediatric exposure in 157 clinical examinations* (manuscript in preparation).

## What it does

The workflow audits a panoramic unit using only its routine DICOM files: no phantom and no dose meter.

| Step | Content |
|---|---|
| Exposure | Reads kV, mA and exposure time from the DICOM header and builds a relative exposure index, I = (kV/72)² · (mA/4) · (t/14.7 s). |
| Dose | Models the dose–area product (DAP), which the unit does not record, from the exposure factors and detector geometry; propagates its uncertainty by Monte-Carlo (100,000 draws); derives an effective dose. |
| Image quality | Computes 9 no-reference metrics on the dento-alveolar region (Immerkær noise, P95–P5 contrast, SNR, CNR, Tenengrad sharpness, entropy, saturation). |
| Noise physics | Noise power spectrum of the high-pass residual, four anatomical zones, log–log slopes compared with the quantum-noise limit. |
| Statistics | Spearman with bootstrap CIs, age- and sex-adjusted regression, Holm correction, 7 sensitivity analyses, dose efficiency (CNR²/DAP), separate kV and mAs effects. |

Main result: noise falls with exposure at about half the quantum rate (log–log slope −0.29 against −0.5), in every anatomical zone, so higher protocols buy little image quality for their extra dose.

## Repository content

| Path | Content |
|---|---|
| `data/Supplementary_Data_S1_anonymised.csv` | One row per image (157 images, 155 patients): exposure factors, exposure index and group, modelled DAP and effective dose, 9 image-quality metrics. Patients coded P001–P155, images S001–S157; only age, sex and exam year are kept. |
| `code/dcmread.py` | DICOM header reader (no pydicom needed). |
| `code/metrics.py` | Pixel loading and image-quality metrics. |
| `code/robust.py` | Main associations, sensitivity analyses, alternative noise estimators, physics slope. |
| `code/phase2.py` | Dose-efficiency figure of merit, kV and mAs effects, Monte-Carlo DAP. |
| `code/phase3.py`, `code/phase3b.py` | Noise power spectrum, regional analysis, effective dose, local typical values, sample size, log–log slopes. |
| `code/figs.py`, `code/fig3.py`, `code/fig4.py`, `code/fig6.py` | Figures 2 to 6 of the article. |
| `figures/` | Figures 2 to 6 (PNG, 300 dpi). |

## How to run

Python 3 with NumPy, SciPy, Matplotlib and scikit-image (`pip install -r requirements.txt`).

The original DICOM files are not shared (patient images). To run the pipeline on your own unit:

1. Put your DICOM files in a folder and set `DICOM_DIR` to it.
2. `python code/dcmread.py <dicom_dir> img_raw.json` reads the headers (files named `1.dcm`, `2.dcm`, …). From it, build `img.json` (one row per file: file, patient code, exam date, age, sex, kV, mA, time) and `files157.json` (list of files in analysis order).
3. `python code/metrics.py <dicom_dir> files157.json met.json` computes the image-quality metrics.
4. From the folder holding these JSON files, run `robust.py`, `phase2.py`, `phase3.py`, `phase3b.py`, then the figure scripts.

All resampling uses the fixed seed 20261121.

## Citation

[To be completed after publication.]

## License

Code: MIT. Data: CC BY 4.0.
