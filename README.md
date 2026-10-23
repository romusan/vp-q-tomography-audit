# Vp–Q tomography parameterisation and bound audit

Research archive for **Auditing parameterisation and bound effects in sequential P-wave velocity and spectral-ratio attenuation tomography of the Middle Magdalena Valley, Colombia**.

Authors: Rómulo Sandoval Flórez, Francisco H. Cabrera Zambrano, Claudia Patricia Parra Medina.

The manuscript is being prepared for Journal of South American Earth Sciences. This repository is distinct from the fast-marching relocation/P-wave manuscript and the depth-window relocation audit.

## Verified scope

The archive contains processed campaign tables, saved grids and model arrays, station terms, retained predictions, manuscript figures, and the existing covariance, clipping and conditional event-bootstrap diagnostics. It supports auditing from saved processed inputs, not complete reproduction from raw waveform acquisition. Raw waveforms are not included. Separate checkerboard and Benioff experiments are represented by their archived summary outputs; full regeneration of those experiments has not been verified here.

The displayed clipped Q fields are provisional. Archived fit is evaluated before clipping. Recompute forward predictions for the displayed bounded model before drawing quantitative interpretations. The event bootstrap conditions on fixed predictions and does not refit every resample.

## Layout

- `tomografia_Q/output_V2/`: saved campaign arrays and processed tables.
- `tomografia_Q/paper_q1/`: current manuscript, original audit scripts, archived diagnostics and figures.
- `environment.json` and `requirements.txt`: environment recorded by the original archive builder.
- `MANIFEST.sha256.json`: file checksums for this release.
- `verify_archive.py`: verifies checksums and independently checks the saved linear-fit component statistics.

## Check the saved audit

Run `python verify_archive.py` from the repository root after installing NumPy and pandas. It verifies checksums and the original event counts, observation counts, RMS and component correlations without rerunning inversions or overwriting archived campaign files.

The original scripts are retained verbatim for provenance. They use relative paths under `tomografia_Q/`. Some make first-run backups of the manuscript and a legacy PDF filename; the release supplies compatibility copies. Running `revision_sugerencias.py` and `segunda_revision_diagnosticos.py` regenerates diagnostics in their corresponding folders and may take longer. Work on a copy. Numerical solver and dependency differences may affect exact reproduction. The verified release check does not claim an independent scientific validation or a complete rerun.

Compile `tomografia_Q/paper_q1/manuscrito.tex` with pdfLaTeX, BibTeX, pdfLaTeX, pdfLaTeX. Sources, macros, style and manuscript figures are included.

## Provenance and rights

Earthquake catalogues, phase picks, waveforms and station information originate from the Servicio Geológico Colombiano. Slab2 originates from the USGS and Hayes et al. (2018), cited in the manuscript. Source-specific terms and attribution remain applicable. The repository does not relicense third-party data or Elsevier's class and bibliography style; those retain their original terms. No blanket software/data licence or DOI is asserted by this archive.

No credentials, raw waveforms or private correspondence are included. The publication status of the manuscript must not be inferred from the public availability of this research archive.

## Cite the archived version

Version 1.0.1: https://github.com/romusan/vp-q-tomography-audit/tree/v1.0.1

Use CITATION.cff for the three authors and version. No DOI is claimed. The immutable commit corresponding to the version tag is recorded by GitHub and can be checked with `git rev-parse v1.0.1`.

Version 1.0.1 adds final caption spacing in the submission PDF. Scientific data, diagnostic scripts and numerical results are unchanged from 1.0.0.
