"""Verify the research archive without rerunning or changing inversions."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / 'MANIFEST.sha256.json').read_text(encoding='utf-8'))
for name, expected in manifest.items():
    path = ROOT / name
    if not path.is_file():
        raise SystemExit(f'Missing file: {name}')
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit(f'Checksum mismatch: {name}')
print(f'PASS: {len(manifest)} archived file checksums')

paper = ROOT / 'tomografia_Q/paper_q1'
expected = json.loads((paper / 'segunda_revision/diagnosticos.json').read_text(encoding='utf-8'))['bootstrap']['campaigns']
campaigns = {
    'Reference': 'sgc_final_q_extra_station_tomography',
    'Expanded': 'sgc_q_ampliado_300km',
    'RayDomain': 'sgc_q_mascara_rayos',
    'Adopted': 'sgc_q_mascara_rayos_damp400',
}
for label, campaign in campaigns.items():
    data = ROOT / 'tomografia_Q/output_V2' / campaign / 'data'
    frame = pd.read_csv(data / 'sgc_1000_q_predictions.csv')
    terms = pd.read_csv(data / 'sgc_1000_q_station_terms.csv').set_index('station').station_term_s
    observed = frame.delta_tstar_s.to_numpy()
    predicted = frame.pred_final_s.to_numpy()
    station = (frame.station_i.map(terms) - frame.station_j.map(terms)).to_numpy()
    structural = predicted - station
    metrics = {
        'E_RMS': float(100 * (1 - np.sqrt(np.sum((observed - predicted)**2) / np.sum(observed**2)))),
        'final_rms_s': float(np.sqrt(np.mean((observed - predicted)**2))),
        'component_corr': float(np.corrcoef(structural, station)[0, 1]),
        'structure_ratio': float(np.std(structural) / np.std(observed)),
        'station_ratio': float(np.std(station) / np.std(observed)),
    }
    assert len(frame) == expected[label]['pairs'], label
    assert frame.event_id.nunique() == expected[label]['events'], label
    for name, value in metrics.items():
        assert np.isclose(value, expected[label][name], rtol=1e-9, atol=1e-10), (label, name, value)
    variance = np.var(structural, ddof=1) + np.var(station, ddof=1) + 2 * np.cov(structural, station)[0, 1]
    assert np.isclose(np.var(predicted, ddof=1), variance, rtol=1e-10, atol=1e-12), label
    print(f'PASS: {label}: {len(frame)} observations; archived RMS and component metrics reproduced')
print('PASS: saved-output audit complete; this is not a rerun of bounded tomography or an independent scientific validation')
