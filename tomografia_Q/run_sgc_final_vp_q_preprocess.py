from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from obspy import UTCDateTime, read, read_inventory
from scipy import stats
from scipy.sparse.linalg import lsmr

from download_sgc_waveforms_for_q_spectral import parse_event_xml, first_domain_p_picks


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output_V2" / "sgc_final_q_extra_station_tomography"
DATA = OUT / "data"
FIG = OUT / "figures"
MODEL = OUT / "models"
CODE = OUT / "codigos_finales"

INVENTORY_CSV = (
    ROOT
    / "output_V2"
    / "data"
    / "sgc_waveforms_q_spectral_final_with_nest_and_extra_stations_inventory.csv"
)
STATIONXML_DIR = ROOT / "station_metadata" / "sgc_response_xml"
CATALOG_CSV = (
    ROOT
    / "output_V2"
    / "data"
    / "sgc_kml_2022_2026_rms_lt1_depth_gt20"
    / "sgc_sismos_area_2022_2026_depth_gt20_rms_lt1_strict.csv"
)
STATIONS_CSV = ROOT / "output_V2" / "data" / "domain_stations_v2_plus_extra_used_for_q.csv"
VP_MODEL_NPZ = ROOT / "modelo_vp_3d_nazca.npz"

LAT0 = 7.0
LON0 = -73.0
LON_MIN, LON_MAX = -74.60, -72.55
LAT_MIN, LAT_MAX = 6.20, 8.05
MAX_MODEL_DEPTH_KM = 180.0
MAX_EVENTS = None

GRID_FACTOR = 4
RAY_STEP_KM = 2.5
VP_MIN_HITS = 10
Q_MIN_HITS = 5
VP_LAPLACE = 3.0
VP_DAMP = 0.25
Q_LAPLACE = 3.0
Q_CELL_DAMP = 0.05
Q_STATION_DAMP = 0.15
LSMR_ITER = 1200

NOISE_START_S = -10.0
NOISE_END_S = -1.0
SIGNAL_START_S = 0.2
SIGNAL_END_S = 4.2
MIN_SIGNAL_LEN_S = 2.0
S_SAFETY_S = 0.5
VP_REF_KM_S = 6.0
VS_REF_KM_S = 3.5
FMIN_HZ = 1.0
FMAX_HZ = 12.0
SNR_FREQ_MIN = 3.0
MIN_FREQ_POINTS = 10
MAX_ABS_DELTA_TSTAR_S = 0.75
REMOVE_INSTRUMENT_RESPONSE = True
RESPONSE_OUTPUT = "VEL"
RESPONSE_PRE_FILT = (0.5, 0.8, 15.0, 20.0)
RESPONSE_WATER_LEVEL = 60.0
QINV_MIN = 1.0 / 300.0
QINV_MAX = 1.0 / 20.0

RESPONSE_INVENTORY_CACHE: dict[tuple[str, str, str, str], object] = {}


def ensure_dirs() -> None:
    for directory in [DATA, FIG, MODEL, CODE]:
        directory.mkdir(parents=True, exist_ok=True)


def latlon_to_xy(lat: float, lon: float) -> tuple[float, float]:
    x = (float(lon) - LON0) * 111.32 * math.cos(math.radians(LAT0))
    y = (float(lat) - LAT0) * 111.32
    return float(x), float(y)


def xy_to_latlon(x: float, y: float) -> tuple[float, float]:
    lat = float(y) / 111.32 + LAT0
    lon = float(x) / (111.32 * math.cos(math.radians(LAT0))) + LON0
    return lat, lon


def depth_km_from_origin(origin) -> float:
    depth = getattr(origin, "depth", None)
    if depth is None:
        return 0.0
    depth = float(depth)
    return depth / 1000.0 if abs(depth) > 1000 else depth


def load_inventory() -> pd.DataFrame:
    inventory = pd.read_csv(INVENTORY_CSV, dtype=str)
    inventory = inventory[inventory["status"].astype(str).str.upper().eq("OK")].copy()
    inventory["origin_time_utc_sort"] = pd.to_datetime(inventory["origin_time_utc"], errors="coerce")
    return inventory.sort_values(["origin_time_utc_sort", "event_id", "station"]).reset_index(drop=True)


def load_event_ids(inventory: pd.DataFrame) -> list[str]:
    by_event = (
        inventory.dropna(subset=["event_id"])
        .groupby("event_id", as_index=False)["origin_time_utc_sort"]
        .min()
        .sort_values(["origin_time_utc_sort", "event_id"])
    )
    event_ids = by_event["event_id"].astype(str).tolist()
    return event_ids if MAX_EVENTS is None else event_ids[:MAX_EVENTS]


def xml_paths_by_event(inventory: pd.DataFrame) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for row in inventory.itertuples(index=False):
        event_id = str(row.event_id)
        path = Path(str(row.xml_file))
        if event_id and path.exists() and event_id not in out:
            out[event_id] = path
    return out


def load_grid_and_vp() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.load(VP_MODEL_NPZ)
    xg = data["xg"][::GRID_FACTOR].astype(float)
    yg = data["yg"][::GRID_FACTOR].astype(float)
    zg = data["zg"][::GRID_FACTOR].astype(float)
    vp0 = data["vp"][::GRID_FACTOR, ::GRID_FACTOR, ::GRID_FACTOR].astype(float)
    return xg, yg, zg, vp0


def idx_to_kji(idx: int, ny: int, nx: int) -> tuple[int, int, int]:
    k = int(idx) // (ny * nx)
    rem = int(idx) - k * ny * nx
    j = rem // nx
    i = rem - j * nx
    return k, j, i


def kji_to_idx(k: int, j: int, i: int, ny: int, nx: int) -> int:
    return int((k * ny + j) * nx + i)


def coord_to_idx(x: float, y: float, z: float, xg, yg, zg) -> int:
    ix = int(np.clip(round((x - xg[0]) / (xg[1] - xg[0])), 0, len(xg) - 1))
    iy = int(np.clip(round((y - yg[0]) / (yg[1] - yg[0])), 0, len(yg) - 1))
    iz = int(np.clip(round((z - zg[0]) / (zg[1] - zg[0])), 0, len(zg) - 1))
    return kji_to_idx(iz, iy, ix, len(yg), len(xg))


def rectangular_domain_mask(xg, yg, zg) -> np.ndarray:
    x_min, y_min = latlon_to_xy(LAT_MIN, LON_MIN)
    x_max, y_max = latlon_to_xy(LAT_MAX, LON_MAX)
    xmin, xmax = sorted([x_min, x_max])
    ymin, ymax = sorted([y_min, y_max])
    zz, yy, xx = np.meshgrid(zg, yg, xg, indexing="ij")
    return (
        (xx >= xmin)
        & (xx <= xmax)
        & (yy >= ymin)
        & (yy <= ymax)
        & (zz <= MAX_MODEL_DEPTH_KM)
    )


def ray_row(src: np.ndarray, rec: np.ndarray, xg, yg, zg, mode="slowness", vp=None) -> dict[int, float]:
    dist = float(np.linalg.norm(rec - src))
    nstep = max(2, int(math.ceil(dist / RAY_STEP_KM)) + 1)
    ds = dist / max(nstep - 1, 1)
    row: dict[int, float] = {}
    for t in np.linspace(0.0, 1.0, nstep):
        x, y, z = src + t * (rec - src)
        if x < xg.min() or x > xg.max() or y < yg.min() or y > yg.max() or z < zg.min() or z > zg.max():
            continue
        idx = coord_to_idx(float(x), float(y), float(z), xg, yg, zg)
        if mode == "slowness":
            w = ds
        else:
            k, j, i = idx_to_kji(idx, len(yg), len(xg))
            w = ds / max(float(vp[k, j, i]), 1.0e-9)
        row[idx] = row.get(idx, 0.0) + w
    return row


def hits_from_rows(rows: list[dict[int, float]], ncell: int) -> np.ndarray:
    hits = np.zeros(ncell, dtype=int)
    for row in rows:
        for idx, val in row.items():
            if abs(val) > 0:
                hits[int(idx)] += 1
    return hits


def rows_to_csr(rows: list[dict[int, float]], active_flat: np.ndarray):
    active = np.flatnonzero(active_flat)
    pos = {int(old): i for i, old in enumerate(active)}
    rr, cc, vv = [], [], []
    kept = []
    for rnum, row in enumerate(rows):
        local = [(pos[int(idx)], float(val)) for idx, val in row.items() if int(idx) in pos]
        if not local:
            continue
        out_r = len(kept)
        kept.append(rnum)
        for col, val in local:
            rr.append(out_r)
            cc.append(col)
            vv.append(val)
    return sp.csr_matrix((vv, (rr, cc)), shape=(len(kept), len(active))), active, kept


def laplacian(mask3: np.ndarray) -> sp.csr_matrix:
    nz, ny, nx = mask3.shape
    active = np.flatnonzero(mask3.ravel())
    pos = {int(old): i for i, old in enumerate(active)}
    rr, cc, vv = [], [], []
    r = 0
    for old in active:
        k, j, i = idx_to_kji(int(old), ny, nx)
        neigh = []
        for dk, dj, di in [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]:
            kk, jj, ii = k + dk, j + dj, i + di
            if 0 <= kk < nz and 0 <= jj < ny and 0 <= ii < nx:
                nb = kji_to_idx(kk, jj, ii, ny, nx)
                if mask3.ravel()[nb]:
                    neigh.append(pos[nb])
        if not neigh:
            continue
        rr.append(r); cc.append(pos[int(old)]); vv.append(float(len(neigh)))
        for nbp in neigh:
            rr.append(r); cc.append(nbp); vv.append(-1.0)
        r += 1
    return sp.csr_matrix((vv, (rr, cc)), shape=(r, len(active)))


def rms(values) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr))) if len(arr) else float("nan")


def build_vp_picks(event_ids: list[str], stations: pd.DataFrame, xml_map: dict[str, Path]) -> pd.DataFrame:
    domain_stations = set(stations["station"].astype(str).str.strip())
    rows = []
    event_rows = []
    for event_id in event_ids:
        xml_path = xml_map.get(event_id)
        if xml_path is None:
            continue
        try:
            event, origin, magnitude = parse_event_xml(xml_path)
        except Exception:
            continue
        depth_km = depth_km_from_origin(origin)
        event_lat = float(origin.latitude)
        event_lon = float(origin.longitude)
        sx, sy = latlon_to_xy(event_lat, event_lon)
        picks = first_domain_p_picks(event, domain_stations)
        event_pick_rows = []
        for pick in picks:
            station = pick["station"]
            srow = stations.loc[stations["station"] == station].iloc[0]
            rx, ry = latlon_to_xy(float(srow["lat"]), float(srow["lon"]))
            travel_time = float(pick["pick_time"] - origin.time)
            if not np.isfinite(travel_time) or travel_time <= 0.2 or travel_time > 300.0:
                continue
            event_pick_rows.append(
                {
                    "event_id": event_id,
                    "origin_time_utc": str(origin.time),
                    "event_lat": event_lat,
                    "event_lon": event_lon,
                    "event_depth_km": depth_km,
                    "source_x_km": sx,
                    "source_y_km": sy,
                    "source_z_km": depth_km,
                    "station": station,
                    "network": pick["network"],
                    "location": pick["location"],
                    "channel": pick["channel"],
                    "station_lat": float(srow["lat"]),
                    "station_lon": float(srow["lon"]),
                    "receiver_x_km": rx,
                    "receiver_y_km": ry,
                    "pick_time_utc": str(pick["pick_time"]),
                    "travel_time_s": travel_time,
                }
            )
        if len(event_pick_rows) >= 3:
            rows.extend(event_pick_rows)
            event_rows.append(
                {
                    "event_id": event_id,
                    "origin_time_utc": str(origin.time),
                    "event_lat": event_lat,
                    "event_lon": event_lon,
                    "event_depth_km": depth_km,
                    "magnitude": getattr(magnitude, "mag", np.nan) if magnitude is not None else np.nan,
                    "picks_used": len(event_pick_rows),
                }
            )
    picks_df = pd.DataFrame(rows)
    events_df = pd.DataFrame(event_rows)
    picks_df.to_csv(DATA / "sgc_1000_vp_picks.csv", index=False)
    events_df.to_csv(DATA / "sgc_1000_events_used_for_vp.csv", index=False)
    return picks_df


def trim_trace(trace, start: UTCDateTime, end: UTCDateTime) -> np.ndarray | None:
    tr = trace.copy().trim(start, end, pad=False)
    data = np.asarray(tr.data, dtype=float)
    if data.size < 32:
        return None
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    return data - float(np.mean(data))


def stationxml_path_for_trace(trace) -> Path:
    loc = str(trace.stats.location) if str(trace.stats.location) else "blank"
    name = f"{trace.stats.network}.{trace.stats.station}.{loc}.{trace.stats.channel}.response.xml"
    return STATIONXML_DIR / name


def load_response_inventory_for_trace(trace):
    key = (
        str(trace.stats.network),
        str(trace.stats.station),
        str(trace.stats.location),
        str(trace.stats.channel),
    )
    if key not in RESPONSE_INVENTORY_CACHE:
        path = stationxml_path_for_trace(trace)
        if not path.exists():
            raise FileNotFoundError(f"Missing StationXML response metadata: {path}")
        RESPONSE_INVENTORY_CACHE[key] = read_inventory(str(path))
    return RESPONSE_INVENTORY_CACHE[key]


def remove_instrument_response(trace):
    if not REMOVE_INSTRUMENT_RESPONSE:
        return trace, False
    inv = load_response_inventory_for_trace(trace)
    corrected = trace.copy()
    corrected.remove_response(
        inventory=inv,
        output=RESPONSE_OUTPUT,
        pre_filt=RESPONSE_PRE_FILT,
        water_level=RESPONSE_WATER_LEVEL,
        zero_mean=True,
        taper=True,
    )
    corrected.detrend("demean")
    return corrected, True


def next_pow2(n: int) -> int:
    return 1 << (max(int(n), 2) - 1).bit_length()


def amplitude_spectrum(y: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    nfft = max(128, next_pow2(len(y)))
    taper = np.hanning(len(y))
    amp = np.abs(np.fft.rfft(y * taper, n=nfft)) + 1.0e-30
    freq = np.fft.rfftfreq(nfft, d=1.0 / float(fs))
    return freq, amp


def hypo_dist_km(event_lat, event_lon, event_depth, sta_lat, sta_lon) -> float:
    ex, ey = latlon_to_xy(event_lat, event_lon)
    sx, sy = latlon_to_xy(sta_lat, sta_lon)
    return float(math.sqrt((ex - sx) ** 2 + (ey - sy) ** 2 + float(event_depth) ** 2))


def build_trace_spectrum(row, stations: pd.DataFrame) -> dict | None:
    try:
        station = str(row.station)
        srow = stations.loc[stations["station"] == station].iloc[0]
        stream = read(str(row.mseed_file))
        stream.merge(method=1, fill_value="interpolate")
        vertical = [tr for tr in stream if str(tr.stats.channel).upper().endswith("Z")]
        if not vertical:
            return None
        trace = sorted(vertical, key=lambda tr: (-float(tr.stats.sampling_rate), str(tr.id)))[0]
        trace.detrend("demean")
        trace, response_removed = remove_instrument_response(trace)
        fs = float(trace.stats.sampling_rate)
        pick = UTCDateTime(str(row.pick_time_utc))
        dist = hypo_dist_km(float(row.event_lat), float(row.event_lon), float(row.event_depth_km), float(srow.lat), float(srow.lon))
        sp_est = dist * (1.0 / VS_REF_KM_S - 1.0 / VP_REF_KM_S)
        sig_end = min(SIGNAL_END_S, sp_est - S_SAFETY_S)
        if sig_end - SIGNAL_START_S < MIN_SIGNAL_LEN_S:
            return None
        noise = trim_trace(trace, pick + NOISE_START_S, pick + NOISE_END_S)
        signal = trim_trace(trace, pick + SIGNAL_START_S, pick + sig_end)
        if noise is None or signal is None:
            return None
        f_sig, a_sig = amplitude_spectrum(signal, fs)
        f_noi, a_noi = amplitude_spectrum(noise, fs)
        noise_interp = np.interp(f_sig, f_noi, a_noi, left=np.nan, right=np.nan)
        snr = a_sig / (noise_interp + 1.0e-30)
        valid = np.isfinite(f_sig) & np.isfinite(a_sig) & np.isfinite(snr) & (a_sig > 0)
        if int(valid.sum()) < MIN_FREQ_POINTS:
            return None
        return {
            "event_id": str(row.event_id),
            "station": station,
            "trace_id": str(trace.id),
            "channel": str(trace.stats.channel),
            "path": str(row.mseed_file),
            "pick_time": str(pick),
            "event_lat": float(row.event_lat),
            "event_lon": float(row.event_lon),
            "event_depth_km": float(row.event_depth_km),
            "station_lat": float(srow.lat),
            "station_lon": float(srow.lon),
            "hypo_dist_km": dist,
            "sp_est_s": sp_est,
            "signal_len_s": float(sig_end - SIGNAL_START_S),
            "fs_hz": fs,
            "instrument_response_removed": bool(response_removed),
            "response_output_unit": RESPONSE_OUTPUT if response_removed else "COUNTS",
            "freq_hz": f_sig[valid],
            "log_signal_amp": np.log(a_sig[valid]),
            "snr": snr[valid],
        }
    except Exception:
        return None


def fit_pair(a: dict, b: dict) -> dict | None:
    fmax = min(FMAX_HZ, 0.45 * min(a["fs_hz"], b["fs_hz"]))
    if fmax <= FMIN_HZ:
        return None
    freq = np.linspace(FMIN_HZ, fmax, 200)
    la = np.interp(freq, a["freq_hz"], a["log_signal_amp"], left=np.nan, right=np.nan)
    lb = np.interp(freq, b["freq_hz"], b["log_signal_amp"], left=np.nan, right=np.nan)
    sa = np.interp(freq, a["freq_hz"], a["snr"], left=np.nan, right=np.nan)
    sb = np.interp(freq, b["freq_hz"], b["snr"], left=np.nan, right=np.nan)
    use = np.isfinite(la) & np.isfinite(lb) & np.isfinite(sa) & np.isfinite(sb) & (sa >= SNR_FREQ_MIN) & (sb >= SNR_FREQ_MIN)
    if int(use.sum()) < MIN_FREQ_POINTS:
        return None
    x = freq[use]
    y = la[use] - lb[use]
    lr = stats.linregress(x, y)
    if not np.isfinite(lr.slope):
        return None
    dt = float(-lr.slope / math.pi)
    return {
        "event_id": a["event_id"],
        "station_i": a["station"],
        "station_j": b["station"],
        "channel_i": a["channel"],
        "channel_j": b["channel"],
        "delta_tstar_s": dt,
        "intercept": float(lr.intercept),
        "slope": float(lr.slope),
        "rvalue": float(lr.rvalue),
        "r2": float(lr.rvalue ** 2),
        "pvalue": float(lr.pvalue),
        "stderr": float(lr.stderr) if lr.stderr is not None else np.nan,
        "n_freq": int(use.sum()),
        "fmin_fit_hz": float(x.min()),
        "fmax_fit_hz": float(x.max()),
        "median_snr_i": float(np.nanmedian(sa[use])),
        "median_snr_j": float(np.nanmedian(sb[use])),
        "hypo_dist_i_km": a["hypo_dist_km"],
        "hypo_dist_j_km": b["hypo_dist_km"],
        "delta_hypo_dist_km": float(a["hypo_dist_km"] - b["hypo_dist_km"]),
        "event_lat": a["event_lat"],
        "event_lon": a["event_lon"],
        "event_depth_km": a["event_depth_km"],
        "station_i_lat": a["station_lat"],
        "station_i_lon": a["station_lon"],
        "station_j_lat": b["station_lat"],
        "station_j_lon": b["station_lon"],
        "path_i": a["path"],
        "path_j": b["path"],
        "pick_i_utc": a["pick_time"],
        "pick_j_utc": b["pick_time"],
        "qc_ok": bool(np.isfinite(dt) and abs(dt) <= MAX_ABS_DELTA_TSTAR_S),
    }


def build_q_pairs(event_ids: list[str], stations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    inv = pd.read_csv(INVENTORY_CSV, dtype=str)
    inv = inv[(inv["status"] == "OK") & (inv["event_id"].isin(event_ids))].copy()
    for col in ["event_lat", "event_lon", "event_depth_km"]:
        inv[col] = pd.to_numeric(inv[col], errors="coerce")
    inv = inv.dropna(subset=["event_lat", "event_lon", "event_depth_km", "mseed_file", "pick_time_utc"])
    specs = []
    for row in inv.itertuples(index=False):
        spec = build_trace_spectrum(row, stations)
        if spec is not None:
            specs.append(spec)
    trace_rows = []
    for spec in specs:
        trace_rows.append({k: v for k, v in spec.items() if k not in {"freq_hz", "log_signal_amp", "snr"}})
        trace_rows[-1]["median_snr_band"] = float(np.nanmedian(spec["snr"][(spec["freq_hz"] >= FMIN_HZ) & (spec["freq_hz"] <= FMAX_HZ)]))
    pairs = []
    by_event: dict[str, list[dict]] = {}
    for spec in specs:
        by_event.setdefault(spec["event_id"], []).append(spec)
    for event_id, items in by_event.items():
        items = sorted(items, key=lambda item: item["station"])
        if len(items) < 2:
            continue
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                pair = fit_pair(items[i], items[j])
                if pair is not None:
                    pairs.append(pair)
    trace_df = pd.DataFrame(trace_rows)
    pair_df = pd.DataFrame(pairs)
    trace_df.to_csv(DATA / "sgc_1000_trace_spectra_inventory.csv", index=False)
    pair_df.to_csv(DATA / "sgc_1000_spectral_tstar_pairs.csv", index=False)
    return trace_df, pair_df


def run_vp_inversion(picks: pd.DataFrame, xg, yg, zg, vp0, domain_mask) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, dict]:
    rows, obs, pred0 = [], [], []
    s0 = (1.0 / np.maximum(vp0, 1.0e-9)).ravel()
    for r in picks.itertuples(index=False):
        src = np.array([r.source_x_km, r.source_y_km, r.source_z_km], dtype=float)
        rec = np.array([r.receiver_x_km, r.receiver_y_km, 0.0], dtype=float)
        row = ray_row(src, rec, xg, yg, zg, "slowness")
        if row:
            rows.append(row)
            obs.append(float(r.travel_time_s))
            pred0.append(sum(s0[idx] * val for idx, val in row.items()))
    hits = hits_from_rows(rows, vp0.size).reshape(vp0.shape)
    mask = domain_mask & (hits >= VP_MIN_HITS)
    A, active, kept = rows_to_csr(rows, mask.ravel())
    obs = np.asarray(obs)[kept]
    pred0 = np.asarray(pred0)[kept]
    b = obs - pred0
    L = laplacian(mask)
    Aaug = sp.vstack([A, VP_LAPLACE * L, VP_DAMP * sp.eye(A.shape[1], format="csr")], format="csr")
    baug = np.concatenate([b, np.zeros(L.shape[0]), np.zeros(A.shape[1])])
    sol = lsmr(Aaug, baug, maxiter=LSMR_ITER)
    ds = sol[0]
    snew = s0.copy()
    snew[active] = np.clip(snew[active] + ds, 1.0 / 15.0, 1.0 / 2.0)
    vp = (1.0 / np.maximum(snew, 1.0e-12)).reshape(vp0.shape)
    pred1 = pred0 + A @ ds
    pred_df = picks.iloc[kept].copy().reset_index(drop=True)
    pred_df["pred_initial_s"] = pred0
    pred_df["pred_final_s"] = pred1
    pred_df["residual_initial_s"] = obs - pred0
    pred_df["residual_final_s"] = obs - pred1
    metrics = {
        "observations": int(len(obs)),
        "input_picks": int(len(picks)),
        "active_cells": int(A.shape[1]),
        "data_parameter_ratio": float(len(obs) / max(A.shape[1], 1)),
        "initial_rms_s": rms(obs - pred0),
        "final_rms_s": rms(obs - pred1),
        "rms_reduction_percent": float(100 * (1 - rms(obs - pred1) / max(rms(obs - pred0), 1.0e-12))),
        "min_hit_threshold": VP_MIN_HITS,
        "lsmr_iterations": int(sol[2]),
        "lsmr_istop": int(sol[1]),
    }
    return vp, mask, pred_df, metrics


def run_q_inversion(pair_df: pd.DataFrame, xg, yg, zg, vp, domain_mask) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame, dict]:
    q0 = np.zeros_like(vp)
    for k, z in enumerate(zg):
        q0[k, :, :] = 40.0 if z < 30.0 else 60.0
    qinv0 = 1.0 / np.maximum(q0, 1.0e-9)
    q0v = qinv0.ravel()
    rows, obs, pred0 = [], [], []
    pairs = pair_df[pair_df["qc_ok"] == True].copy() if "qc_ok" in pair_df.columns else pair_df.copy()
    for r in pairs.itertuples(index=False):
        src = np.array(latlon_to_xy(r.event_lat, r.event_lon) + (float(r.event_depth_km),), dtype=float)
        ri = np.array(latlon_to_xy(r.station_i_lat, r.station_i_lon) + (0.0,), dtype=float)
        rj = np.array(latlon_to_xy(r.station_j_lat, r.station_j_lon) + (0.0,), dtype=float)
        row_i = ray_row(src, ri, xg, yg, zg, "qinv", vp=vp)
        row_j = ray_row(src, rj, xg, yg, zg, "qinv", vp=vp)
        row: dict[int, float] = {}
        for idx, val in row_i.items():
            row[idx] = row.get(idx, 0.0) + val
        for idx, val in row_j.items():
            row[idx] = row.get(idx, 0.0) - val
        if row:
            rows.append(row)
            obs.append(float(r.delta_tstar_s))
            pred0.append(sum(q0v[idx] * val for idx, val in row.items()))
    hits = hits_from_rows(rows, vp.size).reshape(vp.shape)
    mask = domain_mask & (hits >= Q_MIN_HITS)
    A_cells, active, kept = rows_to_csr(rows, mask.ravel())
    obs = np.asarray(obs)[kept]
    pred0 = np.asarray(pred0)[kept]
    pair_kept = pairs.iloc[kept].copy().reset_index(drop=True)
    stations = sorted(set(pair_kept.station_i).union(set(pair_kept.station_j)))
    spos = {s: i for i, s in enumerate(stations)}
    rr, cc, vv = [], [], []
    for rnum, r in enumerate(pair_kept.itertuples(index=False)):
        rr.extend([rnum, rnum])
        cc.extend([A_cells.shape[1] + spos[str(r.station_i)], A_cells.shape[1] + spos[str(r.station_j)]])
        vv.extend([1.0, -1.0])
    A_sta = sp.csr_matrix((vv, (rr, cc)), shape=(len(obs), A_cells.shape[1] + len(stations)))
    A = sp.hstack([A_cells, sp.csr_matrix((A_cells.shape[0], len(stations)))], format="csr") + A_sta
    b = obs - pred0
    L = laplacian(mask)
    blocks = [
        A,
        sp.hstack([Q_LAPLACE * L, sp.csr_matrix((L.shape[0], len(stations)))], format="csr"),
        sp.hstack([Q_CELL_DAMP * sp.eye(A_cells.shape[1], format="csr"), sp.csr_matrix((A_cells.shape[1], len(stations)))], format="csr"),
        sp.hstack([sp.csr_matrix((len(stations), A_cells.shape[1])), Q_STATION_DAMP * sp.eye(len(stations), format="csr")], format="csr"),
    ]
    rhs = [b, np.zeros(L.shape[0]), np.zeros(A_cells.shape[1]), np.zeros(len(stations))]
    sol = lsmr(sp.vstack(blocks, format="csr"), np.concatenate(rhs), maxiter=LSMR_ITER)
    x = sol[0]
    dq = x[: A_cells.shape[1]]
    station_terms = x[A_cells.shape[1] :]
    qinv = qinv0.ravel().copy()
    qinv[active] = np.clip(qinv[active] + dq, QINV_MIN, QINV_MAX)
    pred1 = pred0 + A @ x
    pair_kept["pred_initial_s"] = pred0
    pair_kept["pred_final_s"] = pred1
    pair_kept["residual_initial_s"] = obs - pred0
    pair_kept["residual_final_s"] = obs - pred1
    station_df = pd.DataFrame({"station": stations, "station_term_s": station_terms})
    metrics = {
        "pair_observations_qc": int(len(pairs)),
        "observations_used": int(len(obs)),
        "active_cells": int(A_cells.shape[1]),
        "station_terms": int(len(stations)),
        "data_parameter_ratio_cells_plus_stations": float(len(obs) / max(A.shape[1], 1)),
        "initial_rms_s": rms(obs - pred0),
        "final_rms_s": rms(obs - pred1),
        "rms_reduction_percent": float(100 * (1 - rms(obs - pred1) / max(rms(obs - pred0), 1.0e-12))),
        "median_r2": float(pair_kept["r2"].median()) if len(pair_kept) and "r2" in pair_kept else float("nan"),
        "min_hit_threshold": Q_MIN_HITS,
        "lsmr_iterations": int(sol[2]),
        "lsmr_istop": int(sol[1]),
    }
    return qinv.reshape(vp.shape), q0, mask, pair_kept, station_df, metrics


def plot_area(events_df: pd.DataFrame, stations: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(8, 7), dpi=220)
    ax.scatter(events_df["event_lon"], events_df["event_lat"], s=5, c="black", alpha=0.28, linewidths=0, label=f"Events used (n={len(events_df)})")
    ax.scatter(stations["lon"], stations["lat"], s=70, c="#d7191c", marker="^", edgecolor="white", linewidth=0.8, zorder=4, label=f"Stations (n={len(stations)})")
    ax.plot([LON_MIN, LON_MAX, LON_MAX, LON_MIN, LON_MIN], [LAT_MIN, LAT_MIN, LAT_MAX, LAT_MAX, LAT_MIN], color="#1f78b4", lw=1.8, label="Tomography rectangle")
    ax.set_xlim(LON_MIN - 0.05, LON_MAX + 0.05)
    ax.set_ylim(LAT_MIN - 0.05, LAT_MAX + 0.05)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("SGC-1000 tomography area, stations, and used earthquakes")
    ax.grid(True, color="0.86", lw=0.6)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right", fontsize=8, frameon=True)
    out = FIG / "sgc_1000_area_events_black_stations_red.png"
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_vp(vp0, vp, mask, xg, yg, zg) -> Path:
    dvp = 100 * (vp - vp0) / np.maximum(vp0, 1.0e-9)
    depths = [20, 60, 100, 140]
    vmax = max(5.0, min(35.0, float(np.nanpercentile(np.abs(dvp[mask]), 95)) if np.any(mask) else 20.0))
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=220)
    for ax, depth in zip(axes.ravel(), depths):
        iz = int(np.argmin(np.abs(zg - depth)))
        arr = np.where(mask[iz], dvp[iz], np.nan)
        im = ax.imshow(arr, origin="lower", extent=[xg.min(), xg.max(), yg.min(), yg.max()], cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_title(f"dVp (%) z={zg[iz]:.0f} km")
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label="dVp (%)")
    fig.suptitle("SGC-1000 preliminary P-wave velocity tomography")
    out = FIG / "sgc_1000_vp_slices.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_q(qinv, q0, mask, xg, yg, zg, station_df, pair_df) -> Path:
    dq = 100 * (qinv - 1 / np.maximum(q0, 1.0e-9)) / np.maximum(1 / np.maximum(q0, 1.0e-9), 1.0e-12)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=220)
    if len(pair_df):
        axes[0, 0].scatter(pair_df["residual_initial_s"], pair_df["residual_final_s"], s=12, alpha=0.45, c="black")
        lim = max(abs(pair_df["residual_initial_s"]).max(), abs(pair_df["residual_final_s"]).max())
        axes[0, 0].plot([-lim, lim], [-lim, lim], "r--", lw=1)
    axes[0, 0].set_xlabel("Initial residual (s)")
    axes[0, 0].set_ylabel("Final residual (s)")
    axes[0, 0].set_title("Q residual reduction")
    if len(station_df):
        station_df.sort_values("station_term_s").plot.barh(x="station", y="station_term_s", ax=axes[0, 1], legend=False, color="#4c78a8")
    axes[0, 1].set_title("Station spectral terms")
    axes[0, 1].set_xlabel("Station term (s)")
    vmax = max(25.0, min(100.0, float(np.nanpercentile(np.abs(dq[mask]), 95)) if np.any(mask) else 40.0))
    for ax, depth in zip(axes[1, :], [60, 100]):
        iz = int(np.argmin(np.abs(zg - depth)))
        arr = np.where(mask[iz], dq[iz], np.nan)
        im = ax.imshow(arr, origin="lower", extent=[xg.min(), xg.max(), yg.min(), yg.max()], cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_title(f"pilot dQinv (%) z={zg[iz]:.0f} km")
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
    fig.colorbar(im, ax=axes[1, :].tolist(), shrink=0.85, label="dQinv (%)")
    fig.suptitle("SGC-1000 preliminary spectral-ratio Q inversion")
    out = FIG / "sgc_1000_q_pilot.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_spectral_qc(trace_df: pd.DataFrame, pair_df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=220)
    if len(pair_df):
        qc = pair_df[pair_df["qc_ok"] == True]
        axes[0, 0].hist(pair_df["delta_tstar_s"], bins=45, alpha=0.35, label="all")
        axes[0, 0].hist(qc["delta_tstar_s"], bins=45, alpha=0.75, label="QC")
        axes[0, 0].axvline(0, color="k", lw=1)
        axes[0, 0].legend(fontsize=8)
        axes[0, 1].scatter(qc["delta_hypo_dist_km"], qc["delta_tstar_s"], s=12, alpha=0.45, c="black")
        axes[0, 1].axhline(0, color="0.3", lw=1)
        axes[0, 1].axvline(0, color="0.3", lw=1)
        axes[1, 0].hist(qc["r2"], bins=35, color="#4c78a8", alpha=0.85)
    axes[0, 0].set_title("Differential t* observations")
    axes[0, 0].set_xlabel("Delta t* (s)")
    axes[0, 1].set_title("Distance sanity check")
    axes[0, 1].set_xlabel("Delta hypocentral distance (km)")
    axes[0, 1].set_ylabel("Delta t* (s)")
    axes[1, 0].set_title("Spectral fit R2")
    axes[1, 0].set_xlabel("R2")
    if len(trace_df):
        trace_df["station"].value_counts().sort_values().plot.barh(ax=axes[1, 1], color="#59a14f")
    axes[1, 1].set_title("Accepted spectra by station")
    axes[1, 1].set_xlabel("Trace spectra")
    fig.tight_layout()
    out = FIG / "sgc_1000_spectral_qc.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    ensure_dirs()
    stations = pd.read_csv(STATIONS_CSV)
    inventory = load_inventory()
    event_ids = load_event_ids(inventory)
    xml_map = xml_paths_by_event(inventory)
    if MAX_EVENTS is not None and len(event_ids) < MAX_EVENTS:
        print(f"[warning] only {len(event_ids)} XML events are currently available; continuing with available data.")
    xg, yg, zg, vp0 = load_grid_and_vp()
    domain_mask = rectangular_domain_mask(xg, yg, zg)

    vp_picks = build_vp_picks(event_ids, stations, xml_map)
    trace_df, pair_df = build_q_pairs(event_ids, stations)
    events_df = pd.read_csv(DATA / "sgc_1000_events_used_for_vp.csv")

    vp, vp_mask, vp_pred, vp_metrics = run_vp_inversion(vp_picks, xg, yg, zg, vp0, domain_mask)
    qinv, q0, q_mask, q_pred, station_terms, q_metrics = run_q_inversion(pair_df, xg, yg, zg, vp, domain_mask)

    vp_pred.to_csv(DATA / "sgc_1000_vp_predictions.csv", index=False)
    q_pred.to_csv(DATA / "sgc_1000_q_predictions.csv", index=False)
    station_terms.to_csv(DATA / "sgc_1000_q_station_terms.csv", index=False)
    np.save(MODEL / "sgc_1000_vp_initial.npy", vp0)
    np.save(MODEL / "sgc_1000_vp_final.npy", vp)
    np.save(MODEL / "sgc_1000_vp_active_mask.npy", vp_mask)
    np.save(MODEL / "sgc_1000_qinv_pilot.npy", qinv)
    np.save(MODEL / "sgc_1000_q_pilot.npy", 1.0 / np.maximum(qinv, 1.0e-12))
    np.save(MODEL / "sgc_1000_q_active_mask.npy", q_mask)
    np.savez(MODEL / "sgc_1000_grid.npz", xg=xg, yg=yg, zg=zg)

    area_fig = plot_area(events_df, stations)
    vp_fig = plot_vp(vp0, vp, vp_mask, xg, yg, zg)
    spectral_fig = plot_spectral_qc(trace_df, pair_df)
    q_fig = plot_q(qinv, q0, q_mask, xg, yg, zg, station_terms, q_pred)

    metrics = {
        "run_name": "SGC final Vp and spectral-ratio Q preprocessing with Bucaramanga nest and extra stations",
        "input_catalog": str(CATALOG_CSV),
        "waveform_inventory": str(INVENTORY_CSV),
        "stationxml_response_metadata": str(STATIONXML_DIR),
        "instrument_response_removed_for_q": bool(REMOVE_INSTRUMENT_RESPONSE),
        "instrument_response_output": RESPONSE_OUTPUT if REMOVE_INSTRUMENT_RESPONSE else "COUNTS",
        "instrument_response_pre_filt_hz": list(RESPONSE_PRE_FILT) if REMOVE_INSTRUMENT_RESPONSE else None,
        "requested_max_events": MAX_EVENTS,
        "xml_events_available_for_run": int(len(event_ids)),
        "stations_used": int(stations["station"].nunique()),
        "station_list": sorted(stations["station"].astype(str).unique().tolist()),
        "vp_events_used": int(events_df["event_id"].nunique()) if len(events_df) else 0,
        "vp_picks": int(len(vp_picks)),
        "trace_spectra": int(len(trace_df)),
        "q_pair_observations_total": int(len(pair_df)),
        "q_pair_observations_qc": int(pair_df["qc_ok"].sum()) if len(pair_df) and "qc_ok" in pair_df else 0,
        "domain": {
            "lon_min": LON_MIN,
            "lon_max": LON_MAX,
            "lat_min": LAT_MIN,
            "lat_max": LAT_MAX,
            "max_model_depth_km": MAX_MODEL_DEPTH_KM,
            "rectangular": True,
        },
        "vp": vp_metrics,
        "q_pilot": q_metrics,
        "figures": {
            "area": str(area_fig),
            "vp": str(vp_fig),
            "spectral_qc": str(spectral_fig),
            "q": str(q_fig),
        },
        "limitations": [
            "The inversion is a fast straight-ray diagnostic over a 20-km grid; the publication run should be repeated with the FSM/backtracing operator.",
            "The Q branch is spectral-ratio based and now removes SGC StationXML instrument response before spectral fitting, but it remains pilot-level until stricter spectral QC and independent validation are completed.",
            "Hypocenters are fixed to the SGC catalog; no joint hypocenter-velocity relocation is performed here.",
        ],
    }
    (OUT / "sgc_1000_tomography_metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8")
    shutil.copy2(Path(__file__), CODE / Path(__file__).name)
    print(json.dumps(metrics, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
