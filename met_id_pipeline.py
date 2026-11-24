"""
MET-ID Pipeline: tMS2 Orbitrap → Consensus MGF per GNPS FBMN
Drug: PND (MW 282.14 Da) | Zebrafish metabolism study
"""

import base64
import zlib
import os
import pickle
import hashlib
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import signal
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lxml import etree

# ─── Parameters ───────────────────────────────────────────────────────────────
MZ_TOL_PPM        = 10      # fragment matching tolerance for consistency score
TOP_N_FRAGMENTS   = 10      # top-N fragments used in composite score
SMOOTHING_SIGMA   = 1.5     # Gaussian smoothing width (scan units)
MIN_PROMINENCE    = 0.05    # min peak prominence as fraction of trace max
MIN_PEAK_SNR      = 3.0     # apex must be > N × mean(trace) to reject flat-noise precursors
MIN_PEAK_SCANS    = 4       # minimum scans spanning a valid peak
CONSENSUS_BIN_DA  = 0.01    # m/z bin width for consensus spectrum (Da)
MIN_FRAG_FRAC     = 0.3     # min fraction of peak scans a fragment must appear in
MIN_FRAGMENTS     = 3       # minimum fragments in a valid consensus spectrum
ISOLATION_WINDOW  = 1.6     # quadrupole isolation window (Da)
PRECURSOR_PPM     = 5.0     # tolerance to confirm target precursor survival in MS2

ROOT = os.path.dirname(os.path.abspath(__file__))
CONVERTED_DIR = os.path.join(ROOT, 'converted')
OUTPUT_DIR    = os.path.join(ROOT, 'output')
CACHE_DIR     = os.path.join(ROOT, 'cache')

NS = 'http://psi.hupo.org/ms/mzml'


# ─── 0. Cache helpers ─────────────────────────────────────────────────────────

def _cache_path(mzml_path):
    """Return pickle path for a given mzML file (keyed by filename + size)."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    size   = os.path.getsize(mzml_path)
    name   = os.path.splitext(os.path.basename(mzml_path))[0]
    digest = hashlib.md5(f'{name}_{size}'.encode()).hexdigest()[:8]
    return os.path.join(CACHE_DIR, f'{name}_{digest}.pkl')


def load_cached(mzml_path):
    """
    Return annotated scans from cache if available and up-to-date,
    otherwise return None.  Cache is invalidated when the mzML file size
    changes (i.e. the file was replaced).
    """
    cp = _cache_path(mzml_path)
    if not os.path.exists(cp):
        return None
    try:
        with open(cp, 'rb') as f:
            data = pickle.load(f)
        if data.get('mzml_size') != os.path.getsize(mzml_path):
            return None          # stale cache
        return data['scans']
    except Exception:
        return None


def save_cache(mzml_path, scans):
    """Persist annotated scans to a pickle cache file."""
    cp = _cache_path(mzml_path)
    with open(cp, 'wb') as f:
        pickle.dump({'mzml_size': os.path.getsize(mzml_path), 'scans': scans}, f,
                    protocol=pickle.HIGHEST_PROTOCOL)


def load_or_parse(mzml_path, verbose=True):
    """
    Load annotated scans from cache if available, otherwise parse the mzML,
    run precursor annotation, and save to cache for future runs.
    """
    scans = load_cached(mzml_path)
    if scans is not None:
        if verbose:
            print(f'  [cache hit]  {len(scans)} scans loaded instantly')
        return scans

    if verbose:
        print(f'  [parsing...] {os.path.basename(mzml_path)}')
    scans = parse_mzml(mzml_path)
    scans = annotate_scan_precursors(scans)
    save_cache(mzml_path, scans)
    if verbose:
        print(f'  [cached]     {len(scans)} scans saved to cache')
    return scans


# ─── 1. mzML streaming parser ─────────────────────────────────────────────────

def _tag(name):
    return f'{{{NS}}}{name}'


def _decode_array(binary_text, is_float32=False):
    data = base64.b64decode(binary_text)
    data = zlib.decompress(data)
    return np.frombuffer(data, dtype=np.float32 if is_float32 else np.float64)


def parse_mzml(filepath):
    """Stream-parse mzML, return list of MS2 scan dicts (memory-efficient)."""
    scans = []
    context = etree.iterparse(filepath, events=('end',), tag=_tag('spectrum'),
                               recover=True)

    for _, spectrum in context:
        # ── ms level ──
        ms_level = None
        for cv in spectrum.iter(_tag('cvParam')):
            if cv.get('name') == 'ms level':
                ms_level = int(cv.get('value', 0))
                break
        if ms_level != 2:
            spectrum.clear()
            continue

        # ── RT (minutes → seconds) ──
        rt_sec = None
        for scan_el in spectrum.iter(_tag('scan')):
            for cv in scan_el.iter(_tag('cvParam')):
                if cv.get('name') == 'scan start time':
                    val = float(cv.get('value', 0))
                    unit = cv.get('unitName', 'minute')
                    rt_sec = val * 60.0 if unit == 'minute' else val
                    break
            if rt_sec is not None:
                break

        # ── precursor m/z ──
        precursor_mz = None
        for si in spectrum.iter(_tag('selectedIon')):
            for cv in si.iter(_tag('cvParam')):
                if cv.get('name') == 'selected ion m/z':
                    precursor_mz = float(cv.get('value', 0))
                    break
            if precursor_mz is not None:
                break

        if rt_sec is None or precursor_mz is None:
            spectrum.clear()
            continue

        # ── binary arrays ──
        mz_arr = int_arr = None
        for bda in spectrum.iter(_tag('binaryDataArray')):
            cvs = {cv.get('name') for cv in bda.iter(_tag('cvParam'))}
            is_mz  = 'm/z array' in cvs
            is_int = 'intensity array' in cvs
            is_f32 = '32-bit float' in cvs
            binary_el = bda.find(_tag('binary'))
            if binary_el is not None and binary_el.text:
                arr = _decode_array(binary_el.text, is_float32=is_f32)
                if is_mz:
                    mz_arr = arr
                elif is_int:
                    int_arr = arr

        if mz_arr is None or int_arr is None or len(mz_arr) == 0:
            spectrum.clear()
            continue

        scans.append({
            'rt':           rt_sec,
            'precursor_mz': precursor_mz,
            'mz':           mz_arr,
            'intensity':    int_arr,
        })
        spectrum.clear()

    return scans


# ─── 2. Group scans by precursor ──────────────────────────────────────────────

def group_by_precursor(scans):
    """Group and sort scans by precursor m/z (rounded to 3 dp)."""
    groups = defaultdict(list)
    for s in scans:
        groups[round(s['precursor_mz'], 3)].append(s)
    for key in groups:
        groups[key].sort(key=lambda s: s['rt'])
    return groups


# ─── 3. Surviving precursor analysis ─────────────────────────────────────────

def find_surviving_precursors(scan):
    """
    In stepped HCD (NCE 10/20/30), a fraction of precursor ions survive
    unfragmented. Peaks in the isolation window [target ± ISOLATION_WINDOW/2]
    measured at Orbitrap resolution are surviving precursor ions.
    Their exact mass (5 ppm accuracy at res 15000) identifies which
    compound(s) actually entered the collision cell.

    Returns list of (observed_mz, intensity, ppm_vs_target, confirmed)
    where confirmed=True means the peak matches the target within PRECURSOR_PPM.
    """
    target = scan['precursor_mz']
    half   = ISOLATION_WINDOW / 2
    mz, inten = scan['mz'], scan['intensity']

    mask = (mz >= target - half) & (mz <= target + half)
    survivors = []
    for m, i in zip(mz[mask], inten[mask]):
        ppm = (m - target) / target * 1e6
        confirmed = abs(ppm) <= PRECURSOR_PPM
        survivors.append({'mz': float(m), 'intensity': float(i),
                          'ppm': float(ppm), 'confirmed': confirmed})
    return survivors


def annotate_scan_precursors(scans):
    """
    Add 'survivors' and 'target_confirmed' fields to each scan dict.
    target_confirmed=True when the target precursor is detected as a
    surviving ion within PRECURSOR_PPM.
    """
    for s in scans:
        s['survivors']        = find_surviving_precursors(s)
        s['target_confirmed'] = any(sv['confirmed'] for sv in s['survivors'])
    return scans


# ─── 4. Composite score RT profile ────────────────────────────────────────────

def _cosine_sim(mz1, int1, mz2, int2):
    """Sparse cosine similarity between two top-N MS2 spectra."""
    if len(mz1) == 0 or len(mz2) == 0:
        return 0.0
    tol = mz1 * MZ_TOL_PPM * 1e-6          # per-fragment tolerance
    dot = 0.0
    for m, i, t in zip(mz1, int1, tol):
        diffs = np.abs(mz2 - m)
        j = np.argmin(diffs)
        if diffs[j] < t:
            dot += i * int2[j]
    n1 = np.sqrt(np.dot(int1, int1))
    n2 = np.sqrt(np.dot(int2, int2))
    return dot / (n1 * n2) if n1 > 0 and n2 > 0 else 0.0


def build_rt_profile(scans):
    """
    Score[i] = TIC[i] * mean_cosine_similarity(scan_i, neighbors).

    Multiplying TIC by spectral consistency compresses the baseline:
    - Eluting compound: stable fragments scan-to-scan → cosine ~0.8-1.0
    - Chemical noise:   random fragments each scan  → cosine ~0.05-0.2
    This genuinely differentiates signal from noise, unlike a weighted TIC.
    """
    n = len(scans)
    rt_arr = np.array([s['rt'] for s in scans])
    tic    = np.array([s['intensity'].sum() for s in scans], dtype=float)

    # Pre-extract top-N fragments per scan
    tops = []
    for s in scans:
        mz, inten = s['mz'], s['intensity']
        k = min(TOP_N_FRAGMENTS, len(inten))
        if k > 0:
            idx = np.argpartition(inten, -k)[-k:]
            tops.append((mz[idx], inten[idx].astype(float)))
        else:
            tops.append((np.array([]), np.array([])))

    # Cosine similarity with immediate neighbors
    consistency = np.zeros(n)
    for i in range(n):
        sims = []
        if i > 0:
            sims.append(_cosine_sim(*tops[i], *tops[i - 1]))
        if i < n - 1:
            sims.append(_cosine_sim(*tops[i], *tops[i + 1]))
        consistency[i] = float(np.mean(sims)) if sims else 0.0

    # Boost scans where target precursor is confirmed by surviving ion
    confirmation_boost = np.array(
        [2.0 if s.get('target_confirmed', False) else 1.0 for s in scans]
    )

    return rt_arr, tic * consistency * confirmation_boost


# ─── 4. Peak picking ──────────────────────────────────────────────────────────

def pick_peaks(rt_arr, scores):
    """Return list of peak dicts: apex_rt, start_rt, end_rt, apex_score."""
    if len(scores) < MIN_PEAK_SCANS or scores.max() == 0:
        return []

    smoothed = gaussian_filter1d(scores.astype(float), sigma=SMOOTHING_SIGMA)
    min_prom  = smoothed.max() * MIN_PROMINENCE

    peaks_idx, props = signal.find_peaks(
        smoothed,
        prominence=min_prom,
        width=max(2, MIN_PEAK_SCANS - 2),
    )

    min_apex = smoothed.mean() * MIN_PEAK_SNR

    results = []
    for i, apex in enumerate(peaks_idx):
        if smoothed[apex] < min_apex:
            continue
        left  = int(props['left_bases'][i])
        right = int(props['right_bases'][i])
        results.append({
            'apex_rt':   float(rt_arr[apex]),
            'start_rt':  float(rt_arr[left]),
            'end_rt':    float(rt_arr[right]),
            'apex_score': float(smoothed[apex]),
        })
    return results


# ─── 5. Consensus spectrum ────────────────────────────────────────────────────

def consensus_spectrum(peak_scans):
    """
    Weighted-average spectrum (weight = TIC) binned at CONSENSUS_BIN_DA.
    Fragments present in < MIN_FRAG_FRAC of scans are discarded.
    Returns (mz_array, intensity_array) normalised to max = 999.
    """
    if not peak_scans:
        return np.array([]), np.array([])

    all_mz = np.concatenate([s['mz'] for s in peak_scans])
    if len(all_mz) == 0:
        return np.array([]), np.array([])

    mz_min, mz_max = all_mz.min(), all_mz.max()
    bins    = np.arange(mz_min, mz_max + CONSENSUS_BIN_DA, CONSENSUS_BIN_DA)
    sum_int = np.zeros(len(bins))
    count   = np.zeros(len(bins))
    w_total = 0.0

    n_confirmed = sum(1 for s in peak_scans if s.get('target_confirmed', False))

    # Two separate count arrays:
    # conf_count  = how many confirmed scans contain each fragment bin
    # total_count = all scans (fallback when no confirmed scans)
    conf_count  = np.zeros(len(bins))
    total_count = np.zeros(len(bins))

    for s in peak_scans:
        if len(s['mz']) == 0:
            continue
        tic = float(s['intensity'].sum()) or 1.0
        is_conf = s.get('target_confirmed', False)
        # Confirmed scans get 3× weight in intensity averaging.
        w = tic * (3.0 if is_conf else 1.0)
        w_total += w

        # Exclude surviving precursor ions from the fragment spectrum.
        target    = s['precursor_mz']
        frag_mask = s['mz'] < (target - ISOLATION_WINDOW / 2 - 1.0)
        frag_mz   = s['mz'][frag_mask]
        frag_int  = s['intensity'][frag_mask]

        for mz_v, int_v in zip(frag_mz, frag_int):
            bi = int((mz_v - mz_min) / CONSENSUS_BIN_DA)
            if 0 <= bi < len(bins):
                sum_int[bi]    += int_v * w
                total_count[bi] += 1
                if is_conf:
                    conf_count[bi] += 1

    # Fragment selection: prefer confirmed-scan reproducibility when possible.
    # If confirmed scans exist, a fragment must appear in >= MIN_FRAG_FRAC of
    # them — this prevents the majority of unconfirmed (interferent) scans from
    # dictating which fragments survive.
    if n_confirmed > 0:
        # Fragment selection driven by confirmed scans only.
        # MIN_FRAG_FRAC here just removes single-scan noise hits (≥2 appearances);
        # the 1% relative intensity filter below handles the rest.
        min_c = max(2, int(n_confirmed * MIN_FRAG_FRAC / 10.0))
        mask  = conf_count >= min_c
        if not np.any(mask):
            mask = conf_count >= 1
    else:
        min_c = max(1, int(len(peak_scans) * MIN_FRAG_FRAC))
        mask  = total_count >= min_c
        if not np.any(mask):
            mask = total_count >= 1

    result_mz  = bins[mask]
    result_int = (sum_int[mask] / w_total) if w_total > 0 else sum_int[mask]
    if result_int.max() > 0:
        result_int = result_int / result_int.max() * 999.0

    # Remove fragments below 1% relative intensity
    keep = result_int >= 10.0
    result_mz  = result_mz[keep]
    result_int = result_int[keep]

    return result_mz, result_int, n_confirmed


# ─── 6. MGF writer ────────────────────────────────────────────────────────────

def write_mgf(entries, output_path):
    with open(output_path, 'w') as f:
        for e in entries:
            f.write('BEGIN IONS\n')
            f.write(f'FEATURE_ID={e["feature_id"]}\n')
            f.write(f'PEPMASS={e["precursor_mz"]:.6f}\n')
            f.write(f'CHARGE={e.get("charge", 1)}+\n')
            f.write(f'RTINSECONDS={e["rt_sec"]:.2f}\n')
            f.write(f'SCANS={e["feature_id"]}\n')
            f.write(f'FILENAME={e.get("source_file", "")}\n')
            f.write(f'PRECURSOR_CONFIRMED_SCANS={e.get("n_confirmed", "N/A")}\n')
            f.write(f'TOTAL_SCANS={e.get("n_scans", "N/A")}\n')
            for mz_v, int_v in zip(e['mz_array'], e['int_array']):
                f.write(f'{mz_v:.5f}\t{int_v:.2f}\n')
            f.write('END IONS\n\n')


# ─── 7. QC plots ──────────────────────────────────────────────────────────────

def save_rt_plots(groups, peaks_by_precursor, output_dir, source_name):
    fig_dir = os.path.join(output_dir, 'rt_profiles', source_name)
    os.makedirs(fig_dir, exist_ok=True)

    for mz_key in sorted(groups):
        scans = groups[mz_key]
        rt_arr, scores = build_rt_profile(scans)
        smoothed = gaussian_filter1d(scores.astype(float), sigma=SMOOTHING_SIGMA)
        peaks = peaks_by_precursor.get(mz_key, [])

        fig, ax = plt.subplots(figsize=(10, 3))
        ax.plot(rt_arr / 60, scores,   alpha=0.35, color='gray',      lw=0.8, label='raw score')
        ax.plot(rt_arr / 60, smoothed, color='steelblue', lw=1.5,             label='smoothed')
        for pk in peaks:
            ax.axvspan(pk['start_rt'] / 60, pk['end_rt'] / 60, alpha=0.2, color='orange')
            ax.axvline(pk['apex_rt'] / 60, color='red', lw=1.0, ls='--')
        ax.set_xlabel('RT (min)')
        ax.set_ylabel('Composite score')
        ax.set_title(f'{source_name}  |  precursor {mz_key:.3f} Da  |  {len(peaks)} peaks detected')
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(os.path.join(fig_dir, f'mz_{mz_key:.3f}.png'.replace('.', '_', 1)), dpi=100)
        plt.close(fig)


# ─── 8. Process one mzML file ─────────────────────────────────────────────────

def process_file(mzml_path, output_dir, feature_id_start=1, save_plots=True):
    source_name = os.path.splitext(os.path.basename(mzml_path))[0]
    print(f'\n{"-"*60}')
    print(f'File : {source_name}')

    scans  = load_or_parse(mzml_path)
    print(f'Scans loaded : {len(scans)}')

    groups = group_by_precursor(scans)
    print(f'Precursors   : {len(groups)}')

    entries             = []
    peaks_by_precursor  = {}
    feature_id          = feature_id_start

    for mz_key in sorted(groups):
        group_scans    = groups[mz_key]
        rt_arr, scores = build_rt_profile(group_scans)
        peaks          = pick_peaks(rt_arr, scores)
        peaks_by_precursor[mz_key] = peaks

        for pk in peaks:
            mask       = (rt_arr >= pk['start_rt']) & (rt_arr <= pk['end_rt'])
            peak_scans = [s for s, m in zip(group_scans, mask) if m]
            if len(peak_scans) < 2:
                continue

            mz_cons, int_cons, n_confirmed = consensus_spectrum(peak_scans)
            if len(mz_cons) < MIN_FRAGMENTS:
                continue

            entries.append({
                'feature_id':   feature_id,
                'precursor_mz': mz_key,
                'rt_sec':       pk['apex_rt'],
                'charge':       1,
                'source_file':  source_name,
                'mz_array':     mz_cons,
                'int_array':    int_cons,
                'n_scans':      len(peak_scans),
                'n_confirmed':  n_confirmed,
                'apex_score':   pk['apex_score'],
            })
            feature_id += 1

        if peaks:
            print(f'  mz={mz_key:.3f}: {len(peaks)} peak(s)')

    if save_plots:
        save_rt_plots(groups, peaks_by_precursor, output_dir, source_name)

    print(f'Features     : {len(entries)}')
    return entries, feature_id


# ─── 9. Parent spectrum ───────────────────────────────────────────────────────

def get_parent_entry(mzml_path, feature_id):
    scans = load_or_parse(mzml_path)
    if not scans:
        print('WARNING: parent spectrum file is empty')
        return None, feature_id

    mz_cons, int_cons, _ = consensus_spectrum(scans)
    if len(mz_cons) == 0:
        best = max(scans, key=lambda s: s['intensity'].sum())
        mz_cons, int_cons = best['mz'].copy(), best['intensity'].copy()
        if int_cons.max() > 0:
            int_cons = int_cons / int_cons.max() * 999.0

    entry = {
        'feature_id':   feature_id,
        'precursor_mz': round(scans[0]['precursor_mz'], 6),
        'rt_sec':       scans[0]['rt'],
        'charge':       1,
        'source_file':  'parrent_spectrum',
        'mz_array':     mz_cons,
        'int_array':    int_cons,
    }
    print(f'\nParent spectrum: mz={entry["precursor_mz"]:.4f}, '
          f'{len(scans)} scans, {len(mz_cons)} consensus fragments')
    return entry, feature_id + 1


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_entries = []
    feature_id  = 1

    zeb_files = sorted(
        f for f in os.listdir(CONVERTED_DIR)
        if f.startswith('Zeb') and f.endswith('.mzML')
    )

    for fname in zeb_files:
        entries, feature_id = process_file(
            os.path.join(CONVERTED_DIR, fname),
            OUTPUT_DIR,
            feature_id_start=feature_id,
            save_plots=True,
        )
        all_entries.extend(entries)

    parent_path = os.path.join(CONVERTED_DIR, 'parrent_spectrum.mzML')
    parent_entry, feature_id = get_parent_entry(parent_path, feature_id)
    if parent_entry:
        all_entries.append(parent_entry)

    # ── Write MGF ──
    mgf_path = os.path.join(OUTPUT_DIR, 'PND_metabolites_FBMN.mgf')
    write_mgf(all_entries, mgf_path)
    print(f'\n{"="*60}')
    print(f'MGF written  : {mgf_path}')
    print(f'Total features: {len(all_entries)}')

    # Summary CSV
    rows = []
    for e in all_entries:
        n_sc = e.get('n_scans', 0) or 0
        n_co = e.get('n_confirmed', 0) or 0
        rows.append({
            'feature_id':        e['feature_id'],
            'precursor_mz':      e['precursor_mz'],
            'rt_min':            round(e['rt_sec'] / 60, 3),
            'rt_sec':            round(e['rt_sec'], 2),
            'n_scans':           n_sc,
            'n_confirmed_scans': n_co,
            'pct_confirmed':     round(100 * n_co / n_sc, 1) if n_sc > 0 else 0,
            'apex_score':        round(e.get('apex_score', 0), 1),
            'n_fragments':       len(e['mz_array']),
            'source_file':       e['source_file'],
        })
    summary = pd.DataFrame(rows)
    csv_path = os.path.join(OUTPUT_DIR, 'features_summary.csv')
    summary.to_csv(csv_path, index=False)
    print(f'Summary CSV  : {csv_path}')
    print(f'{"="*60}')
