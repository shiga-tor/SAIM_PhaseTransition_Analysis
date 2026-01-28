"""
calibration_analysis.py

Quality Control / Calibration Analysis Script for Muse EEG Data.
Implements a Hybrid Logic:
1. Validated Manual Timings for S001 (Golden Standard).
2. Automated '59-second Rule' + Sequential Peak Detection for S002+.

Author: Takafumi Shiga
Date: 2026
"""

import numpy as np
import pandas as pd
import glob
import os
from scipy.signal import welch, butter, filtfilt
from scipy.linalg import eigh
from scipy.stats import linregress
from scipy.ndimage import label

# ==========================================
# 1. Automated Detection Engine (For S002+)
# ==========================================
def detect_timings_auto(df, fs=256):
    """
    Automated timing detection logic for subject S002 and later.
    Rule: Events occurring after 59 seconds are assigned sequentially:
          Jaw -> Eye -> Blink.
    """
    # Channel selection
    raw_cols = [c for c in df.columns if 'RAW' in c and ('TP9' in c or 'AF7' in c or 'AF8' in c or 'TP10' in c)]
    if len(raw_cols) < 4:
        alt = ['TP9', 'AF7', 'AF8', 'TP10']
        raw_cols = []
        for t in alt:
             match = [c for c in df.columns if t in c and 'Delta' not in c]
             if match: raw_cols.append(match[0])
    if not raw_cols: return None

    # Data extraction & Interpolation (Handle NaNs)
    data_df = df[raw_cols].interpolate(method='linear', limit_direction='both').fillna(0)
    data = data_df.values
    
    # High-pass Filtering (>30Hz for EMG detection)
    b, a = butter(4, 30 / (fs/2), btype='highpass')
    try:
        filt_data = filtfilt(b, a, data, axis=0)
    except: return None

    # Power Calculation (Log scale to compress dynamic range)
    power = np.mean(filt_data**2, axis=1)
    power = np.maximum(power, 1e-10)
    log_power = np.log10(power)
    
    # Smoothing
    window = int(fs * 0.5)
    smooth_log_power = pd.Series(log_power).rolling(window=window, center=True).mean().fillna(np.min(log_power)).values
    
    # --- 59-Second Rule ---
    scan_start_idx = int(fs * 59.0)
    
    # Threshold calculation (Based only on data after 59s)
    if len(smooth_log_power) > scan_start_idx + (fs * 5):
        ref_data = smooth_log_power[scan_start_idx:]
    else:
        ref_data = smooth_log_power

    floor = np.percentile(ref_data, 5)
    ceiling = np.percentile(ref_data, 99)
    threshold = floor + (ceiling - floor) * 0.30 # Sensitivity: 30%
    
    # Binary activation detection
    is_active = smooth_log_power > threshold
    if len(is_active) > scan_start_idx:
        is_active[:scan_start_idx] = False # Force ignore events before 59s
    
    labeled_array, num_features = label(is_active)
    
    raw_events = []
    for i in range(1, num_features+1):
        idx = np.where(labeled_array == i)[0]
        if len(idx) > fs * 2.0: # Keep events longer than 2s
            raw_events.append((idx[0], idx[-1]))
    
    # Merge events (Fill gaps smaller than 5s)
    merged_events = []
    if raw_events:
        curr_s, curr_e = raw_events[0]
        for next_s, next_e in raw_events[1:]:
            if (next_s - curr_e) < fs * 5.0: 
                curr_e = next_e
            else:
                merged_events.append((curr_s, curr_e))
                curr_s, curr_e = next_s, next_e
        merged_events.append((curr_s, curr_e))
    
    timings = {}
    
    # Rest (Fixed window: 10s to 55s)
    rest_end = 55.0
    if len(df)/fs < 55.0: rest_end = (len(df)/fs) - 5.0
    timings['Rest'] = (10.0, rest_end)
    
    # Tasks Assignment (Sequential)
    task_order = ['Jaw', 'Eye', 'Blink']
    for i, task_name in enumerate(task_order):
        if i < len(merged_events):
            s, e = merged_events[i]
            duration = (e - s) / fs
            margin = min(1.0, duration * 0.1)
            timings[task_name] = ((s/fs)+margin, (e/fs)-margin)
        else:
            timings[task_name] = None
            
    return timings

# ==========================================
# 2. Timing Determination Logic (Hybrid)
# ==========================================
def get_timings_hybrid(df, subject_id, fs=256):
    """
    Returns fixed timings for S001 (Ground Truth), 
    and automated detection for all other subjects.
    """
    if subject_id == "S001":
        # S001 Specific: Manually validated ground truth data
        # Jaw: 69-79s, Eye: 87-93s, Blink: 100s+
        return {
            'Rest': (10.0, 50.0),
            'Jaw':  (71.0, 77.0), # Stable interval
            'Eye':  (88.0, 92.0),
            'Blink':(102.0, 108.0)
        }
    else:
        # S002+: Automated detection
        return detect_timings_auto(df, fs)

# ==========================================
# 3. Common Analysis Class (Robust CCA)
# ==========================================
class RobustCCA:
    def __init__(self, n_components=4):
        self.n_components = n_components
        self.w_ = None 
    def fit(self, X):
        X = np.nan_to_num(X)
        X = X - np.mean(X, axis=1, keepdims=True)
        n_ch, n_times = X.shape
        if n_times < 2: return self
        X_t, X_t1 = X[:, :-1], X[:, 1:]
        try:
            R_xx = np.dot(X_t, X_t.T) / (n_times - 1)
            R_xy = np.dot(X_t, X_t1.T) / (n_times - 1)
            R_xy = (R_xy + R_xy.T) / 2
            eigvals, eigvecs = eigh(R_xy, R_xx)
            idx = np.argsort(eigvals)[::-1]
            self.w_ = eigvecs[:, idx[:self.n_components]]
        except: self.w_ = np.eye(n_ch)[:, :self.n_components]
        return self
    def transform(self, X):
        X = np.nan_to_num(X)
        X = X - np.mean(X, axis=1, keepdims=True)
        if self.w_ is None: return X
        return np.dot(self.w_.T, X)

def calc_simple_exponent(freqs, psd, f_range=[1, 50]):
    idx = np.logical_and(freqs >= f_range[0], freqs <= f_range[1])
    f_sel, p_sel = freqs[idx], psd[idx]
    valid = (f_sel > 0) & (p_sel > 0)
    f_sel, p_sel = f_sel[valid], p_sel[valid]
    if len(f_sel) < 2: return 0.0
    slope, _, _, _, _ = linregress(np.log10(f_sel), np.log10(p_sel))
    return -slope

def verify_segment(raw_df, start_sec, end_sec, fs=256):
    raw_cols = [c for c in raw_df.columns if 'RAW' in c and ('TP9' in c or 'AF7' in c or 'AF8' in c or 'TP10' in c)]
    if len(raw_cols) < 4:
        alt = ['TP9', 'AF7', 'AF8', 'TP10']
        raw_cols = []
        for t in alt:
            match = [c for c in raw_df.columns if t in c and 'Delta' not in c and 'Theta' not in c]
            if match: raw_cols.append(match[0])
    if len(raw_cols) < 4: return None
    start_idx = int(start_sec * fs)
    end_idx = int(end_sec * fs)
    if len(raw_df) < start_idx: return None
    if len(raw_df) < end_idx: end_idx = len(raw_df)
    data = raw_df[raw_cols].iloc[start_idx:end_idx].values.T
    data = np.nan_to_num(data)
    if data.shape[1] < fs: return None
    cca = RobustCCA(n_components=len(raw_cols))
    cca.fit(data)
    comps = cca.transform(data)
    freqs, psd = welch(comps, fs=fs, nperseg=min(fs*2, data.shape[1]))
    if len(freqs) == 0: return 0.0
    avg_psd = np.mean(psd, axis=0)
    exponent = calc_simple_exponent(freqs, avg_psd, f_range=[2, 45])
    return exponent

# ==========================================
# 4. Main Process (Production Mode)
# ==========================================
def batch_process_production(file_pattern="S*_Calibration.csv", fs=256):
    files = sorted(glob.glob(file_pattern))
    if not files:
        print(f"No files found matching: {file_pattern}")
        return

    print(f"Found {len(files)} files. Running PRODUCTION HYBRID ANALYSIS...\n")
    print(f"{'Subject':<8} | {'Rest':<5} | {'Jaw':<5} | {'Eye':<5} | {'Blink':<5} | {'Diff':<6} | {'Status'} | {'Timing Info (s)'}")
    print("-" * 115)

    for file_path in files:
        filename = os.path.basename(file_path)
        subj = filename.split('_')[0]

        try:
            df = pd.read_csv(file_path)
            
            # ★ Call Hybrid Detection Logic ★
            timings = get_timings_hybrid(df, subj, fs)
            
            if timings is None:
                print(f"{subj:<8} | ERROR: Signal processing failed")
                continue
                
            results = {}
            time_str = ""
            for task in ['Rest', 'Jaw', 'Eye', 'Blink']:
                t_range = timings.get(task)
                if t_range:
                    exp = verify_segment(df, t_range[0], t_range[1], fs)
                    results[task] = exp
                    time_str += f"{task[0]}:{int(t_range[0])}-{int(t_range[1])} "
                else:
                    results[task] = None
                    time_str += f"{task[0]}:-- "

            exp_rest = results['Rest']
            exp_jaw = results['Jaw']
            
            if exp_rest is not None and exp_jaw is not None:
                diff = exp_rest - exp_jaw
                # Pass Criteria: Must show significant difference
                if diff > 0.3: status = "OK"
                elif diff > 0.1: status = "WEAK"
                else: status = "CHECK"
                diff_str = f"{diff:.2f}"
            else:
                diff_str = "N/A"
                status = "ERR"

            def fmt(v): return f"{v:.2f}" if v is not None else " - "
            print(f"{subj:<8} | {fmt(results['Rest']):<5} | {fmt(results['Jaw']):<5} | {fmt(results['Eye']):<5} | {fmt(results['Blink']):<5} | {diff_str:<6} | {status:<6} | {time_str}")
            
        except Exception as e:
            print(f"{subj:<8} | ERROR: {str(e)}")

    print("-" * 115)

if __name__ == "__main__":
    # Execute
    batch_process_production("S*_Calibration.csv")
