import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import warnings
from scipy.stats import entropy
from scipy.signal import butter, sosfiltfilt, welch, iirnotch
from sklearn.decomposition import FastICA

# =============================================================================
# SAIM Analysis Pipeline v1.0 
# =============================================================================
# Author: Takafumi Shiga
# Date: January 2026
# =============================================================================

warnings.filterwarnings('ignore')

class SAIMConfig:
    # --- Analysis Parameters ---
    WINDOW_SEC = 10
    STEP_SEC = 2
    FS = 256.0
    EPS = 1e-9

    # --- Frequency Bands ---
    BANDS = {
        'Delta': (1, 4),
        'Theta': (4, 8),
        'Alpha': (8, 13),
        'Beta':  (13, 30),
        'Gamma': (30, 45)
    }

    # --- Metrics to Export ---
    METRIC_COLS = [
        'PE', 'F', 'SII', 'I',
        'NCI', 'NCI_Vol', 'LZC', 'SOM', 'HEMO', 'AUT',
        'Gamma', 'Beta', 'Alpha', 'Theta', 'Delta',
        'HSI_Gamma', 'HSI_Beta', 'HSI_Alpha', 'HSI_Theta', 'HSI_Delta'
    ]

    # --- HEMO Configuration (Muse S Gen2 / Athena MS-03 Mapping) ---
    HEMO_SIGNAL_CHANNELS = [
        'Optics1', 'Optics2', 'Optics3', 'Optics4',
        'Optics5', 'Optics6', 'Optics7', 'Optics8'
    ]
    HEMO_AMBIENT_CHANNELS = ['Optics11', 'Optics12', 'Optics15', 'Optics16']
    
    # --- Visualization Colors ---
    COLORS = {
        'PE': '#DC143C',   'F': '#8B0000',    'SII': '#800080',  'I': '#4682B4',
        'NCI': '#00BFFF',  'NCI_Vol': '#1E90FF', 'LZC': '#DAA520',
        'SOM': '#FFD700',  'HEMO': '#FF4500', 'AUT': '#2E8B57',
        'Gamma': '#FF1493', 'Beta': '#8A2BE2', 'Alpha': '#32CD32', 'Theta': '#FF8C00', 'Delta': '#708090',
        'HSI_Gamma': '#FF69B4', 'HSI_Beta': '#9370DB', 'HSI_Alpha': '#98FB98', 
        'HSI_Theta': '#FFA07A', 'HSI_Delta': '#B0C4DE'
    }

# --- Signal Processing Helper Functions ---
def filter_data(data, fs):
    sos = butter(4, [0.5, 100.0], btype='band', fs=fs, output='sos')
    filt = sosfiltfilt(sos, data, axis=1)
    for f in [50, 60]:
        b, a = iirnotch(f, 30, fs)
        from scipy.signal import filtfilt
        filt = filtfilt(b, a, filt, axis=1)
    return filt

def clean_with_ica(data):
    try:
        ica = FastICA(n_components=4, random_state=42, max_iter=200)
        S = ica.fit_transform(data.T)
        A = ica.mixing_
        vars = np.var(S, axis=0)
        if np.max(vars) > 5 * np.median(vars):
            S[:, np.argmax(vars)] = 0
        return np.dot(S, A.T).T
    except:
        return data

def calc_entropy(s):
    try:
        counts, _ = np.histogram(s.dropna(), bins='fd')
        if np.sum(counts) == 0: return 0.0
        p = counts / np.sum(counts)
        p = p[p > 0]
        ent = -np.sum(p * np.log(p))
        n_bins = len(counts)
        if n_bins > 1: return ent / np.log(n_bins)
        else: return 0.0
    except: return np.nan

def calc_lzc(s):
    try:
        b = (s > s.median()).astype(int)
        seq = "".join(map(str, b)); n = len(seq)
        if n == 0: return np.nan
        i, k, l = 0, 1, 1; k_max = 1
        while i+k+l-1 < n:
            if seq[i+k : i+k+l] in seq[i : i+k+l-1]: l += 1
            else: k_max += 1; i += k; k = 1; l = 1
        return k_max / (n / np.log2(n))
    except: return np.nan

# --- Main Analysis Class ---
class SAIMAnalyzer:
    def __init__(self, subject_id, visit_id, file_map):
        self.subject_id = subject_id
        self.visit_id = visit_id
        self.file_map = file_map 
        self.prefix = f"{self.subject_id}_{self.visit_id}"
        self.results = []
        
    def _process_window(self, df_win):
        metrics = {}
        
        # 1. EEG Processing
        raw_cols = ['RAW_TP9', 'RAW_AF7', 'RAW_AF8', 'RAW_TP10']
        if all(c in df_win.columns for c in raw_cols):
            raw = df_win[raw_cols].values.T
            clean = clean_with_ica(filter_data(raw, SAIMConfig.FS))
            freqs, psd = welch(clean, fs=SAIMConfig.FS, nperseg=int(SAIMConfig.FS))
            
            for band, (l, h) in SAIMConfig.BANDS.items():
                idx = np.logical_and(freqs >= l, freqs <= h)
                pwr = np.sum(psd[:, idx], axis=1)
                metrics[band] = np.mean(pwr)
                P_L = np.mean(pwr[[0, 1]]); P_R = np.mean(pwr[[2, 3]])
                metrics[f'HSI_{band}'] = 1.0 - abs(P_L - P_R)/(P_L + P_R + SAIMConfig.EPS)
                
            metrics['PE'] = metrics['Gamma'] / metrics['Delta'] if metrics['Delta'] > 0 else 0
            metrics['NCI'] = np.mean([calc_entropy(df_win[c]) for c in raw_cols])
            metrics['LZC'] = np.mean([calc_lzc(df_win[c]) for c in raw_cols])

        # 2. Physiological
        acc = [c for c in df_win.columns if 'Accelerometer' in c]
        if acc: metrics['SOM'] = 1.0 / (1.0 + np.std(np.sqrt(np.sum(df_win[acc]**2, axis=1))))
        
        # 3. HEMO (v1.0 Logic)
        sig_cols = [c for c in SAIMConfig.HEMO_SIGNAL_CHANNELS if c in df_win.columns]
        amb_cols = [c for c in SAIMConfig.HEMO_AMBIENT_CHANNELS if c in df_win.columns]
        
        if sig_cols:
            sig_var = df_win[sig_cols].var().mean()
            amb_var = df_win[amb_cols].var().mean() if amb_cols else 0.0
            metrics['HEMO'] = 1.0 / (1.0 + (sig_var / (amb_var + SAIMConfig.EPS)))
        else:
            metrics['HEMO'] = np.nan
        
        if 'Heart_Rate' in df_win.columns: metrics['AUT'] = calc_entropy(df_win['Heart_Rate'])

        # 4. Systemic Integration
        comps = [metrics.get(k, np.nan) for k in ['HSI_Gamma', 'SOM', 'HEMO', 'AUT']]
        metrics['I'] = np.nanmean(comps)
        if not np.isnan(metrics['I']) and 'PE' in metrics:
            metrics['F'] = (1 - metrics['I']) + 0.5 * metrics['PE']
            metrics['SII'] = 1.0 / (1.0 + np.exp(-5 * (metrics['I'] - metrics['PE'])))
            
        return metrics

    def process(self):
        print(f"Processing Session: {self.subject_id} - {self.visit_id}")
        
        # Define strict phase order for concatenation
        phase_order = [
            'Pre_BL1', 'Pre_Stress', 'Pre_BL2',
            'Post_BL1', 'Post_Stress', 'Post_BL2'
        ]
        
        global_time_offset = 0.0
        
        for phase_name in phase_order:
            if phase_name not in self.file_map: continue
            
            fpath = self.file_map[phase_name]
            print(f"  -> Analyzing Phase: {phase_name} ({os.path.basename(fpath)})")
            
            try:
                df = pd.read_csv(fpath)
                if 'TimeStamp' in df.columns: df['TimeStamp'] = pd.to_datetime(df['TimeStamp'], errors='coerce')
                df = df.dropna(subset=['RAW_TP9']).sort_values('TimeStamp').reset_index(drop=True)
                
                if len(df) < SAIMConfig.FS * 5: continue
                
                start_t = df['TimeStamp'].iloc[0]
                end_t = df['TimeStamp'].iloc[-1]
                curr_t = start_t
                
                condition = "Stress" if "Stress" in phase_name else "Rest"
                phase_group = phase_name.split('_')[0] # 'Pre' or 'Post'
                
                while curr_t + pd.Timedelta(seconds=SAIMConfig.WINDOW_SEC) <= end_t:
                    end_win = curr_t + pd.Timedelta(seconds=SAIMConfig.WINDOW_SEC)
                    win = df[(df['TimeStamp'] >= curr_t) & (df['TimeStamp'] < end_win)]
                    
                    if len(win) > SAIMConfig.FS:
                        m = self._process_window(win)
                        
                        # Stitching & Metadata
                        local_time = (curr_t - start_t).total_seconds()
                        m['Time'] = global_time_offset + local_time
                        m['Phase'] = phase_name
                        m['Phase_Group'] = phase_group # Crucial for Overall Stats
                        m['Session'] = self.visit_id
                        m['Condition'] = condition
                        self.results.append(m)
                        
                    curr_t += pd.Timedelta(seconds=SAIMConfig.STEP_SEC)
                
                if not df.empty:
                    duration = (end_t - start_t).total_seconds()
                    global_time_offset += duration + 10.0 
                    
            except Exception as e:
                print(f"Error in {phase_name}: {e}")

        self._save_outputs()

    def _plot_continuous_dynamics(self, df):
        prefix = f"{self.subject_id}_{self.visit_id}"
        print(f"   -> Generating Continuous Dynamics Plot for {prefix}...")

        target_metrics = ['F', 'PE', 'SII', 'I', 'NCI', 'LZC', 'SOM', 'HEMO', 'AUT']
        valid_metrics = [m for m in target_metrics if m in df.columns]

        if not valid_metrics: return

        n_metrics = len(valid_metrics)
        fig, axes = plt.subplots(n_metrics, 1, figsize=(15, 2.5 * n_metrics), sharex=True)
        if n_metrics == 1: axes = [axes]

        phases = df['Phase'].unique()
        phase_changes = []
        for p in phases:
            t = df[df['Phase'] == p]['Time'].min()
            phase_changes.append((p, t))

        for i, metric in enumerate(valid_metrics):
            ax = axes[i]
            color = SAIMConfig.COLORS.get(metric, 'black')

            series = df[metric].interpolate(method='linear', limit_direction='both')
            smooth_data = series.rolling(window=5, center=True, min_periods=1).mean()

            ax.plot(df['Time'], smooth_data, color=color, linewidth=2.0, label=metric)

            for p_name, p_time in phase_changes:
                ax.axvline(x=p_time, color='gray', linestyle='--', alpha=0.5)
                if i == 0:
                    ax.text(p_time, ax.get_ylim()[1], p_name, rotation=45, va='bottom', fontsize=8)

            ylabel = metric
            if metric == 'HSI_Gamma': ylabel = 'HSI (Gamma)'
            ax.set_ylabel(ylabel, fontweight='bold', fontsize=9)
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel('Time (s)')
        plt.suptitle(f'SAIM Continuous Dynamics: {self.subject_id} ({self.visit_id})', fontweight='bold', y=1.01)
        plt.tight_layout()
        plt.savefig(f"{prefix}_Continuous_Dynamics.png", dpi=300, bbox_inches='tight')
        plt.close()

    def _save_outputs(self):
        if not self.results: return
        df = pd.DataFrame(self.results)
        prefix = f"{self.subject_id}_{self.visit_id}"
        
        if 'NCI' in df.columns:
            df['NCI_Vol'] = df['NCI'].rolling(5, center=True).std().fillna(0)

        # 1. Save Full Time Series CSV
        df.to_csv(f"{prefix}_TimeSeries.csv", index=False)

        # 2. Save Overall Stats CSV (Aggregated by Phase_Group: Pre vs Post)
        metric_cols = [c for c in SAIMConfig.METRIC_COLS if c in df.columns]
        df.groupby(['Session', 'Phase_Group'])[metric_cols].mean().to_csv(f"{prefix}_Overall_Stats.csv")
        
        # 3. Save SubPhase Stats CSV
        df.groupby(['Session', 'Phase'])[metric_cols].agg(['mean', 'std']).to_csv(f"{prefix}_SubPhase_Stats.csv")

        # 4. Plots
        self._plot_continuous_dynamics(df)
        self.plot_omni(df)

        print(f"Generated outputs for {prefix}")

    def plot_omni(self, df):
        fig, axes = plt.subplots(5, 4, figsize=(20, 20))
        axes = axes.flatten()
        plot_list = [
            'F', 'PE', 'SII', 'SOM',
            'I', 'LZC', 'NCI', 'HEMO',
            'Gamma', 'Beta', 'Alpha', 'Theta',
            'Delta', 'HSI_Gamma', 'HSI_Beta', 'HSI_Alpha',
            'HSI_Theta', 'HSI_Delta', 'AUT', 'NCI_Vol'
        ]
        stats = df.groupby('Phase')[plot_list].mean().reset_index()
        phase_order = ['Pre_BL1', 'Pre_Stress', 'Pre_BL2', 'Post_BL1', 'Post_Stress', 'Post_BL2']
        stats['Phase'] = pd.Categorical(stats['Phase'], categories=phase_order, ordered=True)
        stats = stats.sort_values('Phase')

        for i, metric in enumerate(plot_list):
            if i < len(axes) and metric in stats.columns:
                sns.barplot(data=stats, x='Phase', y=metric, ax=axes[i], 
                            color=SAIMConfig.COLORS.get(metric, 'grey'))
                axes[i].set_title(metric)
                axes[i].tick_params(axis='x', rotation=45)
        plt.tight_layout()
        out_name = f"{self.prefix}_OmniPanel.png"
        plt.savefig(out_name, dpi=300)
        plt.close()

def main():
    files = glob.glob("*.csv")
    dataset = {}
    
    print("Scanning files...")
    for f in files:
        if "Stats" in f or "TimeSeries" in f: continue
        
        parts = f.replace('.csv', '').split('_')
        if len(parts) >= 5:
            subj = parts[0]
            visit = parts[2]
            phase_name = f"{parts[3]}_{parts[4]}"
            key = (subj, visit)
            if key not in dataset: dataset[key] = {}
            dataset[key][phase_name] = f

    print(f"Found {len(dataset)} sessions to process.")
    for (subj, visit), file_map in dataset.items():
        analyzer = SAIMAnalyzer(subj, visit, file_map)
        analyzer.process()

if __name__ == "__main__":
    main()
