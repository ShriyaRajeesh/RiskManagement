# plot_results.py
# Usage: python plot_results.py
# Reads returns_<scenario>.csv files and scenario_summary.csv (created by stress_test.cpp)
# Produces:
#  - Overlaid histogram of scenario sample returns
#  - Bar chart of VaR/CVaR/Mean per scenario

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# find sample files created by the C++ program
sample_files = sorted(glob.glob("returns_*.csv"))
if len(sample_files) == 0:
    print("No returns_*.csv files found in current directory. Run the C++ program first.")
    # Optional: create a fake demo to show code works
    demo = True
else:
    demo = False

# Load summary CSV if present
summary_file = "scenario_summary.csv"
if os.path.exists(summary_file):
    summary_df = pd.read_csv(summary_file)
else:
    summary_df = None

# Helper to load a sample file
def load_sample(fname):
    try:
        df = pd.read_csv(fname)
        if "return" in df.columns:
            return df["return"].values
        else:
            # assume single column of numeric values
            return df.iloc[:,0].values
    except Exception as e:
        print(f"Error reading {fname}: {e}")
        return np.array([])

# If no real files, produce demo random data
if demo:
    print("Creating demo synthetic data (no real CSVs present).")
    scenario_names = ["Baseline", "Rate-Hike+200bp", "Market-Crash-30pct", "Rate+Crash"]
    samples = {
        "Baseline": np.random.normal(loc=0.02, scale=0.06, size=50000),
        "Rate-Hike+200bp": np.random.normal(loc=0.0, scale=0.07, size=50000),
        "Market-Crash-30pct": np.random.normal(loc=-0.12, scale=0.15, size=50000),
        "Rate+Crash": np.random.normal(loc=-0.10, scale=0.16, size=50000)
    }
else:
    # load all sample files
    samples = {}
    scenario_names = []
    for f in sample_files:
        # turn filename into scenario name
        # returns_Rate-Hike+200bp.csv -> Rate-Hike+200bp
        name = os.path.basename(f)[len("returns_"):-len(".csv")]
        arr = load_sample(f)
        if arr.size == 0: continue
        samples[name] = arr
        scenario_names.append(name)

# Plot overlaid histograms
plt.figure(figsize=(11,7))
bins = 200
for name in scenario_names:
    r = samples[name]
    plt.hist(r, bins=bins, density=True, alpha=0.45, label=name, histtype='stepfilled')
    # compute VaR95
    var95 = np.percentile(r, 5)
    ymax = plt.gca().get_ylim()[1]
    plt.axvline(var95, linestyle='--', linewidth=1)
    plt.text(var95, ymax*0.6, f"{name} VaR95={var95:.3f}", rotation=90, fontsize=8)

plt.title("Portfolio return distributions under scenarios")
plt.xlabel("Return")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()

# If summary exists, plot VaR/CVaR/Mean bar chart
if summary_df is not None:
    # Order bars to match summary_df order
    names = list(summary_df['scenario'])
    vars_ = list(summary_df['VaR95'])
    cvars_ = list(summary_df['CVaR95'])
    means_ = list(summary_df['mean'])

    x = np.arange(len(names))
    width = 0.25
    plt.figure(figsize=(10,5))
    plt.bar(x - width, vars_, width, label='VaR95')
    plt.bar(x, cvars_, width, label='CVaR95')
    plt.bar(x + width, means_, width, label='Mean')
    plt.xticks(x, names, rotation=30)
    plt.ylabel("Return")
    plt.title("Risk metrics by scenario (from scenario_summary.csv)")
    plt.legend()
    plt.tight_layout()
    plt.show()

# Also plot a simple violin-like overlay using boxplots for quick comparison
plt.figure(figsize=(10,6))
data_for_box = [samples[n] for n in scenario_names]
plt.boxplot(data_for_box, labels=scenario_names, showfliers=False)
plt.title("Return distribution summary (boxplot) by scenario")
plt.ylabel("Return")
plt.xticks(rotation=20)
plt.tight_layout()
plt.show()
