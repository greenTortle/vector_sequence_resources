"""
plot_k_zgc.py  —  plot the CSV produced by k_zgc_scanner
=========================================================
"""

import sys
import pathlib
import csv
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector


def load_csv(path: str):
    metadata = {}
    k_vals, z_vals, g_vals, c_vals = [], [], [], []

    with open(path, newline="", encoding="utf-8") as f:
        lines = f.readlines()

    # Extract metadata
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#") and ":" in stripped:
            try:
                key_part, value_part = stripped[1:].split(":", 1)
                metadata[key_part.strip()] = value_part.strip()
            except:
                continue

    # Read data with more robust parsing
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header_found = False

        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue

            # Detect header
            if not header_found:
                if "k_decimal" in row[0]:
                    header_found = True
                continue

            # Process data row
            if len(row) >= 4:
                try:
                    k_vals.append(float(row[0].strip()))
                    z_vals.append(int(row[1].strip()))
                    g_vals.append(int(row[2].strip()))
                    c_vals.append(int(row[3].strip()))
                except (ValueError, IndexError):
                    continue

    return k_vals, z_vals, g_vals, c_vals, metadata


def plot(k, z, g, c, metadata):
    rule = metadata.get("Rule", "Unknown Rule")
    k_range = metadata.get("k range", "Unknown")
    x0_range = metadata.get("x0 range", "Unknown")
    check_n = metadata.get("check_n (max v_n)", "Unknown")

    fig, ax = plt.subplots(figsize=(13, 7))
    
    ax.plot(k, z, label="Z: zeroed", color="steelblue", linewidth=1.5)
    ax.plot(k, g, label="G: grown",  color="darkorange", linewidth=1.5)
    ax.plot(k, c, label="C: cycled", color="seagreen", linewidth=1.1, alpha=0.85)

    full_title = f"{rule}\n"
    full_title += f"k range: {k_range}   |   x₀: {x0_range}   |   max v_n: {check_n}\n"
    full_title += f"{len(k)} data points"

    ax.set_title(full_title, fontsize=14, pad=20)
    ax.set_xlabel("k")
    ax.set_ylabel("Number of sequences")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    def on_select(xmin, xmax):
        if abs(xmax - xmin) < 1e-12:
            return
        lo, hi = sorted((xmin, xmax))
        print(f"\nSelected k range: {lo:.10g} .. {hi:.10g}")
        print(f"Re-run: ./k_zgc_scanner {lo:.10g} {hi:.10g} <points> <x0s> <x0e> <chkn> <rule>")

    _span = SpanSelector(ax, on_select, "horizontal",
                         useblit=True, props=dict(alpha=0.25, facecolor="gray"),
                         interactive=False)

    plt.show()


if __name__ == "__main__":
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "k_zgc_results.csv"

    if not pathlib.Path(csv_path).exists():
        print(f"File not found: {csv_path}")
        sys.exit(1)

    print(f"Loading {csv_path} ...")
    k, z, g, c, meta = load_csv(csv_path)
    print(f"  {len(k)} data points loaded.")

    if len(k) == 0:
        print("No data points found.")
    else:
        plot(k, z, g, c, meta)