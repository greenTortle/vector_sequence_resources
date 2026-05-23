"""
ZGC_animate.py  —  Animate k vs Z/G/C as q (hailstone coefficient) evolves
===========================================================================

Requires the compiled C scanner (k_zgc_scanner.exe / ZGC_grapher.exe) to be
in the same folder as this script. Build it with:

    gcc -O3 -fopenmp -o k_zgc_scanner k_zgc_scanner.c -lm          (Linux/Mac)
    gcc -O3 -fopenmp -o ZGC_grapher.exe ZGC_grapher.c -lm           (Windows)

USAGE:
    python ZGC_animate.py <q_start> <q_end> <k_points> <x0_start> <x0_end> <check_n> [options]

POSITIONAL ARGUMENTS (required, in this order):
    q_start    First hailstone coefficient to animate  (integer >= 1)
    q_end      Last  hailstone coefficient to animate  (integer >= q_start)
    k_points   Number of k sample points per frame
    x0_start   First x0 value
    x0_end     Last  x0 value
    check_n    Max iterations (max v_n) per sequence

    For each q the k range is automatically set to: 1.0 .. q + 1.001

OPTIONAL FLAGS:
    --step N       Only compute every Nth q value (default 1 = every q)
                   e.g. --step 10 with q 3..100 gives q = 3, 13, 23, ...
    --interval MS  Milliseconds each frame is displayed in the animation
                   (default 800 = 0.8 seconds per frame)
    --workers N    Number of q-values computed in parallel (default = all CPU cores)
                   Each worker also uses OpenMP threads internally. On a machine
                   with e.g. 8 cores, try --workers 2 with OMP_NUM_THREADS=4.
    --save FILE    Save the animation to a file instead of opening a window.
                   Supported formats:
                     .gif  — no extra dependencies needed
                     .mp4  — requires ffmpeg (see below)
    --dpi N        Resolution of saved file in dots-per-inch (default 100)
    --colors N     GIF palette size 2-256 (default 32, ignored for MP4)
    --width N      Figure width  in inches (default 10)
    --height N     Figure height in inches (default 5)

EXAMPLES:
    # Open interactive animation window, all q from 3 to 20
    python ZGC_animate.py 3 20 500 2 1000 200

    # Save as GIF, every 5th q, slow animation
    python ZGC_animate.py 3 100 500 2 1000 200 --step 5 --interval 1200 --save out.gif

    # Save as MP4 (much smaller file than GIF)
    python ZGC_animate.py 3 100 500 2 1000 200 --step 5 --save out.mp4

    # Limit parallelism (Windows PowerShell)
    $env:OMP_NUM_THREADS=4; python ZGC_animate.py 3 50 500 2 1000 200 --workers 2 --save out.gif

SAVING AS MP4 (10-20x smaller than GIF at the same quality):
    Install ffmpeg, then just use --save out.mp4. To install ffmpeg on Windows:
      winget install Gyan.FFmpeg        (recommended, built into Windows 10/11)
      choco install ffmpeg              (if you use Chocolatey)
    After installing, close and reopen PowerShell, then verify with: ffmpeg -version

CACHING:
    CSVs are saved to ./zgc_cache/ with filenames that encode all parameters.
    Re-running with the same parameters skips the C scanner entirely — only
    new q values are computed. Delete ./zgc_cache/ to force a full recompute.
"""

import sys
import os
import subprocess
import pathlib
import csv
import argparse
import matplotlib
# Backend is chosen after args are parsed (Agg for --save, TkAgg for interactive)

# ── locate the C binary ────────────────────────────────────────────────────────
def find_scanner():
    import shutil
    script_dir = pathlib.Path(__file__).parent.resolve()
    cwd = pathlib.Path.cwd().resolve()
    names = ["ZGC_grapher", "k_zgc_scanner"]
    exts  = ["", ".exe"]
    search_dirs = [cwd, script_dir]
    for d in search_dirs:
        for name in names:
            for ext in exts:
                p = d / (name + ext)
                if p.exists():
                    return str(p)
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None

# ── CSV helpers ────────────────────────────────────────────────────────────────
def load_csv(path):
    k_vals, z_vals, g_vals, c_vals = [], [], [], []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header_found = False
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if not header_found:
                if "k_decimal" in row[0]:
                    header_found = True
                continue
            if len(row) >= 4:
                try:
                    k_vals.append(float(row[0].strip()))
                    z_vals.append(int(row[1].strip()))
                    g_vals.append(int(row[2].strip()))
                    c_vals.append(int(row[3].strip()))
                except (ValueError, IndexError):
                    continue
    return k_vals, z_vals, g_vals, c_vals

def csv_path_for(cache_dir, q, k_pts, x0s, x0e, chkn):
    return cache_dir / f"zgc_q{q:04d}_kpts{k_pts}_x{x0s}-{x0e}_n{chkn}.csv"

# ── run scanner for one q (in its own temp working dir to avoid CSV collisions) ─
def run_scanner(args_tuple):
    """
    Runs the C scanner in an isolated temp directory so parallel jobs don't
    clobber each other's k_zgc_results.csv output.
    Returns (q, success).
    """
    scanner, q, k_pts, x0s, x0e, chkn, out_csv = args_tuple
    import tempfile, shutil

    if pathlib.Path(out_csv).exists():
        return (q, True)   # already cached — nothing to do

    k_start = "1.0"
    k_end   = f"{q + 1.001:.4f}"
    cmd = [scanner, k_start, k_end, str(k_pts),
           str(x0s), str(x0e), str(chkn), "H", str(q)]

    with tempfile.TemporaryDirectory() as tmp_dir:
        result = subprocess.run(
            cmd,
            cwd=tmp_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        tmp_csv = pathlib.Path(tmp_dir) / "k_zgc_results.csv"
        if result.returncode != 0 or not tmp_csv.exists():
            return (q, False)
        shutil.move(str(tmp_csv), str(out_csv))

    return (q, True)

# ── main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Animate k vs Z/G/C as hailstone coefficient q evolves."
    )
    parser.add_argument("q_start",  type=int)
    parser.add_argument("q_end",    type=int)
    parser.add_argument("k_points", type=int)
    parser.add_argument("x0_start", type=int)
    parser.add_argument("x0_end",   type=int)
    parser.add_argument("check_n",  type=int)
    parser.add_argument("--interval", type=int, default=800,
                        help="Milliseconds per frame (default 800)")
    parser.add_argument("--save", type=str, default=None,
                        help="Save animation to this file (e.g. out.gif or out.mp4)")
    parser.add_argument("--step", type=int, default=1,
                        help="Step size between q values (default 1, e.g. --step 10 gives q=3,13,23,...)")
    parser.add_argument("--dpi", type=int, default=100,
                        help="Resolution in dots-per-inch for saved file (default 100)")
    parser.add_argument("--colors", type=int, default=32,
                        help="GIF color palette size: 2-256 (default 32). "
                             "Ignored for MP4.")
    parser.add_argument("--width", type=float, default=10.0,
                        help="Figure width in inches (default 10)")
    parser.add_argument("--height", type=float, default=5.0,
                        help="Figure height in inches (default 5)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel q-jobs (default: number of CPU cores). "
                             "Each job already uses OpenMP threads internally, so "
                             "consider setting --workers to half your core count if "
                             "the machine feels overloaded.")
    args = parser.parse_args()

    q_start  = args.q_start
    q_end    = args.q_end
    k_pts    = args.k_points
    x0s      = args.x0_start
    x0e      = args.x0_end
    chkn     = args.check_n
    interval = args.interval
    dpi      = max(36, min(args.dpi, 300))
    colors   = max(2,  min(args.colors, 256))
    fig_w    = args.width
    fig_h    = args.height

    if q_start < 1 or q_end < q_start:
        print("Error: need 1 <= q_start <= q_end"); sys.exit(1)

    # Choose backend before importing pyplot:
    # Agg is headless (for --save); TkAgg needs a display (interactive)
    if args.save:
        matplotlib.use("Agg")
    else:
        matplotlib.use("TkAgg")

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from concurrent.futures import ProcessPoolExecutor, as_completed

    try:
        from tqdm import tqdm
    except ImportError:
        # Graceful fallback: tqdm is optional
        class tqdm:
            def __init__(self, it=None, total=None, desc="", unit="", **kwargs):
                self._it    = it
                self._total = total
                self._desc  = desc
                self._n     = 0
            def __iter__(self):
                for item in self._it:
                    yield item
            def update(self, n=1):
                self._n += n
                t = self._total or "?"
                print(f"\r{self._desc}: {self._n}/{t}", end="", flush=True)
            def close(self):
                print()

    scanner = find_scanner()
    if scanner is None:
        print("Error: cannot find ZGC_grapher or k_zgc_scanner binary in current directory or PATH.")
        print("Build it first:  gcc -O3 -fopenmp -o k_zgc_scanner k_zgc_scanner.c -lm")
        sys.exit(1)
    print(f"Using scanner: {scanner}")

    cache_dir = pathlib.Path("zgc_cache")
    cache_dir.mkdir(exist_ok=True)

    step = args.step
    if step < 1:
        print("Error: --step must be >= 1"); sys.exit(1)

    q_values = list(range(q_start, q_end + 1, step))

    # ── phase 1: generate CSVs in parallel, then load ─────────────────────────
    all_data = {}   # q -> (k, z, g, c)
    n_q = len(q_values)

    # Separate already-cached from needs-computing
    needs_run = []
    cached    = []
    for q in q_values:
        csv_file = csv_path_for(cache_dir, q, k_pts, x0s, x0e, chkn)
        if csv_file.exists():
            cached.append(q)
        else:
            needs_run.append(q)

    if cached:
        print(f"\n{len(cached)} q-value(s) already cached, "
              f"{len(needs_run)} to compute.\n")
    else:
        print(f"\nGenerating data for q = {q_start} .. {q_end}  ({n_q} frames)\n")

    # Decide worker count: default to os.cpu_count(), but never more than jobs
    import os
    max_workers = args.workers or os.cpu_count() or 1
    max_workers = min(max_workers, max(len(needs_run), 1))

    if needs_run:
        print(f"Running {len(needs_run)} scanner job(s) with up to {max_workers} parallel worker(s).")
        print("(Each worker also uses OpenMP threads internally.)\n")

        job_args = [
            (scanner, q, k_pts, x0s, x0e, chkn,
             csv_path_for(cache_dir, q, k_pts, x0s, x0e, chkn))
            for q in needs_run
        ]

        failed = []
        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            futures = {pool.submit(run_scanner, a): a[1] for a in job_args}
            pbar = tqdm(total=len(futures), desc="Computing q frames", unit="q")
            for fut in as_completed(futures):
                q_done = futures[fut]
                ok = False
                try:
                    _, ok = fut.result()
                except Exception as exc:
                    print(f"\n  [q={q_done}] Exception: {exc}")
                if not ok:
                    failed.append(q_done)
                pbar.update(1)
            pbar.close()

        if failed:
            print(f"\nWarning: {len(failed)} q-value(s) failed: {failed}")

    # Load all CSVs (cached + freshly computed)
    print("\nLoading CSVs...")
    for q in tqdm(q_values, desc="Loading", unit="q"):
        csv_file = csv_path_for(cache_dir, q, k_pts, x0s, x0e, chkn)
        if not csv_file.exists():
            continue
        k, z, g, c = load_csv(csv_file)
        if k:
            all_data[q] = (k, z, g, c)

    loaded_qs = sorted(all_data.keys())
    if not loaded_qs:
        print("No data loaded. Exiting."); sys.exit(1)
    print(f"\nLoaded {len(loaded_qs)} frames. Building animation...\n")

    # ── phase 2: animate ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    plt.subplots_adjust(top=0.82)

    # Compute a stable y-axis max across all frames so the plot doesn't jump
    global_ymax = max(
        max(max(z), max(g), max(c))
        for (k, z, g, c) in all_data.values()
        if z and g and c
    )

    # Initial frame
    q0 = loaded_qs[0]
    k0, z0, g0, c0 = all_data[q0]

    line_z, = ax.plot(k0, z0, color="steelblue",  linewidth=1.5, label="Z: zeroed")
    line_g, = ax.plot(k0, g0, color="darkorange", linewidth=1.5, label="G: grown")
    line_c, = ax.plot(k0, c0, color="seagreen",   linewidth=1.1, alpha=0.85, label="C: cycled")

    ax.set_xlabel("k", fontsize=12)
    ax.set_ylabel("Number of sequences", fontsize=12)
    ax.set_ylim(0, global_ymax * 1.08)
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    title = ax.set_title("", fontsize=14)

    def make_title(q):
        k_end_val = q + 1.001
        return (
            f"Hailstone ({q}x+1)\n"
            f"k range: 1.0 .. {k_end_val:.3f}   |   "
            f"x₀: {x0s} .. {x0e}   |   max v_n: {chkn}\n"
            f"{k_pts} k-points   |   frame {loaded_qs.index(q)+1}/{len(loaded_qs)}"
        )

    title.set_text(make_title(q0))

    def update(frame_idx):
        q = loaded_qs[frame_idx]
        k, z, g, c = all_data[q]
        line_z.set_data(k, z)
        line_g.set_data(k, g)
        line_c.set_data(k, c)
        ax.set_xlim(min(k), max(k))
        title.set_text(make_title(q))
        return line_z, line_g, line_c, title

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(loaded_qs),
        interval=interval,
        blit=False,
        repeat=True
    )

    if args.save:
        save_path = args.save
        ext = pathlib.Path(save_path).suffix.lower()
        fps = max(1, 1000 // interval)
        print(f"Saving animation to {save_path}  (dpi={dpi}) ...")
        if ext == ".gif":
            # PillowWriter accepts savefig_kwargs for dpi; color quantization
            # is applied frame-by-frame via a post-save Pillow reopen+quantize.
            writer = animation.PillowWriter(fps=fps)
            ani.save(save_path, writer=writer, dpi=dpi)
            # Post-process: re-quantize to reduce palette and optimize
            try:
                from PIL import Image
                print(f"Optimizing GIF palette ({colors} colors)...")
                src = Image.open(save_path)
                frames_out = []
                durations  = []
                try:
                    while True:
                        frame = src.convert("RGB").quantize(
                            colors=colors,
                            method=Image.Quantize.MEDIANCUT,
                            dither=Image.Dither.NONE
                        )
                        frames_out.append(frame)
                        durations.append(src.info.get("duration", 1000 // fps))
                        src.seek(src.tell() + 1)
                except EOFError:
                    pass
                if frames_out:
                    frames_out[0].save(
                        save_path,
                        save_all=True,
                        append_images=frames_out[1:],
                        loop=0,
                        duration=durations,
                        optimize=True,
                    )
                size_mb = pathlib.Path(save_path).stat().st_size / 1e6
                print(f"Saved: {save_path}  ({size_mb:.1f} MB, "
                      f"{len(frames_out)} frames, {colors} colors, {dpi} dpi)")
            except ImportError:
                print("(Pillow not available for palette optimization — "
                      "install it with: pip install Pillow)")
                size_mb = pathlib.Path(save_path).stat().st_size / 1e6
                print(f"Saved: {save_path}  ({size_mb:.1f} MB)")
        elif ext in (".mp4", ".mov"):
            writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
            ani.save(save_path, writer=writer, dpi=dpi)
            size_mb = pathlib.Path(save_path).stat().st_size / 1e6
            print(f"Saved: {save_path}  ({size_mb:.1f} MB)")
        else:
            print(f"Unknown extension '{ext}', defaulting to GIF.")
            writer = animation.PillowWriter(fps=fps)
            ani.save(save_path, writer=writer, dpi=dpi)
            size_mb = pathlib.Path(save_path).stat().st_size / 1e6
            print(f"Saved: {save_path}  ({size_mb:.1f} MB)")
    else:
        plt.show()


if __name__ == "__main__":
    main()
