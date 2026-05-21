"""
Vector Sequence Multi-Graph Visualizer
==============================================
Visualizes 1 to 4 growing-vector Collatz-type recurrences at once:

    v_{n+1} = (floor(x_0/k), ..., floor(x_n/k), C(sum of floors))

where S(m) is selected by the user: Hailstone, Aliquot, Euler Totient, or Sum of Divisors.

This version lets you specify:
    - the number of graphs to show: 1, 2, 3, or 4
    - a separate k value for each graph
    - a separate starting value x0 for each graph, including negative x0
    - the sequence rule S(.) used for the appended term

Requirements:
    pip install matplotlib

Run:
    python collatz_vector_multigraph_simulator.py
"""

import math
from decimal import Decimal, InvalidOperation, ROUND_FLOOR, getcontext
import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.ticker as ticker


# Keep extra precision for floor(x / k), especially when k has many digits.
getcontext().prec = 60


# ── Sequence/vector logic ─────────────────────────────────────────────────────

SEQUENCE_RULES = ("Hailstone", "Aliquot", "Euler Totient", "Sum of Divisors")


def hailstone_S(m: int) -> int:
    """Hailstone map extended directly to any integer m."""
    return m // 2 if m % 2 == 0 else 3 * m + 1


def _factorization(n: int) -> dict[int, int]:
    """Return prime factorization of positive n as {prime: exponent}."""
    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def sum_of_divisors_S(m: int) -> int:
    """sigma(m): sum of all positive divisors. Returns 0 for m <= 0."""
    if m <= 0:
        return 0
    total = 1
    for p, a in _factorization(m).items():
        total *= (p ** (a + 1) - 1) // (p - 1)
    return total


def aliquot_S(m: int) -> int:
    """s(m): sum of proper positive divisors. Returns 0 for m <= 1."""
    if m <= 1:
        return 0
    return sum_of_divisors_S(m) - m


def euler_totient_S(m: int) -> int:
    """phi(m): Euler totient. Returns 0 for m <= 0."""
    if m <= 0:
        return 0
    result = m
    for p in _factorization(m):
        result -= result // p
    return result


def apply_sequence_rule(m: int, rule_name: str) -> int:
    if rule_name == "Hailstone":
        return hailstone_S(m)
    if rule_name == "Aliquot":
        return aliquot_S(m)
    if rule_name == "Euler Totient":
        return euler_totient_S(m)
    if rule_name == "Sum of Divisors":
        return sum_of_divisors_S(m)
    raise ValueError(f"Unknown sequence rule: {rule_name}")


def floor_divide_by_k(x: int, k: Decimal) -> int:
    """Return floor(x / k) using Decimal arithmetic."""
    return int((Decimal(x) / k).to_integral_value(rounding=ROUND_FLOOR))


def build_sequence(x0: int, k: Decimal, rule_name: str, max_steps: int = 100) -> list[dict]:
    """Return list of dicts {vec, Hn} for steps 0..max_steps."""
    steps = [{"vec": [x0], "Hn": None}]
    vec = [x0]

    for _ in range(max_steps):
        reduced = [floor_divide_by_k(x, k) for x in vec]
        Hn = sum(reduced)
        new_vec = reduced + [apply_sequence_rule(Hn, rule_name)]
        steps.append({"vec": new_vec, "Hn": Hn})
        vec = new_vec

        if all(v == 0 for v in vec):
            break

    return steps


def parse_decimal(text: str, field_name: str) -> Decimal:
    try:
        value = Decimal(text.strip())
    except (InvalidOperation, ValueError):
        raise ValueError(f"{field_name} must be a valid number.")

    if not value.is_finite():
        raise ValueError(f"{field_name} must be finite.")
    if value == 0:
        raise ValueError(f"{field_name} cannot be 0.")

    return value


# ── App ───────────────────────────────────────────────────────────────────────

class VectorSequenceApp:
    MAX_STEPS = 100
    PLAY_INTERVAL_MS = 1200
    MAX_GRAPHS = 4

    DOT_COLOR = "#378ADD"
    NEW_COLOR = "#E24B4A"
    TRAIL_BASE = (55 / 255, 138 / 255, 221 / 255)

    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Vector Sequence Multi-Graph Visualizer")
        self.root.resizable(True, True)

        self.current_n = 0
        self.playing = False
        self._after_id = None

        self.graph_count_var = tk.IntVar(value=1)
        self.rule_var = tk.StringVar(value="Hailstone")
        self.graph_rows = []
        self.configs = []
        self.sequences = []

        self._build_ui()
        self._recompute()

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        outer = ttk.Frame(self.root, padding=8)
        outer.pack(fill=tk.X)

        top = ttk.Frame(outer)
        top.pack(fill=tk.X)

        ttk.Label(top, text="Graphs:").grid(row=0, column=0, sticky=tk.W)
        graph_count = ttk.Spinbox(
            top,
            from_=1,
            to=self.MAX_GRAPHS,
            width=4,
            textvariable=self.graph_count_var,
            command=self._on_graph_count_change,
        )
        graph_count.grid(row=0, column=1, padx=(4, 14), sticky=tk.W)
        graph_count.bind("<Return>", self._on_graph_count_change)
        graph_count.bind("<FocusOut>", self._on_graph_count_change)

        ttk.Label(top, text="S rule:").grid(row=0, column=2, sticky=tk.W)
        rule_combo = ttk.Combobox(
            top,
            width=16,
            textvariable=self.rule_var,
            values=SEQUENCE_RULES,
            state="readonly",
        )
        rule_combo.grid(row=0, column=3, padx=(4, 14), sticky=tk.W)
        rule_combo.bind("<<ComboboxSelected>>", self._on_apply)

        ttk.Label(top, text="Step n:").grid(row=0, column=4, sticky=tk.W)
        self.n_var = tk.IntVar(value=0)
        self.n_slider = ttk.Scale(
            top,
            from_=0,
            to=self.MAX_STEPS,
            orient=tk.HORIZONTAL,
            variable=self.n_var,
            length=260,
            command=self._on_n_slider,
        )
        self.n_slider.grid(row=0, column=5, padx=6)
        self.n_label = ttk.Label(top, text="0", width=5)
        self.n_label.grid(row=0, column=6, sticky=tk.W)

        self.play_btn = ttk.Button(top, text="▶ Play", command=self._toggle_play)
        self.play_btn.grid(row=0, column=7, padx=4)
        self.reset_btn = ttk.Button(top, text="↺ Reset", command=self._reset)
        self.reset_btn.grid(row=0, column=8, padx=4)
        self.apply_btn = ttk.Button(top, text="Apply / Recompute", command=self._on_apply)
        self.apply_btn.grid(row=0, column=9, padx=(12, 4))

        entry_frame = ttk.LabelFrame(outer, text="Graph settings", padding=6)
        entry_frame.pack(fill=tk.X, pady=(8, 0))

        ttk.Label(entry_frame, text="Graph").grid(row=0, column=0, padx=4, sticky=tk.W)
        ttk.Label(entry_frame, text="k value").grid(row=0, column=1, padx=4, sticky=tk.W)
        ttk.Label(entry_frame, text="x₀ value").grid(row=0, column=2, padx=4, sticky=tk.W)

        defaults = [("2", "7"), ("2.5", "7"), ("2.666", "7"), ("2.718281828", "7")]
        for idx in range(self.MAX_GRAPHS):
            k_default, x_default = defaults[idx]
            row = idx + 1
            label = ttk.Label(entry_frame, text=f"{idx + 1}")
            label.grid(row=row, column=0, padx=4, pady=2, sticky=tk.W)

            k_var = tk.StringVar(value=k_default)
            k_entry = ttk.Entry(entry_frame, width=18, textvariable=k_var)
            k_entry.grid(row=row, column=1, padx=4, pady=2, sticky=tk.W)
            k_entry.bind("<Return>", self._on_apply)

            x_var = tk.StringVar(value=x_default)
            x_entry = ttk.Entry(entry_frame, width=14, textvariable=x_var)
            x_entry.grid(row=row, column=2, padx=4, pady=2, sticky=tk.W)
            x_entry.bind("<Return>", self._on_apply)

            self.graph_rows.append({
                "label": label,
                "k_var": k_var,
                "k_entry": k_entry,
                "x_var": x_var,
                "x_entry": x_entry,
            })

        self._update_graph_row_visibility()

        self.fig = plt.Figure(figsize=(10, 7))
        self.axes = []
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self.status_var = tk.StringVar(value="")
        ttk.Label(
            self.root,
            textvariable=self.status_var,
            font=("Courier", 10),
            anchor=tk.W,
            padding=(8, 2),
        ).pack(fill=tk.X)

        self.vec_var = tk.StringVar(value="")
        ttk.Label(
            self.root,
            textvariable=self.vec_var,
            font=("Courier", 9),
            anchor=tk.W,
            padding=(8, 0, 8, 4),
        ).pack(fill=tk.X)

    def _update_graph_row_visibility(self):
        count = self._get_graph_count_safely()
        for idx, row in enumerate(self.graph_rows):
            state = "normal" if idx < count else "disabled"
            row["k_entry"].configure(state=state)
            row["x_entry"].configure(state=state)

    # ── parsing/recompute ─────────────────────────────────────────────────────

    def _get_graph_count_safely(self) -> int:
        try:
            count = int(self.graph_count_var.get())
        except (tk.TclError, ValueError):
            count = 1
        count = max(1, min(self.MAX_GRAPHS, count))
        self.graph_count_var.set(count)
        return count

    def _read_configs(self):
        count = self._get_graph_count_safely()
        rule_name = self.rule_var.get()
        if rule_name not in SEQUENCE_RULES:
            raise ValueError("Select a valid S rule.")

        configs = []
        for idx in range(count):
            row = self.graph_rows[idx]
            try:
                k = parse_decimal(row["k_var"].get(), f"Graph {idx + 1} k")
                x0 = int(row["x_var"].get().strip())
            except ValueError as exc:
                raise ValueError(str(exc))
            configs.append({"k": k, "x0": x0, "rule_name": rule_name})
        return configs

    def _recompute(self):
        try:
            configs = self._read_configs()
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return

        self.configs = configs
        self.sequences = [build_sequence(c["x0"], c["k"], c["rule_name"], self.MAX_STEPS) for c in configs]

        max_n = max(len(seq) - 1 for seq in self.sequences)
        self.n_slider.configure(to=max_n)
        self.current_n = min(self.current_n, max_n)
        self.n_var.set(self.current_n)
        self.n_label.configure(text=str(self.current_n))
        self._draw()

    # ── drawing ───────────────────────────────────────────────────────────────

    @staticmethod
    def _fmt(v: float) -> str:
        av = abs(v)
        sign = "-" if v < 0 else ""
        if av >= 1e9:
            return f"{sign}{av/1e9:.3g}B"
        if av >= 1e6:
            return f"{sign}{av/1e6:.3g}M"
        if av >= 1e3:
            return f"{sign}{av/1e3:.3g}k"
        return str(int(v))

    @staticmethod
    def _active_step(sequence: list[dict], n: int) -> tuple[int, dict]:
        actual_n = min(n, len(sequence) - 1)
        return actual_n, sequence[actual_n]

    @staticmethod
    def _max_abs_up_to(sequence: list[dict], n: int) -> int:
        m = 1
        stop = min(n, len(sequence) - 1)
        for step_idx in range(stop + 1):
            for v in sequence[step_idx]["vec"]:
                m = max(m, abs(v))
        return m

    @staticmethod
    def _symlog_limit(max_abs: int) -> int:
        if max_abs <= 1:
            return 10
        return 10 ** math.ceil(math.log10(max_abs))

    def _make_axes(self, count: int):
        self.fig.clear()
        if count == 1:
            axes = [self.fig.add_subplot(1, 1, 1)]
        elif count == 2:
            axes = [self.fig.add_subplot(1, 2, i + 1) for i in range(2)]
        else:
            axes = [self.fig.add_subplot(2, 2, i + 1) for i in range(count)]
        self.axes = axes
        return axes

    def _draw(self):
        if not self.sequences:
            return

        count = len(self.sequences)
        axes = self._make_axes(count)
        n = self.current_n
        status_parts = []
        vec_lines = []

        for graph_idx, (ax, sequence, config) in enumerate(zip(axes, self.sequences, self.configs), start=1):
            actual_n, step = self._active_step(sequence, n)
            vec = step["vec"]
            Hn = step["Hn"]

            max_abs = self._max_abs_up_to(sequence, actual_n)
            y_limit = self._symlog_limit(max_abs)

            ax.set_yscale("symlog", linthresh=1)
            ax.set_ylim(-y_limit, y_limit)

            # Keep every subplot on the same vector-index scale.
            # Without these two lines, matplotlib autoscaling/compression can make
            # later subplots show -0.5..2.5 when the current vector is short.
            ax.set_xlim(0, self.MAX_STEPS)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

            ax.set_xlabel("Vector index i")
            ax.set_ylabel("Value (signed symlog)")
            ax.set_title(
                f"Graph {graph_idx}: S={config['rule_name']}, k={config['k']}, x₀={config['x0']}, n={actual_n}",
                fontsize=10,
            )
            ax.grid(True, which="major", linestyle="-", linewidth=0.4, alpha=0.3)
            ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.15)
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: self._fmt(v)))

            for step_idx in range(actual_n):
                age = actual_n - step_idx
                alpha = max(0.06, 0.22 * (0.82 ** (age - 1)))
                sv = sequence[step_idx]["vec"]
                xs = list(range(len(sv)))
                ys = sv
                if xs:
                    r, g, b = self.TRAIL_BASE
                    ax.scatter(xs, ys, s=22, color=(r, g, b, alpha), zorder=2, linewidths=0)

            xs = list(range(len(vec)))
            ys = vec
            if len(xs) > 1:
                ax.plot(xs, ys, color=self.DOT_COLOR, linewidth=1.2, alpha=0.45, zorder=3)

            ax.scatter(xs, ys, s=40, color=self.DOT_COLOR, edgecolors="white", linewidths=1.2, zorder=4)
            ax.axhline(0, linewidth=0.8, alpha=0.35)

            if actual_n > 0 and vec:
                last_i = len(vec) - 1
                last_v = vec[last_i]
                ax.scatter(
                    [last_i],
                    [last_v],
                    s=120,
                    facecolors="none",
                    edgecolors=self.NEW_COLOR,
                    linewidths=2,
                    zorder=5,
                )

            cur_max = max(vec) if vec else 0
            cur_min = min(vec) if vec else 0
            Hn_str = str(Hn) if Hn is not None else "—"
            stopped = " stopped" if actual_n < n else ""
            status_parts.append(
                f"G{graph_idx}: n={actual_n}{stopped}, len={len(vec)}, Hₙ={Hn_str}, min={cur_min:,}, max={cur_max:,}"
            )
            display = ", ".join(map(str, vec[:16])) + (" …" if len(vec) > 16 else "")
            vec_lines.append(f"G{graph_idx} v{actual_n} = [{display}]")

        self.fig.tight_layout(pad=2.0)
        self.canvas.draw()
        self.status_var.set("  |  ".join(status_parts))
        self.vec_var.set("    ".join(vec_lines))

    # ── controls ──────────────────────────────────────────────────────────────

    def _on_graph_count_change(self, _=None):
        self._stop_play()
        self._update_graph_row_visibility()
        self.current_n = 0
        self._recompute()

    def _on_apply(self, _=None):
        self._stop_play()
        self.current_n = 0
        self._recompute()

    def _on_n_slider(self, _=None):
        self._stop_play()
        self.current_n = int(self.n_var.get())
        self.n_label.configure(text=str(self.current_n))
        self._draw()

    def _toggle_play(self):
        if self.playing:
            self._stop_play()
        else:
            max_n = max(len(seq) - 1 for seq in self.sequences)
            if self.current_n >= max_n:
                self.current_n = 0
            self.playing = True
            self.play_btn.configure(text="⏸ Pause")
            self._step_play()

    def _step_play(self):
        if not self.playing:
            return

        self.current_n += 1
        self.n_var.set(self.current_n)
        self.n_label.configure(text=str(self.current_n))
        self._draw()

        max_n = max(len(seq) - 1 for seq in self.sequences)
        if self.current_n >= max_n:
            self._stop_play()
        else:
            self._after_id = self.root.after(self.PLAY_INTERVAL_MS, self._step_play)

    def _stop_play(self):
        self.playing = False
        self.play_btn.configure(text="▶ Play")
        if self._after_id is not None:
            self.root.after_cancel(self._after_id)
            self._after_id = None

    def _reset(self):
        self._stop_play()
        self.current_n = 0
        self.n_var.set(0)
        self.n_label.configure(text="0")
        self._draw()


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    root = tk.Tk()
    app = VectorSequenceApp(root)
    root.mainloop()
