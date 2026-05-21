import tkinter as tk
from tkinter import ttk, messagebox
from math import floor
from collections import Counter, defaultdict


DEFAULT_CHECK_N = 1000
DEFAULT_PRINT_N = 50


SEQUENCE_RULES = ("Hailstone", "Aliquot", "Euler Totient", "Sum of Divisors")


def hailstone_S(m):
    if m == 0:
        return 0
    return m // 2 if m % 2 == 0 else 3 * m + 1


def _factorization(n):
    """Return prime factorization of positive n as {prime: exponent}."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def sum_of_divisors_S(m):
    """sigma(m): sum of all positive divisors. Returns 0 for m <= 0."""
    if m <= 0:
        return 0
    total = 1
    for p, a in _factorization(m).items():
        total *= (p ** (a + 1) - 1) // (p - 1)
    return total


def aliquot_S(m):
    """s(m): sum of proper positive divisors. Returns 0 for m <= 1."""
    if m <= 1:
        return 0
    return sum_of_divisors_S(m) - m


def euler_totient_S(m):
    """phi(m): Euler totient. Returns 0 for m <= 0."""
    if m <= 0:
        return 0
    result = m
    for p in _factorization(m):
        result -= result // p
    return result


def apply_sequence_rule(m, rule_name):
    if rule_name == "Hailstone":
        return hailstone_S(m)
    if rule_name == "Aliquot":
        return aliquot_S(m)
    if rule_name == "Euler Totient":
        return euler_totient_S(m)
    if rule_name == "Sum of Divisors":
        return sum_of_divisors_S(m)
    raise ValueError(f"Unknown sequence rule: {rule_name}")


def next_vector(v, k, rule_name="Hailstone"):
    decayed = [floor(x / k) for x in v]
    return tuple(decayed + [apply_sequence_rule(sum(decayed), rule_name)])


def tail_signature(v):
    i = 0
    while i < len(v) and v[i] == 0:
        i += 1
    return v[i:]


def compress_leading_zeros(v):
    tail = tail_signature(v)

    if tail == ():
        return "(0...)"

    if len(tail) == len(v):
        return str(v)

    return f"(0..., {', '.join(map(str, tail))})"


def classify_sequence(x0, k, check_n=DEFAULT_CHECK_N, print_n=DEFAULT_PRINT_N, rule_name="Hailstone"):
    v = (x0,)
    printed_seq = [v]
    full_seq = [v]

    seen_tails = {tail_signature(v): 0}

    max_value_seen = x0
    max_sum_seen = x0

    for n in range(1, check_n + 1):
        v = next_vector(v, k, rule_name)
        full_seq.append(v)

        max_value_seen = max(max_value_seen, max(v))
        max_sum_seen = max(max_sum_seen, sum(v))

        if n <= print_n:
            printed_seq.append(v)

        tail = tail_signature(v)

        if tail == ():
            if n <= print_n:
                printed_seq = printed_seq[:-1]

            return {
                "x0": x0,
                "status": "zeroed",
                "stop_n": n,
                "printed_sequence": printed_seq,
                "event_message": f"v_{n} zeroes",
                "max_value_seen": max_value_seen,
                "max_sum_seen": max_sum_seen,
                "cycle": None,
            }

        if tail in seen_tails:
            cycle_start = seen_tails[tail]
            cycle = [tail_signature(vec) for vec in full_seq[cycle_start:n]]

            if n <= print_n:
                printed_seq = printed_seq[:-1]

            return {
                "x0": x0,
                "status": "cycled",
                "stop_n": n,
                "printed_sequence": printed_seq,
                "event_message": (
                    f"v_{n} has the same nonzero tail as v_{cycle_start}: {tail}"
                ),
                "cycle_start": cycle_start,
                "cycle_length": n - cycle_start,
                "cycle": cycle,
                "max_value_seen": max_value_seen,
                "max_sum_seen": max_sum_seen,
            }

        seen_tails[tail] = n

    return {
        "x0": x0,
        "status": "grown",
        "stop_n": check_n,
        "printed_sequence": printed_seq,
        "event_message": f"No zero/cycle found by v_{check_n}; classified as grown",
        "max_value_seen": max_value_seen,
        "max_sum_seen": max_sum_seen,
        "cycle": None,
    }


class VectorRecurrenceGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Vector Recurrence S-Rule Explorer")
        self.root.geometry("1150x750")

        self.results = {}
        self.cycle_groups = defaultdict(list)

        self.build_controls()
        self.build_summary()
        self.build_grid_area()

    def build_controls(self):
        frame = ttk.Frame(self.root, padding=10)
        frame.pack(fill="x")

        ttk.Label(frame, text="k:").grid(row=0, column=0, padx=5, sticky="e")
        self.k_entry = ttk.Entry(frame, width=10)
        self.k_entry.insert(0, "2.5")
        self.k_entry.grid(row=0, column=1, padx=5)

        ttk.Label(frame, text="S rule:").grid(row=0, column=2, padx=5, sticky="e")
        self.rule_var = tk.StringVar(value="Hailstone")
        self.rule_combo = ttk.Combobox(
            frame,
            width=16,
            textvariable=self.rule_var,
            values=SEQUENCE_RULES,
            state="readonly",
        )
        self.rule_combo.grid(row=0, column=3, padx=5)

        ttk.Label(frame, text="start x₀:").grid(row=0, column=4, padx=5, sticky="e")
        self.start_entry = ttk.Entry(frame, width=10)
        self.start_entry.insert(0, "2")
        self.start_entry.grid(row=0, column=5, padx=5)

        ttk.Label(frame, text="end x₀:").grid(row=0, column=6, padx=5, sticky="e")
        self.end_entry = ttk.Entry(frame, width=10)
        self.end_entry.insert(0, "1000")
        self.end_entry.grid(row=0, column=7, padx=5)

        ttk.Label(frame, text="check to vₙ:").grid(row=0, column=8, padx=5, sticky="e")
        self.check_entry = ttk.Entry(frame, width=10)
        self.check_entry.insert(0, str(DEFAULT_CHECK_N))
        self.check_entry.grid(row=0, column=9, padx=5)

        ttk.Label(frame, text="print to vₙ:").grid(row=0, column=10, padx=5, sticky="e")
        self.print_entry = ttk.Entry(frame, width=10)
        self.print_entry.insert(0, str(DEFAULT_PRINT_N))
        self.print_entry.grid(row=0, column=11, padx=5)

        run_button = ttk.Button(frame, text="Run", command=self.run_experiment)
        run_button.grid(row=0, column=12, padx=10)

        cycles_button = ttk.Button(frame, text="Show Cycles", command=self.show_cycles)
        cycles_button.grid(row=0, column=13, padx=10)

    def build_summary(self):
        frame = ttk.Frame(self.root, padding=10)
        frame.pack(fill="x")

        self.summary_label = ttk.Label(
            frame,
            text="Run an experiment to see grown / zeroed / cycled tallies.",
            font=("Arial", 12, "bold")
        )
        self.summary_label.pack(anchor="w")

    def build_grid_area(self):
        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True, padx=10, pady=10)

        self.canvas = tk.Canvas(outer)
        scrollbar = ttk.Scrollbar(outer, orient="vertical", command=self.canvas.yview)

        self.grid_frame = ttk.Frame(self.canvas)

        self.grid_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        self.canvas.create_window((0, 0), window=self.grid_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def add_cycle_to_groups(self, new_cycle, x0):
        """
        Groups cycles that are the same up to rotation/overlap.

        A new cycle is considered the same as a previous cycle if the first
        tail of the new cycle appears anywhere inside a previous cycle.
        """
        if not new_cycle:
            return

        first_line = new_cycle[0]

        for existing_cycle in self.cycle_groups:
            if first_line in existing_cycle:
                self.cycle_groups[existing_cycle].append(x0)
                return

        self.cycle_groups[tuple(new_cycle)].append(x0)

    def run_experiment(self):
        try:
            k = float(self.k_entry.get())
            start = int(self.start_entry.get())
            end = int(self.end_entry.get())
            check_n = int(self.check_entry.get())
            print_n = int(self.print_entry.get())
            rule_name = self.rule_var.get()
            if rule_name not in SEQUENCE_RULES:
                raise ValueError("Select a valid S rule.")

            if k <= 0:
                raise ValueError("k must be positive.")
            if start > end:
                raise ValueError("start x₀ must be less than or equal to end x₀.")
            if check_n <= 0:
                raise ValueError("check to vₙ must be positive.")
            if print_n <= 0:
                raise ValueError("print to vₙ must be positive.")
            if print_n > check_n:
                raise ValueError("print to vₙ cannot be greater than check to vₙ.")

        except ValueError as e:
            messagebox.showerror("Invalid input", str(e))
            return

        for widget in self.grid_frame.winfo_children():
            widget.destroy()

        self.results.clear()
        self.cycle_groups.clear()

        for x0 in range(start, end + 1):
            result = classify_sequence(
                x0,
                k,
                check_n=check_n,
                print_n=print_n,
                rule_name=rule_name
            )

            self.results[x0] = result

            if result["status"] == "cycled":
    	        self.add_cycle_to_groups(result["cycle"], x0)

        counts = Counter(r["status"] for r in self.results.values())

        self.summary_label.config(
            text=(
                f"S = {rule_name} | "
                f"k = {k} | "
                f"x₀ = {start} to {end} | "
                f"checked to v_{check_n} | "
                f"printed to v_{print_n} | "
                f"grown: {counts['grown']} | "
                f"zeroed: {counts['zeroed']} | "
                f"cycled: {counts['cycled']}"
            )
        )

        self.populate_x0_grid()

    def populate_x0_grid(self):
        columns = 20

        for idx, x0 in enumerate(self.results):
            result = self.results[x0]
            status = result["status"]

            label = f"{x0}\n{status}"

            button = tk.Button(
                self.grid_frame,
                text=label,
                width=9,
                height=3,
                command=lambda value=x0: self.show_sequence(value)
            )

            if status == "grown":
                button.configure(bg="#ffd6d6")
            elif status == "zeroed":
                button.configure(bg="#d6ffd6")
            elif status == "cycled":
                button.configure(bg="#d6e4ff")

            row = idx // columns
            col = idx % columns
            button.grid(row=row, column=col, padx=3, pady=3)

    def show_sequence(self, x0):
        result = self.results[x0]

        win = tk.Toplevel(self.root)
        win.title(f"Sequence for x₀ = {x0}")
        win.geometry("950x600")

        text = tk.Text(win, wrap="none")
        text.pack(fill="both", expand=True)

        text.insert("end", f"x₀ = {x0}\n")
        text.insert("end", f"status = {result['status']}\n")
        text.insert("end", f"stop_n = {result['stop_n']}\n")
        text.insert("end", f"max_value_seen = {result['max_value_seen']}\n")
        text.insert("end", f"max_sum_seen = {result['max_sum_seen']}\n\n")

        for i, v in enumerate(result["printed_sequence"]):
            text.insert("end", f"v_{i} = {compress_leading_zeros(v)}\n")

        text.insert("end", f"\n{result['event_message']}\n")

        text.config(state="disabled")

    def show_cycles(self):
        win = tk.Toplevel(self.root)
        win.title("Cycles")
        win.geometry("1000x650")

        text = tk.Text(win, wrap="none")
        text.pack(fill="both", expand=True)

        if not self.cycle_groups:
            text.insert("end", "No cycles found.\n")
            text.config(state="disabled")
            return

        text.insert("end", f"Number of distinct cycles: {len(self.cycle_groups)}\n\n")

        for idx, (cycle, x0_values) in enumerate(self.cycle_groups.items(), start=1):
            text.insert("end", "=" * 80 + "\n")
            text.insert("end", f"Cycle {idx}\n")
            text.insert("end", f"cycle length: {len(cycle)}\n")
            text.insert("end", f"number of x₀ values ending in this cycle: {len(x0_values)}\n")
            text.insert("end", f"x₀ values: {x0_values}\n\n")

            for j, tail in enumerate(cycle):
                text.insert("end", f"cycle[{j}] tail = {tail}\n")

            text.insert("end", "\n")

        text.config(state="disabled")


if __name__ == "__main__":
    root = tk.Tk()
    app = VectorRecurrenceGUI(root)
    root.mainloop()