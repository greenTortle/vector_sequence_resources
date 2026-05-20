import matplotlib.pyplot as plt
from fractions import Fraction
from collections import Counter
import csv


# ============================================================
# User-adjustable defaults
# ============================================================

DEFAULT_START_K = "2.0"
DEFAULT_END_K = "2.5"
DEFAULT_NUM_K_POINTS = 500

DEFAULT_START_X0 = 2
DEFAULT_END_X0 = 1000

DEFAULT_CHECK_N = 1000


# ============================================================
# Exact rational k handling
# ============================================================

def decimal_string_to_fraction(s):
    """
    Converts a decimal string into an exact Fraction.

    Examples:
    "2.5"       -> Fraction(5, 2)
    "2.6667"    -> Fraction(26667, 10000)
    "2.665781"  -> exact rational version of that decimal
    """
    return Fraction(s.strip())


def fraction_to_decimal_string(frac, digits=20):
    """
    Converts a Fraction to a decimal string for saving/display.
    Computation remains exact.
    """
    value = frac.numerator / frac.denominator
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def exact_floor_div_by_k(x, k_frac):
    """
    Computes floor(x / k) exactly.

    If k = numerator / denominator, then:

        floor(x / k) = floor(x * denominator / numerator)
    """
    return (x * k_frac.denominator) // k_frac.numerator


def exact_linspace_fraction(start_k, end_k, num_points):
    """
    Creates exactly equally spaced Fraction k-values.
    """
    if num_points <= 0:
        raise ValueError("Number of k points must be positive.")

    if num_points == 1:
        return [start_k]

    step = (end_k - start_k) / (num_points - 1)
    return [start_k + i * step for i in range(num_points)]


# ============================================================
# Recurrence logic
# ============================================================

def hailstone_S(m):
    """
    One Hailstone / Collatz step.
    """
    if m == 0:
        return 0
    return m // 2 if m % 2 == 0 else 3 * m + 1


def next_vector(v, k_frac):
    """
    Applies the vector recurrence.
    """
    decayed = [exact_floor_div_by_k(x, k_frac) for x in v]
    return tuple(decayed + [hailstone_S(sum(decayed))])


def tail_signature(v):
    """
    Removes leading zeros and returns the nonzero tail.
    """
    i = 0
    while i < len(v) and v[i] == 0:
        i += 1
    return v[i:]


def classify_sequence(x0, k_frac, check_n=DEFAULT_CHECK_N):
    """
    zeroed:
        The nonzero tail becomes empty.

    cycled:
        The same nonzero tail appears again.

    grown:
        The sequence has not zeroed or cycled by v_check_n.
    """
    v = (x0,)
    seen_tails = {tail_signature(v): 0}

    for n in range(1, check_n + 1):
        v = next_vector(v, k_frac)
        tail = tail_signature(v)

        if tail == ():
            return "zeroed"

        if tail in seen_tails:
            return "cycled"

        seen_tails[tail] = n

    return "grown"


def compute_Z_G_for_k(k_frac, start_x0, end_x0, check_n):
    """
    Computes Z, G, and C for one fixed k.
    """
    counts = Counter()

    for x0 in range(start_x0, end_x0 + 1):
        status = classify_sequence(x0, k_frac, check_n=check_n)
        counts[status] += 1

    return counts["zeroed"], counts["grown"], counts["cycled"]


# ============================================================
# Graphing logic
# ============================================================

def make_k_vs_Z_G_graph(
    start_k_str=DEFAULT_START_K,
    end_k_str=DEFAULT_END_K,
    num_k_points=DEFAULT_NUM_K_POINTS,
    start_x0=DEFAULT_START_X0,
    end_x0=DEFAULT_END_X0,
    check_n=DEFAULT_CHECK_N,
    save_csv=True,
    save_plot=True
):
    """
    Computes and plots k vs Z and G.

    Z = number of zeroed sequences.
    G = number of grown sequences.
    C = number of cycled sequences.

    Only the graph is shown. Individual k results are not printed.
    """

    start_k = decimal_string_to_fraction(start_k_str)
    end_k = decimal_string_to_fraction(end_k_str)

    k_values = exact_linspace_fraction(start_k, end_k, num_k_points)

    Z_values = []
    G_values = []
    C_values = []

    for k_frac in k_values:
        Z, G, C = compute_Z_G_for_k(
            k_frac=k_frac,
            start_x0=start_x0,
            end_x0=end_x0,
            check_n=check_n
        )

        Z_values.append(Z)
        G_values.append(G)
        C_values.append(C)

    safe_start = start_k_str.replace(".", "p")
    safe_end = end_k_str.replace(".", "p")

    if save_csv:
        csv_name = f"k_Z_G_{safe_start}_to_{safe_end}_{num_k_points}_points_exact.csv"

        with open(csv_name, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["k_decimal", "k_fraction", "zeroed", "grown", "cycled"])

            for k_frac, Z, G, C in zip(k_values, Z_values, G_values, C_values):
                writer.writerow([
                    fraction_to_decimal_string(k_frac, digits=40),
                    f"{k_frac.numerator}/{k_frac.denominator}",
                    Z,
                    G,
                    C
                ])

    k_plot_values = [float(k) for k in k_values]

    plt.figure(figsize=(12, 6))

    plt.plot(
        k_plot_values,
        Z_values,
        label="Z: zeroed",
        color="blue",
        linewidth=1.5
    )

    plt.plot(
        k_plot_values,
        G_values,
        label="G: grown",
        color="orange",
        linewidth=1.5
    )

    plt.title(
        f"k vs. Z and G "
        f"(x₀ = {start_x0}..{end_x0}, checked to v_{check_n})"
    )

    plt.xlabel("k")
    plt.ylabel("Number of sequences")
    plt.grid(True, alpha=0.35)
    plt.legend()
    plt.tight_layout()

    if save_plot:
        plot_name = f"k_vs_Z_G_{safe_start}_to_{safe_end}_{num_k_points}_points_exact.png"
        plt.savefig(plot_name, dpi=300)

    plt.show()

    return k_values, Z_values, G_values, C_values


# ============================================================
# Main program
# ============================================================

if __name__ == "__main__":
    print("k vs Z and G graph generator")
    print("All recurrence calculations use exact rational arithmetic.")
    print()

    start_k_str = input(f"Starting k [{DEFAULT_START_K}]: ").strip() or DEFAULT_START_K
    end_k_str = input(f"Ending k [{DEFAULT_END_K}]: ").strip() or DEFAULT_END_K

    num_k_points = int(
        input(f"Number of equally spaced k points [{DEFAULT_NUM_K_POINTS}]: ")
        or DEFAULT_NUM_K_POINTS
    )

    start_x0 = int(input(f"Starting x0 [{DEFAULT_START_X0}]: ") or DEFAULT_START_X0)
    end_x0 = int(input(f"Ending x0 [{DEFAULT_END_X0}]: ") or DEFAULT_END_X0)

    check_n = int(input(f"Check through v_n [{DEFAULT_CHECK_N}]: ") or DEFAULT_CHECK_N)

    make_k_vs_Z_G_graph(
        start_k_str=start_k_str,
        end_k_str=end_k_str,
        num_k_points=num_k_points,
        start_x0=start_x0,
        end_x0=end_x0,
        check_n=check_n,
        save_csv=True,
        save_plot=True
    )