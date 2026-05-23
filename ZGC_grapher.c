/* 
 * k vs Z/G/C Scanner — pure C, single-file, OpenMP parallelised
 * =============================================================
 *
 * Scans k-values for the vector recurrence:
 *
 * v_{n+1} = (floor(x_0/k), ..., floor(x_n/k), S(sum of floors))
 *
 * where S(m) is selected at runtime: hailstone, aliquot, totient, or sigma.
 *
 * For each k, x0 values are classified as:
 * Z (zeroed) : the nonzero tail becomes empty
 * G (grown)  : neither zeroed nor cycled within check_n steps
 * C (cycled) : the same nonzero tail appears again
 *
 * Build (Mac / Linux)  — parallel (recommended):
 * gcc -O3 -march=native -fopenmp -o k_zgc_scanner k_zgc_scanner.c -lm
 *
 * Build without OpenMP (single-threaded fallback):
 * gcc -O3 -march=native -o ZGC_grapher ZGC_grapher.c -lm
 *
 * Build (Windows, MinGW / MSYS2):
 * gcc -O3 -fopenmp -o ZGC_grapher.exe ZGC_grapher.c -lm
 *
 * Control thread count (default = all cores):
 * $env:OMP_NUM_THREADS=4
 *
 * RUN:
 * ./ZGC_grapher 2.0 2.5 500 2 100 200 H
 * (start end npts x0s x0e chkn rule [coeff])
 * rule: H [coeff], A, T, S
 *   H [coeff] = generalized hailstone: odd→coeff*x+1, even→x/2
 *               coeff defaults to 3 (classic Collatz) if omitted
 *   A = aliquot, T = Euler totient, S = sum of divisors
 *
 * Examples:
 *   ./ZGC_grapher 2.0 2.5 500 2 100 200 H        (classic 3x+1)
 *   ./ZGC_grapher 2.0 2.5 500 2 100 200 H 5      (5x+1 hailstone)
 *   ./ZGC_grapher 0.5 2.0 500 2 100 200 H 3      (k < 1 allowed)
 *
 * NOTE: For k <= 1 the sequence values can grow, so runs may be slow.
 *       For k ~ 1..2.5 the decay rate is slow so the program slows WAY down.
 *
 * Output:
 * Prints a progress bar + summary table to stdout and writes
 * k_zgc_results.csv in the current directory.
 *
 * PRINT PLOT:
 * python ZGC_grapher.py
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

/* ============================================================
 * Configuration / defaults
 * ============================================================ */
#define DEFAULT_START_K_NUM 2LL
#define DEFAULT_START_K_DEN 1LL
#define DEFAULT_END_K_NUM 5LL
#define DEFAULT_END_K_DEN 2LL
#define DEFAULT_NUM_K_POINTS 500
#define DEFAULT_START_X0 2
#define DEFAULT_END_X0 100
#define DEFAULT_CHECK_N 200
#define DEFAULT_HAILSTONE_COEFF 3LL

#define OUTPUT_CSV "k_zgc_results.csv"

/* Maximum tail length */
#define MAX_TAIL 2048

/* Cache for expensive S() functions */
#define S_CACHE_SIZE 4194304

/* ============================================================
 * Cycle detection
 * ============================================================ */
#define SEEN_SLOTS (1 << 13)
#define SEEN_MASK (SEEN_SLOTS - 1)
#define SEEN_EMPTY 0ULL

typedef struct {
    uint64_t slots[SEEN_SLOTS];
} SeenSet;

static inline void seen_clear(SeenSet *s) {
    memset(s->slots, 0, sizeof(s->slots));
}

static inline uint64_t fnv1a(const int64_t *arr, int len) {
    uint64_t h = 14695981039346656037ULL;
    const uint8_t *p = (const uint8_t *)arr;
    int bytes = len * (int)sizeof(int64_t);
    for (int i = 0; i < bytes; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h ? h : 1ULL;
}

static inline int seen_insert(SeenSet *s, uint64_t h) {
    uint32_t slot = (uint32_t)(h & SEEN_MASK);
    for (int probe = 0; probe < SEEN_SLOTS; probe++) {
        uint64_t v = s->slots[slot];
        if (v == SEEN_EMPTY) { s->slots[slot] = h; return 0; }
        if (v == h) return 1;
        slot = (slot + 1) & SEEN_MASK;
    }
    return 0;
}

/* ============================================================
 * Sequence rules + caching
 * ============================================================ */
typedef enum {
    RULE_HAILSTONE = 0,
    RULE_ALIQUOT,
    RULE_TOTIENT,
    RULE_SIGMA
} SequenceRule;

/* Hailstone coefficient (generalized: odd -> coeff*x + 1, even -> x/2) */
static int64_t g_hailstone_coeff = DEFAULT_HAILSTONE_COEFF;

static const char *rule_display_name(SequenceRule rule) {
    switch (rule) {
        case RULE_HAILSTONE: return "Hailstone";
        case RULE_ALIQUOT:   return "Aliquot";
        case RULE_TOTIENT:   return "Euler Totient";
        case RULE_SIGMA:     return "Sum of Divisors";
        default:             return "Hailstone";
    }
}

static SequenceRule parse_rule(const char *s) {
    if (!s) return RULE_HAILSTONE;
    if (!strcmp(s, "hailstone") || !strcmp(s, "Hailstone") || !strcmp(s, "collatz") || !strcmp(s, "Collatz") || !strcmp(s, "H"))
        return RULE_HAILSTONE;
    if (!strcmp(s, "aliquot") || !strcmp(s, "Aliquot") || !strcmp(s, "A"))
        return RULE_ALIQUOT;
    if (!strcmp(s, "totient") || !strcmp(s, "Totient") || !strcmp(s, "euler") || !strcmp(s, "phi") || !strcmp(s, "T"))
        return RULE_TOTIENT;
    if (!strcmp(s, "sigma") || !strcmp(s, "Sigma") || !strcmp(s, "sumdiv") || !strcmp(s, "sum-of-divisors") || !strcmp(s, "S"))
        return RULE_SIGMA;

    fprintf(stderr, "Warning: unknown rule '%s'; using hailstone.\n", s);
    return RULE_HAILSTONE;
}

static inline int64_t sat_from_i128(__int128 v) {
    if (v > (__int128)INT64_MAX) return INT64_MAX;
    if (v < (__int128)INT64_MIN) return INT64_MIN;
    return (int64_t)v;
}

static int64_t sigma_S(int64_t m) {
    if (m <= 0) return 0;
    int64_t n = m;
    __int128 total = 1;
    for (int64_t p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
        if (n % p != 0) continue;
        __int128 term = 1;
        __int128 powp = 1;
        while (n % p == 0) {
            n /= p;
            powp *= p;
            term += powp;
        }
        total *= term;
        if (total > (__int128)INT64_MAX) return INT64_MAX;
    }
    if (n > 1) {
        total *= ((__int128)n + 1);
        if (total > (__int128)INT64_MAX) return INT64_MAX;
    }
    return (int64_t)total;
}

static int64_t aliquot_S(int64_t m) {
    if (m <= 1) return 0;
    return sat_from_i128((__int128)sigma_S(m) - m);
}

static int64_t totient_S(int64_t m) {
    if (m <= 0) return 0;
    int64_t n = m;
    int64_t result = m;
    for (int64_t p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
        if (n % p != 0) continue;
        while (n % p == 0) n /= p;
        result -= result / p;
    }
    if (n > 1) result -= result / n;
    return result;
}

/* Generalized hailstone: odd -> coeff*x + 1, even -> x/2
 * Uses the global g_hailstone_coeff. */
static int64_t hailstone_S(int64_t m) {
    if (m & 1)
        return sat_from_i128((__int128)g_hailstone_coeff * m + 1);
    else
        return m >> 1;
}

static int64_t sigma_cache[S_CACHE_SIZE];
static int64_t totient_cache[S_CACHE_SIZE];
static uint8_t cache_init = 0;

static void init_s_caches(void) {
    if (cache_init) return;
    memset(sigma_cache, -1, sizeof(sigma_cache));
    memset(totient_cache, -1, sizeof(totient_cache));
    cache_init = 1;
}

static int64_t cached_sigma(int64_t m) {
    if (m <= 0) return 0;
    if (m < S_CACHE_SIZE) {
        if (sigma_cache[m] != -1) return sigma_cache[m];
        return sigma_cache[m] = sigma_S(m);
    }
    return sigma_S(m);
}

static int64_t cached_totient(int64_t m) {
    if (m <= 0) return 0;
    if (m < S_CACHE_SIZE) {
        if (totient_cache[m] != -1) return totient_cache[m];
        return totient_cache[m] = totient_S(m);
    }
    return totient_S(m);
}

static int64_t apply_sequence_rule(int64_t m, SequenceRule rule) {
    switch (rule) {
        case RULE_HAILSTONE: return hailstone_S(m);
        case RULE_ALIQUOT:   return aliquot_S(m);
        case RULE_TOTIENT:   return cached_totient(m);
        case RULE_SIGMA:     return cached_sigma(m);
        default:             return hailstone_S(m);
    }
}

/* ============================================================
 * Core classification
 * ============================================================ */
static int classify_sequence(int64_t k_num, int64_t k_den,
                             int64_t x0, int check_n,
                             SequenceRule rule,
                             int64_t *tail, int64_t *next_tail,
                             SeenSet *seen) {
    tail[0] = x0;
    int tail_len = 1;

    seen_clear(seen);
    seen_insert(seen, fnv1a(tail, tail_len));

    for (int step = 0; step < check_n; step++) {
        int64_t h = 0;
        int nz = 0;
        int next_len = 0;

        for (int i = 0; i < tail_len; i++) {
            /* floor(x / k) = floor(x * k_den / k_num)
             * Works for any positive rational k, including k < 1 and k = 1. */
            int64_t y = (tail[i] * k_den) / k_num;
            h += y;
            if (y || nz) {
                nz = 1;
                if (next_len < MAX_TAIL)
                    next_tail[next_len++] = y;
            }
        }

        if (next_len == 0) return 0;

        h = apply_sequence_rule(h, rule);

        if (next_len < MAX_TAIL)
            next_tail[next_len++] = h;

        memcpy(tail, next_tail, (size_t)next_len * sizeof(int64_t));
        tail_len = next_len;

        if (tail_len >= MAX_TAIL) return 1;

        uint64_t h_tail = fnv1a(tail, tail_len);
        if (seen_insert(seen, h_tail)) return 2;
    }
    return 1;
}

/* ============================================================
 * Rational arithmetic
 * ============================================================ */
static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = b; b = a % b; a = t; }
    return a ? a : 1;
}

typedef struct { int64_t num, den; } Frac;

static Frac frac(int64_t n, int64_t d) {
    int64_t g = gcd64(n < 0 ? -n : n, d);
    Frac f = { n / g, d / g };
    if (f.den < 0) { f.num = -f.num; f.den = -f.den; }
    return f;
}

static Frac linspace_frac(Frac start, Frac end, int i, int n) {
    if (n <= 1 || i == 0) return start;
    if (i == n - 1) return end;

    double s = (double)start.num / start.den;
    double e = (double)end.num / end.den;
    double v = s + (double)i * (e - s) / (double)(n - 1);

    int64_t best_num = (int64_t)round(v);
    int64_t best_den = 1;
    double best_err = fabs(v - (double)best_num);

    for (int64_t d = 2; d <= 20000 && best_err > 1e-14; d++) {
        int64_t nn = (int64_t)round(v * d);
        double err = fabs(v - (double)nn / d);
        if (err < best_err) {
            best_err = err;
            best_num = nn;
            best_den = d;
        }
    }
    return frac(best_num, best_den);
}

static Frac parse_frac(const char *s) {
    const char *slash = strchr(s, '/');
    if (slash) return frac(atoll(s), atoll(slash + 1));

    const char *dot = strchr(s, '.');
    if (!dot) return frac(atoll(s), 1);

    int places = 0;
    for (const char *p = dot + 1; *p >= '0' && *p <= '9'; p++) places++;
    int64_t den = 1;
    for (int i = 0; i < places; i++) den *= 10;

    int64_t int_part = atoll(s);
    int negative = (s[0] == '-');
    char fbuf[32] = {0};
    strncpy(fbuf, dot + 1, (size_t)places < sizeof(fbuf)-1 ? (size_t)places : sizeof(fbuf)-1);
    int64_t frac_part = atoll(fbuf);
    if (negative) frac_part = -frac_part;

    return frac(int_part * den + frac_part, den);
}

/* ============================================================
 * Progress bar
 * ============================================================ */
static void print_progress(int done, int total, double elapsed) {
    int w = 40;
    int filled = (total > 0) ? (int)((double)done / total * w) : 0;
    double rate = (done > 0 && elapsed > 0) ? done / elapsed : 0;
    double eta = (rate > 0 && done < total) ? (total - done) / rate : 0;

    printf("\r [");
    for (int i = 0; i < w; i++) putchar(i < filled ? '=' : ' ');
    printf("] %d/%d %.1fs", done, total, elapsed);
    if (eta > 0.5) printf(" ETA %.0fs", eta);
    printf(" ");
    fflush(stdout);
}

/* ============================================================
 * Main
 * ============================================================ */
int main(int argc, char **argv) {
    Frac start_k = frac(DEFAULT_START_K_NUM, DEFAULT_START_K_DEN);
    Frac end_k = frac(DEFAULT_END_K_NUM, DEFAULT_END_K_DEN);
    int num_points = DEFAULT_NUM_K_POINTS;
    int64_t start_x0 = DEFAULT_START_X0;
    int64_t end_x0 = DEFAULT_END_X0;
    int check_n = DEFAULT_CHECK_N;
    SequenceRule rule = RULE_HAILSTONE;

    if (argc >= 3) { start_k = parse_frac(argv[1]); end_k = parse_frac(argv[2]); }
    if (argc >= 4) num_points = atoi(argv[3]);
    if (argc >= 5) start_x0 = atoll(argv[4]);
    if (argc >= 6) end_x0 = atoll(argv[5]);
    if (argc >= 7) check_n = atoi(argv[6]);
    if (argc >= 8) rule = parse_rule(argv[7]);

    /* Optional hailstone coefficient: argv[8], only meaningful for H rule */
    if (argc >= 9) {
        if (rule == RULE_HAILSTONE) {
            int64_t c = atoll(argv[8]);
            if (c < 1) {
                fprintf(stderr, "Error: hailstone coefficient must be >= 1.\n");
                return 1;
            }
            g_hailstone_coeff = c;
        } else {
            fprintf(stderr, "Warning: coefficient argument ignored for non-hailstone rule.\n");
        }
    }

    /* Validate k: must be positive (k > 0, not necessarily > 1) */
    double sk = (double)start_k.num / start_k.den;
    double ek = (double)end_k.num / end_k.den;
    if (sk <= 0.0 || ek <= 0.0) {
        fprintf(stderr, "Error: k must be > 0.\n"); return 1;
    }
    if (num_points <= 0 || check_n <= 0 || end_x0 < start_x0) {
        fprintf(stderr, "Error: invalid parameters.\n"); return 1;
    }

    int64_t x0_count = end_x0 - start_x0 + 1;

    /* Build display string for the rule, including coeff for hailstone */
    char rule_str[64];
    if (rule == RULE_HAILSTONE)
        snprintf(rule_str, sizeof(rule_str), "Hailstone (%" PRId64 "x+1)",
                 g_hailstone_coeff);
    else
        snprintf(rule_str, sizeof(rule_str), "%s", rule_display_name(rule));

    printf("\n k vs Z/G/C Scanner\n");
    printf(" ==================\n");
    printf(" k range  : %.10g .. %.10g\n", sk, ek);
    printf(" k points : %d\n", num_points);
    printf(" x0 range : %" PRId64 " .. %" PRId64 " (%" PRId64 " values)\n", start_x0, end_x0, x0_count);
    printf(" check_n  : %d\n", check_n);
    printf(" S rule   : %s\n", rule_str);
    printf(" Est. ops : %.2fM\n", (double)num_points * x0_count * check_n / 1e6);
    printf(" Output CSV : %s\n\n", OUTPUT_CSV);

    /* Pre-compute all k fractions */
    Frac *k_fracs = (Frac *)malloc((size_t)num_points * sizeof(Frac));
    double *k_arr = (double *)malloc((size_t)num_points * sizeof(double));
    char **k_str = (char **)malloc((size_t)num_points * sizeof(char *));

    for (int ki = 0; ki < num_points; ki++) {
        k_fracs[ki] = linspace_frac(start_k, end_k, ki, num_points);
        double kv = (double)k_fracs[ki].num / k_fracs[ki].den;
        k_arr[ki] = kv;
        k_str[ki] = (char *)malloc(48);
        snprintf(k_str[ki], 48, "%.15g", kv);
    }

    init_s_caches();

    /* Result arrays */
    int *z_arr = (int *)calloc((size_t)num_points, sizeof(int));
    int *g_arr = (int *)calloc((size_t)num_points, sizeof(int));
    int *c_arr = (int *)calloc((size_t)num_points, sizeof(int));

    if (!z_arr || !g_arr || !c_arr || !k_arr || !k_str || !k_fracs) {
        fprintf(stderr, "Memory allocation failed.\n"); return 1;
    }

    /* Per-thread working buffers — each thread needs its own tail/next_tail/seen
     * so they never alias.  Allocate nthreads sets up front. */
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    int64_t **tails      = (int64_t **)malloc((size_t)nthreads * sizeof(int64_t *));
    int64_t **next_tails = (int64_t **)malloc((size_t)nthreads * sizeof(int64_t *));
    SeenSet **seens      = (SeenSet **)malloc((size_t)nthreads * sizeof(SeenSet *));
    if (!tails || !next_tails || !seens) {
        fprintf(stderr, "Memory allocation failed.\n"); return 1;
    }
    for (int t = 0; t < nthreads; t++) {
        tails[t]      = (int64_t *)malloc((size_t)MAX_TAIL * sizeof(int64_t));
        next_tails[t] = (int64_t *)malloc((size_t)MAX_TAIL * sizeof(int64_t));
        seens[t]      = (SeenSet *)malloc(sizeof(SeenSet));
        if (!tails[t] || !next_tails[t] || !seens[t]) {
            fprintf(stderr, "Memory allocation failed.\n"); return 1;
        }
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

#ifdef _OPENMP
    printf(" Threads     : %d\n", nthreads);
#endif
    printf(" Starting calculations...\n");
    print_progress(0, num_points, 0.0);

    /* Shared progress counter (updated atomically via critical section) */
    int done = 0;

#pragma omp parallel for schedule(dynamic, 1) default(none) \
    shared(k_fracs, z_arr, g_arr, c_arr, tails, next_tails, seens, \
           num_points, start_x0, end_x0, check_n, rule, t0, done)
    for (int ki = 0; ki < num_points; ki++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        Frac k = k_fracs[ki];
        int64_t *tail      = tails[tid];
        int64_t *next_tail = next_tails[tid];
        SeenSet *seen      = seens[tid];

        int zeroed = 0, grown = 0, cycled = 0;

        for (int64_t x0 = start_x0; x0 <= end_x0; x0++) {
            int cls = classify_sequence(k.num, k.den, x0, check_n, rule,
                                        tail, next_tail, seen);
            if (cls == 0) zeroed++;
            else if (cls == 2) cycled++;
            else grown++;
        }

        z_arr[ki] = zeroed;
        g_arr[ki] = grown;
        c_arr[ki] = cycled;

        /* Progress bar: only update from one thread at a time */
#pragma omp critical
        {
            done++;
            struct timespec tnow;
            clock_gettime(CLOCK_MONOTONIC, &tnow);
            double el = (tnow.tv_sec - t0.tv_sec) + (tnow.tv_nsec - t0.tv_nsec) * 1e-9;
            print_progress(done, num_points, el);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("\n\n Done in %.3f seconds.\n\n",
           (double)(t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9);

    /* === Write CSV with metadata === */
    FILE *f = fopen(OUTPUT_CSV, "w");
    if (f) {
        fprintf(f, "# Rule: %s\n", rule_str);
        fprintf(f, "# k range: %.10g .. %.10g  (%d points)\n", sk, ek, num_points);
        fprintf(f, "# x0 range: %" PRId64 " .. %" PRId64 "\n", start_x0, end_x0);
        fprintf(f, "# check_n (max v_n): %d\n", check_n);
        fprintf(f, "# Generated on: %s", ctime((time_t*)&t1.tv_sec));
        fprintf(f, "k_decimal,zeroed,grown,cycled\n");

        for (int ki = 0; ki < num_points; ki++)
            fprintf(f, "%s,%d,%d,%d\n", k_str[ki], z_arr[ki], g_arr[ki], c_arr[ki]);
        fclose(f);
        printf(" Saved: %s\n\n", OUTPUT_CSV);
        printf(" Run 'python ZGC_grapher.py' to print plot\n\n");
    }

    /* Cleanup */
    for (int t = 0; t < nthreads; t++) {
        free(tails[t]); free(next_tails[t]); free(seens[t]);
    }
    free(tails); free(next_tails); free(seens);
    free(z_arr); free(g_arr); free(c_arr);
    free(k_arr); free(k_fracs);
    for (int ki = 0; ki < num_points; ki++) free(k_str[ki]);
    free(k_str);

    return 0;
}
