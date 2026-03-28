"""
collatz_inverse_pairs.py
========================
Analysis of loop structure and inverse pairs for Collatz-type functions 3x+n.

Three main analyses:
  1. LOOP FINDER      -- finds all nontrivial loops for n=3..N_MAX, x=1..SEARCH_RANGE
  2. INVERSE PAIR CHECK -- verifies every loop touches an inverse pair value
  3. MODULAR ANALYSIS -- classifies which n produce loops (even / powers of 3 / other odd)

Conjugacy theorem (proved in the paper):
  For n = 3^k, f_n(3^k * x) = 3^k * f_1(x) for all odd x.
  Hence any loop in 3x+3^k would imply a loop in 3x+1.

Usage:
  python collatz_inverse_pairs.py
  Adjust the constants below to change search ranges.
"""

from collections import Counter
import math

# ── Parameters ────────────────────────────────────────────────────────────────
N_MAX        = 500    # scan n = 3 .. N_MAX
SEARCH_RANGE = 300    # starting x values: odd integers 1 .. SEARCH_RANGE
IP_SCAN      = 10000  # initial inverse pair scan range
IP_SCAN_EXT  = 30000  # extended scan for resolving apparent counterexamples
MAX_ITER     = 5000   # iteration cap per trajectory (increase if loops are missed)

# ── Core map ──────────────────────────────────────────────────────────────────

def odd_only_next(x, n):
    """
    One step of the odd-only reduced map for 3x+n.
    Applies 3x+n then divides out all factors of 2.
    """
    v = 3 * x + n
    while v % 2 == 0:
        v //= 2
    return v


def find_loop(start, n):
    """
    Follow the odd-only map from `start`.
    Returns frozenset of cycle members if a cycle is found, else None.
    Only odd starting values make sense; even starting values return None.
    """
    if start % 2 == 0:
        return None
    path, seen = [], {}
    x = start
    for _ in range(MAX_ITER):
        if x in seen:
            return frozenset(path[seen[x]:])
        seen[x] = len(path)
        path.append(x)
        x = odd_only_next(x, n)
    return None  # trajectory did not cycle within MAX_ITER steps


def get_inverse_pair_values(n, scan_range):
    """
    Scan odd integers 1..scan_range under the odd-only map for n.
    An inverse pair is a pair (a, b) with d(a) = -d(b), where d(x) = f(x) - x.
    Returns the set of all x values that belong to any inverse pair.
    """
    diffs = {}
    for x in range(1, scan_range, 2):
        d = odd_only_next(x, n) - x
        diffs.setdefault(d, []).append(x)
    ip_values = set()
    for d, xs in diffs.items():
        if d > 0 and -d in diffs:
            ip_values.update(xs)
            ip_values.update(diffs[-d])
    return ip_values


def is_power_of_3(n):
    """Returns True if n is a positive power of 3."""
    if n < 3:
        return False
    while n % 3 == 0:
        n //= 3
    return n == 1

# ── Analysis 1: Loop finder ────────────────────────────────────────────────────

def run_loop_finder():
    print("=" * 65)
    print("ANALYSIS 1: Loop finder")
    print(f"n = 3..{N_MAX},  x = 1..{SEARCH_RANGE} (odd),  MAX_ITER = {MAX_ITER}")
    print("=" * 65)

    results = {}  # n -> list of distinct loops (as frozensets)

    for n in range(3, N_MAX + 1):
        seen_cycles = set()
        for start in range(1, SEARCH_RANGE + 1, 2):
            c = find_loop(start, n)
            if c and len(c) > 1 and c not in seen_cycles:
                seen_cycles.add(c)
        results[n] = list(seen_cycles)

    # Print loops grouped by n (only n with loops)
    for n in range(3, N_MAX + 1):
        loops = results[n]
        if not loops:
            continue
        print(f"\nn = {n}:")
        for c in sorted(loops, key=min):
            ordered = sorted(c)
            print(f"  length {len(c):3d}  min={min(c):6d}  max={max(c):6d}  "
                  f"members={ordered[:8]}{'...' if len(ordered) > 8 else ''}")

    total = sum(len(v) for v in results.values())
    n_with = sum(1 for v in results.values() if v)
    print(f"\nTotal distinct nontrivial loops: {total}")
    print(f"n values with at least one loop: {n_with} / {N_MAX - 2}")
    return results

# ── Analysis 2: Inverse pair necessity check ──────────────────────────────────

def run_ip_check(results):
    print("\n" + "=" * 65)
    print("ANALYSIS 2: Inverse pair necessity check")
    print(f"IP initial scan: 1..{IP_SCAN},  extended scan: 1..{IP_SCAN_EXT}")
    print("=" * 65)

    total = touch = 0
    counterexamples = []

    for n, loops in results.items():
        if not loops:
            continue
        ip = get_inverse_pair_values(n, IP_SCAN)
        for c in loops:
            total += 1
            if c & ip:
                touch += 1
            else:
                # Retry with extended scan
                ip_ext = get_inverse_pair_values(n, IP_SCAN_EXT)
                if c & ip_ext:
                    touch += 1
                    print(f"  n={n}: resolved with extended scan "
                          f"(max loop member={max(c)})")
                else:
                    counterexamples.append((n, sorted(c)))

    print(f"\nNontrivial loops checked:       {total}")
    print(f"Touching an inverse pair:       {touch}")
    print(f"Genuine counterexamples:        {len(counterexamples)}")
    if counterexamples:
        print("COUNTEREXAMPLES:")
        for n, c in counterexamples:
            print(f"  n={n}: {c[:10]}{'...' if len(c) > 10 else ''}")
    else:
        print("RESULT: Every nontrivial loop passes through at least one "
              "inverse pair value.")

# ── Analysis 3: Modular classification ────────────────────────────────────────

def run_modular_analysis(results):
    print("\n" + "=" * 65)
    print("ANALYSIS 3: Classification of n by loop behaviour")
    print("=" * 65)

    loop_ns    = [n for n, v in results.items() if v]
    no_loop_ns = [n for n, v in results.items() if not v]

    even_loop     = [n for n in loop_ns if n % 2 == 0]
    even_no_loop  = [n for n in no_loop_ns if n % 2 == 0]
    odd_loop      = [n for n in loop_ns if n % 2 == 1]
    odd_no_loop   = [n for n in no_loop_ns if n % 2 == 1]

    pow3_loop    = [n for n in odd_loop if is_power_of_3(n)]
    pow3_no_loop = [n for n in odd_no_loop if is_power_of_3(n)]
    other_loop   = [n for n in odd_loop if not is_power_of_3(n)]
    other_no_loop = [n for n in odd_no_loop if not is_power_of_3(n)]

    print(f"\nTotal n values in range:        {N_MAX - 2}")
    print(f"  With loops:                   {len(loop_ns)}")
    print(f"  Without loops:                {len(no_loop_ns)}")
    print()
    print(f"Even n with loops:              {len(even_loop)}")
    print(f"Even n without loops:           {len(even_no_loop)}")
    print()
    print(f"Odd n (powers of 3) with loops:    {len(pow3_loop)}")
    print(f"Odd n (powers of 3) without loops: {len(pow3_no_loop)}"
          f"  {pow3_no_loop}")
    print()
    print(f"Other odd n with loops:         {len(other_loop)}")
    print(f"Other odd n without loops:      {len(other_no_loop)}"
          + (f"  {other_no_loop}" if other_no_loop else ""))

    print()
    conj_holds = (len(even_loop) == 0 and
                  len(pow3_loop) == 0 and
                  len(other_no_loop) == 0)
    print("Conjecture holds in search range "
          "(even=no loop, 3^k=no loop, other odd=loop)?", conj_holds)

# ── Conjugacy theorem verification ────────────────────────────────────────────

def verify_conjugacy():
    print("\n" + "=" * 65)
    print("CONJUGACY VERIFICATION: f_{3^k}(3^k * x) = 3^k * f_1(x)")
    print("=" * 65)
    for k in range(1, 6):
        n = 3 ** k
        if n > N_MAX:
            break
        factor = 3 ** k
        ok = all(
            odd_only_next(factor * x, n) == factor * odd_only_next(x, 1)
            for x in range(1, 10001, 2)
        )
        print(f"  k={k}, n=3^{k}={n:4d}: holds for x=1..9999 (odd)? {ok}")

# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    results = run_loop_finder()
    run_ip_check(results)
    run_modular_analysis(results)
    verify_conjugacy()
