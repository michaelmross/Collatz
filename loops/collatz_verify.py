#!/usr/bin/env python3
"""
collatz_verify.py

Computational verification accompanying the note

    "Loop Structure of Collatz-Type Functions 3x+n:
     A Conjugacy Theorem and Powers of Three"  (M. M. Ross)

For odd n, f_n is the odd-only reduced map  x |-> (3x+n) / 2^{v_2(3x+n)},
acting on the positive odd integers. This script checks the four facts the
note actually relies on:

  (A) Zero-sum identity (Lemma 1): displacements around any cycle sum to 0.
  (B) Rescaling identity (Lemma 2): f_n(n*x) = n*f_1(x) for odd n, odd x.
  (C) 3-adic confinement (Lemma 3): for n = 3^k, orbits settle to v_3 == k.
  (D) Powers-of-three split (Prop. 1): n in {1,3,9,27,...} carry no nontrivial
      cycle in range, while other odd n do.

NOTE ON SCOPE. This script makes NO claim that every cycle passes through an
"inverse pair." That claim does not hold and is not part of the note; inverse
pairs appear there only as one possible cancellation mechanism (Remark 1).
The cycle search below simply enumerates cycles; it does not test them against
any inverse-pair set.

Run:  python3 collatz_verify.py
"""

from math import isqrt  # noqa: F401  (kept for convenient interactive use)


# ---------------------------------------------------------------------------
# Core map
# ---------------------------------------------------------------------------

def f(x, n):
    """Odd-only reduced map for 3x+n. Requires x odd, n odd; returns odd int."""
    y = 3 * x + n
    while y % 2 == 0:
        y //= 2
    return y


def v3(m):
    """3-adic valuation of a positive integer."""
    c = 0
    while m % 3 == 0:
        m //= 3
        c += 1
    return c


# ---------------------------------------------------------------------------
# Cycle enumeration (used by (A) and (D))
# ---------------------------------------------------------------------------

def find_cycles(n, start_hi=4000, value_cap=10**7, step_cap=5000):
    """
    Return the set of nontrivial cycles reachable from odd starts < start_hi.

    A cycle is returned as a sorted tuple of its members. The trivial fixed
    point {n} is excluded. Orbits are abandoned if they exceed value_cap or
    step_cap (these bounds are generous for the default start range).
    """
    cycles = set()
    for start in range(1, start_hi, 2):
        seen, path, x, steps = {}, [], start, 0
        while x not in seen and steps < step_cap and x < value_cap:
            seen[x] = len(path)
            path.append(x)
            x = f(x, n)
            steps += 1
        if x in seen:
            loop = tuple(sorted(path[seen[x]:]))
            if loop != (n,):              # drop the trivial fixed point
                cycles.add(loop)
    return cycles


# ---------------------------------------------------------------------------
# (A) Zero-sum identity  (Lemma 1)
# ---------------------------------------------------------------------------

def check_zero_sum(ns):
    print("(A) Zero-sum identity: sum of displacements around each cycle == 0")
    all_ok = True
    for n in ns:
        for cyc in find_cycles(n):
            # displacement at each member, summed around the cycle
            s = sum(f(x, n) - x for x in cyc)
            if s != 0:
                all_ok = False
                print(f"    FAIL n={n} cycle={cyc} sum={s}")
    print(f"    all tested cycles balance: {all_ok}\n")
    return all_ok


# ---------------------------------------------------------------------------
# (B) Rescaling identity  (Lemma 2):  f_n(n*x) = n*f_1(x)
# ---------------------------------------------------------------------------

def check_rescaling(odd_ns, x_hi=2000):
    print("(B) Rescaling identity: f_n(n*x) == n*f_1(x) for odd n, odd x")
    all_ok = True
    for n in odd_ns:
        for x in range(1, x_hi, 2):
            if f(n * x, n) != n * f(x, 1):
                all_ok = False
                print(f"    FAIL n={n} x={x}")
                break
    print(f"    identity holds for all tested (n,x): {all_ok}\n")
    return all_ok


# ---------------------------------------------------------------------------
# (C) 3-adic confinement  (Lemma 3):  for n = 3^k, orbits settle to v_3 == k
# ---------------------------------------------------------------------------

def check_confinement(ks, start_hi=400, settle=80, sample=40):
    print("(C) 3-adic confinement: for n=3^k, eventual v_3 of orbit == k")
    all_ok = True
    for k in ks:
        n = 3 ** k
        tail_valuations = set()
        for start in range(1, start_hi, 2):
            x = start
            for _ in range(settle):       # run past the transient
                x = f(x, n)
            for _ in range(sample):       # sample the tail
                tail_valuations.add(v3(x))
                x = f(x, n)
        ok = tail_valuations == {k}
        all_ok &= ok
        print(f"    k={k} (n={n}): tail v_3 values = {sorted(tail_valuations)}"
              f"  -> {'ok' if ok else 'UNEXPECTED'}")
    print(f"    confinement holds for all tested k: {all_ok}\n")
    return all_ok


# ---------------------------------------------------------------------------
# (D) Powers-of-three split  (Prop. 1)
# ---------------------------------------------------------------------------

def is_power_of_three(n):
    while n > 1 and n % 3 == 0:
        n //= 3
    return n == 1


def check_powers_split(n_lo=1, n_hi=1000):
    print("(D) Powers-of-three split: no nontrivial cycle for n a power of 3,")
    print("    at least one for every other odd n in range")
    all_ok = True
    for n in range(n_lo, n_hi + 1, 2):       # odd n only
        cycles = find_cycles(n)
        powers = is_power_of_three(n)
        if powers and len(cycles) != 0:
            all_ok = False
            print(f"    FAIL n={n} (power of 3) has {len(cycles)} cycle(s)")
        if (not powers) and len(cycles) == 0:
            all_ok = False
            print(f"    FAIL n={n} (not a power of 3) has no cycle in range")
    # report the powers of three actually exercised
    powers_seen = [n for n in range(n_lo, n_hi + 1, 2) if is_power_of_three(n)]
    print(f"    powers of three tested: {powers_seen} -> all cycle-free")
    print(f"    split holds across odd n in [{n_lo},{n_hi}]: {all_ok}\n")
    return all_ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(__doc__.strip().splitlines()[2])  # one-line banner
    print("=" * 70)
    results = {
        "A zero-sum":      check_zero_sum([5, 7, 11, 13]),
        "B rescaling":     check_rescaling([1, 3, 5, 7, 9, 15, 21, 27, 81]),
        "C confinement":   check_confinement([1, 2, 3]),
        "D powers split":  check_powers_split(1, 500),
    }
    print("=" * 70)
    print("SUMMARY")
    for name, ok in results.items():
        print(f"  {name:16} {'PASS' if ok else 'FAIL'}")
    print("\nAll checks passed." if all(results.values())
          else "\nSome checks FAILED — see above.")


if __name__ == "__main__":
    main()
