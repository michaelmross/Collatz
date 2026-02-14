# Read https://github.com/michaelmross/tasktracker/blob/collatz/spike_certify_readme.md

import math
import itertools
import time
import bisect
import os
import multiprocessing as mp

# =========================
# CLEAN PARAMETERS
# =========================
Lmin = 50
Lmax = 200

R = 2
Amax = 20
max_a = 63
T_bits = 32
b_fixed = 3

PRIMES = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73]

NWORKERS = max(1, (os.cpu_count() or 1))


# =========================
# Helpers
# =========================
def v2(n: int) -> int:
    return (n & -n).bit_length() - 1

def verify_cycle_primitive(x1: int, L: int, max_a: int) -> bool:
    x = x1
    for _ in range(L):
        t = 3*x + 1
        a = v2(t)
        if a > max_a:
            return False
        x = t >> a
    if x != x1:
        return False

    divs = set()
    for d in range(1, int(L**0.5) + 1):
        if L % d == 0:
            divs.add(d); divs.add(L//d)
    divs.discard(L)
    for d in divs:
        y = x1
        for _ in range(d):
            t = 3*y + 1
            y = t >> v2(t)
        if y == x1:
            return False
    return True

def necessary_S(L: int) -> int:
    return int(math.floor(L * math.log2(3))) + 1

def compositions(n, k):
    for cuts in itertools.combinations(range(1, n), k-1):
        prev = 0
        parts = []
        for c in cuts:
            parts.append(c-prev)
            prev = c
        parts.append(n-prev)
        yield parts

def count_twos_from_runs(start_bit: int, runs):
    twos = 0
    bit = start_bit
    for Lr in runs:
        if bit == 2:
            twos += Lr
        bit = 1 if bit == 2 else 2
    return twos

def precompute_run_structures(m: int, R: int):
    items = []
    for r in range(1, min(R, m) + 1):
        for runs in compositions(m, r):
            for start_bit in (1, 2):
                n2 = count_twos_from_runs(start_bit, runs)
                items.append((n2, start_bit, runs))
    items.sort(key=lambda t: t[0])
    n2_list = [t[0] for t in items]
    return items, n2_list

def build_a_seq_from_runs(L: int, big_pos, big_vals, start_bit: int, runs):
    a = [0] * L
    pos_set = set(big_pos)
    for p, v in zip(big_pos, big_vals):
        a[p] = v
    bit = start_bit
    run_i = 0
    left = runs[0]
    for i in range(L):
        if i in pos_set:
            continue
        a[i] = bit
        left -= 1
        if left == 0:
            run_i += 1
            if run_i < len(runs):
                bit = 1 if bit == 2 else 2
                left = runs[run_i]
    return a


# =========================
# 2-adic fast tables
# =========================
def precompute_mod2T_tables(Lmax: int, T_bits: int):
    mod = 1 << T_bits
    pow3 = [1] * (Lmax + 1)
    for k in range(1, Lmax + 1):
        pow3[k] = (pow3[k-1] * 3) % mod
    inv_minus_3L = [0] * (Lmax + 1)
    for L in range(1, Lmax + 1):
        D = (-pow3[L]) % mod
        inv_minus_3L[L] = pow(D, -1, mod)
    return mod, pow3, inv_minus_3L

def two_adic_full_cycle_sieve_fast(a_seq, L: int, mod: int, pow3: list[int], inv_minus_3L: list[int], T_bits: int) -> bool:
    N = 0
    A = 0
    for j in range(L):
        if j > 0:
            A += a_seq[j-1]
        twoA = 0 if A >= T_bits else (1 << A)
        N = (N + (pow3[L-1-j] * twoA)) % mod

    x1 = (N * inv_minus_3L[L]) % mod
    if x1 % 2 == 0:
        return False

    x = x1
    for a in a_seq:
        val = (3 * x + 1) % mod
        if val == 0:
            v = T_bits
        else:
            v = (val & -val).bit_length() - 1
        if v < min(a, T_bits):
            return False
        if a < T_bits and v != a:
            return False
        x = (val >> a) % mod
        if x % 2 == 0:
            return False

    return x == x1


# =========================
# Prime sieve tables
# =========================
def precompute_prime_tables(Lmax: int, Smax: int):
    pow2 = {}
    pow3 = {}
    for p in PRIMES:
        a2 = [1] * (Smax + 1)
        a3 = [1] * (Lmax + 1)
        for e in range(1, Smax + 1):
            a2[e] = (a2[e-1] * 2) % p
        for j in range(1, Lmax + 1):
            a3[j] = (a3[j-1] * 3) % p
        pow2[p] = a2
        pow3[p] = a3
    return pow2, pow3

def prime_sieve_pass(a_seq, pow2, pow3) -> bool:
    L = len(a_seq)
    S = sum(a_seq)
    for p in PRIMES:
        Dp = (pow2[p][S] - pow3[p][L]) % p
        if Dp != 0:
            continue
        Np = 0
        A = 0
        for j in range(L):
            if j > 0:
                A += a_seq[j-1]
            Np += (pow3[p][L-1-j] * pow2[p][A]) % p
        if (Np % p) != 0:
            return False
    return True

def cycle_equation_x1(a_seq):
    L = len(a_seq)
    S = sum(a_seq)
    D = (1 << S) - pow(3, L)
    if D <= 0:
        return None
    N = 0
    A = 0
    for j in range(L):
        if j > 0:
            A += a_seq[j-1]
        N += pow(3, L-1-j) * (1 << A)
    if N % D != 0:
        return None
    x1 = N // D
    if x1 <= 0 or x1 % 2 == 0 or x1 == 1:
        return None
    return x1


# =========================
# Big triples
# =========================
def precompute_big_triples_sorted(Amax: int):
    vals = range(3, Amax + 1)
    triples = []
    for a in vals:
        for b in vals:
            for c in vals:
                triples.append((a+b+c, (a, b, c)))
    triples.sort(key=lambda t: t[0])
    sums = [t[0] for t in triples]
    return triples, sums


# =========================
# Worker (NO local functions)
# =========================
def worker(args):
    (wid, L, m, Smin, pairs, eligible_triples, run_items, run_n2s,
     mod2T, pow3_2T, inv_minus_3L, pow2p, pow3p) = args

    candidates = 0

    for (p2, p3) in pairs:
        big_pos = (0, p2, p3)

        for big_sum, big_vals in eligible_triples:
            n2_min = max(0, Smin - big_sum - m)
            if n2_min > m:
                continue

            j0 = bisect.bisect_left(run_n2s, n2_min)
            for (_, start_bit, runs) in run_items[j0:]:
                a_seq = build_a_seq_from_runs(L, big_pos, big_vals, start_bit, runs)
                candidates += 1

                if not two_adic_full_cycle_sieve_fast(a_seq, L, mod2T, pow3_2T, inv_minus_3L, T_bits):
                    continue

                # Only if anything ever survives:
                if not prime_sieve_pass(a_seq, pow2p, pow3p):
                    continue
                x1 = cycle_equation_x1(a_seq)
                if x1 is None:
                    continue
                if verify_cycle_primitive(x1, L, max_a):
                    return ("FOUND", wid, L, candidates, big_pos, big_vals, start_bit, runs, x1)

    return ("DONE", wid, L, candidates)


# =========================
# Main
# =========================
def main():
    triples, triple_sums = precompute_big_triples_sorted(Amax)
    mod2T, pow3_2T, inv_minus_3L = precompute_mod2T_tables(Lmax, T_bits)
    pow2p, pow3p = precompute_prime_tables(Lmax=Lmax, Smax=Lmax*Amax)

    run_cache = {}

    for L in range(Lmin, Lmax + 1):
        Smin = necessary_S(L)
        m = L - b_fixed
        if m <= 0:
            continue

        if m not in run_cache:
            run_cache[m] = precompute_run_structures(m, R)
        run_items, run_n2s = run_cache[m]

        # forced lower bound for big_sum
        big_sum_min = Smin - 2*m
        idx = bisect.bisect_left(triple_sums, max(9, big_sum_min))
        eligible_triples = triples[idx:]

        all_pairs = list(itertools.combinations(range(1, L), 2))
        shards = [[] for _ in range(NWORKERS)]
        for i, pr in enumerate(all_pairs):
            shards[i % NWORKERS].append(pr)

        print(f"\n=== L={L} | pairs={len(all_pairs)} | eligible_triples={len(eligible_triples)} | run_items={len(run_items)} | workers={NWORKERS} ===")
        tL0 = time.time()

        with mp.Pool(processes=NWORKERS) as pool:
            tasks = [
                (wid, L, m, Smin, shards[wid], eligible_triples, run_items, run_n2s,
                 mod2T, pow3_2T, inv_minus_3L, pow2p, pow3p)
                for wid in range(NWORKERS)
            ]

            total_candidates = 0
            for res in pool.imap_unordered(worker, tasks):
                if res[0] == "FOUND":
                    _, wid, L, cand, big_pos, big_vals, start_bit, runs, x1 = res
                    print("\n✅ FOUND cycle!")
                    print(f"L={L} x1={x1} worker={wid} candidates_checked_by_worker={cand}")
                    print(f"big_pos={big_pos} big_vals={big_vals} start={start_bit} runs={runs}")
                    pool.terminate()
                    pool.join()
                    return

                _, wid, L, cand = res
                total_candidates += cand

        dt = time.time() - tL0
        rate = total_candidates / dt if dt else 0
        print(f"=== L={L} DONE | candidates={total_candidates:,} | time={dt:.1f}s | rate={rate:,.0f}/s ===")

    print("\n✅ DONE: No cycles found in the requested L-range under constraints.")


if __name__ == "__main__":
    mp.freeze_support()
    main()
