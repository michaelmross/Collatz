#!/usr/bin/env python3
"""
spike_certify_fast.py  —  Optimized Collatz R=2 cycle search
=============================================================

Requirements:  Python >= 3.9, numpy
Runtime: ~10-25 minutes
See README.md for details.
"""

import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
           "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import time, math, itertools, multiprocessing as mp
from itertools import combinations
import numpy as np

# ── Parameters ────────────────────────────────────────────────────────────────
Lmin     = 120
Lmax     = 128
Amax     = 20
max_a    = 63          # used only in final exact verification
T_BITS   = 32
MOD      = 1 << T_BITS
MOD64    = np.uint64(MOD)
B_VALUES = [1, 2, 3]
PRIMES   = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73]
NWORKERS = max(1, os.cpu_count() or 1)

# ── Basic helpers ─────────────────────────────────────────────────────────────
def v2(n: int) -> int:
    return (n & -n).bit_length() - 1

def verify_cycle_primitive(x1: int, L: int) -> bool:
    x = x1
    for _ in range(L):
        t = 3 * x + 1
        a = v2(t)
        if a > max_a: return False
        x = t >> a
    if x != x1: return False
    for d in range(2, math.isqrt(L) + 2):
        if L % d == 0:
            for dd in {d, L // d}:
                if 1 < dd < L:
                    y = x1
                    for _ in range(dd):
                        t = 3 * y + 1; y = t >> v2(t)
                    if y == x1: return False
    return True

def cycle_x1_exact(a_seq: list):
    L, S = len(a_seq), sum(a_seq)
    D = (1 << S) - pow(3, L)
    if D <= 0: return None
    N = A = 0
    for j in range(L):
        N += pow(3, L - 1 - j) * (1 << A); A += a_seq[j]
    if N % D != 0: return None
    x1 = N // D
    return x1 if (x1 > 1 and x1 % 2 == 1) else None

def prime_sieve(S: int, L: int, a_seq: list, pw2: dict, pw3: dict) -> bool:
    for p in PRIMES:
        if (pw2[p][S] - pw3[p][L]) % p != 0: continue
        A = Np = 0
        for j in range(L):
            Np = (Np + pw3[p][L-1-j] * pw2[p][A]) % p; A += a_seq[j]
        if Np % p != 0: return False
    return True

def rebuild_seq(L: int, b: int, pos: tuple, vals: tuple, n2: int) -> list:
    m = L - b; a = [0] * L
    for p, v in zip(pos, vals): a[p] = v
    fill = [2]*n2 + [1]*(m-n2); fi = 0
    for j in range(L):
        if a[j] == 0: a[j] = fill[fi]; fi += 1
    return a

def make_prime_tables(Lmax: int, Smax: int) -> tuple:
    pw2, pw3 = {}, {}
    for p in PRIMES:
        pw2[p] = [pow(2, e, p) for e in range(Smax + 1)]
        pw3[p] = [pow(3, j, p) for j in range(Lmax + 1)]
    return pw2, pw3

# ── Segment weight prefix sums ────────────────────────────────────────────────
def make_Q_all(L: int, b: int, n2_values) -> dict:
    """
    For each n2 in n2_values, compute prefix-sum arrays Q_1 .. Q_b where
        Q_k[j] = sum_{i=1}^{j} w_k[i],
        w_k[i] = 3^{L-1-i} * 2^{(i-k) + min(i-k, n2)}   for i >= k.

    Returns {n2: (q0, Q_1, ..., Q_b)}.
    """
    p3 = np.array([pow(3, L-1-j, MOD) for j in range(L)], dtype=np.uint64)
    q0 = int(p3[0])          # j=0: ns_count=0, 2^0=1
    result = {}
    for n2 in n2_values:
        Qs = []
        for k in range(1, b + 1):
            q = np.zeros(L, dtype=np.uint64)
            for j in range(k, L):
                ns = j - k; e = ns + min(ns, n2)
                q[j] = p3[j] * (np.uint64(0) if e >= T_BITS
                                 else np.uint64(1 << e)) % MOD64
            Qs.append(np.cumsum(q, dtype=np.uint64) % MOD64)
        result[n2] = (q0, *Qs)
    return result

# ── Segment matrices (specialised for each b) ─────────────────────────────────
def seg_b3(pairs_np: np.ndarray, Qs: tuple) -> np.ndarray:
    q0, Q1, Q2, Q3 = Qs
    p1 = pairs_np[:, 0].astype(np.int64)
    p2 = pairs_np[:, 1].astype(np.int64)
    seg = np.empty((len(pairs_np), 4), dtype=np.uint64)
    seg[:, 0] = np.uint64(q0)
    seg[:, 1] = Q1[p1]
    seg[:, 2] = (Q2[p2].astype(np.int64) - Q2[p1].astype(np.int64)) % MOD
    seg[:, 3] = (int(Q3[-1]) - Q3[p2].astype(np.int64)) % MOD
    return seg

def seg_b2(pairs_np: np.ndarray, Qs: tuple) -> np.ndarray:
    q0, Q1, Q2 = Qs
    p1 = pairs_np[:, 0].astype(np.int64)    # (P,)
    seg = np.empty((len(pairs_np), 3), dtype=np.uint64)
    seg[:, 0] = np.uint64(q0)
    seg[:, 1] = Q1[p1]                       # (P,) indexing into Q1
    seg[:, 2] = (int(Q2[-1]) - Q2[p1].astype(np.int64)) % MOD
    return seg

def seg_b1(Qs: tuple) -> np.ndarray:
    q0, Q1 = Qs
    return np.array([[q0, int(Q1[-1])]], dtype=np.uint64)

# ── Power matrices ────────────────────────────────────────────────────────────
def pm_b3(spike_list: list) -> np.ndarray:
    T = len(spike_list); pm = np.zeros((4, T), dtype=np.uint64); pm[0] = 1
    for t, (sv0, sv1, sv2) in enumerate(spike_list):
        pm[1,t]=pow(2,sv0,MOD); pm[2,t]=pow(2,sv0+sv1,MOD); pm[3,t]=pow(2,sv0+sv1+sv2,MOD)
    return pm

def pm_b2(spike_list: list) -> np.ndarray:
    T = len(spike_list); pm = np.zeros((3, T), dtype=np.uint64); pm[0] = 1
    for t, (sv0, sv1) in enumerate(spike_list):
        pm[1,t]=pow(2,sv0,MOD); pm[2,t]=pow(2,sv0+sv1,MOD)
    return pm

def pm_b1(spike_list: list) -> np.ndarray:
    T = len(spike_list); pm = np.zeros((2, T), dtype=np.uint64); pm[0] = 1
    for t, (sv0,) in enumerate(spike_list): pm[1,t]=pow(2,sv0,MOD)
    return pm

# ── Modular matrix multiply (16-bit-split float64 trick) ─────────────────────
def matmul_mod(seg: np.ndarray, pm: np.ndarray) -> np.ndarray:
    """
    (seg @ pm) mod 2^32.  Both inputs uint64 with values < 2^32.

    16-bit split: seg = seg_hi*2^16 + seg_lo, both halves < 2^16.
    Each matmul product  ≤ K * (2^16-1) * (2^32-1) < 4 * 2^48 < 2^53  (exact in float64).
    We take mod before the *65536 scale so that (hi%MOD)*65536 < 2^48 < 2^53.
    """
    seg_lo = (seg & np.uint64(0xFFFF)).view(np.int64).astype(np.float64)
    seg_hi = (seg >> np.uint64(16)).view(np.int64).astype(np.float64)
    pm_f   = pm.astype(np.float64)
    fMOD   = float(MOD)
    lo     = seg_lo @ pm_f                     # each entry < 2^50 < 2^53 ✓
    hi     = seg_hi @ pm_f                     # each entry < 2^50 < 2^53 ✓
    N_f    = (lo % fMOD + (hi % fMOD) * 65536.0) % fMOD
    return N_f.astype(np.uint64)

# ── Triple catalogue ──────────────────────────────────────────────────────────
def build_catalogue(b: int, Amax: int) -> dict:
    cat: dict = {}
    for vals in itertools.product(range(3, Amax + 1), repeat=b):
        cat.setdefault(sum(vals), []).append(vals)
    return cat

# ── Worker ────────────────────────────────────────────────────────────────────
def worker(args: tuple) -> tuple:
    wid, L, b, pair_shard, cat_b, pw2, pw3 = args
    m     = L - b
    D     = (pow(2, 2*L, MOD) - pow(3, L, MOD)) % MOD
    inv_D = np.uint64(pow(int(D), -1, MOD))

    valid = {big_sum: L+b-big_sum
             for big_sum in range(b*3, b*Amax+1)
             if 0 <= L+b-big_sum <= m and big_sum in cat_b}
    if not valid:
        return ("DONE", wid, L, b, 0)

    Q_cache = make_Q_all(L, b, sorted(set(valid.values())))

    if b == 1:
        pairs_np = np.zeros((1, 0), dtype=np.int64)
    else:
        pairs_np = np.array(pair_shard, dtype=np.int64)  # (P, b-1)
    if len(pairs_np) == 0:
        return ("DONE", wid, L, b, 0)

    total = 0
    for big_sum, n2 in valid.items():
        spike_list = cat_b[big_sum]
        Qs = Q_cache[n2]

        if   b == 3: seg = seg_b3(pairs_np, Qs); pm = pm_b3(spike_list)
        elif b == 2: seg = seg_b2(pairs_np, Qs); pm = pm_b2(spike_list)
        else:        seg = seg_b1(Qs);            pm = pm_b1(spike_list)

        N_mat  = matmul_mod(seg, pm)
        x1_mat = N_mat * inv_D % MOD64
        total += len(pairs_np) * len(spike_list)

        mask = (x1_mat % np.uint64(2) == np.uint64(1)) & (x1_mat > np.uint64(1))
        if not mask.any():
            continue

        # ── Filter 2: vectorised N mod p cascade ─────────────────────────
        # N mod p == 0 is a necessary condition (cycle equation).
        # Each prime independently rejects ~(1-1/p) of survivors.
        # Chain 3,5,7,11,13,17,19 reduces 3.4M candidates to ~0 for L=50.
        seg_f = seg.astype(np.float64)
        pm_f2 = pm.astype(np.float64)
        for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
            if not mask.any():
                break
            N_p = (seg_f % p) @ (pm_f2 % p) % p
            mask &= (N_p.astype(np.int32) == 0)

        if not mask.any():
            continue

        # ── Python fallback: exact verification (extremely rare) ─────────
        for pi, ti in np.argwhere(mask):
            sv  = spike_list[ti]
            if b == 1: pos = (0,)
            elif b == 2: pos = (0, int(pairs_np[pi, 0]))
            else: pos = (0, int(pairs_np[pi, 0]), int(pairs_np[pi, 1]))

            seq = rebuild_seq(L, b, pos, sv, n2)
            x1m = int(x1_mat[pi, ti])

            x, ok = x1m, True
            for a in seq:
                val = (3*x+1) % MOD
                vv  = T_BITS if val == 0 else (val & -val).bit_length() - 1
                if a < T_BITS and vv != a: ok = False; break
                x = (val >> a) % MOD
                if x % 2 == 0: ok = False; break
            if not ok or x != x1m: continue

            if not prime_sieve(2*L, L, seq, pw2, pw3): continue
            x1e = cycle_x1_exact(seq)
            if x1e is None: continue
            if verify_cycle_primitive(x1e, L):
                return ("FOUND", wid, L, b, total, pos, sv, int(x1e))

    return ("DONE", wid, L, b, total)

# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    print(f"R=2 Collatz cycle search  L∈[{Lmin},{Lmax}]  "
          f"Amax={Amax}  b∈{B_VALUES}  workers={NWORKERS}\n")

    pw2, pw3 = make_prime_tables(Lmax, 2*Lmax + 1)
    cat      = {b: build_catalogue(b, Amax) for b in B_VALUES}
    t_start  = time.time()
    total_all = 0

    with mp.Pool(NWORKERS) as pool:
        for L in range(Lmin, Lmax + 1):
            for b in B_VALUES:
                if L - b <= 0: continue

                if b == 1:   all_pos = [()]
                elif b == 2: all_pos = [(p,) for p in range(1, L)]
                else:        all_pos = list(combinations(range(1, L), b-1))

                shards = [all_pos[i::NWORKERS] for i in range(NWORKERS)]
                tasks  = [(wid, L, b, shards[wid], cat[b], pw2, pw3)
                          for wid in range(NWORKERS) if shards[wid]]
                t0 = time.time()
                Lc = 0
                for res in pool.imap_unordered(worker, tasks):
                    if res[0] == "FOUND":
                        _, wid, Lf, bf, cand, pos, sv, x1 = res
                        print(f"\n✅ CYCLE FOUND  L={Lf}  b={bf}  x1={x1}")
                        print(f"   positions={pos}  spike_vals={sv}")
                        print(f"   Elapsed: {time.time()-t_start:.1f}s")
                        pool.terminate(); pool.join(); return
                    Lc += res[4]

                dt = time.time() - t0; total_all += Lc
                rate = Lc/dt if dt else 0
                print(f"L={L:3d} b={b}  cands={Lc:>12,}  "
                      f"t={dt:6.2f}s  {rate:>10,.0f}/s  "
                      f"total={time.time()-t_start:.0f}s")

    print(f"\n✅ DONE – no R=2 cycles found.")
    print(f"Total candidates: {total_all:,}  Time: {time.time()-t_start:.1f}s")

if __name__ == "__main__":
    mp.freeze_support(); main()
