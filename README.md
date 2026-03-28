# Collatz R=2 Cycle Search

Computational support for:

**"Spike Structures, Finite-State Exclusion, and Bounded-Exponent Collatz Cycles"**
Michael M. Ross

*Zenodo Record:* [19039376](https://zenodo.org/records/19039376)

## What it does

Exhaustively searches for R=2 Collatz cycles — cycles whose average 2-adic division exponent equals 2 — across cycle lengths L ∈ [50, 200] with maximum division exponent bounded by Amax = 20.

**Result:** No cycles found. This computationally verifies the predictions of Theorem 6 across 7.56 × 10⁹ candidate exponent sequences.

## Key concepts

| Term | Meaning |
|------|---------|
| **L** | Cycle length — number of odd integers in the cycle |
| **b** | Number of spike positions (exponents ≥ 3) in the sequence |
| **n2** | Number of exponent-2 steps among the non-spike positions |
| **candidate** | A (position, spike-value) combination tested by the sieve |

For R=2 and b spikes, the total exponent sum is fixed at S = 2L, which uniquely determines n2 = L + b − big\_sum for each spike-value tuple. This eliminates the combinatorial explosion present in naive enumeration.

Each value of L produces three output rows, one per b ∈ {1, 2, 3}. The b=1 row always shows 18 candidates — one for each valid spike value in [3, 20], independent of L.

## Algorithm

1. **Segment decomposition** — For fixed (L, n2), the cycle equation residue N mod 2³² factors as a matrix product `seg_matrix @ pow_matrix`, where rows index spike positions and columns index spike values. Computed once per (L, n2) bucket.

2. **Vectorized prime sieve** — N ≡ 0 (mod p) is a necessary condition for a cycle. Applying primes 3, 5, 7, 11, 13, 17, 19 in sequence as vectorized matrix operations eliminates ~100% of candidates before any Python-level loop is entered.

3. **Exact verification** — The rare survivors (none found in practice) are checked by solving the cycle equation exactly and confirming primitivity by iteration.

## Requirements

```
Python >= 3.9
numpy
```

## Usage

```bash
python spike_certify_fast.py
```

Parameters are set at the top of the file:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Lmin`, `Lmax` | 50, 200 | Cycle length range |
| `Amax` | 20 | Maximum spike exponent |
| `B_VALUES` | [1, 2, 3] | Spike counts to search |
| `T_BITS` | 32 | Bit width for modular sieve |
| `NWORKERS` | all cores | Parallel worker count |

## Performance

Full L ∈ [50, 200] range on an 8-core machine: approximately 22 minutes.
Runtime scales roughly as O(L³) and linearly with available cores.

## Output

For each (L, b) pair:
```
L= 50 b=3  cands=6,759,648  t=0.49s  13,818,733/s  total=1s
```

A cycle found would be reported immediately with x₁, position tuple, and spike values. The final line confirms completion:
```
✅ DONE – no R=2 cycles found.
Total candidates: 7,558,109,808  Time: 1344.4s
```

## Citation

Michael M. Ross, "Spike Structures, Finite-State Exclusion, and
Bounded-Exponent Collatz Cycles" *[journal/preprint details]*.
