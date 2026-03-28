# collatz_inverse_pairs.py

Computational companion to:

> **Inverse Pairs and the Loop Structure of Collatz-Type Functions**  
> Michael M. Ross, 2026

---

## Overview

This script investigates the loop structure of Collatz-type functions of the
form `3x + n` using an **odd-only reduced map** — the standard Collatz map
with even intermediates removed. It performs three analyses:

1. **Loop finder** — finds all distinct nontrivial cycles for each `n` in a
   specified range, starting from odd values of `x` up to a search bound.

2. **Inverse pair necessity check** — verifies that every nontrivial loop
   passes through at least one *inverse pair value*. An inverse pair is a
   pair of inputs `(a, b)` whose signed differences satisfy `d(a) = −d(b)`,
   providing the cancellation required by the zero-sum loop condition.

3. **Modular classification** — classifies each `n` by loop behaviour,
   confirming the pattern:
   - **Even `n`**: no loops
   - **Odd `n` that are powers of 3**: no loops (proved by conjugacy to
     `3x + 1` in the paper)
   - **All other odd `n`**: at least one loop found

A fourth section **verifies the conjugacy theorem** algebraically:
`f_{3^k}(3^k · x) = 3^k · f_1(x)` for all odd `x`.

---

## Key result

Over `n = 3` to `500` and starting values `x = 1` to `300` (odd only),
**856 distinct nontrivial loops** were found. Every single one passes through
at least one inverse pair value. Zero counterexamples were found after
extending the inverse pair scan range to `30,000` to resolve apparent
exceptions at shorter scan bounds.

---

## Requirements

- Python 3.7 or later
- No external dependencies (standard library only)

---

## Usage

```bash
python collatz_inverse_pairs.py
```

Output is printed to stdout. For the full `n = 3` to `500` run, expect a
runtime of several minutes. Redirect to a file if needed:

```bash
python collatz_inverse_pairs.py > results.txt
```

---

## Parameters

All parameters are defined as constants near the top of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_MAX` | 500 | Upper bound for `n` (scans `n = 3` to `N_MAX`) |
| `SEARCH_RANGE` | 300 | Odd starting values `x = 1` to `SEARCH_RANGE` |
| `IP_SCAN` | 10000 | Initial scan range for inverse pair detection |
| `IP_SCAN_EXT` | 30000 | Extended scan range used when a loop appears to miss all inverse pairs at the shorter range |
| `MAX_ITER` | 5000 | Iteration cap per trajectory; increase if loops are missed for large `n` |

For a quick test run, set `N_MAX = 100` and `SEARCH_RANGE = 100`.

---

## Background

### The odd-only reduced map

For a Collatz-type function `3x + n`, the odd-only map `f` applies `3x + n`
and then divides out all factors of 2:

```
f(x) = (3x + n) / 2^v
```

where `v` is the 2-adic valuation of `3x + n`. The three congruence cases are:

- If `x + n ≡ 0 (mod 4)`: `f(x) = (3x + n) / 2`, `d(x) = (x + n) / 2`
- If `x − n ≡ 0 (mod 8)`: `f(x) = (3x + n) / 4`, `d(x) = (n − x) / 4`
- Else: further halvings

### Zero-sum condition

Any loop `x₁ → x₂ → ⋯ → xₖ → x₁` must satisfy:

```
Σ d(xᵢ) = 0
```

This is a pure arithmetic identity (telescoping sum). Inverse pairs are the
mechanism by which this cancellation is achieved.

### Inverse pairs for n > 7

For odd `n > 7`, two inverse pairs always exist at explicit inputs:

| Input | Class | Signed difference |
|-------|-------|------------------|
| `x = n − 6` | 1 | `+(n − 3)` |
| `x = 5n − 12` | 2 | `−(n − 3)` |
| `x = n + 6` | 1 | `+(n + 3)` |
| `x = 5n + 12` | 2 | `−(n + 3)` |

This is proved algebraically in the paper (Proposition 6).

### The conjugacy theorem

For `n = 3^k`:

```
f_{3^k}(3^k · x) = 3^k · f_1(x)   for all odd x
```

**Proof:** `3·(3^k·x) + 3^k = 3^k·(3x + 1)`. Since `3^k` is odd,
`v₂(3^k·(3x+1)) = v₂(3x+1)`, so dividing out powers of 2 gives
`3^k · f_1(x)`.

This means any loop in `3x + 3^k` would imply a loop in `3x + 1` — and
since `3x + 1` terminates at the floor of the positive integers via its
inverse pair, no such loop is possible.

---

## Related work

- Ross, M. M. *Modular Spike Structures and Finite-State Exclusion of
  Bounded-Exponent Collatz Cycles*. Zenodo, 2025.
