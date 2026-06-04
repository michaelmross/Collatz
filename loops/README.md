# collatz_verify.py

Computational companion to:

> **Loop Structure of Collatz-Type Functions `3x+n`: A Conjugacy Theorem and Powers of Three**
> Michael M. Ross, 2026

---

## Overview

This script studies the loop structure of Collatz-type functions `3x + n`
(for odd `n`) using the **odd-only reduced map** — the standard Collatz map
with even intermediates divided out. It verifies the four facts the note
relies on, each labeled with the result it supports:

1. **(A) Zero-sum identity** *(Lemma 1)* — the signed displacements around any
   cycle sum to zero. A pure telescoping identity, checked here on the cycles
   the search finds.

2. **(B) Rescaling identity** *(Lemma 2)* — `f_n(n·x) = n·f_1(x)` for every
   odd `n` and odd `x`. This is the algebraic engine of the note and holds
   unconditionally.

3. **(C) 3-adic confinement** *(Lemma 3)* — for `n = 3^k`, every orbit
   eventually has 3-adic valuation exactly `k`, i.e. settles onto the
   sublattice on which the rescaling identity is an exact copy of `3x + 1`.

4. **(D) Powers-of-three split** *(Proposition 1)* — over a range of odd `n`,
   the powers of three carry no nontrivial cycle (matching `n = 1`), while
   every other odd `n` carries at least one.

Together, (B) and (C) give the conjugacy theorem: `3x + 3^k` is a faithful
rescaled copy of `3x + 1` and therefore shares its cycle structure exactly.
This is the structural core of the note and is unconditional. Statements about
the **absence** of cycles for a power of three are conditional on the same
statement for `3x + 1`.

---

## A note on inverse pairs

Inverse pairs — odd inputs whose displacements are equal in magnitude and
opposite in sign — appear in the note only as **one possible mechanism** by
which a cycle can meet the zero-sum condition of Lemma 1 (see Remark 1). They
are a common but **not necessary** feature of cycles: the zero-sum balance can
also be met by three or more displacements with no two cancelling, and there
are nontrivial cycles that contain no inverse pair.

Accordingly, **this script does not test cycles against any inverse-pair
set**. The cycle search simply enumerates cycles; it makes no
inverse-pair-universality claim, and no such claim is part of the note.

---

## Requirements

- Python 3.7 or later
- No external dependencies (standard library only)

---

## Usage

```bash
python3 collatz_verify.py
```

All four checks run in sequence and print a `PASS`/`FAIL` summary at the end.
Runtime is well under a minute at the default ranges. Redirect to a file if
you want a record:

```bash
python3 collatz_verify.py > results.txt
```

---

## Parameters

The script has no global constants; ranges are passed as arguments in `main()`
and as function defaults. The main knobs:

| Where | Argument | Default | Description |
|-------|----------|---------|-------------|
| `check_powers_split` | `n_hi` | 500 | Upper bound for the odd-`n` scan in check (D) |
| `find_cycles` | `start_hi` | 4000 | Odd starting values `x = 1 … start_hi` for cycle search |
| `find_cycles` | `value_cap` | `10**7` | Orbit values above this are abandoned |
| `find_cycles` | `step_cap` | 5000 | Iteration cap per trajectory |
| `check_rescaling` | `x_hi` | 2000 | Odd `x` range for the identity check (B) |
| `check_confinement` | `ks` | `[1,2,3]` | Exponents `k` tested for `n = 3^k` in (C) |

To increase `n` change the value `n_hi` (e.g. 1000) in the `check_powers_split` call in
`main()`.  If cycles are missed for larger `n`, raise `step_cap` and
`value_cap` in `find_cycles`.

---

## Background

### The odd-only reduced map

For `3x + n`, the odd-only map applies `3x + n` and divides out all factors
of 2:

```
f(x) = (3x + n) / 2^v ,    v = v_2(3x + n)
```

Because `n` and `x` are odd, `3x + n` is even, so `f(x)` is again a positive
odd integer. The fixed point `x = n` (the trivial cycle `{n}`) is always
present, since `3n + n = 4n = 2^2·n`.

### Zero-sum condition

With displacement `d(x) = f(x) − x`, any cycle `x₁ → x₂ → ⋯ → x_L → x₁`
satisfies

```
Σ d(xᵢ) = 0
```

by telescoping. This is what check (A) confirms on the enumerated cycles.

### The conjugacy theorem

For every odd `n`, `f_n(n·x) = n·f_1(x)`, because `3(n·x) + n = n(3x + 1)`
and `n` odd gives `v_2(n(3x+1)) = v_2(3x+1)`. For `n = 3^k` a 3-adic
confinement drives every orbit onto the sublattice `3^k·(odd, 3 ∤ ·)`, where
this identity becomes an exact rescaling of `3x + 1`:

```
f_{3^k}(3^k · x) = 3^k · f_1(x)   for all odd x.
```

Hence `3x + 3^k` has a nontrivial cycle if and only if `3x + 1` does. This is
why the powers of three pattern with `n = 1`. The reduction fails for odd `n`
that are not powers of three (including multiples of 3 such as 15 or 21, which
reach 3-adic valuation `k` but spread across a sublattice strictly larger than
the conjugacy copy), and such `n` do carry nontrivial cycles — as check (D)
confirms.
