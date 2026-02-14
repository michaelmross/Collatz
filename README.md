# Collatz Cycle Search - Computational Verification

This code implements the computational verification described in the paper "Modular Spike Structures and Finite-State Exclusion of Bounded-Exponent Collatz Cycles" by Michael M. Ross.

## Overview

The program searches for R=2 Collatz cycles (cycles with average division exponent equal to 2) across cycle lengths L ∈ [50, 200] with maximum division exponent bounded by Amax = 20.

## Algorithm

The search uses a multi-stage filtering approach:

1. **Exponent sequence enumeration**: Generate candidate sequences respecting R=2 constraints and structural restrictions (no 1→2 transitions, run structures)

2. **2-adic sieve**: Test cycle equation consistency modulo 2^32

3. **Prime sieve**: Test consistency modulo primes {3, 5, 7, ..., 73}

4. **Exact verification**: For surviving candidates, compute x₁ exactly and verify the cycle

## Requirements

- Python 3.7 or higher
- Standard library only (no external dependencies)
- Multi-core CPU recommended (uses all available cores by default)

## Usage

```bash
python certify_mp_fixed.py
```

The program will search through all cycle lengths from L=50 to L=200 and report:
- Progress for each L value
- Candidate counts and search rates
- Any cycles found (none expected based on theoretical results)

## Parameters

Key parameters can be modified at the top of the file:

- `Lmin`, `Lmax`: Range of cycle lengths to search
- `R`: Maximum number of runs in the {1,2} portion
- `Amax`: Maximum spike value to consider
- `max_a`: Overall maximum division exponent
- `T_bits`: Bit width for 2-adic sieve (32 bits standard)
- `b_fixed`: Number of spike positions (fixed at 3)

## Expected Runtime

On a modern multi-core CPU (16+ cores):
- Each individual L value: 1-10 minutes
- Full L ∈ [50, 200] range: ~10-20 hours

Runtime scales approximately linearly with:
- Number of CPU cores (parallel speedup)
- Cycle length L (quadratically due to combination explosion)

## Output

For each L, the program reports:
- Number of exponent sequences examined
- Search time and candidate processing rate
- Final status: "DONE" (no cycles found) or details if a cycle is found

## Verification

The code includes internal consistency checks:
- `verify_cycle_primitive()`: Confirms x₁ actually forms a cycle under iteration
- Divisor check: Ensures the cycle is primitive (not a multiple of a shorter cycle)
- Exponent verification: Confirms the division exponents match the claimed sequence

## Code Structure

- `v2()`: Fast 2-adic valuation using bit operations
- `two_adic_full_cycle_sieve_fast()`: 2-adic modular consistency check
- `prime_sieve_pass()`: Prime modular filtering
- `cycle_equation_x1()`: Exact solution of cycle equation
- `verify_cycle_primitive()`: Final verification with primitivity check
- `worker()`: Parallel worker function for multi-core execution
- `precompute_*()`: Optimization tables for modular arithmetic

## Theoretical Justification

This computational search complements Theorem 12.6 in the paper, which proves that no bounded-exponent R=2 cycles can exist. The computation:

1. Verifies the theorem's predictions for concrete parameter ranges
2. Demonstrates the effectiveness of the structural constraints (no 1→2 transitions, run structures)
3. Provides empirical confidence in the theoretical framework

## Citation

If you use or modify this code, please cite:

Michael M. Ross, "Modular Spike Structures and Finite-State Exclusion of Bounded-Exponent Collatz Cycles," [journal/preprint details]

## Contact

[Your email/contact information]

## License

[Specify license - e.g., MIT, GPL, or "provided as supplementary material for the above paper"]
