# WITHDRAWN: spike_certify_fast.py

This directory preserves the balanced-regime
(R = 2) cycle-word search from an earlier version of *Spike structures*
(`spike_certify_fast.py` and its output `spike_certify_fast_results.txt`).

**The result that search reported — "no R = 2 cycle for L ∈ [50, 200]" — is
withdrawn.** Its results are invalid because:

1. **The enumeration was incomplete.** Candidate words were generated under
   two structural restrictions — that the {1,2}-part consist of alternating
   runs, and that no 1 → 2 adjacency occur. Both were consequences of the
   "no 1 → 2 transition" claim of the first draft, which is false: the
   transition a = 1 → a = 2 occurs precisely on x ≡ 11 (mod 16), a quarter
   of the family D (see the correction remark in the current paper). The
   search therefore omitted a large class of admissible words, and its
   negative report cannot be concluded from it.

2. **The regime it searched is empty for classical reasons.** For any
   nontrivial cycle, multiplying 2^(a_i) = (3x_i + 1)/x_(i+1) around the
   cycle gives 2^S = ∏(3 + 1/x_i) < (10/3)^L, so the average exponent
   satisfies S/L < log2(10/3) ≈ 1.737 < 2. There was never an R = 2 regime
   to search; the exclusion is a one-line consequence of the cycle
   equation (see the average-exponent proposition in the current paper).

3. **The length range was already excluded unconditionally.** Cycles of
   length L ≤ 200 — indeed of any length below roughly 10^10 odd terms —
   are ruled out by known results (Hercher 2023, building on Simons–de
   Weger and the computational verification bound), for every admissible
   average exponent, not only R = 2.

The files are retained rather than deleted to keep the correction record
legible, matching the versioned-correction practice of the paper itself:
see the Version note and the Computational Status section of *Spike
structures and 2-adic transition laws in the accelerated Collatz map*
(https://doi.org/10.5281/zenodo.21049703).
<!-- TODO: update to the DOI of the revised version once deposited. -->

The repository's current computational artifact is `residue_graph.py`
(repository root), which reproduces the strongly connected-component
computations reported in the current paper.
