# Collatz

Code supporting Michael M. Ross's papers on the structure of the accelerated
Collatz map.

## Contents

| Path | Purpose |
|---|---|
| `residue_graph.py` | Constructs the nondeterministic residue graph G(K,M) and reports its strongly connected components. The repository's primary computational artifact. |
| `loops/` | Search code for the companion paper on Collatz-type maps 3x+n (cited there; do not relocate — the paper links this path). |
| `archive/` | Superseded code retained for the record. **The result it reported is withdrawn** — see `archive/README.md`. |

## residue_graph.py

Python 3, standard library only.

The graph G(K,M) has a node for each odd residue r mod 2^M and an edge
r → s whenever some integer x ≡ r with v2(3x+1) ≤ max_a = K+3 has
U(x) ≡ s, where U(x) = (3x+1)/2^(v2(3x+1)) is the accelerated map. Knowing
x mod 2^M determines U(x) only mod 2^(M−a), so each node has a full coset
of 2^a successors — the graph is deliberately nondeterministic. Any
bounded-exponent cycle of integers would project to a directed cycle here;
the computation shows the graph is nowhere near acyclic.

### Usage

Single graph (default K=7, M=13):

    $ python3 residue_graph.py 7 13
    G_{K=7,M=13}: 4096 nodes, 40960 edges
    SCC sizes (top 5): [4088, 1, 1, 1, 1];  singletons: 8

Sweep over moduli at fixed K:

    $ python3 residue_graph.py --sweep 7 11 20
      M     nodes       edges  giant SCC  outside
     11     1,024      10,240      1,023        1
     12     2,048      20,480      2,045        3
     13     4,096      40,960      4,088        8
     14     8,192      81,920      8,172       20
     15    16,384     163,840     16,336       48
     16    32,768     327,680     32,656      112
     17    65,536     655,360     65,280      256
     18   131,072   1,310,720    130,496      576
     19   262,144   2,621,440    260,864    1,280
     20   524,288   5,242,880    521,472    2,816

At every modulus there is a single giant strongly connected component;
the nodes outside it number (M−9)·2^(M−12) for K = 7 — the residue
classes whose exponent exceeds max_a, plus the nodes stranded behind
them — a fraction that vanishes as M grows. This reproduces the
finite-state section of the paper below: the sound residue graph is
recurrent, so acyclicity of any finite quotient cannot exclude cycles.
(The deterministic variant that *is* acyclic omits the very edges a
genuine cycle would use; see the paper.)

## Papers and provenance

- M. M. Ross, *Spike structures and 2-adic transition laws in the
  accelerated Collatz map*, Zenodo.
  https://doi.org/10.5281/zenodo.21049703
  <!-- TODO: update to the DOI of the revised version once deposited,
       or switch to the record's concept DOI. -->
- M. M. Ross, *Loop structure of Collatz-type functions 3x+n: a conjugacy
  theorem and powers of three*, Zenodo.
  https://doi.org/10.5281/zenodo.20545164
- Earlier structural notes (2023), Figshare:
  https://doi.org/10.6084/m9.figshare.24778641 and
  https://doi.org/10.6084/m9.figshare.24778536
- *A "perfect union" of the nth terms of convergence* (Quora note, 2023), archival copy at Figshare:
  https://doi.org/10.6084/m9.figshare.33389983

An interactive companion to `residue_graph.py` — the graph drawn in
bit-reversed (2-adic) order, plus the exponent spike train — is at:
https://michaelmross.github.io/residue-graph.html
