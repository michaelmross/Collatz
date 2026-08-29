#!/usr/bin/env python3
"""
residue_graph.py -- Construct the nondeterministic residue graph G_{K,M}
of "Spike Structures and 2-adic Transition Laws in the Accelerated Collatz
Map" (Definition in the Finite-State section) and report its strongly
connected components.

Usage:
    python3 residue_graph.py [K M]              single graph (default K=7, M=13)
    python3 residue_graph.py --sweep K M1 M2    table for M = M1..M2 at fixed K

For K=7 (max_a=10), M=13: 4096 odd residues mod 2^13; the paper reports a
single nontrivial SCC of size 4088. The sweep reported in the paper is
    python3 residue_graph.py --sweep 7 11 20
whose "outside" column follows (M-9)*2^(M-12) for M >= 13.
"""
import sys
sys.setrecursionlimit(10**6)

def v2(n): return (n & -n).bit_length() - 1

def build(K=7, M=13):
    max_a, mod = K + 3, 1 << M
    nodes = list(range(1, mod, 2))
    idx = {r: i for i, r in enumerate(nodes)}
    adj = [[] for _ in nodes]
    for r in nodes:
        t = (3 * r + 1) % mod
        if t == 0:          # v2(3x+1) >= M for the whole class; all exceed max_a
            continue
        a = v2(t)
        if a > max_a:
            continue
        base = ((3 * r + 1) >> a) % (mod >> a)   # target determined mod 2^(M-a)
        step = mod >> a
        for k in range(1 << a):                  # top a bits free
            adj[idx[r]].append(idx[(base + k * step) % mod])
    return nodes, adj

def sccs(adj):
    n = len(adj); order = []; seen = [False]*n
    for s in range(n):                            # iterative Kosaraju pass 1
        if seen[s]: continue
        stack = [(s, 0)]; seen[s] = True
        while stack:
            u, i = stack.pop()
            if i < len(adj[u]):
                stack.append((u, i+1))
                w = adj[u][i]
                if not seen[w]: seen[w] = True; stack.append((w, 0))
            else: order.append(u)
    radj = [[] for _ in range(n)]
    for u in range(n):
        for w in adj[u]: radj[w].append(u)
    comp = [-1]*n; c = 0
    for s in reversed(order):                     # pass 2
        if comp[s] != -1: continue
        stack = [s]; comp[s] = c
        while stack:
            u = stack.pop()
            for w in radj[u]:
                if comp[w] == -1: comp[w] = c; stack.append(w)
        c += 1
    from collections import Counter
    return Counter(comp)

def report(K, M):
    nodes, adj = build(K, M)
    sizes = sorted(sccs(adj).values(), reverse=True)
    return len(nodes), sum(map(len, adj)), sizes

def main(argv):
    if argv and argv[0] == "--sweep":
        if len(argv) != 4:
            sys.exit("usage: residue_graph.py --sweep K M1 M2")
        K, M1, M2 = map(int, argv[1:4])
        if M1 <= K + 3:
            sys.exit(f"need M > max_a = K+3 = {K+3}")
        print(f"{'M':>3} {'nodes':>9} {'edges':>11} {'giant SCC':>10} {'outside':>8}")
        for M in range(M1, M2 + 1):
            n, e, sizes = report(K, M)
            print(f"{M:>3} {n:>9,} {e:>11,} {sizes[0]:>10,} {n - sizes[0]:>8,}")
        return
    if len(argv) not in (0, 2):
        sys.exit("usage: residue_graph.py [K M]  or  residue_graph.py --sweep K M1 M2")
    K, M = (int(argv[0]), int(argv[1])) if argv else (7, 13)
    if M <= K + 3:
        sys.exit(f"need M > max_a = K+3 = {K+3}")
    n, e, sizes = report(K, M)
    print(f"G_{{K={K},M={M}}}: {n} nodes, {e} edges")
    print(f"SCC sizes (top 5): {sizes[:5]};  singletons: {sizes.count(1)}")

if __name__ == "__main__":
    main(sys.argv[1:])
