# Reconstructing Evolutionary Relationships

> Module 15: How sequences are turned into evolutionary trees — the principles of molecular evolution, tree terminology, distance-based (neighbor-joining), character-based (parsimony), maximum likelihood and Bayesian methods, tree visualization, molecular clocks, phylogenomics, tracking pathogen outbreaks, bootstrap support, and a complete worked case study.

---

## Table of Contents

1. [Principles of Molecular Evolution](#1-principles-of-molecular-evolution)
2. [Building Phylogenetic Trees — Concepts and Terminology](#2-building-phylogenetic-trees--concepts-and-terminology)
3. [Distance-Based Methods — Neighbor-Joining](#3-distance-based-methods--neighbor-joining)
4. [Character-Based Methods — Maximum Parsimony](#4-character-based-methods--maximum-parsimony)
5. [Maximum Likelihood and Bayesian Approaches](#5-maximum-likelihood-and-bayesian-approaches)
6. [Tree Visualization and Interpretation Tools](#6-tree-visualization-and-interpretation-tools)
7. [Molecular Clocks and Divergence Time Estimation](#7-molecular-clocks-and-divergence-time-estimation)
8. [Phylogenomics — Whole-Genome Evolutionary Analysis](#8-phylogenomics--whole-genome-evolutionary-analysis)
9. [Applications in Tracking Pathogen Evolution and Outbreaks](#9-applications-in-tracking-pathogen-evolution-and-outbreaks)
10. [Evaluating Tree Reliability — Bootstrap and Confidence Metrics](#10-evaluating-tree-reliability--bootstrap-and-confidence-metrics)
11. [Case Study — Reconstructing an Evolutionary Tree](#11-case-study--reconstructing-an-evolutionary-tree)
12. [Key Takeaways](#12-key-takeaways)
13. [What's Next](#13-whats-next)

---

## 1. Principles of Molecular Evolution

**Molecular evolution tracks changes in DNA, RNA, and proteins across generations.** Mutations accumulate over time, and they leave **patterns of similarity and difference** between sequences. Sequences therefore act like **historical records**: they preserve clues about shared ancestry.

The central rule of this module is simple:

- **Greater sequence similarity usually suggests a more recent common ancestor.**
- **More sequence differences often indicate a more distant evolutionary relationship.**

The slide's picture shows two pairs of sequences. The first pair, `ACCTGACCT` and `ACCTGACCG`, differ at only one position, so their common ancestor was recent. The second pair is very different, so their common ancestor lies far back. Everything else in the module is a method for turning that idea into a tree.

The whole workflow of the module in one picture:

```
homologous sequences
   -> multiple sequence alignment        (Module 12)
   -> choose a method and a model        (distance, parsimony, likelihood, Bayesian)
   -> tree
   -> support values                      (bootstrap, posterior probability)
   -> interpretation                      (clades, dates, outbreaks)
```

### Slow and fast evolving regions

Not every part of a genome changes at the same speed. The slide's closing line: **slow and fast evolving regions reveal deep versus recent evolutionary history.**

| Region | Rate of change | What it is good for |
|---|---|---|
| Strongly constrained regions, such as the core of essential genes | Slow | Relationships that go far back in time, because the signal has not been erased |
| Weakly constrained regions, such as many non-coding sites | Fast | Recent relationships, such as strains of one species, because they change often enough to tell close relatives apart |

A region that changes too slowly has no differences between close relatives, and a region that changes too fast has been overwritten many times between distant ones. The choice of region has to match the question.

### Counting differences honestly: p-distance and the Jukes-Cantor correction

The simplest measure of difference is the **p-distance**, the fraction of aligned positions that differ. It has a flaw. Over long times, the same position can mutate more than once (A to G, and later G to T, or even back to A). Only the final state is visible, so p-distance **undercounts** the true amount of change. With four possible bases, two completely unrelated random sequences still match at about a quarter of positions, so p cannot rise above about 0.75.

The **Jukes-Cantor (JC) correction** estimates the real number of substitutions per site, assuming all substitutions are equally likely: d = −(3/4) × ln(1 − 4p/3).

```python
import math

def p_distance(a, b):
    diffs = sum(1 for x, y in zip(a, b) if x != y)
    return diffs / len(a)

def jukes_cantor(p):
    """Expected substitutions per site, correcting for multiple hits."""
    return -0.75 * math.log(1 - 4 * p / 3)

s1, s2 = "ACCTGACCT", "ACCTGACCG"
p = p_distance(s1, s2)
print(round(p, 3), round(jukes_cantor(p), 3))      # 0.111 0.12

for p in (0.05, 0.10, 0.30, 0.50, 0.70):
    print(p, round(jukes_cantor(p), 3))
# 0.05 0.052
# 0.1 0.107
# 0.3 0.383
# 0.5 0.824
# 0.7 2.031
```

For closely related sequences the correction is tiny (p = 0.05 gives 0.052). For distant ones it is large (p = 0.5 gives 0.824, and p = 0.7 gives more than 2 substitutions per site). As p approaches 0.75 the corrected distance grows without limit, which is the signal that the sequences are **saturated**: so many changes have piled up that the history can no longer be recovered.

> **Key idea:** sequence differences are a record of shared history. More similarity means a more recent ancestor, and the right region and the right correction are needed to read that record.

---

## 2. Building Phylogenetic Trees — Concepts and Terminology

A **phylogenetic tree** shows evolutionary relationships among genes, proteins, species, or sequences. The slide's checklist is a good way to read any tree: **tips, branches, and whether it is rooted.**

| Term | Meaning |
|---|---|
| **Tips (leaves)** | The observed organisms or sequences being compared |
| **Internal nodes** | **Inferred** common ancestors, estimated from the data. Nobody observed them |
| **Branches** | Evolutionary lineages. **Branch length may indicate the amount of evolutionary change** |
| **Root** | The oldest point of the tree: the common ancestor of everything in it |
| **Rooted tree** | Shows **direction through time** and a common ancestor |
| **Unrooted tree** | Shows relationships among the tips but no direction in time |
| **Clade** | An ancestor together with **all** of its descendants |
| **Sister taxa** | Two tips (or groups) that share an ancestor not shared with anything else |
| **Outgroup** | A lineage known to branch off before the group of interest, used to place the root |

Two points are easy to get wrong:

- **Only the branching pattern matters, not the left-to-right order of the tips.** Every internal node can be rotated like a mobile without changing the tree. `((A,B),C)` and `(C,(B,A))` are the same tree.
- **Closeness on the page is not closeness in evolution.** Two tips side by side but joined by a node far back are not close relatives. Always trace back to the node that joins them.

The slide's final line: **read trees by checking tips, branches, and whether the tree is rooted or unrooted.**

### The Newick format

Trees are stored as text in **Newick** format. Parentheses show nesting and the numbers after colons are branch lengths. The tree below has tips A, B, C, D:

`((A:0.1,B:0.2):0.05,(C:0.15,D:0.1):0.07);`

A small parser reads this and answers structural questions:

```python
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0
    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [branch()]
            while s[pos] == ",":
                pos += 1
                kids.append(branch())
            pos += 1                                # skip ")"
            return {"name": None, "kids": kids}
        start = pos
        while s[pos] not in ":,)":
            pos += 1
        return {"name": s[start:pos], "kids": []}
    def branch():
        nonlocal pos
        n = node()
        length = 0.0
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",)":
                pos += 1
            length = float(s[start:pos])
        n["length"] = length
        return n
    return branch()

def leaves(n):
    return [n["name"]] if not n["kids"] else [x for k in n["kids"] for x in leaves(k)]

def internal_count(n):
    return 0 if not n["kids"] else 1 + sum(internal_count(k) for k in n["kids"])

def total_length(n):
    return n["length"] + sum(total_length(k) for k in n["kids"])

def path_length(n, a, b):
    """Distance between tips a and b through the tree."""
    def find(n, target, path):
        if not n["kids"]:
            return path + [n] if n["name"] == target else None
        for k in n["kids"]:
            r = find(k, target, path + [n])
            if r: return r
        return None
    pa, pb = find(n, a, []), find(n, b, [])
    i = 0
    while i < min(len(pa), len(pb)) and pa[i] is pb[i]:
        i += 1
    return sum(x["length"] for x in pa[i:]) + sum(x["length"] for x in pb[i:])

tree = parse_newick("((A:0.1,B:0.2):0.05,(C:0.15,D:0.1):0.07);")
print(leaves(tree), internal_count(tree), round(total_length(tree), 2))
# ['A', 'B', 'C', 'D'] 3 0.67
print(round(path_length(tree, "A", "B"), 2), round(path_length(tree, "A", "C"), 2))
# 0.3 0.37
```

Reading the output: the tree has 4 tips, 3 internal nodes (two ancestors of the pairs, plus the root), and a total of 0.67 substitutions per site of branch length. A and B are 0.1 + 0.2 = 0.3 apart, while A and C must travel up through two internal branches, 0.1 + 0.05 + 0.07 + 0.15 = 0.37. This tip-to-tip distance along the tree is what branch lengths mean.

### Why the number of possible trees matters

For n tips, the number of possible unrooted trees is (2n − 5)!! and rooted trees (2n − 3)!!, where `!!` multiplies every second number down to 1:

```python
def double_factorial(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r

for n in (4, 5, 10, 20):
    print(n, double_factorial(2 * n - 5), double_factorial(2 * n - 3))
# 4 3 15
# 5 15 105
# 10 2027025 34459425
# 20 221643095476699771875 8200794532637891559375
```

| Tips | Unrooted trees | Rooted trees |
|---|---|---|
| 4 | 3 | 15 |
| 5 | 15 | 105 |
| 10 | 2,027,025 | 34,459,425 |
| 20 | about 2.2 × 10²⁰ | about 8.2 × 10²¹ |

With only 20 sequences there are far too many trees to test one by one. This is the same kind of problem as in Module 12 (too many possible alignments) and Module 13 (too many possible assemblies), and the answer is again the same: clever algorithms and heuristic searches.

> **Key idea:** a tree is read from its branching pattern. Tips are observed, internal nodes are inferred, branch lengths show change, and the root gives direction.

---

## 3. Distance-Based Methods — Neighbor-Joining

**Neighbor-Joining (NJ)** is a fast, practical distance-based method. The slide's description of how it works:

1. It starts with a **matrix of pairwise evolutionary distances** (for example, the Jukes-Cantor distances of §1).
2. It asks **which taxa are closest overall**, not every character-by-character pattern.
3. It **repeatedly joins the best neighbouring taxa or clusters**.
4. **After each join, the distance matrix is recalculated** until the tree is complete.

### The steps

For n taxa and distances d(i, j):

1. For each taxon i, compute r(i), its total distance to all the others.
2. For every pair, compute **Q(i, j) = (n − 2) × d(i, j) − r(i) − r(j)**.
3. Join the pair with the **smallest Q**. This is not simply the closest pair: subtracting r corrects for taxa that are far from everything and would otherwise be joined wrongly.
4. Give the new node its branch lengths: d(i, u) = d(i, j)/2 + (r(i) − r(j)) / (2(n − 2)), and d(j, u) = d(i, j) − d(i, u).
5. Replace i and j by the new node u, with d(u, k) = (d(i, k) + d(j, k) − d(i, j)) / 2.
6. Repeat until three nodes are left, and join them at a central point.

### A worked example

A distance matrix for five taxa:

| | a | b | c | d | e |
|---|---|---|---|---|---|
| **a** | 0 | 5 | 9 | 9 | 8 |
| **b** | 5 | 0 | 10 | 10 | 9 |
| **c** | 9 | 10 | 0 | 8 | 7 |
| **d** | 9 | 10 | 8 | 0 | 3 |
| **e** | 8 | 9 | 7 | 3 | 0 |

The row sums are r(a) = 31, r(b) = 34, r(c) = 34, r(d) = 30, r(e) = 27. With n = 5, the first-step Q values are:

| Pair | Q |
|---|---|
| a, b | **−50** |
| d, e | −48 |
| c, d | −40 |
| c, e | −40 |
| a, c | −38 |
| b, c | −38 |
| a, d / a, e / b, d / b, e | −34 |

The smallest Q is (a, b), so they are joined first. The code carries out all the steps:

```python
import itertools

def neighbor_joining(names, dist):
    """names: list of taxa; dist: dict {(a, b): distance} for every pair (either order)."""
    d = {}
    for (a, b), v in dist.items():
        d[(a, b)] = d[(b, a)] = v
    label = {n: n for n in names}               # node -> Newick fragment
    active = list(names)
    count = 0
    while len(active) > 3:
        n = len(active)
        r = {i: sum(d[(i, j)] for j in active if j != i) for i in active}
        q = {(i, j): (n - 2) * d[(i, j)] - r[i] - r[j]
             for i, j in itertools.combinations(active, 2)}
        i, j = min(q, key=q.get)                # pair with the smallest Q
        li = d[(i, j)] / 2 + (r[i] - r[j]) / (2 * (n - 2))
        lj = d[(i, j)] - li
        count += 1
        u = f"u{count}"
        label[u] = f"({label[i]}:{li:.3f},{label[j]}:{lj:.3f})"
        for k in active:
            if k not in (i, j):
                d[(u, k)] = d[(k, u)] = (d[(i, k)] + d[(j, k)] - d[(i, j)]) / 2
        active = [k for k in active if k not in (i, j)] + [u]
    x, y, z = active                            # three nodes left: join at a centre
    lx = (d[(x, y)] + d[(x, z)] - d[(y, z)]) / 2
    ly = (d[(x, y)] + d[(y, z)] - d[(x, z)]) / 2
    lz = (d[(x, z)] + d[(y, z)] - d[(x, y)]) / 2
    return f"({label[x]}:{lx:.3f},{label[y]}:{ly:.3f},{label[z]}:{lz:.3f});"

names = ["a", "b", "c", "d", "e"]
dist = {("a","b"): 5, ("a","c"): 9, ("a","d"): 9, ("a","e"): 8, ("b","c"): 10,
        ("b","d"): 10, ("b","e"): 9, ("c","d"): 8, ("c","e"): 7, ("d","e"): 3}
print(neighbor_joining(names, dist))
# (d:2.000,e:1.000,(c:4.000,(a:2.000,b:3.000):3.000):2.000);
```

The result: a and b join with branches of 2 and 3, d and e join with branches of 2 and 1, and c sits between the two groups on a branch of length 4. Rebuild the distances from the tree to check: a to b is 2 + 3 = 5, and a to d is 2 + 3 + 2 + 2 = 9. Both match the matrix, which is expected here because this matrix is exactly tree-like.

### Speed against detail

The slide's picture is a scale with **speed** on one side and **detail** on the other, and its summary: **NJ is fast but may hide evolutionary detail, so it works well as a first pass.**

- **Why it is fast:** it needs only one number per pair of sequences, and it builds a single tree directly, with no search through the huge number of possible trees.
- **Why it hides detail:** each pair of sequences is squeezed into one number, so which sites differ and how is thrown away. Two very different patterns of change can give the same distance.

Unlike a simpler clustering method (UPGMA), NJ does not assume that all lineages evolve at the same rate, so it is more reliable on real data.

> **Key idea:** NJ builds a tree quickly from a distance matrix by repeatedly joining the best neighbours. It is an excellent first look, and it is not the final word.

---

## 4. Character-Based Methods — Maximum Parsimony

**Maximum parsimony** is a **character-based** alternative to distance-based methods. It **evaluates each character directly**, where a character is, for example, each nucleotide position.

The key question: **which tree explains the data with the fewest changes?** Parsimony **favours the tree needing the fewest substitutions overall**, and candidate trees are scored by **counting changes across taxa and characters**.

### Counting changes on a tree: Fitch's method

For each site, work from the tips toward the root. At each internal node, combine the sets of states from its two children:

- if the two sets share a state, keep the shared states (no change needed);
- if they share nothing, take the union and count **one change**.

### A worked example

Four taxa and eight sites. For four taxa there are three possible unrooted trees: `((A,B),(C,D))`, `((A,C),(B,D))`, and `((A,D),(B,C))`.

```python
def fitch_site(tree, states):
    """tree: nested tuples of taxon names. Returns (possible states at root, changes)."""
    if isinstance(tree, str):
        return {states[tree]}, 0
    (s1, c1), (s2, c2) = fitch_site(tree[0], states), fitch_site(tree[1], states)
    common = s1 & s2
    if common:
        return common, c1 + c2
    return s1 | s2, c1 + c2 + 1                  # no overlap: one change is needed

def parsimony_score(tree, alignment):
    n_sites = len(next(iter(alignment.values())))
    return sum(fitch_site(tree, {t: seq[i] for t, seq in alignment.items()})[1]
               for i in range(n_sites))

alignment = {
    "A": "ACATGAAC",
    "B": "ACATAAAC",
    "C": "GTACGAAC",
    "D": "GTGCAAAC",
}
trees = {
    "((A,B),(C,D))": (("A", "B"), ("C", "D")),
    "((A,C),(B,D))": (("A", "C"), ("B", "D")),
    "((A,D),(B,C))": (("A", "D"), ("B", "C")),
}
for name, t in trees.items():
    print(name, parsimony_score(t, alignment))
# ((A,B),(C,D)) 6
# ((A,C),(B,D)) 8
# ((A,D),(B,C)) 9
```

Site by site, for the winning tree `((A,B),(C,D))`:

| Site | A B C D | Changes | Comment |
|---|---|---|---|
| 1 | A A G G | 1 | Supports AB with CD |
| 2 | C C T T | 1 | Supports AB with CD |
| 3 | A A A G | 1 | Only D differs: fits every tree equally |
| 4 | T T C C | 1 | Supports AB with CD |
| 5 | G A G A | 2 | Supports AC with BD, against the winning tree |
| 6 to 8 | constant | 0 | No information |

The totals are 6, 8, and 9 changes, so parsimony selects `((A,B),(C,D))`. Only some sites actually help choose between trees: a **parsimony-informative site** needs at least two states, each occurring in at least two taxa. Here those are sites 1, 2, 4, and 5. Site 3 (one odd taxon) and sites 6-8 (no variation) cannot favour any tree.

The slide's conclusion: **character patterns may reveal relationships hidden by distance summaries**, because the pattern at each site is used, not just a total difference.

### Limits

- Parsimony has **no explicit model**: it treats every change as equally likely and counts it once, even on a very long branch where many hidden changes probably happened.
- It can be misled by **long-branch attraction**: two rapidly evolving lineages share many chance similarities and get grouped together, though they are not close relatives.
- Several trees can tie on the same score, so a single "best" tree may not exist.

These limits motivate the methods of the next section.

> **Key idea:** parsimony looks at each site and chooses the tree that needs the fewest evolutionary changes. It uses more information than distances, but it assumes change is rare and simple.

---

## 5. Maximum Likelihood and Bayesian Approaches

Maximum likelihood (ML) and Bayesian inference **improve on maximum parsimony**.

### Maximum likelihood

ML asks: **which tree makes the observed sequence data most probable?** It uses **an explicit evolutionary model, not just a count of changes.** A model states how likely each kind of substitution is over a given amount of time. The simplest model is Jukes-Cantor (§1). Richer ones separate common substitution types from rare ones, or allow some sites to evolve faster than others, and the model is chosen by comparing how well each fits the data.

This is why the slide notes that **trees with similar change counts may differ once substitution patterns are modelled**: two trees that need the same number of changes can still have different probabilities.

**ML estimates a branch length, too.** For two sequences of length n that differ at k sites, the likelihood of a branch length t under Jukes-Cantor can be computed directly. Searching for the t with the highest likelihood gives the same answer as the Jukes-Cantor formula of §1:

```python
import math

def p_same(t): return 0.25 + 0.75 * math.exp(-4 * t / 3)
def p_diff(t): return 0.25 - 0.25 * math.exp(-4 * t / 3)       # one specific other base

def loglik_pair(n, k, t):
    """Two sequences of length n that differ at k sites, separated by distance t (JC69)."""
    return n * math.log(0.25) + (n - k) * math.log(p_same(t)) + k * math.log(p_diff(t))

n, k = 100, 10
best_t = max((i / 10000 for i in range(1, 20000)), key=lambda t: loglik_pair(n, k, t))
print(best_t, round(-0.75 * math.log(1 - 4 * (k / n) / 3), 4))     # 0.1073 0.1073
```

With 10 differences in 100 sites, the most likely distance is 0.1073, the same as the corrected distance, so the distance of §1 is the maximum-likelihood estimate under that model.

### Comparing trees by likelihood

For a whole tree, the likelihood of each site is found by moving from the tips to the root and summing over the unknown ancestral states (Felsenstein's pruning algorithm). Using the eight-site alignment of §4, with every branch fixed at 0.1 substitutions per site:

```python
BASES = "ACGT"

def node_like(node, site):
    """node: a taxon name, or ((child, length), (child, length))."""
    if isinstance(node, str):
        return {b: 1.0 if b == site[node] else 0.0 for b in BASES}
    like = {b: 1.0 for b in BASES}
    for child, length in node:
        cl = node_like(child, site)
        for b in BASES:
            like[b] *= sum((p_same(length) if b == c else p_diff(length)) * cl[c] for c in BASES)
    return like

def with_lengths(topology, t):
    """Give every external branch length t, and split t across the middle branch."""
    (x, y), (z, w) = topology
    return ((((x, t), (y, t)), t / 2), (((z, t), (w, t)), t / 2))

def tree_loglik(topology, alignment, t):
    root = with_lengths(topology, t)
    total = 0.0
    for i in range(len(next(iter(alignment.values())))):
        site = {tax: seq[i] for tax, seq in alignment.items()}
        like = node_like(root, site)
        total += math.log(sum(0.25 * like[b] for b in BASES))
    return total

logl = {name: tree_loglik(top, alignment, 0.1) for name, top in trees.items()}
for name, v in logl.items():
    print(name, round(v, 2))
# ((A,B),(C,D)) -34.25
# ((A,C),(B,D)) -39.67
# ((A,D),(B,C)) -42.37
```

The tree with the highest log-likelihood wins: `((A,B),(C,D))` at −34.25, which here agrees with parsimony. The gap to the next tree, 5.4 log units, means the data are more than 200 times more probable on the winning tree (e^5.4 is about 226). (This is a toy: real programs estimate the branch lengths and model parameters too, instead of fixing them.)

### Bayesian inference

Bayesian inference **combines the data with prior assumptions to estimate tree probabilities**: posterior ∝ prior × likelihood. The slide's picture shows three candidate trees with probabilities (0.2, 0.5, 0.3) and a "Prior + Data" box leading to a "Posterior" distribution.

```python
def posterior(loglik, prior):
    w = {name: prior[name] * math.exp(v) for name, v in loglik.items()}
    s = sum(w.values())
    return {name: round(x / s, 3) for name, x in w.items()}

flat = {name: 1 / 3 for name in trees}
print(posterior(logl, flat))
# {'((A,B),(C,D))': 0.995, '((A,C),(B,D))': 0.004, '((A,D),(B,C))': 0.0}

skeptical = {"((A,B),(C,D))": 0.01, "((A,C),(B,D))": 0.495, "((A,D),(B,C))": 0.495}
print(posterior(logl, skeptical))
# {'((A,B),(C,D))': 0.81, '((A,C),(B,D))': 0.178, '((A,D),(B,C))': 0.012}
```

With no preference among the trees (the "flat" prior), the posterior probability of the winning tree is 0.995. With a prior that starts out strongly against it (0.01), the data still lift it to 0.81. A prior matters most when the data are weak, and a strong signal can overcome a doubtful prior.

The slide's closing line: **Bayesian methods show uncertainty through a distribution of likely trees.** Real programs explore the enormous space of trees with a random-walk procedure (MCMC) and report how often each clade appears across the sampled trees. A reader then sees not one tree but a range of plausible ones, and which parts of the tree are stable.

### Summary of the three families

| | Distance (NJ) | Parsimony | ML and Bayesian |
|---|---|---|---|
| **Uses** | One distance per pair | Each site, counting changes | Each site, with an explicit probability model |
| **Output** | One tree | One or several best trees | A best tree (ML), or a distribution of trees (Bayesian) |
| **Speed** | Fastest | Moderate | Slowest |
| **Weak point** | Hides site-level detail | No model, long-branch attraction | Depends on the model, computationally heavy |

> **Key idea:** ML scores trees by how probable the data are under an explicit model, and Bayesian inference adds prior assumptions and reports uncertainty as a distribution of trees.

---

## 6. Tree Visualization and Interpretation Tools

**Visualization turns raw trees into readable biological evidence.** A tree file is just nested text, and a viewer turns it into a figure. Common viewers include FigTree, iTOL, and MEGA.

### What a good tree figure shows

- **Branch lengths**, with a scale bar so that distance means something.
- **Support values** at the nodes (§10), so strong and weak groupings are told apart.
- **Labels** that carry useful information (taxon, host, place, date).
- **A useful layout**: rectangular for reading relationships, circular for large trees, with clades collapsed or coloured where that makes the pattern clearer.

The slide's rule: **a good tree should reveal signal, not hide it behind complexity**, and **clear visualization exposes patterns; poor design can bury them.** As the lecture's narration puts it, the goal is not just to make the tree look neat enough for a publication figure.

### Interpretation: reading the clades

The slide's board shows clades A, B, and C and an **unusual group**. **Interpretation tools compare clades and flag unexpected groupings.** An unexpected grouping might be a real discovery, but it has three different possible explanations, shown at the bottom of the slide:

| Possible explanation | What to check |
|---|---|
| **Biology** | Real evolutionary events, such as gene transfer between distant species, or very similar traits evolving independently |
| **Sampling** | Missing taxa, uneven sampling of lineages, or a mislabelled sample |
| **Analysis error** | A poor alignment, a wrong model, contamination, or long-branch attraction |

The last line of the slide: **strongly supported clustering increases confidence in relationships.** An unexpected grouping with strong support from several methods deserves attention as biology. An unexpected grouping with weak support is more likely to be an artefact.

### A reading checklist

1. Is the tree rooted, and where is the root?
2. What does the scale bar mean (substitutions per site, or time)?
3. Which clades are strongly supported, and which are not?
4. Are the groupings consistent with known biology?
5. Could sampling, alignment, or method explain anything surprising?

> **Key idea:** visualization is the step where a tree becomes evidence. Show the branch lengths and the support, and then check any surprising pattern before believing it.

---

## 7. Molecular Clocks and Divergence Time Estimation

A tree shows who is related to whom. **Molecular clocks ask when lineages split, not just who is related.** The slide's image is a stopwatch: **sequence differences can act like a biological stopwatch**, because **more sequence differences often mean an older common ancestor.**

### The calculation

If two lineages separated T years ago and each accumulated substitutions at a rate r (substitutions per site per year), the distance between them is d = 2 × r × T, since both lineages have been changing. So **T = d / (2r)**.

Distance on its own gives only relative time. **Calibration points convert relative differences into real time.** A calibration is a split whose date is known independently, for example from a fossil. It lets us solve for the rate, and then date other splits:

```python
import math

# calibrate with a fossil-dated split, then date another split
d_cal, t_cal = 0.08, 20e6                        # distance 0.08, split dated to 20 million years
rate = d_cal / (2 * t_cal)                       # substitutions per site per year, per lineage
print(f"{rate:.1e}")                             # 2.0e-09
d_new = 0.05
print(round(d_new / (2 * rate) / 1e6, 1))        # 12.5 million years
```

The calibration gives a rate of 2 × 10⁻⁹ substitutions per site per year, and a pair with distance 0.05 would then have split about 12.5 million years ago.

### The same idea on a fast-evolving virus

For a virus, the clock runs in days rather than millions of years. Suppose two isolates of a virus with a genome of 29,903 bases differ at 3 sites, and the virus accumulates about 8 × 10⁻⁴ substitutions per site per year:

```python
genome, snps, rate_virus = 29903, 3, 8e-4        # rate: substitutions per site per year
d = snps / genome
t_years = d / (2 * rate_virus)
print(round(t_years * 365), "days")              # 23 days
print(round(genome * rate_virus, 1))             # 23.9 substitutions per lineage per year
```

That rate corresponds to about 24 changes per lineage per year, so 3 differences between two isolates suggest a common ancestor roughly 23 days earlier. The estimate is rough: a count of 3 could easily have been 1 or 6 by chance, which would put the answer anywhere between about 8 and 46 days.

### The big caveat: rates are not constant

The slide's warning: **evolutionary rates vary across genes and lineages.** A **strict clock** assumes one rate for everything, and a **relaxed clock** lets the rate differ between branches. Because real rates do vary, relaxed clocks are generally safer, and calibration is the weakest link: a wrong calibration shifts every date.

The slide's last line: **divergence dating links evolution, pathogens, and Earth history.** The same principle dates the origin of a group of species, the start of an outbreak, and the order of events on a geological timeline.

> **Key idea:** a molecular clock turns genetic distance into time. It needs a calibration, and its dates are only as good as the assumption about the rate.

---

## 8. Phylogenomics — Whole-Genome Evolutionary Analysis

**Phylogenomics expands analysis from a few genes to entire genomes.** It **uses genome-wide data instead of a single marker for richer evidence.** The slide's analogy: judging a full series, not one episode.

### What it adds

**Comparing many genes reveals species relationships and genome change.** Beyond sequence substitutions, whole genomes record large-scale events, which the slide's picture shows:

| Event | What happens |
|---|---|
| **Gene loss** | A gene present in an ancestor is absent in a descendant |
| **Gene gain and duplication** | A new gene appears, or an existing gene is copied |
| **Rearrangement** | Blocks of genes change their order or orientation |

### Why one gene can be misleading

**It is useful when single genes conflict; whole genomes give a broader view.** Different genes can tell different stories because:

- different genes have different histories, as ancestral variation sorts differently into descendants;
- genes can move between species (horizontal transfer);
- genes may be paralogs (copies made by duplication) rather than orthologs (versions of the same gene in different species), and a tree of paralogs mixes duplication with speciation;
- or one gene may simply have a weak signal or an error.

### From many gene trees to one species tree

A simple way to see the principle is to build a tree for each gene, and count the answers:

```python
from collections import Counter

gene_trees = ["AB|CD"] * 6 + ["AC|BD"] * 3 + ["AD|BC"]
counts = Counter(gene_trees)
print(counts.most_common())                      # [('AB|CD', 6), ('AC|BD', 3), ('AD|BC', 1)]
top, n = counts.most_common(1)[0]
print(top, f"{100 * n / len(gene_trees):.0f}%")  # AB|CD 60%
```

Of ten genes, six support `AB|CD`, three support `AC|BD`, and one supports `AD|BC`. Any single gene would have given a confident-looking tree, but if one gene were picked at random there would be a 40% chance of getting a tree that differs from the majority. The majority view is only a starting point, because real methods either combine all genes into one large alignment, or combine the separate gene trees with models that expect some disagreement. And disagreement between genes can be real biology, not just noise.

> **Key idea:** a single gene is one episode, and a genome is the whole series. Phylogenomics gives a broader and more reliable picture, and shows events such as gene loss, duplication, and rearrangement.

---

## 9. Applications in Tracking Pathogen Evolution and Outbreaks

**Phylogenetic trees compare pathogen genomes to show how strains are related over time.** The slide's central point: **outbreaks are evolutionary processes unfolding in real time, not just medical events.** A pathogen mutates as it spreads, and those mutations are the record.

### How to read the tree

| Pattern | Interpretation |
|---|---|
| **Few mutations between isolates** | May suggest **linked transmission** or a **shared source** |
| **Distinct branches** | Can indicate **multiple introductions** into the outbreak |
| **New mutations in a lineage** | Act as **biological clues** to trace spread and **detect emerging variants** |

### An example with SNP distances

Six isolates of 20 bases each (a toy), with the number of single-base differences (SNPs) between each pair, and a rule that isolates within a chosen number of SNPs of each other belong to the same cluster:

```python
import itertools

isolates = {
    "P1": "ACGTACGTACGTACGTACGT",
    "P2": "ACGTACGTACGTACGAACGT",
    "P3": "ACGTACGTACCTACGAACGT",
    "P4": "ACGTTCGTACGTACGTACGT",
    "P5": "TCGAACGTTCGTACCTACGT",
    "P6": "TCGAACGTTCGTACCTACGA",
}

def snp_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

def clusters(isolates, threshold):
    names = sorted(isolates)
    parent = {n: n for n in names}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in itertools.combinations(names, 2):
        if snp_distance(isolates[a], isolates[b]) <= threshold:
            parent[find(a)] = find(b)
    groups = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return sorted(groups.values())

print(clusters(isolates, 2))    # [['P1', 'P2', 'P3', 'P4'], ['P5', 'P6']]
print(clusters(isolates, 4))    # [['P1', 'P2', 'P3', 'P4', 'P5', 'P6']]
```

At a threshold of 2 SNPs there are two clusters: P1 to P4, which differ by only 1 to 3 SNPs, and P5 and P6, which differ from P1 by 4 to 5 and from each other by 1. That pattern fits **two separate introductions** into the population, followed by local spread within each group. (P3 and P4 are 3 SNPs apart, so they are in the same cluster only because both are within 2 SNPs of P1. This chaining is called single linkage.) At a threshold of 4, P5 is linked to P1 and the clusters merge. So the threshold is an assumption, and it has to be matched to how fast the pathogen mutates and how much time has passed.

### What the tree cannot say

- A phylogeny shows relatedness, not who infected whom. Two cases with identical genomes may be linked directly, or through unsampled people.
- Sampling is rarely complete, and missing samples can create false gaps.
- A pathogen's rate of mutation sets the scale: a handful of changes means days for one virus and years for another (§7).

### Why it matters

The slide's last line: **phylogenetics supports surveillance, containment, and resource decisions in epidemics.** Knowing whether cases come from one source or many, and whether a new variant is spreading, decides where to test, what to contain, and where to send resources.

> **Key idea:** mutations are a trace of spread. Few differences point to linked cases, distinct branches point to separate introductions, and both depend on honest assumptions.

---

## 10. Evaluating Tree Reliability — Bootstrap and Confidence Metrics

The slide begins with a reminder: **phylogenetic trees are estimates, not fixed facts.** And: **small data changes can alter branches, especially with weak signal.** So for every branch the question is: **which branches are truly supported by the data?**

### The bootstrap

The **bootstrap resamples the data and rebuilds trees many times**:

1. From an alignment with N sites, draw N sites **with replacement**: some sites appear twice or more, others not at all.
2. Build a tree from this resampled alignment.
3. Repeat many times (100 to 1,000 or more).
4. For each branch of the original tree, count the fraction of resampled trees that contain it.

**Branches seen repeatedly get higher support values.** The slide's picture shows a branch with support 98 and one with support 55.

### A worked example

A 40-site alignment of four taxa (the sequences were generated by simulation on a tree where A and B are sisters), analysed with parsimony. Ties among trees share the credit equally:

```python
import random

alignment40 = {
    "A": "ATAGAATAGCCCTTAGATCTCACTGGGGCGCGTGAATCTT",
    "B": "ATAGAATAGCCCATAAACCTCACGGGCGCGCGTGAATCTT",
    "C": "ATAGAATAGCCCATAGACACCACGTGCTCGCGTGAATCAT",
    "D": "ATAGCATAACCCATAAACATCACTTGTTCGCGTGAATCAT",
}

def bootstrap_support(alignment, replicates, seed):
    rng = random.Random(seed)
    n = len(next(iter(alignment.values())))
    support = {name: 0.0 for name in trees}
    for _ in range(replicates):
        cols = [int(rng.random() * n) for _ in range(n)]            # sample sites with replacement
        resampled = {t: "".join(seq[c] for c in cols) for t, seq in alignment.items()}
        scores = {name: parsimony_score(top, resampled) for name, top in trees.items()}
        best = min(scores.values())
        winners = [name for name, s in scores.items() if s == best]
        for name in winners:
            support[name] += 1 / len(winners)                        # ties share the credit
    return {name: round(100 * v / replicates) for name, v in support.items()}

print({k: parsimony_score(v, alignment40) for k, v in trees.items()})
# {'((A,B),(C,D))': 15, '((A,C),(B,D))': 18, '((A,D),(B,C))': 18}
print(bootstrap_support(alignment40, 100, 1))
# {'((A,B),(C,D))': 92, '((A,C),(B,D))': 5, '((A,D),(B,C))': 3}
```

The best tree needs 15 changes against 18 for the other two, and the grouping of A with B appears in 92% of resampled datasets. Now the same experiment with only the first part of the alignment:

| Sites used | Changes on the trees that group A with B, A with C, A with D | Bootstrap support (%) for A+B, A+C, A+D |
|---|---|---|
| 12 | 2, 2, 2 | 33, 33, 33 |
| 24 | 10, 10, 10 | 26, 35, 39 |
| 30 | 14, 16, 16 | 77, 9, 13 |
| 40 | 15, 18, 18 | 92, 5, 3 |

With 12 sites the three trees tie exactly, and support is spread evenly, so the data say nothing. With 24 sites the grouping the data truly came from is not even the top-supported one. Support tends to rise as the amount of data grows (it does not rise smoothly, as the dip between 12 and 24 sites shows), which is what **small data changes can alter branches, especially with weak signal** looks like in practice.

### Reading support values

- **High values mean stronger support; low values need caution.** A common working guide: 70% or more is moderate support and 90 to 95% or more is strong, but there is no universal threshold.
- A bootstrap value is **not the probability that the branch is true.** It measures how consistently the data support the branch when they are resampled.
- The numbers themselves are noisy: the same clade above gets 86% instead of 92% when 1,000 replicates are used (the other two trees get 6% and 7%), so differences of a few percent are not meaningful.
- **High support is not proof of correctness.** A wrong model, a poor alignment, or long-branch attraction can produce a consistently wrong tree with high support. Support measures sampling consistency, not systematic error.

### Other confidence measures

For Bayesian trees, the **posterior probability** of a clade is its frequency among the sampled trees (§5). Likelihood-based programs offer faster alternatives to the standard bootstrap. The principle is the same in every case: **even evolutionary trees need receipts.**

> **Key idea:** a tree is an estimate. Bootstrap and similar measures show which branches the data support consistently, and low support is a warning, not a result.

---

## 11. Case Study — Reconstructing an Evolutionary Tree

The slide's case study walks through the whole process. **Start with sequences from organisms or genes you want to compare. The goal is to infer relatedness and represent it as a branching tree.** And the central warning: **poor sequence choice or alignment can produce misleading trees.**

### The steps

| Step | What to do | Where it was covered |
|---|---|---|
| 1. **Choose sequences** | **Choose sequences that fit the biological question being asked** | §1 |
| 2. **Collect them** | Retrieve the sequences from databases, with accessions recorded | Module 11 |
| 3. **Align** | **Align positions to reflect shared evolutionary history** | Module 12 |
| 4. **Clean the alignment** | Remove ambiguous or badly aligned regions | Module 12 |
| 5. **Build the tree** | NJ for a quick first look; ML or Bayesian for the final tree | §3 to §5 |
| 6. **Test support** | Bootstrap or posterior probabilities | §10 |
| 7. **Root and display** | Root with an outgroup, then draw with branch lengths and support | §2, §6 |
| 8. **Interpret** | Read clades, check unexpected groupings, estimate dates if needed | §6, §7 |

### Choosing sequences

- All sequences must be **homologous** and ideally **orthologous**, the same gene in each organism.
- The rate of change must suit the question: slow-evolving regions for deep relationships, fast-evolving ones for close relatives (§1).
- Include an **outgroup**, a lineage that branched off before the others, to place the root.
- Check for contamination and for mislabelled samples before anything else.

### A complete miniature example

Six sequences of 30 bases (a toy), where `Out` is the outgroup. The code computes Jukes-Cantor distances for every pair, builds the tree with the neighbor-joining function of §3, and prints it:

```python
alignment = {
    "Out": "GCATTGCGAGTTTCAAGCTAACTGGGACCT",
    "S1": "ATGTCGCAAACTTTGAACTAGTTAGGATCC",
    "S2": "ATGTCGTAAACTTTGAACTAGCTAGGATCC",
    "S3": "ATGTCGCAAACCTCGAACTAGCTAGGATCC",
    "S4": "ACGTTACAGGCTCCGAGCCAGCTAGAATCC",
    "S5": "ACGTTACAAGCTCCGAGCCAGCTAGAATTC",
}
names = list(alignment)
dist = {}
for a, b in itertools.combinations(names, 2):
    dist[(a, b)] = jukes_cantor(p_distance(alignment[a], alignment[b]))

print(neighbor_joining(names, dist))
# (S3:0.015,(Out:0.520,(S4:0.035,S5:0.035):0.175):0.206,(S1:0.035,S2:0.035):0.058);
```

A few of the distances that drive it:

| Pair | JC distance | Reading |
|---|---|---|
| S1 and S2 | 0.070 | Very close |
| S1 and S3 | 0.107 | Close |
| S4 and S5 | 0.070 | Very close |
| S1 and S4 | 0.503 | Far |
| Out and every S | 0.73 to 0.82 | Furthest of all |

### Reading the result

The printed tree is unrooted. It joins three groups at the centre: `S3`, the pair `(S1,S2)`, and a branch that leads to `Out` and the pair `(S4,S5)`. To root it, use the outgroup: place the root on the branch that leads to `Out`. The rooted tree reads `(Out, ((S4,S5), (S3,(S1,S2))))`, and its relationships are:

- S1, S2, and S3 form one clade, with S1 and S2 as closest relatives within it.
- S4 and S5 form a second clade, a sister to the first.
- `Out` is the outgroup, which branched off before both.

The outgroup has the longest branch (0.52), because it differs most from everything else. This is consistent with the distances in the table.

### What the tree is, and is not

The slide's final line: **trees aid species, gene, and pathogen studies, but remain testable hypotheses.** The tree above is one reasonable hypothesis from one method and 30 sites. A fuller study would repeat it with likelihood-based inference, check support (§10), add more genes (§8), and test the result against other evidence. In the lecture's words, the tree is a hypothesis: useful and powerful, and still something to question thoughtfully.

> **Key idea:** the quality of a tree starts with the choice of sequences and the alignment. Build it, test its support, and treat the result as a hypothesis.

---

## 12. Key Takeaways

- Molecular evolution leaves a record in sequences: more similarity suggests a more recent ancestor, more difference a more distant one. Slow and fast regions inform deep and recent history. Raw differences must be corrected for multiple hits (Jukes-Cantor), and they saturate.
- In a tree, tips are observed, internal nodes are inferred, branch lengths show change, and the root gives direction in time. Read trees by their branching pattern, not the order of the tips.
- The number of possible trees explodes with the number of tips (about 2.0 million unrooted trees for 10 tips), so algorithms must search cleverly.
- Neighbor-joining builds a tree from a distance matrix by repeatedly joining the best neighbours. It is fast and a good first pass, but it hides site-level detail.
- Maximum parsimony chooses the tree with the fewest changes across all sites. It uses site patterns but has no explicit model, and can be misled by long branches.
- Maximum likelihood scores trees by the probability of the data under an explicit model, and Bayesian inference adds priors and reports a distribution of trees.
- A good tree figure reveals signal, shows branch lengths and support, and prompts checks of unexpected groupings for biology, sampling, or analysis error.
- Molecular clocks turn distance into time, T = d / (2r), but need calibration and depend on the assumption about rate.
- Phylogenomics uses whole genomes, which gives a broader view when single genes conflict and reveals gene loss, duplication, and rearrangement.
- In outbreaks, few mutations suggest linked transmission and distinct branches suggest multiple introductions, and trees support surveillance and containment decisions.
- Bootstrap and posterior probabilities show which branches the data support. Support grows with data, high values mean stronger support, and low values need caution. Neither proves a branch is correct.
- A tree is a testable hypothesis. Good sequence choice and a good alignment come first.

---

## 13. What's Next

With sequences, alignments, genomes, structures, and evolutionary trees in hand, each of the tools so far describes one molecule or one history. The next step is to combine these layers and ask how they work together to explain what happens inside a living cell.
