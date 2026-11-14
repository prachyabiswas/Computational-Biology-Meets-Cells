# Understanding Protein Structure and Function

> Module 14: How a protein's amino acid sequence becomes a three-dimensional shape and a biological job — the four levels of structure, X-ray crystallography, NMR and cryo-EM, structure visualization, sequence-structure-function relationships, domains and motifs, homology modeling, protein-protein interaction networks, molecular docking and drug targets, AI structure prediction with AlphaFold, structural bioinformatics in drug discovery, and how to judge a predicted structure.

---

## Table of Contents

1. [Levels of Protein Structure — Primary to Quaternary](#1-levels-of-protein-structure--primary-to-quaternary)
2. [Methods for Determining Protein Structure — X-ray, NMR, Cryo-EM](#2-methods-for-determining-protein-structure--x-ray-nmr-cryo-em)
3. [Structure Visualization Tools — PyMOL and Chimera](#3-structure-visualization-tools--pymol-and-chimera)
4. [Sequence-Structure-Function Relationships](#4-sequence-structure-function-relationships)
5. [Protein Domain and Motif Analysis](#5-protein-domain-and-motif-analysis)
6. [Homology Modeling Principles](#6-homology-modeling-principles)
7. [Protein-Protein Interaction Networks](#7-protein-protein-interaction-networks)
8. [Molecular Docking and Drug Target Analysis](#8-molecular-docking-and-drug-target-analysis)
9. [AI-Driven Structure Prediction — Understanding Tools Like AlphaFold](#9-ai-driven-structure-prediction--understanding-tools-like-alphafold)
10. [Structural Bioinformatics in Drug Discovery](#10-structural-bioinformatics-in-drug-discovery)
11. [Interpreting and Validating Predicted Structures](#11-interpreting-and-validating-predicted-structures)
12. [Key Takeaways](#12-key-takeaways)
13. [What's Next](#13-whats-next)

---

## 1. Levels of Protein Structure — Primary to Quaternary

Earlier modules worked with sequences. This module asks what a sequence does once it becomes a molecule. The idea that runs through all eleven sections is one chain of reasoning:

```
amino acid sequence
   -> structure          (how the chain folds in three dimensions)
   -> function           (what that shape lets the protein do)
   -> interactions       (which other proteins and molecules it works with)
   -> application        (drug targets and drug candidates)
```

A protein's structure is described at four levels, and **each level builds on the previous one**.

| Level | What it is | What holds it together |
|---|---|---|
| **Primary** | The amino acid sequence, read from the N-terminus to the C-terminus | Peptide bonds between neighbouring amino acids |
| **Secondary** | Local repeating shapes: **alpha helices** and **beta sheets** | Hydrogen bonds between backbone atoms |
| **Tertiary** | The complete three-dimensional shape of **one polypeptide chain** | Interactions between side chains: a water-avoiding hydrophobic core, hydrogen bonds, salt bridges, and in some proteins disulfide bonds |
| **Quaternary** | Several chains assembled into one functional protein | The same kinds of interactions, now between different chains |

### Secondary structure in numbers

- In an **alpha helix** the backbone coils with about 3.6 residues per turn, and each residue adds about 1.5 Å to the length. A helix of 20 residues is therefore about 20 × 1.5 = 30 Å long. Each backbone C=O forms a hydrogen bond with the N-H four residues further along.
- In a **beta sheet**, stretches of chain (strands) lie side by side and are joined by hydrogen bonds between neighbouring strands. Strands can run in the same direction (parallel) or opposite directions (antiparallel).
- Stretches that are neither helix nor strand (loops and turns) are often flexible, and often sit on the surface. §11 explains why that matters for judging models.

### Quaternary structure and hemoglobin

The slide's example of quaternary structure is **hemoglobin**, which is built from four chains (two alpha and two beta), each holding a heme group. One chain on its own is a tertiary structure; four chains together form the quaternary structure that carries oxygen.

### A change at the bottom propagates upward

The slide's picture shows a red marker at the primary level that reappears at every level above it, labelled "change propagates". The classic example is sickle-cell disease. In the beta chain of hemoglobin, residue 6 is normally **glutamate** (negatively charged, water-loving). In the sickle allele a single base changes the codon from `GAG` to `GTG`, and residue 6 becomes **valine** (hydrophobic). That one substitution creates a sticky patch on the surface of the folded chain, and under low oxygen the hemoglobin molecules stick together and form fibres that distort the red blood cell.

The code below makes the substitution on a short coding sequence and compares the two amino acids using the Kyte-Doolittle hydropathy scale (positive values are hydrophobic, negative values are hydrophilic):

```python
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(BASES)
               for j, b in enumerate(BASES)
               for k, c in enumerate(BASES)}

HYDROPATHY = {   # Kyte-Doolittle scale: positive = hydrophobic, negative = hydrophilic
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

def translate(dna):
    return "".join(CODON_TABLE[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3))

normal = "ATGGTGCACCTGACTCCTGAGGAGAAG"        # opening of the beta-globin coding sequence
mutant = normal[:19] + "T" + normal[20:]      # codon 6: GAG -> GTG (one base changed)

p_normal, p_mutant = translate(normal), translate(mutant)
print(p_normal, p_mutant)                     # MVHLTPEEK MVHLTPVEK
for pos, (a, b) in enumerate(zip(p_normal, p_mutant)):
    if a != b:
        print(f"residue {pos}: {a} -> {b}, hydropathy {HYDROPATHY[a]} -> {HYDROPATHY[b]}")
        # residue 6: E -> V, hydropathy -3.5 -> 4.2
```

The initial methionine is counted as position 0 in the code, so the changed residue appears as position 6, matching the usual numbering of the mature chain.

One changed base, one changed amino acid, a change in hydropathy from −3.5 to +4.2, and then a change in how the whole molecule behaves. That is the slide's point: **small changes at the primary level can alter outcomes at every level above.**

> **Key idea:** structure is a hierarchy. Primary sequence determines secondary and tertiary structure, assembly of chains gives quaternary structure, and each level helps determine function.

---

## 2. Methods for Determining Protein Structure — X-ray, NMR, Cryo-EM

Everything above is only useful if we can actually see the shape. Three main experimental methods do this. They aim at the same goal, **each using a different experimental approach**.

| | X-ray crystallography | NMR spectroscopy | Cryo-electron microscopy (cryo-EM) |
|---|---|---|---|
| **Sample** | Protein **crystals** | Protein dissolved in solution | Protein flash-frozen in a thin layer of ice |
| **What is measured** | The pattern of X-rays diffracted by the crystal | Magnetic responses of atomic nuclei, which report which atoms are close together | Thousands of electron-microscope images of individual particles |
| **How the structure is obtained** | The diffraction pattern is converted into an electron-density map, and an atomic model is built into it | Distance information between nearby atoms is used to calculate models consistent with it | Images from many orientations are combined computationally into a 3D map |
| **Strength** | The foundational method, able to give very detailed atomic structures | Studies smaller proteins in solution, **closer to natural conditions**, and can show flexibility | **Excels for large complexes**, and needs no crystal |
| **Main limitation** | **Getting a good crystal is often difficult**, and the crystal may hold the protein in just one form | Mostly limited to smaller proteins | Very small proteins are harder, and the equipment is specialised |

The slide's picture is a fork in the road with three signposts and the answer "It depends". The summary line says the same thing: **method choice depends on the question.**

Three practical examples:

- A small, stable protein that crystallises well, and you want fine atomic detail: X-ray.
- A small protein that will not crystallise, or whose motion in solution matters: NMR.
- A large assembly such as a protein-RNA machine or a membrane complex that cannot be crystallised: cryo-EM.

### Resolution

Experimental structures are reported with a **resolution** in ångströms. A smaller number means finer detail. At high resolution, individual side chains and even ordered water molecules can be placed confidently; at lower resolution only the overall shape and the backbone path are clear. When you use a structure, its resolution tells you how much detail to trust.

### Where the results go

Solved structures are deposited as coordinate files in the Protein Data Bank (PDB, from Module 11). The next section shows what is inside such a file.

> **Key idea:** X-ray, NMR, and cryo-EM answer the same question in different ways. The protein, its size, and the biological question decide which one to use.

---

## 3. Structure Visualization Tools — PyMOL and Chimera

A coordinate file is a long list of atoms and numbers. Visualization tools **convert raw coordinate files into interpretable 3D molecular models**. The slide stresses that they are used after the structure has been determined by X-ray, NMR, or cryo-EM.

| Tool | What it is good for |
|---|---|
| **PyMOL** | Producing clear images, and viewing helices, sheets, ligands, and binding sites |
| **Chimera** | **Interactive exploration** and analysis, including large and complex molecular assemblies (its newer version is called ChimeraX) |

Both can show the same structure in different styles: a cartoon that emphasises helices and sheets, sticks for individual residues and ligands, or a surface that shows the outer shape.

### What is inside a PDB file

Each atom is one line, in fixed columns. Protein atoms begin with `ATOM`, and everything that is not part of the protein chain (a bound drug, for example) begins with `HETATM`. A line gives the atom name, the residue name and number, the chain, and the x, y, z coordinates in ångströms.

### Why visualization matters: where is the mutation?

The slide's picture labels a residue (residue 145), an active site, and a stabilizing region. The summary: **visualization shows whether a mutation is in a pocket, on the surface, or stabilizing the fold.** The three cases have different meanings:

| Where the residue sits | What a change there may do |
|---|---|
| **In a pocket or active site** | Can block binding or catalysis directly |
| **On the surface** | May change how the protein meets other proteins, or may have little effect |
| **Buried in a stabilizing region** | May weaken folding and make the protein unstable |

A viewer shows this by eye. The same question can be asked in code. Below is a small made-up coordinate block (invented numbers, with a ligand as residue 301) and a function that finds which residues come within 4 Å of the ligand, which is the kind of "residues near the ligand" selection a viewer performs:

```python
import math

PDB_TEXT = """\
ATOM      1  N   GLY A 143      10.000  12.000   5.000  1.00 20.00           N
ATOM      2  CA  GLY A 143      11.200  12.500   5.500  1.00 20.00           C
ATOM      3  N   SER A 145      15.000  10.000   8.000  1.00 18.00           N
ATOM      4  CA  SER A 145      16.100  10.400   8.900  1.00 18.00           C
ATOM      5  OG  SER A 145      17.300  11.000   8.300  1.00 18.00           O
ATOM      6  CA  LEU A 146      20.000  13.000  15.000  1.00 22.00           C
ATOM      7  CA  ALA A 190      30.000  25.000  20.000  1.00 25.00           C
HETATM    8  C1  LIG A 301      19.000  13.000  10.000  1.00 30.00           C
HETATM    9  O1  LIG A 301      20.200  13.400  10.600  1.00 30.00           O
"""

def parse_pdb(text):
    atoms = []
    for line in text.splitlines():
        if line[:6].strip() in ("ATOM", "HETATM"):
            atoms.append({
                "hetero": line[:6].strip() == "HETATM",
                "name": line[12:16].strip(),
                "res": line[17:20].strip(),
                "chain": line[21],
                "num": int(line[22:26]),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            })
    return atoms

atoms = parse_pdb(PDB_TEXT)
protein = [a for a in atoms if not a["hetero"]]
ligand = [a for a in atoms if a["hetero"]]

closest = {}                                    # residue -> smallest distance to the ligand
for a in protein:
    d = min(math.dist(a["xyz"], l["xyz"]) for l in ligand)
    key = (a["res"], a["num"])
    closest[key] = min(d, closest.get(key, d))

for (res, num), d in sorted(closest.items(), key=lambda kv: kv[0][1]):
    flag = "contact" if d <= 4.0 else ""
    print(f"{res}{num}: {d:.2f} A  {flag}")
# GLY143: 9.02 A
# SER145: 3.13 A  contact
# LEU146: 4.42 A
# ALA190: 17.86 A
```

Reading the result: residue 145 (a serine) sits 3.13 Å from the ligand, so it is in the pocket and a change there could affect binding. Residue 146 is 4.42 Å away, just outside the cut-off, so it is close but not touching. Residue 190 is far from the ligand, and a change there is more likely to matter through stability or interactions elsewhere than through direct binding.

The slide's closing line states the purpose of all of this: **these tools support hypothesis generation, communication, and biological interpretation.** A picture suggests a hypothesis, a figure communicates it, and the interpretation connects it back to the biology.

> **Key idea:** visualization turns coordinates into understanding. It shows where a residue lies, and that location is the first clue to what a change in it might do.

---

## 4. Sequence-Structure-Function Relationships

The three words in the title form one chain, and **they are tightly linked in biology.**

1. **Sequence guides folding.** The amino acid sequence contains the information that guides the chain into its 3D structure. This was shown by Anfinsen's classic experiments, in which a denatured enzyme refolded to its working shape when conditions were restored.
2. **Shape determines what the protein can do.** A protein works by fitting other molecules.
3. **A sequence change can affect folding, stability, or binding.** Even one change can.

The slide's picture is a **lock and key**: an active site is a pocket, one key matches ("Match") and another does not ("No Match"). The summary: **active site changes can block proper target molecule interaction.**

### The two directions of bioinformatics

The slide's closing line describes how bioinformatics uses this chain: **sequence is used to predict structure, and structure is used to infer function.**

| Direction | Question | Methods in this module |
|---|---|---|
| Sequence → structure | What shape will this sequence adopt? | Homology modeling (§6), AI prediction (§9) |
| Structure → function | What can a protein with this shape do? | Domains and motifs (§5), docking (§8), interaction networks (§7) |

A useful habit to build, in the words of the lecture: look at a sequence and ask *what shape could this encode, and what job might that shape perform?*

### A caution

The chain is strong but not absolute. Different sequences can fold into similar shapes, and the same shape can serve different jobs in different proteins. So a structure is evidence for a function, not a guarantee, exactly as with sequence similarity in Module 12.

> **Key idea:** sequence determines structure, and structure determines function. Bioinformatics moves along that chain in both directions.

---

## 5. Protein Domain and Motif Analysis

Proteins are built from smaller functional pieces.

| | Domain | Motif |
|---|---|---|
| **What it is** | A larger, functional unit of the protein | A shorter, recurring functional pattern |
| **Typical size** | Tens to a couple of hundred residues | A handful to a few dozen residues |
| **Behaviour** | **Often folds independently** | May not act alone |
| **Role** | Linked to DNA binding, interactions, or catalysis | Serves as a recognisable **biological fingerprint** |

The slide's picture shows a domain as a ribbon, three possible roles (DNA binding, interaction, catalysis) above it, and a motif as a magnified signature inside it.

### Why they predict function

Known domains and motifs **help predict protein function before lab validation.** The slide's two examples:

- **Kinase domains suggest phosphorylation**: the protein probably transfers phosphate groups onto other molecules.
- **DNA-binding motifs suggest gene regulation**: the protein probably controls which genes are switched on.

The slide's bridge picture sums it up: sequence data crosses a bridge of domains and motifs to functional insight. In its words, domain and motif analysis gives **fast, evidence-based insight.**

### Where to look

Domain and motif libraries record families that have been characterised before. Pfam, InterPro, and PROSITE were introduced in Module 12 §10. A query sequence is scanned against them, and a match suggests the domain or motif is present.

### A motif scan

Some motifs can be written as patterns. The **P-loop**, a motif that helps bind ATP or GTP, is written in PROSITE style as `[AG]-x(4)-G-K-[ST]`: an A or G, then any four residues, then G, K, and an S or T. As a regular expression this is `[AG].{4}GK[ST]`. It can be tested on the start of the human HRAS protein:

```python
import re

# P-loop (ATP/GTP-binding): [AG]-x(4)-G-K-[ST]
P_LOOP = re.compile(r"[AG].{4}GK[ST]")
ras_start = "MTEYKLVVVGAGGVGKSALTIQLIQ"
for m in P_LOOP.finditer(ras_start):
    print(m.start() + 1, m.group())             # 10 GAGGVGKS

# chance matches: probability for one random position (uniform amino acid usage)
p = (2/20) * (1/20) * (1/20) * (2/20)
print(f"{p:.1e}", round(10_000_000 * p))        # 2.5e-05 250
```

The motif is found at position 10. The second half of the code is a warning: a pattern this short has a chance of 2.5 × 10⁻⁵ at any random position, so in a database of 10 million residues it would match about **250 times by chance alone**. A motif hit is therefore a clue to follow up (is it in a folded region? is the protein from a plausible family?), not a proof.

> **Key idea:** domains and motifs bridge sequence to function quickly. Treat them as strong clues that still need support from context.

---

## 6. Homology Modeling Principles

Most proteins have no experimentally solved structure. **Homology modeling predicts the 3D structure from evolutionary sequence similarity.**

The logic: **similar sequences often fold similarly, because evolution tends to conserve functionally useful shapes.** So a **known protein structure serves as the template** for an unknown but related protein (the target).

### The workflow

The slide's workflow is to choose a template, align the target to it, and build the model from the matched regions. In practice a fourth step, checking the model, completes it:

1. **Choose a template.** Search a structure database for a related protein with a solved structure. This is the BLAST search of Module 12, aimed at the PDB.
2. **Align the target to the template.** The quality of this alignment is the most important factor; errors here become errors in the model.
3. **Build the model from the matched regions.** The backbone of aligned regions is copied from the template. Insertions, deletions, and gaps (loops) have to be modelled separately, and the side chains are placed.
4. **Evaluate the model** (see §11).

### Where the model is trustworthy

The slide's picture shows an alignment coloured green where the sequences match and orange where they do not, with a confidence bar chart above it: **model confidence is highest in well-aligned regions and lower in weak or gapped segments.**

As a rule of thumb, the higher the sequence identity between target and template, the more reliable the model. Above roughly 30% identity a template is usually good enough for a useful model, while at 20-30% (sometimes called the "twilight zone") the alignment itself becomes uncertain and the model needs much more care.

### Identity and confidence in code

```python
template = "MKTAYIAKQRQISFVKSHFSRQ"
target   = "MKSAYLAKQR-ISFVKTHFARQ"

def compare(template, target):
    ident = mism = gaps = 0
    marks = []
    for a, b in zip(template, target):
        if a == "-" or b == "-":
            gaps += 1; marks.append(".")
        elif a == b:
            ident += 1; marks.append("|")
        else:
            mism += 1; marks.append(".")
    aligned = ident + mism
    return ident, mism, gaps, "".join(marks), 100 * ident / aligned

ident, mism, gaps, marks, pct = compare(template, target)
print(template); print(marks); print(target)
print(ident, mism, gaps, round(pct, 1))
# MKTAYIAKQRQISFVKSHFSRQ
# ||.||.||||.|||||.||.||
# MKSAYLAKQR-ISFVKTHFARQ
# 17 4 1 81.0
```

Across the 21 aligned positions there are 17 identical residues (81.0% identity), 4 mismatches, and 1 gap. Each `|` is a position where the template's coordinates can be used with confidence, and each `.` marks a mismatch or gap that is less reliable, like the orange stretches in the slide's picture. Real modeling programs weigh many more factors, but identity along the alignment is the starting point.

### Why it is used

**Homology modeling is faster and cheaper than experiments** and **guides further analysis**: a model suggests where the active site is, which residues to mutate, and what to test, even before any experiment.

> **Key idea:** a known structure is the template for a related protein. Models are most reliable where the alignment is good, and they are hypotheses to be tested.

---

## 7. Protein-Protein Interaction Networks

**Proteins usually function through interactions, not in isolation.** Signalling, metabolism, and gene regulation all depend on groups of proteins working together.

### The network picture

In a protein-protein interaction (PPI) network:

- each **node** is a protein,
- each **edge** is an interaction between two proteins.

The slide's picture is a transit map, with stations for proteins and coloured lines for interactions. It makes the central point: **networks reveal coordinated cellular behaviour, not just single-protein roles.**

Interaction data comes from experiments that detect which proteins bind each other, from curated databases (STRING, BioGRID, and IntAct are examples), and from computational prediction.

### What to read from a network

| Feature | Meaning |
|---|---|
| **Degree** | The number of interaction partners a protein has |
| **Hub** | A protein with many partners, connecting many parts of the network |
| **Module** | A densely connected group of proteins that work together |
| **Pathway** | A chain of interactions carrying a process, such as a signal |

The slide's summary: **PPI networks reveal hubs, modules, functions, pathways, and drug targets.**

### A small network in code

Eleven interactions among nine proteins (the names are placeholders). The code finds each protein's degree, picks out the hub, and then removes the hub to see what happens to the network:

```python
from collections import defaultdict

edges = [("P1","P2"),("P1","P3"),("P1","P4"),("P1","P5"),("P1","P6"),
         ("P2","P3"),("P5","P6"),("P6","P7"),("P7","P8"),("P7","P9"),("P8","P9")]
graph = defaultdict(set)
for a, b in edges:
    graph[a].add(b); graph[b].add(a)

degree = {p: len(n) for p, n in graph.items()}
hub = max(degree, key=degree.get)
print(sorted(degree.items(), key=lambda kv: -kv[1])[:3], hub)
# [('P1', 5), ('P6', 3), ('P7', 3)] P1

def components(graph, removed=()):
    seen, groups = set(removed), []
    for start in graph:
        if start in seen: continue
        stack, group = [start], set()
        while stack:
            node = stack.pop()
            if node in seen: continue
            seen.add(node); group.add(node)
            stack.extend(n for n in graph[node] if n not in seen)
        groups.append(sorted(group))
    return groups

print(components(graph))
# [['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9']]
print(components(graph, removed={hub}))
# [['P2', 'P3'], ['P4'], ['P5', 'P6', 'P7', 'P8', 'P9']]
```

P1 is the hub with five partners. With P1 present the network is one connected whole. With P1 removed it breaks into three separate pieces, and the protein P4 is left completely isolated. This is the slide's last point in miniature: **disrupted interactions can affect the wider system and cell function.** The tightly connected triangles (P2-P3 with P1, and P7-P8-P9) are what a module looks like.

### Cautions

- Edges differ in quality: some come from direct experiments and some from weak or indirect evidence, so networks contain false positives and gaps.
- An interaction may happen only in certain cells or conditions.
- A hub is not automatically a good drug target. A hub that every cell needs may be dangerous to inhibit.

> **Key idea:** a protein's role is defined partly by its neighbours. PPI networks show hubs, modules, and pathways, and point to proteins whose disruption matters.

---

## 8. Molecular Docking and Drug Target Analysis

Two different questions come up in drug work, and the slide's picture puts them side by side: **"Docking: What binds?"** and **"Target analysis: Worthwhile?"**

### Docking

**Docking predicts how a ligand may bind to a protein target.** A ligand is a small molecule (a candidate drug, for example), and the protein is the target. Docking does three things:

1. It estimates the **binding pose**: where the ligand sits in the protein, and in what orientation.
2. It **scores and ranks** the interactions by likely binding strength.
3. It does this allowing for flexibility: **the ligand and the protein can both be flexible, unlike a rigid lock-and-key.**

A typical docking run prepares the protein, defines the region to search (usually a pocket), generates many candidate poses, scores them, and ranks them. AutoDock Vina is a widely used program of this kind.

### Reading a docking score

Scores are usually reported as an estimated binding energy in kcal/mol, where more negative means stronger predicted binding. A binding free energy ΔG corresponds to a dissociation constant Kd = e^(ΔG / RT), where a smaller Kd means tighter binding:

```python
import math

R = 0.001987          # gas constant, kcal/(mol*K)
T = 298.0             # kelvin

def kd_from_dg(dg_kcal):
    """Dissociation constant (molar) from a binding free energy: Kd = exp(dG / RT)."""
    return math.exp(dg_kcal / (R * T))

for dg in (-6.0, -7.4, -9.0):
    kd = kd_from_dg(dg)
    print(dg, f"{kd:.2e} M", f"{kd * 1e9:.0f} nM")
# -6.0 3.97e-05 M 39747 nM
# -7.4 3.74e-06 M 3737 nM
# -9.0 2.51e-07 M 251 nM

print(round(R * T * math.log(10), 2))      # 1.36 kcal/mol per 10-fold change in Kd
```

At room temperature, every **1.36 kcal/mol** more negative means about a 10-fold tighter binder. A score of −9.0 corresponds to a Kd of roughly 250 nM, while −6.0 corresponds to about 40,000 nM. Docking scores are only rough estimates, though, and can be wrong by much more than these tidy numbers suggest. They are good for ranking and screening, and poor as exact predictions.

### Checking a pose: RMSD

If an experimental structure of the ligand bound to the protein is known, a docked pose can be compared with it using the **root-mean-square deviation (RMSD)**, the average distance between matching atoms. A common rule of thumb is that a pose within about 2 Å RMSD of the experimental one counts as correctly predicted:

```python
import math

def rmsd(pose_a, pose_b):
    """Root-mean-square deviation between two poses with matching atom order."""
    sq = [sum((p - q) ** 2 for p, q in zip(a, b)) for a, b in zip(pose_a, pose_b)]
    return math.sqrt(sum(sq) / len(sq))

crystal    = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0), (3.0, 1.5, 0.0)]
docked_ok  = [(0.3, 0.2, 0.0), (1.7, 0.3, 0.1), (3.2, -0.2, 0.0), (3.3, 1.4, 0.2)]
docked_bad = [(0.0, 4.0, 0.0), (1.5, 4.0, 0.0), (3.0, 4.0, 0.0), (3.0, 5.5, 0.0)]
print(round(rmsd(crystal, docked_ok), 2), round(rmsd(crystal, docked_bad), 2))
# 0.35 4.0
```

The first pose is almost on top of the reference (0.35 Å) and the second is the same shape in the wrong place (4.0 Å).

### Target analysis

The slide's warning is explicit: **a strong docking result alone does not prove a useful drug target.** A molecule can fit a pocket well on a computer and still be pointless. Target analysis asks the second question, whether the target is worth pursuing, and the slide lists three checks:

| Check | Question |
|---|---|
| **Disease relevance** | Is this protein actually involved in the disease? |
| **Accessibility** | Can a drug reach it, and does it have a pocket a drug can bind? |
| **Likely effect** | If its activity is changed, is the biological effect the one we want? |

> **Key idea:** docking asks what binds; target analysis asks whether it matters. A good drug effort needs both answers.

---

## 9. AI-Driven Structure Prediction — Understanding Tools Like AlphaFold

For decades, solving a structure was **slow, costly, and often incomplete**. AI changed this: **it made protein structure prediction faster and more accessible.**

### How it works, in outline

- The AI **learns folding patterns from known structures** and then predicts new ones.
- **AlphaFold predicts an accurate 3D structure from an amino acid sequence.** Beyond the sequence, it makes use of evolutionary information from related sequences, the multiple sequence alignments of Module 12.
- Predicted structures are available in public databases for very large numbers of proteins, so for many sequences a model can be looked up rather than calculated.
- As the slide says, **structure helps explain function, interactions, and biological role.**

### Confidence comes with the prediction

AlphaFold reports a per-residue confidence score called **pLDDT**, from 0 to 100. The usual bands are:

| pLDDT | Band | Meaning |
|---|---|---|
| Above 90 | Very high | Backbone and often side chains are expected to be accurate |
| 70 to 90 | Confident | The backbone is generally reliable |
| 50 to 70 | Low | Treat with caution |
| Below 50 | Very low | Often a flexible or disordered region, and may be meaningless as a shape |

These values are stored in the B-factor column of the predicted coordinate file, so a viewer from §3 can colour the structure by confidence. A second output, the predicted aligned error (PAE), tells you how confidently two parts of the protein are placed relative to each other, which matters for multi-domain proteins.

```python
def band(plddt):
    if plddt > 90:  return "very high"
    if plddt > 70:  return "confident"
    if plddt > 50:  return "low"
    return "very low"

plddt = [45, 52, 60, 72, 85, 91, 95, 96, 94, 93, 90, 88, 92, 95, 97, 96, 74, 68, 55, 48, 42, 40]
counts = {}
for v in plddt:
    counts[band(v)] = counts.get(band(v), 0) + 1
print(counts)
# {'very low': 4, 'low': 4, 'confident': 5, 'very high': 9}

def low_segments(scores, cutoff=70, min_len=3):
    """Runs of at least min_len consecutive residues (numbered from 1) below the cutoff."""
    segments, start = [], None
    for i, v in enumerate(scores + [100], start=1):    # sentinel closes a final run
        if v < cutoff and start is None:
            start = i
        elif v >= cutoff and start is not None:
            if i - start >= min_len:
                segments.append((start, i - 1))
            start = None
    return segments

print(low_segments(plddt))        # [(1, 3), (18, 22)]
print(round(sum(plddt) / len(plddt), 1))   # 75.8
```

For this 22-residue example, the average pLDDT is 75.8, which sounds acceptable, but the details matter more: the confident core is in the middle, and residues 1-3 and 18-22 (the two ends) are low-confidence. A single average hides that.

### Not magic

The slide's picture crosses out a crystal ball labelled "Not Magic". **AI predictions speed research but still need biological interpretation.** Typical limits:

- A prediction is a single shape, while real proteins move and can adopt several.
- Flexible or disordered regions are poorly predicted, as the low scores show.
- It does not by itself tell you how a single point mutation changes stability or function.

> **Key idea:** AlphaFold gives a fast, often accurate model and an honest map of its own confidence. Use both.

---

## 10. Structural Bioinformatics in Drug Discovery

Structure predictions become valuable when they are applied: **AI-predicted structures become useful when applied to real drug discovery questions.** Structural bioinformatics **uses 3D molecular shape to predict drug interactions.**

The slide's central picture: **protein shape works like a lock; better shape knowledge improves how well a molecule fits.**

### What structures contribute

| Use | How |
|---|---|
| **Finding binding sites** | Pockets on the surface or inside the protein where a drug could bind |
| **Comparing proteins** | Checking related proteins for differences in the pocket, to design drugs that hit one target and spare its relatives |
| **Guiding design** | Seeing how a candidate fits, and which changes would improve it |

The slide notes that **many drugs work by binding abnormal proteins and changing their activity.** Kinase inhibitors such as imatinib are an example: they bind the working part of an abnormal kinase and reduce its activity.

### A structure-based workflow

```
disease protein identified            (target analysis, §8)
   -> structure obtained               (experiment, §2, or prediction, §6 and §9)
   -> binding site located             (visualization, §3)
   -> candidate molecules docked       (docking, §8)
   -> best candidates tested           (laboratory experiments)
   -> candidates improved using the structure
   -> selectivity checked against related proteins
```

Structure does not replace the laboratory: docking narrows a huge set of possible molecules down to a few worth testing, and experiments decide. The benefit is the one the slide closes with: **structural insight turns molecular models into targeted drug discovery decisions**, making the search more targeted, efficient, and informed.

> **Key idea:** shape lets you design for fit. Structure guides which molecules to test, and experiments confirm.

---

## 11. Interpreting and Validating Predicted Structures

A predicted structure is useful only when it is **interpreted carefully and validated sensibly.** The slide's pointed warning: **a visually impressive model is not automatically scientifically reliable.** Predicted structures look equally smooth and professional whether they are right or wrong.

### How to read a prediction

1. **Check confidence by separating well-predicted regions from uncertain ones.** Use pLDDT (§9) region by region, not just the average.
2. **Know which parts to trust.** **Core regions may be reliable; loops and flexible tails often need caution.** The ends and loops of a chain are where models are weakest.
3. **Compare the fold and key features with expected biology and function.** Does the model contain the domains the sequence predicts? Are the known active-site residues gathered together in the same pocket? Are conserved residues where conservation says they should be?

### Ways to validate

| Check | What it asks |
|---|---|
| **Confidence scores** | Is the region you care about high-confidence? (pLDDT, PAE) |
| **Comparison with a known structure** | If a related protein or the same protein has an experimental structure, how closely do they agree? RMSD (§8) measures this |
| **Geometry** | Are bond lengths, backbone angles (the Ramachandran plot), and atom contacts physically possible, with no overlapping atoms? Programs such as MolProbity check this |
| **Biological consistency** | Does the model agree with known mutations, binding data, and conservation from Module 12? |
| **Experiment** | Does a direct test (for example, mutating a predicted contact residue) behave as the model predicts? |

The slide's conclusion: **validate with known structures or experiments; use informed judgment, not blind trust.**

### A sensible reading of the earlier example

For the 22-residue pLDDT list in §9: the central stretch (residues 4-17) is confident to very high, so questions about that core are reasonable. Residues 1-3 and 18-22 are low confidence, so any claim that depends on how those ends sit (a binding site, a contact with another chain) should be treated as unsupported until something else backs it.

> **Key idea:** a model is evidence, not an answer. Judge it region by region, compare it with biology, and confirm what matters experimentally.

---

## 12. Key Takeaways

- Protein structure has four levels: primary (sequence), secondary (helices and sheets from backbone hydrogen bonds), tertiary (the 3D shape of one chain, held by side-chain interactions), and quaternary (several chains together, as in hemoglobin). Each level builds on the one before.
- A small primary change can propagate upward: in sickle-cell disease, one base change (GAG to GTG) swaps glutamate for valine and changes how hemoglobin behaves.
- X-ray crystallography, NMR, and cryo-EM reach the same goal differently: X-ray needs good crystals, NMR suits smaller proteins in solution, and cryo-EM suits large complexes. The question decides the method.
- PyMOL and Chimera turn coordinate files into models, and show whether a mutation sits in a pocket, on the surface, or in a stabilizing region.
- Sequence guides folding, shape determines function, and bioinformatics predicts structure from sequence and infers function from structure.
- Domains are larger functional units, motifs are short recognisable fingerprints, and both predict function. Short motifs also match by chance, so they are clues, not proof.
- Homology modeling builds a model from a related template. Models are most reliable where the alignment is good.
- PPI networks put proteins in context, revealing hubs, modules, and pathways, and showing how a disrupted interaction affects the system.
- Docking predicts the pose and strength of ligand binding, but a good score does not prove a useful drug target. Target analysis checks disease relevance, accessibility, and likely effect.
- AlphaFold predicts structure from sequence and reports per-residue confidence. It is fast and often accurate, and still needs interpretation.
- Structural bioinformatics guides drug discovery by locating binding sites, comparing proteins, and guiding design.
- A pretty model is not a reliable model. Check confidence region by region, compare with biology, and validate with known structures or experiments.

---

## 13. What's Next

With sequences, genomes, and structures understood, the next questions look beyond a single molecule: how genes and proteins vary across species and populations, and how they act together in whole cells and systems.
