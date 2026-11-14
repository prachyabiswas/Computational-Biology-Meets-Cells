# From Raw Reads to Annotated Genomes

> Module 13: How a genome is read, assembled, and interpreted — Sanger versus next-generation sequencing, short and long reads, de novo versus reference-guided assembly, overlap-layout-consensus and de Bruijn graph assemblers, N50 and other quality metrics, gene prediction, structural and functional annotation, regulatory elements, comparative annotation, and why repeats make genomes hard to assemble.

---

## Table of Contents

1. [Genome Sequencing Strategies — Sanger vs Next-Generation Sequencing](#1-genome-sequencing-strategies--sanger-vs-next-generation-sequencing)
2. [Read Types and Sequencing Platforms Overview](#2-read-types-and-sequencing-platforms-overview)
3. [De Novo vs Reference-Guided Assembly](#3-de-novo-vs-reference-guided-assembly)
4. [Assembly Algorithms — Overlap-Layout-Consensus and De Bruijn Graphs](#4-assembly-algorithms--overlap-layout-consensus-and-de-bruijn-graphs)
5. [Assessing Assembly Quality — N50 and Other Metrics](#5-assessing-assembly-quality--n50-and-other-metrics)
6. [Genome Annotation Fundamentals — Gene Prediction](#6-genome-annotation-fundamentals--gene-prediction)
7. [Structural vs Functional Annotation](#7-structural-vs-functional-annotation)
8. [Annotation Pipelines and Tools](#8-annotation-pipelines-and-tools)
9. [Identifying Regulatory Elements and Non-Coding Regions](#9-identifying-regulatory-elements-and-non-coding-regions)
10. [Comparative Genome Annotation](#10-comparative-genome-annotation)
11. [Challenges in Assembling Complex and Repetitive Genomes](#11-challenges-in-assembling-complex-and-repetitive-genomes)
12. [Key Takeaways](#12-key-takeaways)
13. [What's Next](#13-whats-next)

---

## 1. Genome Sequencing Strategies — Sanger vs Next-Generation Sequencing

A sequencer does not read a genome from one end to the other. It reads millions of short pieces, and software puts them back together. This whole module follows that route:

```
DNA sample
   -> sequencing                 (reads, stored as FASTQ with quality scores)
   -> assembly                   (reads become contigs and scaffolds, stored as FASTA)
   -> assembly quality check     (N50, completeness, misassemblies)
   -> structural annotation      (where are the genes?)
   -> functional annotation      (what do the genes do?)
   -> annotated genome           (stored as GFF/GTF plus protein FASTA)
```

The file formats at each step are the ones from Module 11. The comparisons that power annotation are the alignments from Module 12.

### The two main strategies

| | Sanger sequencing | Next-generation sequencing (NGS) |
|---|---|---|
| **How it reads** | One DNA fragment at a time | Millions to billions of fragments **in parallel** |
| **Read length** | Roughly 500 to 1,000 bases per read | Short reads of about 50 to 300 bases, or long reads of many kilobases (see §2) |
| **Accuracy** | Very high for each read | Each read is less reliable, but errors are outvoted when every position is read many times |
| **Scale** | Small | Very large |
| **Cost per base** | High | Low |
| **Typical use** | Confirming a mutation, checking one gene region, verifying a cloned insert | Whole genomes, many samples at once |

The slide's picture is a balance: Sanger on one side (accuracy, small scale) and NGS on the other (speed, large scale). Neither side is "better". The lesson of the slide is to **choose based on the biological question, scale, time, budget, and workflow needs**.

Two quick decisions to make this concrete:

- You suspect one specific mutation in one gene of one patient. Sanger is the fit: a single, accurate answer for a single region.
- You have a newly isolated bacterium and want its complete genome. NGS is the fit: no one would sequence millions of bases one fragment at a time.

### How many reads are enough? Coverage

Because NGS reads are short, a genome is covered by many overlapping reads. **Coverage** (also called depth) is the average number of reads that cover a position:

```
Coverage C = (number of reads N x read length L) / genome size G
```

Worked example: a bacterial genome of 5 million bases, read at 30x coverage with 150-base reads:

```
N = C x G / L = 30 x 5,000,000 / 150 = 1,000,000 reads
```

If reads land randomly, the fraction of the genome that no read touches is about e^(−C). The idealised numbers:

| Coverage | Fraction of genome with no read | Reading |
|---|---|---|
| 1x | about 37% | Large gaps |
| 5x | about 0.67% | Still thousands of gaps in a large genome |
| 10x | about 0.0045% | Few gaps |
| 30x | about 9 × 10⁻¹⁴ | Effectively none from chance alone |

Real coverage is never perfectly random (some regions are harder to sequence), so practical targets are higher than this table suggests. The point is why NGS works at all: **redundancy** is what turns many unreliable short reads into one reliable sequence.

> **Key idea:** Sanger gives accuracy on a small scale, NGS gives throughput on a large scale. The biological question decides which one you use.

---

## 2. Read Types and Sequencing Platforms Overview

A **read** is the sequence a machine reports for one fragment of DNA. Reads fall into two broad groups, and the group you pick shapes every downstream step.

| | Short reads | Long reads |
|---|---|---|
| **Length** | Roughly 50 to 300 bases | Thousands to tens of thousands of bases, and in some platforms far longer |
| **Strengths** | Fast, high-throughput, cost-effective for routine genomics | Span large regions, so they help resolve complex assemblies and variants |
| **Weakness** | Struggle in repetitive regions, because similar fragments are hard to place | Higher cost, and historically a higher error rate per base (improving) |
| **Example platform families** | Illumina | PacBio, Oxford Nanopore |

The reason short reads struggle with repeats is simple. Imagine a read of 150 bases that lies entirely inside a 5,000-base repeat that occurs twice in the genome. That read fits equally well in both copies, so nothing tells the assembler where it belongs. A read of 10,000 bases would reach out of the repeat into the unique sequence on each side, and that anchors it. §11 returns to this.

The slide's balance shows four things competing: **read length, throughput, accuracy, and cost**. No platform wins on all four, which is why platform choice depends on the project. A common compromise is to combine them: long reads to get the large-scale structure right, and short, accurate reads to correct small errors.

Every read comes with a quality value for each base (the Phred scores of Module 11). Assemblers and variant callers use these to decide how far to trust a base.

> **Key idea:** platform choice balances read length, throughput, accuracy, and cost. Read length decides how well you can reconstruct and interpret a genome.

---

## 3. De Novo vs Reference-Guided Assembly

**Assembly** turns sequencing reads into a genome sequence. There are two main strategies.

| | De novo assembly | Reference-guided assembly |
|---|---|---|
| **Idea** | Use the overlaps between reads to build the sequence from scratch | Map (align) each read to a related, already-known genome |
| **Analogy** | A jigsaw puzzle without the picture on the box (the slide's analogy) | The same puzzle, with the picture on the box to check against |
| **Needs a reference?** | No | Yes, and a good one |
| **Best when** | No good reference exists, or you are looking for novel sequences and structural changes | A close reference exists and you want speed and simplicity |
| **Cost** | Heavy computation and enough coverage; repeats make it hard | Faster and simpler, since each read only has to find its place |
| **Main risk** | Fragmented or wrongly joined assemblies | **Reference bias**: sequence that is absent or very different in the reference is missed or mis-placed |

The slide's picture shows the risk of reference-guided work with a single question: "New variation?". If a stretch of your sample's genome does not exist in the reference, there is nowhere for those reads to map, and the new sequence simply disappears from the result.

Two concrete cases:

- A newly isolated bacterium with no close relative in the databases needs de novo assembly.
- Re-sequencing a human sample to find small differences from the human reference is reference-guided work, because a high-quality reference exists and the question is about differences.

The slide's conclusion is that **the best choice depends on the organism, the reference quality, and whether discovery or efficiency matters more**.

Mapping reads to a reference uses the same idea as the local alignments of Module 12, but organised to handle millions of reads quickly. Specialised read mappers do this; the details are beyond this lecture.

> **Key idea:** de novo assembly can find what no reference contains; reference-guided assembly is faster when a good reference exists.

---

## 4. Assembly Algorithms — Overlap-Layout-Consensus and De Bruijn Graphs

Both classic algorithms try to answer the same question: *which reads belong next to which?* They differ in how they represent the reads.

| | Overlap-Layout-Consensus (OLC) | De Bruijn graph |
|---|---|---|
| **Works with** | Whole reads | **k-mers**, which are substrings of length k cut from the reads |
| **Best with** | Long reads, because long overlaps are easy to detect reliably | Very large short-read datasets |
| **Main cost** | Comparing reads against each other | Storing and walking the k-mer graph |
| **Weak point** | Slow when there are very many reads | Sensitive to the choice of k, sequencing errors, and repeats |

The slide's rule: **algorithm choice depends on how sequence information is represented and on read length.**

### Overlap-Layout-Consensus (OLC)

Three stages, as the name says:

1. **Overlap.** Compare every read with every other read and find pairs where the end of one matches the start of the other.
2. **Layout.** Use these overlaps to place the reads in order along the genome.
3. **Consensus.** At every position, take a vote across the reads that cover it. A random sequencing error in one read is outvoted by the other reads.

Worked example. The genome is `ATGGCTACCGTTAGCAAGT`, and the four reads are:

| Read | Sequence |
|---|---|
| 1 | `ATGGCTAC` |
| 2 | `GCTACCGT` |
| 3 | `CCGTTAGC` |
| 4 | `TAGCAAGT` |

The overlaps between neighbours (end of one equals start of the next):

| Pair | Overlap | Length |
|---|---|---|
| Read 1 → Read 2 | `GCTAC` | 5 |
| Read 2 → Read 3 | `CCGT` | 4 |
| Read 3 → Read 4 | `TAGC` | 4 |

Chain them, writing each overlap only once, and the result is `ATGGCTACCGTTAGCAAGT`, the original genome.

Why OLC works best with long reads: a long overlap is unlikely to occur by chance, so the "overlap" stage rarely produces false pairs. Why it struggles with huge datasets: the overlap stage compares every pair. For 1,000,000 reads that is about 5 × 10¹¹ pairs.

A minimal version in plain Python (a greedy version that keeps extending one contig; real assemblers are far more careful):

```python
def overlap(a, b, min_len=3):
    """Longest suffix of a that equals a prefix of b."""
    best = 0
    for k in range(min_len, min(len(a), len(b)) + 1):
        if a[-k:] == b[:k]:
            best = k
    return best

def olc(reads):
    contig = reads[0]
    remaining = reads[1:]
    while remaining:
        best_len, best_read = 0, None
        for r in remaining:
            o = overlap(contig, r)
            if o > best_len:
                best_len, best_read = o, r
        if best_read is None:
            break
        contig += best_read[best_len:]
        remaining.remove(best_read)
    return contig

reads = ["ATGGCTAC", "GCTACCGT", "CCGTTAGC", "TAGCAAGT"]
print(olc(reads))        # ATGGCTACCGTTAGCAAGT
```

### De Bruijn graph assembly

Instead of comparing whole reads, this method **cuts every read into k-mers** and builds a graph:

- Each **edge** is a k-mer.
- Each **node** is a (k−1)-mer, the first or last k−1 letters of a k-mer.
- An edge runs from the first k−1 letters of a k-mer to its last k−1 letters.

Assembling the genome then means finding a path that uses every edge exactly once. The slide's picture shows this as a branching graph of k-mers fed by many short reads.

Worked example with k = 4. The reads are `ATGGCGT`, `GGCGTGC`, and `GCGTGCA`. Cutting them into 4-mers gives seven distinct k-mers:

`ATGG`, `TGGC`, `GGCG`, `GCGT`, `CGTG`, `GTGC`, `TGCA`

The graph (each line is one edge, written as node → node):

| K-mer | Edge |
|---|---|
| `ATGG` | `ATG` → `TGG` |
| `TGGC` | `TGG` → `GGC` |
| `GGCG` | `GGC` → `GCG` |
| `GCGT` | `GCG` → `CGT` |
| `CGTG` | `CGT` → `GTG` |
| `GTGC` | `GTG` → `TGC` |
| `TGCA` | `TGC` → `GCA` |

The graph is a single chain, and walking it from `ATG` gives the path `ATG → TGG → GGC → GCG → CGT → GTG → TGC → GCA`. Writing the first node and then the last letter of each later node gives `ATGGCGTGCA`, the genome.

Why this scales: each read of length L gives L − k + 1 k-mers (a 150-base read with k = 31 gives 120), and repeated k-mers collapse into one edge. There is no all-against-all comparison.

The same idea in code, including a function that lists every way of walking the graph:

```python
from collections import defaultdict

def kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

def build_graph(reads, k):
    graph = defaultdict(set)
    for r in reads:
        for km in kmers(r, k):
            graph[km[:-1]].add(km[1:])       # edge: (k-1)-mer -> (k-1)-mer
    return graph

def all_paths(graph):
    """Every path that uses each edge exactly once (fine for tiny graphs)."""
    edges = {(u, v) for u, vs in graph.items() for v in vs}
    indeg = defaultdict(int)
    for u, v in edges:
        indeg[v] += 1
    nodes = {u for u, _ in edges} | {v for _, v in edges}
    starts = [n for n in nodes
              if sum(1 for u, _ in edges if u == n) - indeg[n] == 1] or sorted(nodes)
    results = []
    def walk(node, left, path):
        if not left:
            results.append(path)
            return
        for (u, v) in sorted(left):
            if u == node:
                walk(v, left - {(u, v)}, path + [v])
    for s in starts:
        walk(s, edges, [s])
    return results

def to_sequence(path):
    return path[0] + "".join(node[-1] for node in path[1:])

g = build_graph(["ATGGCGT", "GGCGTGC", "GCGTGCA"], k=4)
print([to_sequence(p) for p in all_paths(g)])        # ['ATGGCGTGCA']
```

### Choosing k, and why errors hurt

- A **small k** produces k-mers that occur by chance in many places, so the graph tangles.
- A **large k** gives more distinctive k-mers, but needs more coverage, because a k-mer is only useful if some read contains it completely.
- A single sequencing error in the middle of a read corrupts **up to k k-mers at once**. These false k-mers create short dead-end branches and bubbles in the graph that the assembler has to clean away.

§11 shows what repeats do to this graph with a worked example.

> **Key idea:** OLC compares whole reads and suits long reads; de Bruijn graphs break reads into k-mers and scale to huge short-read datasets. Neither is universally better.

---

## 5. Assessing Assembly Quality — N50 and Other Metrics

An assembly is a hypothesis about the genome, and it needs checking. The slide's first point is that **quality metrics help you judge whether an assembly is useful and reliable.**

The pieces an assembler produces are called **contigs** (continuous stretches with no gaps) and, after joining them using extra information, **scaffolds** (contigs ordered and oriented with gaps of estimated size between them).

### N50

N50 measures **contiguity**: half of the assembly is contained in contigs of length N50 or longer.

How to compute it:

1. Sort the contigs from longest to shortest.
2. Add up lengths from the longest down, until the running total reaches half of the total assembly length.
3. The length of the contig that crosses the halfway point is the **N50**. The number of contigs used to get there is the **L50**.

Worked example. Both assemblies below total 300 kb:

| Assembly | Contig lengths (kb) | Half of total | Running total | N50 | L50 |
|---|---|---|---|---|---|
| A | 80, 70, 50, 40, 30, 20, 10 | 150 | 80, then 150 (reached at the 2nd contig) | **70 kb** | 2 |
| C | ten contigs of 30 each | 150 | 30, 60, 90, 120, 150 (reached at the 5th contig) | **30 kb** | 5 |

Assembly A is more contiguous: the same total sequence sits in fewer, longer pieces.

The same calculation in code:

```python
def n50(lengths):
    total = sum(lengths)
    running = 0
    for count, length in enumerate(sorted(lengths, reverse=True), start=1):
        running += length
        if running >= total / 2:
            return length, count          # (N50, L50)

assembly_a = [80, 70, 50, 40, 30, 20, 10]     # kb, total 300
assembly_c = [30] * 10                        # kb, total 300
print(n50(assembly_a))   # (70, 2)
print(n50(assembly_c))   # (30, 5)
print(n50([150, 150]))   # (150, 1)
```

### A high N50 does not guarantee a correct assembly

The last line of the code shows the trap. An assembly of two contigs of 150 kb has the best possible N50 for 300 kb of sequence, yet if the assembler wrongly joined unrelated pieces to make those two contigs, the assembly is **long but wrong**. N50 measures only length, never correctness.

So N50 must be read with other evidence:

| Check | What it tells you |
|---|---|
| **Total assembly length** | Compare with the expected genome size. Far too short suggests missing regions; far too long suggests duplicated or contaminating sequence |
| **Number of contigs** | Many short contigs suggest fragmentation |
| **Largest contig** | Shows the best-case contiguity, but is a single number |
| **Completeness** | How many expected genes, or how many of the original reads, are present. Tools such as BUSCO (expected single-copy genes) and mapping the reads back to the assembly are used for this |
| **Accuracy / misassemblies** | Whether contigs contain false joins. Tools such as QUAST check this, especially when a reference is available |

The slide's own warning: **unexpected length or many short contigs can signal problems like fragmentation.** The picture's podium makes the same point: length, count, and completeness stand together, and no single number should be trusted alone.

> **Key idea:** N50 measures contiguity only. Read it with total length, contig count, and completeness or accuracy evidence.

---

## 6. Genome Annotation Fundamentals — Gene Prediction

An assembly is just a long string of A, C, G, and T. **Annotation** is the work of marking what is in it. The first step is **gene prediction**: identifying the DNA regions that are likely to encode genes.

The slide's analogy: finding real sentences in a book written with no spaces or punctuation. The letters are all there, but the boundaries are not marked.

### The signals

| Signal | What it is |
|---|---|
| **ORF (open reading frame)** | A stretch of codons that starts at a start codon, runs without interruption, and ends at a stop codon |
| **Start site** | The start codon, usually `ATG` |
| **Stop site** | One of the three stop codons: `TAA`, `TAG`, `TGA` |
| **Other coding features** | The statistical "feel" of coding DNA, such as codon usage and the three-base periodicity of coding sequence |

Because DNA is double-stranded and a gene can lie on either strand, a program has to look in **six reading frames**: three on the forward strand and three on the reverse complement.

### Why a long ORF is meaningful

In random sequence, a stop codon turns up quickly. Three of the 64 codons are stops, so on average a stop appears after about 64 / 3 ≈ 21 codons. The chance that a random stretch of N codons contains no stop is (61/64)^N:

| ORF length (codons) | Chance of this length with no stop, in random sequence |
|---|---|
| 20 | about 38% |
| 50 | about 9% |
| 100 | about 0.8% |
| 300 | about 6 × 10⁻⁷ |

So a short ORF means little, but a long ORF is unlikely to be an accident and is probably a real gene.

### A working ORF finder

The following reports ORFs on both strands. Positions are given on the forward strand, and the length counts the codons before the stop. A minimum of three codons is used here only to keep the demonstration short:

```python
STOPS = {"TAA", "TAG", "TGA"}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def revcomp(s):
    return "".join(COMP[b] for b in reversed(s))

def find_orfs(seq, strand="+", min_codons=3):
    orfs = []
    L = len(seq)
    for frame in range(3):
        i = frame
        while i + 3 <= L:
            if seq[i:i + 3] == "ATG":
                j = i
                while j + 3 <= L and seq[j:j + 3] not in STOPS:
                    j += 3
                if j + 3 <= L:                   # a stop codon was found
                    n = (j - i) // 3             # codons before the stop
                    if n >= min_codons:
                        start, end = i + 1, j + 3
                        if strand == "-":        # report on the forward strand
                            start, end = L - end + 1, L - start + 1
                        orfs.append((strand, start, end, n))
                    i = j + 3
                    continue
                break
            i += 3
    return orfs

dna = "CCATGAAAGGGTTTCCCTAAGGCTATTTGGGCCCAAATTTCATGG"
for orf in find_orfs(dna, "+") + find_orfs(revcomp(dna), "-"):
    print(orf)
# ('+', 3, 20, 5)
# ('-', 23, 43, 6)
```

The two results:

| Strand | Positions | Sequence read as a gene | Codons before stop |
|---|---|---|---|
| Forward | 3 to 20 | `ATGAAAGGGTTTCCCTAA` | 5 |
| Reverse | 23 to 43 | `ATGAAATTTGGGCCCAAATAG` (on the reverse strand) | 6 |

The reverse-strand gene is invisible if you only read the forward strand, which is why both strands are always searched.

### Why ORFs alone are not enough

- **Bacterial genes** are continuous and tightly packed, so ORF-based prediction works well, though not every gene starts with `ATG`.
- **Eukaryotic genes** are interrupted by introns, so the coding sequence is split across several exons. Such a gene does not appear as one long ORF, so the predictor must also find the splice boundaries. (Introns typically begin with `GT` and end with `AG`.)

### Where the evidence comes from

The slide says methods **may use DNA alone or add known genes and transcript evidence**:

| Approach | Uses | Strength |
|---|---|---|
| **Ab initio** | The DNA sequence only, with statistical models of what genes look like | Works on any genome, even a new one |
| **Evidence-based** | Alignment of known proteins, known genes, and transcripts (RNA-seq) to the genome | Anchored in real data |
| **Combined** | Both | Usually the most accurate |

Common predictors include Prodigal and Glimmer for bacteria, and AUGUSTUS for eukaryotes.

The output is an initial, organised map of likely genes. The slide's final point is that **predictions guide all downstream analysis**: a gene that is missed here is missed everywhere later, and a false gene is carried along with it.

> **Key idea:** gene prediction turns a raw sequence into a first gene map. Everything that follows depends on it being right.

---

## 7. Structural vs Functional Annotation

Annotation has **two layers**, and they answer different questions.

| | Structural annotation | Functional annotation |
|---|---|---|
| **Question** | Where are the features? | What do they do? |
| **Marks** | Genes, exons, introns, start and stop sites, and other sequence elements | Predicted biological role: enzyme, transporter, regulatory protein, and so on |
| **Slide's picture** | A street map with a "?" over every building: the layout is there, the function is not | Labels saying what each place is for |
| **Result** | The genome's blueprint before meaning is assigned | A reason each feature matters |
| **Typical output** | Coordinates and gene models (GFF/GTF) | Names, descriptions, domains, and functional terms attached to each gene |

The slide's summary: **structure shows where features are; function explains why they matter biologically.** A map with streets but no labels (the slide's "No Function" side) shows the layout but not the story.

### How function is inferred

Most functional annotation is **inference by similarity**, which is why Module 12 matters here:

1. Translate each predicted gene into a protein.
2. Compare the protein with databases using BLAST-style searches (UniProt, from Module 11, is a common target).
3. Search for known domains and motifs (Module 12 §10).
4. Assign the function of the best well-supported match.

When nothing matches, the honest annotation is "hypothetical protein": a real predicted gene whose role is not known.

### The risk

Function by similarity is evidence, not proof. If a database entry was itself assigned by similarity to another entry, an early mistake can be copied from genome to genome. This is the same warning as Module 12's pitfalls: a high score does not guarantee the same function.

> **Key idea:** structural annotation draws the map; functional annotation says what each place is for.

---

## 8. Annotation Pipelines and Tools

Annotating by hand is impossible for a genome with thousands of genes. An **annotation pipeline** is a chain of tools that does the work automatically and at scale.

The slide says pipelines identify genes, predict their structures, and assign likely functions. A typical pipeline runs in the order below:

```
assembled contigs (FASTA)
   -> mask repeats                         (so repeats are not mistaken for genes)
   -> predict genes                        (ORFs, exon-intron boundaries)
   -> predict non-coding genes             (tRNA, rRNA)
   -> compare with databases               (BLAST, domain searches)
   -> assign names and functions
   -> annotated genome (GFF3 plus protein FASTA)
```

| Tool type | What it does | Examples |
|---|---|---|
| **Structural tools** | Find ORFs, exon-intron boundaries, and tRNA genes | Prodigal and Glimmer (bacterial genes), AUGUSTUS (eukaryotic genes), tRNAscan-SE (tRNA genes) |
| **Functional tools** | Compare sequences with databases to infer roles | BLAST, UniProt, domain databases such as Pfam and InterPro |
| **Combined pipelines** | Run structural and functional steps together | Prokka and Bakta (bacteria); MAKER and BRAKER (eukaryotes) |

Pipelines improve **consistency**: every gene in every genome goes through the same steps, so results can be compared. The cost is that they inherit the weaknesses of what feeds them.

The slide's warning is exact: **results depend on assembly quality, databases, and algorithms.** In practice:

- A fragmented assembly can split one gene across two contigs, so the pipeline sees two broken half-genes.
- A missing gene in the assembly can never be annotated.
- A gene whose relatives are absent from the databases will stay "hypothetical".
- Different predictors, run on the same genome, give somewhat different gene sets.

> **Key idea:** pipelines turn raw sequence into biologically meaningful information at scale, but their output is only as good as the assembly, the databases, and the algorithms behind it.

---

## 9. Identifying Regulatory Elements and Non-Coding Regions

Not every meaningful stretch of a genome codes for a protein. In the human genome, protein-coding exons make up only a small percentage of the sequence, while much of the rest has other roles. The slide compares non-coding regions to **the backstage crew in a theatre**: the audience sees the actors (proteins), but the crew decides when and how the show runs.

So annotation is **not only about finding proteins, but also about finding gene-control elements.**

### Regulatory elements

| Element | What it does |
|---|---|
| **Promoter** | The region near the start of a gene where transcription begins |
| **Enhancer** | A region that increases transcription of a gene, and can act from some distance |
| **Silencer** | A region that reduces or switches off transcription |
| **Insulator** | A boundary element that prevents regulatory signals from leaking into neighbouring genes |

### Why they matter

The slide's two panels make the point. In the first, the protein sequence is in the spotlight. In the second, a **regulatory mutation** leaves the protein itself unchanged but alters **when, where, or how much** protein is made. A change in a non-coding region can therefore change development or cause disease without touching any coding sequence.

### Other non-coding features

Some non-coding regions are transcribed into RNA that is functional in its own right: tRNA, rRNA, and many regulatory RNAs.

### How they are found

Regulatory regions have no start codon or ORF to look for, so they are harder to identify. Common clues are:

- **Conservation**: regions that stay similar across species despite not coding for protein are likely to be doing something (an alignment idea from Module 12).
- **Known sequence motifs**, such as binding sites for regulatory proteins.
- **Experimental evidence**, such as which regions are bound by regulatory proteins or are open and accessible in cells.

The slide's conclusion: **finding these regions improves interpretation of development, disease, and variation.**

> **Key idea:** annotation means finding the controls as well as the proteins. A variant outside any gene can still matter if it lands in a regulatory element.

---

## 10. Comparative Genome Annotation

When a genome is new, you do not have to annotate it from nothing. **Comparative annotation** compares it with a **trusted reference genome** that has already been carefully annotated. The slide's analogy: using a labelled city map to understand a similar new city.

### How it works

1. Align the new genome (or its predicted genes) to the reference. This is the alignment of Module 12, used at genome scale.
2. Look for **conserved features**: genes, exon-intron structures, and functional elements that appear in both.
3. Where a close match to a known reference gene exists, transfer its annotation to the new genome.

The slide's picture shows this clearly: two genomes with the same arrangement of gene, exons, and intron, so the structure is **conserved**.

### Why it helps

| Benefit | Why |
|---|---|
| **Faster interpretation** | Known genes are recognised instead of being predicted from scratch |
| **Strong evidence** | A close match to a known reference gene is good support for a gene model |
| **Useful when validation is limited** | When time, budget, or resources do not allow experimental checks |
| **Consistency** | Annotation follows the same conventions across related genomes |
| **Separating biology from noise** | A feature found in both genomes is less likely to be a technical artefact than one seen in a single assembly |

### Limits

- The annotation is only as good as the reference. An error in the reference is copied into every genome annotated from it.
- Genes that are specific to the new organism have no counterpart, so comparison cannot find them.
- The more distant the relationship, the weaker the evidence. A close relative is a good guide; a distant one is only a hint.

> **Key idea:** comparing with a trusted genome makes annotation faster and more consistent, and helps separate real biology from technical noise.

---

## 11. Challenges in Assembling Complex and Repetitive Genomes

Not all genomes are equally easy to assemble. **Repeats and complexity raise the difficulty.**

### Why repeats are a problem

**Repeated DNA removes context.** If the same sequence occurs in many places, a read from inside it matches many genomic locations, and nothing says which one it came from. The slide's analogy is a shredded book with repeated sentences: all the pieces fit, but placement is unclear.

A repeat can be resolved only by a read (or k-mer) **longer than the repeat**, because only then does the piece reach into the unique sequence on each side and anchor itself.

### A worked example in the de Bruijn graph

Take this 22-base genome, in which the repeats `GCG` and `CAC` each occur twice, interleaved:

```
TA GCG TA CAC AT GCG AT CAC TA
   ^^^    ^^^    ^^^    ^^^    <- the repeats: GCG, CAC, GCG, CAC
```

Build the de Bruijn graph (using `build_graph`, `all_paths`, and `to_sequence` from §4) and list every walk that uses every edge once:

```python
# interleaved repeats: GCG and CAC each occur twice, in the order X..Y..X..Y
genome = "TAGCGTACACATGCGATCACTA"
for k in (4, 5):
    graph = build_graph([genome], k)
    candidates = sorted(to_sequence(p) for p in all_paths(graph))
    print(k, len(candidates), candidates)
```

The output:

| k | Number of genomes consistent with the k-mers | Candidates |
|---|---|---|
| 4 | **3** | `TAGCGATCACATGCGTACACTA`, `TAGCGTACACATGCGATCACTA` (the true genome), `TAGCGTACATGCGATCACACTA` |
| 5 | **1** | `TAGCGTACACATGCGATCACTA` |

With k = 4, the k-mers cannot tell which flank belongs to which copy of each repeat, so three different genomes explain the data equally well. With k = 5, a k-mer is long enough to include one unique base on each side of a three-base repeat, so it bridges from unique sequence through the repeat into unique sequence on the other side. Only the true genome remains.

Real repeats are far longer than three bases (transposable elements and tandem repeats of thousands of bases are common), so reads and k-mers must be correspondingly long, or extra long-range information is needed.

### Complexity beyond repeats

Other features also make assembly harder, such as regions with extreme base composition, a genome with high variation between its two copies of a chromosome, or one with more than two copies of each chromosome. All of them create several versions of "the same" sequence.

### Assembly errors

The slide names three:

| Error | What happens | Typical cause |
|---|---|---|
| **Fragmentation** | The genome is broken into many short contigs | Assembler stops at repeats it cannot resolve |
| **Missing regions** | Part of the genome is absent from the assembly | Collapsed repeats (two copies merged into one) and unsequenced difficult regions |
| **False joins** | Pieces that are not neighbours are joined | Repeats misleading the path through the graph |

### Consequences

These errors **distort downstream work, including gene discovery and evolution studies**. A fragmented assembly splits genes. A collapsed repeat removes gene copies, so copy number is wrong. A false join creates a rearrangement that is not real, and a study of genome evolution would report it as biology.

### What helps

- **Long reads** that span repeats.
- Several values of k, or longer k-mers, when there is enough coverage.
- Extra long-range information to order and orient contigs.
- Checking every assembly with N50, completeness, and misassembly tests (§5) rather than trusting any one number.

The slide's conclusion: **good assembly means spotting when repeats and complexity mislead reconstruction.**

> **Key idea:** repeats cannot be resolved by reads shorter than the repeat. Knowing where an assembly may be misleading is part of using it.

---

## 12. Key Takeaways

- Sanger gives accuracy on a small scale; NGS gives throughput on a large scale. Choose by the biological question, scale, time, and budget.
- NGS works because each position is read many times: coverage = N × L / G.
- Short reads are cheap and fast but fail in repeats; long reads span repeats but cost more. Platform choice balances length, throughput, accuracy, and cost.
- De novo assembly builds from overlaps and can find what no reference holds; reference-guided assembly is faster but suffers from reference bias.
- OLC compares whole reads and suits long reads; de Bruijn graphs use k-mers and scale to huge short-read data. A sequencing error can corrupt up to k k-mers.
- N50 measures contiguity only. A high N50 does not guarantee correctness, so read it with total length, contig count, and completeness or accuracy evidence.
- Gene prediction finds likely genes from ORFs, start and stop signals, and other coding features, on both strands; it can use DNA alone or add known genes and transcript evidence.
- Structural annotation says where the features are; functional annotation says what they do, usually by similarity to known sequences.
- Pipelines give consistent annotation at scale but depend on assembly quality, databases, and algorithms.
- Regulatory and other non-coding regions control when, where, and how much a gene is expressed, and changes there can matter even when no protein changes.
- Comparative annotation uses a trusted reference to annotate faster and more consistently, but copies the reference's errors.
- Repeats and complexity cause fragmentation, missing regions, and false joins, and these distort gene discovery and evolution studies.

---

## 13. What's Next

With an assembled and annotated genome in hand, the genes and regulatory regions of an organism are mapped. The next questions are comparative: how do these genes and genomes relate to those of other organisms, and what can sequence differences tell us about function, disease, and evolution?
