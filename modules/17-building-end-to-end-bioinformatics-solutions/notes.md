# Building End-to-End Bioinformatics Solutions

> Module 17: How individual analysis steps become a complete, trustworthy solution — reproducible pipelines, workflow tools (Snakemake and Nextflow), cloud and high-performance computing, Python and R, statistical foundations, machine learning, AI-assisted literature mining, and clinical and translational applications.

---

## Table of Contents

1. [Building Reproducible Bioinformatics Pipelines](#1-building-reproducible-bioinformatics-pipelines)
2. [Workflow Management Tools — Snakemake and Nextflow](#2-workflow-management-tools--snakemake-and-nextflow)
3. [Cloud Computing and High-Performance Computing](#3-cloud-computing-and-high-performance-computing)
4. [Programming Foundations — Python and R](#4-programming-foundations--python-and-r)
5. [Statistical Foundations for Biological Data Analysis](#5-statistical-foundations-for-biological-data-analysis)
6. [Machine Learning Applications in Bioinformatics](#6-machine-learning-applications-in-bioinformatics)
7. [AI-Assisted Literature Mining and Hypothesis Generation](#7-ai-assisted-literature-mining-and-hypothesis-generation)
8. [Clinical and Translational Bioinformatics Applications](#8-clinical-and-translational-bioinformatics-applications)
9. [Key Takeaways](#9-key-takeaways)
10. [What's Next](#10-whats-next)

---

## How to read this module

Modules 10 to 16 taught individual techniques: sequences, databases, alignment, assembly, structure, phylogenetics, and omics analysis. A real project needs all of them to work together, repeatedly, on large data, with results that others can trust. That is an **end-to-end solution**, and this module covers the eight pieces that make one possible.

```
biological question
   -> reproducible pipeline        (§1: inputs, tools, parameters written down)
   -> run by a workflow tool       (§2: Snakemake or Nextflow run the steps in order)
   -> on suitable computers        (§3: laptop, HPC cluster, or cloud)
   -> analysed with code           (§4: Python for data preparation, R for analysis)
   -> checked with statistics      (§5: is it signal or noise?)
   -> mined with machine learning  (§6: patterns, predictions)
   -> placed in the literature     (§7: AI-assisted reading, hypotheses)
   -> applied in the clinic        (§8: biomarkers, precision medicine)
```

The code in this module is Python (using `numpy`, `pandas`, `scipy`, and `scikit-learn` where natural) and was run to produce every printed output. Examples of **Snakemake, Nextflow, Slurm, and R are shown as text for reading and have not been run**, because those tools were not installed on the machine used to check the code. The datasets are small and invented, so that each result can be checked by hand.

---

## 1. Building Reproducible Bioinformatics Pipelines

**A pipeline organizes many analysis steps so that bioinformatics work stays clear and consistent.** The lecture compares it with cooking. A cook who relies on **memory** forgets the amount of salt, the oven time, or which pan was used. A cook who follows a **written recipe** gets the same dish each time, and someone else can cook it too. A reproducible pipeline is the written recipe for an analysis.

### What reproducibility means

A **reproducible** pipeline **uses the same inputs, tools, and settings to get consistent results.** The slide's picture is a staircase of four things that must be fixed:

| Layer | What must be recorded |
|---|---|
| **Input** | Exactly which data files, with a checksum (a fingerprint) to prove they did not change |
| **Tools** | Exactly which programs and which versions |
| **Parameters** | Every setting: thresholds, random seeds, reference genome version |
| **Result** | Follows from the three above, and can therefore be checked |

This matters, as the slide says, **because bioinformatics includes large datasets and many easy-to-forget choices.** A minimum-quality threshold typed into a terminal six months ago is exactly the kind of choice that nobody remembers.

### Why it is worth the effort

**Reproducibility improves accuracy, collaboration, method tracking, and trust in results.** And: **it lets teams rerun workflows with new data instead of rebuilding analyses from scratch.** When new samples arrive, you run the same recipe again.

### The idea in code: a run record and a fingerprint

The code below writes down input, parameters, and tool versions for a small analysis, and computes a short fingerprint of that record. It then shows three things: identical settings give an identical fingerprint, changing one parameter changes both the fingerprint and the result, and a pinned random seed makes a random step repeatable.

```python
import hashlib, json, os, random, tempfile

def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()

def run_analysis(input_path, params):
    """A tiny analysis: keep reads that are long enough and have enough GC, then subsample."""
    reads = [ln.strip() for ln in open(input_path) if ln.strip()]
    kept = [r for r in reads
            if len(r) >= params["min_length"]
            and (r.count("G") + r.count("C")) / len(r) >= params["min_gc"]]
    rng = random.Random(params["seed"])             # a pinned random seed
    return sorted(rng.sample(kept, k=min(params["subsample"], len(kept))))

def run_record(input_path, params, tool_versions):
    """Everything needed to repeat the run, plus a fingerprint of it."""
    record = {"input_sha256": sha256_of(input_path), "params": params, "tools": tool_versions}
    blob = json.dumps(record, sort_keys=True).encode()
    record["fingerprint"] = hashlib.sha256(blob).hexdigest()[:12]
    return record

with tempfile.TemporaryDirectory() as d:
    path = os.path.join(d, "reads.txt")
    toy_reads = ["ATGCGC", "ATATAT", "GGCCGGAA", "GCGCGCGC", "ATGC",
                 "CCGGTTAA", "GGGCCCAT", "TTTTTTTT", "GCATGCAT", "CGCGATAT"]
    open(path, "w").write("\n".join(toy_reads) + "\n")
    tools = {"python": "3.x", "analysis_script": "v1"}
    base = {"min_length": 6, "min_gc": 0.5, "subsample": 3, "seed": 42}

    r1 = run_record(path, base, tools)
    r2 = run_record(path, dict(base), tools)                       # identical settings
    r3 = run_record(path, {**base, "min_gc": 0.6}, tools)          # one parameter changed
    print("run 1 vs run 2 (same settings):", r1["fingerprint"], r2["fingerprint"],
          "-> same" if r1["fingerprint"] == r2["fingerprint"] else "-> different")
    print("run 1 vs run 3 (min_gc 0.5 -> 0.6):", r1["fingerprint"], r3["fingerprint"])

    print("result of run 1:", run_analysis(path, base))
    print("repeat of run 1:", run_analysis(path, base))
    print("result of run 3:", run_analysis(path, {**base, "min_gc": 0.6}))
    print("different seed :", run_analysis(path, {**base, "seed": 7}))

    open(path, "a").write("AAAAGGGC\n")                            # the input file changes
    print("after changing the input file:", run_record(path, base, tools)["fingerprint"])

# Output:
# run 1 vs run 2 (same settings): 520eb7f28ed6 520eb7f28ed6 -> same
# run 1 vs run 3 (min_gc 0.5 -> 0.6): 520eb7f28ed6 be8fd856bc94
# result of run 1: ['ATGCGC', 'CGCGATAT', 'GCATGCAT']
# repeat of run 1: ['ATGCGC', 'CGCGATAT', 'GCATGCAT']
# result of run 3: ['ATGCGC', 'GGCCGGAA', 'GGGCCCAT']
# different seed : ['CCGGTTAA', 'GCGCGCGC', 'GGCCGGAA']
# after changing the input file: 60838ec21407
```

Reading the output:

- Runs 1 and 2 have the **same** fingerprint (`520eb7f28ed6`): the same recipe.
- Changing only `min_gc` from 0.5 to 0.6 gives a different fingerprint, and the analysis returns different reads. A result is only meaningful together with the settings that produced it.
- With the **seed pinned** (42), repeating run 1 returns the same random subsample. A different seed (7) gives a different subsample, so the seed is a parameter too, and belongs in the record.
- Appending one read to the input file changes the fingerprint, which is how you would notice that "the same analysis" was quietly run on different data.

### Practical tools for reproducibility

- **Version control** (Git) for all code and configuration files.
- **Environments and containers** (conda environments, Docker, Singularity/Apptainer) that pin the exact tool versions.
- **Workflow tools** (§2) so that the order of steps and their settings live in files and not in someone's memory.
- **Plain-text records**: a README with the commands, the genome and database versions used, and the date.

> **Key idea:** a reproducible pipeline is a written recipe: the same inputs, tools, and parameters give the same result, and anyone can rerun it when new data arrive.

---

## 2. Workflow Management Tools — Snakemake and Nextflow

A pipeline written as a pile of shell commands still depends on a person to run them in the right order. **Workflow tools organize pipelines into clear, executable steps.** The slide lists what they provide:

- **Each step runs in the right order with correct inputs and outputs.**
- They **prevent manual reruns and reduce chaotic, unsupervised analysis.**
- **Snakemake uses rules linking input files, outputs, and commands.**
- **Nextflow builds portable workflows across laptops, clusters, and cloud.**
- **Both improve automation, reproducibility, reliability, and save debugging time.**

### How a workflow engine thinks

The central idea is a **dependency graph**. Each step says "I need these files and I produce this file." The engine then works backwards from the final result you ask for:

1. To make the report, I need the counts file.
2. To make the counts file, I need the cleaned file.
3. To make the cleaned file, I need the raw file, which already exists.

It runs steps in that order, and **skips any step whose output is already up to date**, which means a rerun after new data repeats only the steps affected. The toy engine below implements that logic in about 25 lines, using files and modification times. It is a teaching model of what Snakemake does, not Snakemake itself.

```python
import os, tempfile, time

class Workflow:
    """A very small make/Snakemake-style engine: rules link input files to output files."""
    def __init__(self):
        self.rules = {}                                  # output file -> (name, inputs, function)

    def rule(self, name, inputs, output, func):
        self.rules[output] = (name, inputs, func)

    def needs_run(self, output):
        if not os.path.exists(output):
            return True
        newest_input = max(os.path.getmtime(i) for i in self.rules[output][1])
        return os.path.getmtime(output) < newest_input    # an input is newer than the output

    def build(self, target, log):
        if target not in self.rules:                      # a raw input file: nothing to build
            return
        name, inputs, func = self.rules[target]
        for i in inputs:                                  # first build what this step needs
            self.build(i, log)
        if self.needs_run(target):
            func(inputs, target)
            log.append(f"ran  {name}")
        else:
            log.append(f"skip {name} (up to date)")

with tempfile.TemporaryDirectory() as d:
    raw, clean, counts, report = (os.path.join(d, f) for f in
                                  ("raw.txt", "clean.txt", "counts.txt", "report.txt"))
    open(raw, "w").write("# header\nGENE1 10\nGENE2 0\n# comment\nGENE3 25\n")

    def clean_step(ins, out):
        lines = [ln for ln in open(ins[0]) if not ln.startswith("#")]
        open(out, "w").writelines(lines)

    def count_step(ins, out):
        total = sum(int(ln.split()[1]) for ln in open(ins[0]))
        open(out, "w").write(f"total_reads {total}\n")

    def report_step(ins, out):
        open(out, "w").write("REPORT: " + open(ins[0]).read())

    wf = Workflow()
    wf.rule("clean", [raw], clean, clean_step)
    wf.rule("count", [clean], counts, count_step)
    wf.rule("report", [counts], report, report_step)

    def run(label):
        log = []
        wf.build(report, log)
        print(label, "->", "; ".join(log))

    run("first run          ")
    run("run again          ")
    earlier = time.time() - 10                            # pretend the results are 10 s old
    for f in (clean, counts, report):
        os.utime(f, (earlier, earlier))
    open(raw, "a").write("GENE4 5\n")                     # new data arrive now
    run("after new data     ")
    os.remove(counts)                                     # an intermediate file is lost
    run("after deleting file")
    print(open(report).read().strip())

# Output:
# first run           -> ran  clean; ran  count; ran  report
# run again           -> skip clean (up to date); skip count (up to date); skip report (up
#   to date)
# after new data      -> ran  clean; ran  count; ran  report
# after deleting file -> skip clean (up to date); ran  count; ran  report
# REPORT: total_reads 40
```

- First run: all three steps execute.
- Run again: all three are skipped, because nothing changed.
- After new data are appended to the raw file: everything downstream is rerun, because every step depends, directly or not, on the raw file.
- After one intermediate file is deleted: only the steps from that file onwards are rerun (the cleaning step is skipped).

This is the "no manual reruns" benefit: the engine decides what needs doing.

### Snakemake: rules linking files

In **Snakemake**, a workflow is a file called `Snakefile` made of **rules**. Each rule has `input`, `output`, and a `shell` command (or script). Names in braces such as `{sample}` are **wildcards**, filled in from the file names requested. The final rule (`all`) lists what you want, and Snakemake works out the rest. *(Shown for reading; not run here.)*

```
# Snakefile
SAMPLES = ["s1", "s2", "s3"]

rule all:                                   # the final targets
    input:
        expand("results/{sample}.summary.txt", sample=SAMPLES)

rule clean:
    input:  "data/{sample}.raw.txt"
    output: "results/{sample}.clean.txt"
    shell:  "grep -v '^#' {input} > {output}"

rule summarise:
    input:  "results/{sample}.clean.txt"
    output: "results/{sample}.summary.txt"
    shell:  "wc -l {input} > {output}"
```

Typical commands: `snakemake -n` (a dry run that only lists what would be done) and `snakemake --cores 4` (run, using up to 4 CPU cores in parallel). The three samples are independent, so Snakemake can run them at the same time.

### Nextflow: processes and channels

**Nextflow** describes a workflow as **processes** (steps) connected by **channels** (streams of files or values). Its strength, as the slide says, is **portability**: the same workflow can run on a laptop, a cluster, or the cloud, with the choice of computer kept in a configuration file and not in the workflow. *(Shown for reading; not run here, and the syntax details depend on the Nextflow version.)*

```
// main.nf
params.reads = "data/*.raw.txt"

process CLEAN {
    input:
    path raw
    output:
    path "${raw.baseName}.clean.txt"
    script:
    """
    grep -v '^#' ${raw} > ${raw.baseName}.clean.txt
    """
}

process SUMMARISE {
    input:
    path clean
    output:
    path "${clean.baseName}.summary.txt"
    script:
    """
    wc -l ${clean} > ${clean.baseName}.summary.txt
    """
}

workflow {
    Channel.fromPath(params.reads) | CLEAN | SUMMARISE
}
```

```
// nextflow.config: choose WHERE the same workflow runs
profiles {
    local { process.executor = 'local' }
    slurm { process.executor = 'slurm' }
    cloud { process.executor = 'awsbatch' }
}
```

`nextflow run main.nf -profile slurm -resume` runs the workflow on a Slurm cluster, and `-resume` reuses the results of steps that have already finished.

### Choosing between them

| | Snakemake | Nextflow |
|---|---|---|
| **Style** | Rules linking files (Python-based syntax) | Processes and channels (Groovy-based syntax) |
| **Strength** | Easy to read, natural for Python users | Very portable across laptops, clusters, and cloud |
| **Both** | Automation, reproducibility, restart after failure, container and environment support | |

Both are widely used, and the choice often depends on the team and on existing pipelines.

> **Key idea:** a workflow tool turns a list of commands into a graph of steps with declared inputs and outputs, runs it in the right order, and reruns only what changed.

---

## 3. Cloud Computing and High-Performance Computing

**Workflow tools need scalable environments to run pipelines at research scale.** A genome project may have hundreds of samples and each sample may need hours of computing. The slide gives two answers:

- **Cloud computing** provides **on-demand storage and compute over the internet.** It lets you **scale without owning hardware full-time**: you rent the computers when you need them.
- **High-performance computing (HPC)** uses **clusters or supercomputers for fast parallel computation**: many computers connected by a fast network and shared among many users through a scheduler.

The lecture's image is a set of pipes (cloud and HPC) feeding the research workflow, and a third pipe leaking, labelled **failure**. That is a fair warning: at scale, some jobs will fail, run out of memory, or hit a time limit, so a pipeline must be able to detect failure and resume (§2).

| | HPC cluster | Cloud |
|---|---|---|
| **Who owns it** | An institution (a university cluster) | A provider, rented by the hour |
| **Cost model** | Usually allocated or free to members | Pay for what you use |
| **Capacity** | Fixed size, queues when busy | Very large, started on demand |
| **Best for** | Steady, heavy research computing | Bursts of demand, sharing, and collaboration |

The lecture's conclusion: **at real-world scale, cloud and HPC are essential infrastructure**, because modern bioinformatics needs **speed, capacity, reproducibility, and practical timelines.**

### What scaling does and does not buy you

There are two different kinds of parallelism:

1. **Many independent jobs** (for example, one alignment per sample). These parallelise perfectly: 100 samples can run on 100 computers at the same time.
2. **One job made faster by more threads.** This is limited by the part of the job that cannot be split up (**Amdahl's law**): if a fraction *f* of the work can be parallelised, the speed-up on *p* threads is 1 / ((1 − *f*) + *f* / *p*).

```python
import math

def wall_time(n_jobs, hours_per_job, cores):
    """Jobs are independent (one sample each): they run in waves of `cores` jobs."""
    return math.ceil(n_jobs / cores) * hours_per_job

def amdahl_speedup(parallel_fraction, threads):
    return 1 / ((1 - parallel_fraction) + parallel_fraction / threads)

# 100 samples, 2 core-hours each, price in $ per core-hour (invented)
n_jobs, hours_per_job, price = 100, 2.0, 0.05
print("cores  wall-clock hours  core-hours billed  cost")
for cores in (1, 10, 25, 50, 100):
    t = wall_time(n_jobs, hours_per_job, cores)
    core_hours = n_jobs * hours_per_job
    print(f"{cores:5d}  {t:16.0f}  {core_hours:17.0f}  ${core_hours * price:.2f}")

print("one tool, 90% of its work parallelisable (Amdahl's law):")
for threads in (1, 2, 4, 8, 16, 32, 64):
    print(f"  {threads:2d} threads -> speed-up {amdahl_speedup(0.90, threads):4.1f}x")
print("upper limit as threads -> infinity:", round(1 / (1 - 0.90), 1), "x")

# Output:
# cores  wall-clock hours  core-hours billed  cost
#     1               200                200  $10.00
#    10                20                200  $10.00
#    25                 8                200  $10.00
#    50                 4                200  $10.00
#   100                 2                200  $10.00
# one tool, 90% of its work parallelisable (Amdahl's law):
#    1 threads -> speed-up  1.0x
#    2 threads -> speed-up  1.8x
#    4 threads -> speed-up  3.1x
#    8 threads -> speed-up  4.7x
#   16 threads -> speed-up  6.4x
#   32 threads -> speed-up  7.8x
#   64 threads -> speed-up  8.8x
# upper limit as threads -> infinity: 10.0 x
```

Two lessons stand out:

- Because the jobs are independent, 100 cores finish in 2 hours instead of 200, and **the number of core-hours billed is the same** (200), so in this simple model the cost is the same too. This is the real attraction of the cloud: you can buy time instead of waiting for it. (The price is invented, and real bills also include storage, data transfer, and idle time.)
- For one tool with 10% serial work, 64 threads give only 8.8 times the speed, and the ceiling is 10 times however many threads you add. More hardware is not a cure for every slow step.

### An HPC job in practice

On a cluster, you do not run heavy jobs directly: you describe the resources you need and submit the job to the **scheduler** (commonly Slurm), which starts it when resources are free. A **job array** runs the same script for many samples. *(Shown for reading; not run here.)*

```
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-100                 # 100 jobs, numbered 1 to 100

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)     # the Nth sample name
bwa mem -t $SLURM_CPUS_PER_TASK ref.fa data/${sample}_R1.fq.gz data/${sample}_R2.fq.gz \
    > results/${sample}.sam
```

You submit it with `sbatch align.sh` and watch it with `squeue`. A workflow tool such as Snakemake or Nextflow can do this submission for you, which is why the two pieces are taught together.

> **Key idea:** research-scale analysis needs infrastructure. Cloud rents capacity on demand and HPC shares a cluster, and a workflow tool lets the same pipeline run on either.

---

## 4. Programming Foundations — Python and R

**Python and R are the two main languages used in bioinformatics.** They are strong in different places, which the lecture sums up as **Python for data preparation and R for analysis**.

| | Python | R |
|---|---|---|
| **Strong at** | **Automation, data cleaning, pipelines, and file handling** | **Statistics, visualization, and result exploration** |
| **Typical libraries** | pandas, NumPy, SciPy, scikit-learn, Biopython | ggplot2 and the Bioconductor packages (such as DESeq2 and limma) |
| **Role in a project** | **Organizes workflows** | **Helps interpret data and communicate findings** |

Bioinformatics **often uses both languages for different stages of analysis**, and **the key skill is knowing where each language adds the most value.** The usual hand-over is a plain file (such as CSV) that Python writes and R reads.

### Python for data preparation

Real data have missing values. The example below prepares a small expression table (invented numbers) with missing entries. It is written as a **function with a docstring and error checks**, which is what moves code from a one-off experiment toward a reusable pipeline step.

```python
import numpy as np
import pandas as pd

# invented expression table with missing values; SampleID is text, the gene columns are numbers
raw = pd.DataFrame({
    "SampleID": ["S1", "S2", "S3", "S4", "S5", "S6"],
    "GeneA": [100, 150, np.nan, 200, 120, 180],
    "GeneB": [np.nan, 300, 250, 400, 280, np.nan],
    "GeneC": [500, 600, 700, np.nan, 650, 550],
})

def clean_sequencing_data(df, strategy="mean"):
    """Return a cleaned copy of the table. Strategies: 'mean', 'median', 'drop'."""
    if "SampleID" not in df.columns:
        raise ValueError("table needs a SampleID column")
    genes = df.select_dtypes(include="number").columns       # work on the numeric columns only
    out = df.copy()
    if strategy == "mean":
        out[genes] = out[genes].fillna(out[genes].mean())
    elif strategy == "median":
        out[genes] = out[genes].fillna(out[genes].median())
    elif strategy == "drop":
        out = out.dropna(subset=genes)
    else:
        raise ValueError(f"unknown strategy: {strategy}")
    return out

print("missing values per column:\n", raw.isna().sum().to_string())
for strategy in ("mean", "median", "drop"):
    cleaned = clean_sequencing_data(raw, strategy)
    genes = ["GeneA", "GeneB", "GeneC"]
    print(f"\nstrategy = {strategy}: {len(cleaned)} samples kept")
    print(cleaned[genes].agg(["mean", "std"]).round(1).to_string())

print("\nstd of GeneB using only the measured values:", round(raw["GeneB"].std(), 1))

cleaned = clean_sequencing_data(raw, "median")
cleaned.to_csv("cleaned_sequencing_data.csv", index=False)      # hand-over file for R
print(open("cleaned_sequencing_data.csv").read().splitlines()[:3])

# Output:
# missing values per column:
#  SampleID    0
# GeneA       1
# GeneB       2
# GeneC       1
# 
# strategy = mean: 6 samples kept
#       GeneA  GeneB  GeneC
# mean  150.0  307.5  600.0
# std    36.9   50.3   70.7
# 
# strategy = median: 6 samples kept
#       GeneA  GeneB  GeneC
# mean  150.0  301.7  600.0
# std    36.9   51.2   70.7
# 
# strategy = drop: 2 samples kept
#       GeneA  GeneB  GeneC
# mean  135.0  290.0  625.0
# std    21.2   14.1   35.4
# 
# std of GeneB using only the measured values: 65.0
# ['SampleID,GeneA,GeneB,GeneC', 'S1,100.0,290.0,500.0', 'S2,150.0,300.0,600.0']
```

Four points are worth noticing:

- **The text column matters.** A common first attempt is `df.fillna(df.mean())`. On the current pandas used here (3.0.5), this fails with a `TypeError` because the `SampleID` column holds text and cannot be averaged. The fix is to work on the numeric columns only, as the function does with `select_dtypes`, or to pass `numeric_only=True`. (As far as I know, older pandas versions skipped text columns silently, so code written for them can break after an upgrade.)
- **Mean filling hides uncertainty.** For `GeneB` the standard deviation of the measured values is 65.0, but after mean filling it shrinks to 50.3, because the invented values sit exactly in the middle. Downstream tests then see less variation than really exists and become over-confident.
- **Dropping rows loses data.** Removing every sample with any missing value leaves only 2 of 6 samples.
- **There is no universally correct strategy.** For real counts, a missing value and a zero are different things, and each assay has its own accepted handling, so the choice should be made deliberately and recorded (§1).

The slide's suggestions for improving such code are exactly the three habits shown: use **NumPy** (vectorised operations, instead of Python loops) for large numerical data, **modularise** the code into functions, and **add error handling** for file reading, writing, and missing data.

### R for analysis

Once the cleaned table is saved to CSV, R takes over for statistics and plots. *(Shown for reading; R was not available here, so this was not run.)*

```
df <- read.csv("cleaned_sequencing_data.csv")
df$group <- c("control", "control", "control", "treated", "treated", "treated")

t.test(GeneC ~ group, data = df)                    # compare the two groups
boxplot(GeneC ~ group, data = df, ylab = "GeneC expression")
```

> **Key idea:** Python prepares and organises, R explores and communicates. Write each step as a function with error handling, record the choices you make, and know which language fits which stage.

---

## 5. Statistical Foundations for Biological Data Analysis

**Statistics helps interpret the outputs of Python and R.** The slide's warning is that **biological data may look persuasive but still reflect noise, not signal.** Its picture shows two bar charts, Group A and Group B, with one taller, and a shield labelled "Statistical Test" blocking a banner that reads "Breakthrough!". The rule is that **differences in expression, variants, or proteins must be tested before trust.**

**Statistical reasoning compares groups, measures variability, and estimates confidence.** The three demonstrations below show what that means in practice. (Multiple-testing correction for thousands of genes is covered in Module 16, §8.)

### 1. Is a difference bigger than shuffling would produce?

A **permutation test** answers this directly. Shuffle the group labels many times, and count how often the shuffled groups differ at least as much as the real ones.

### 2. How likely is a real difference to be detected? (power)

The **power** of an experiment is the probability that it detects a difference that really exists. It depends on the size of the difference, the variability, and the **number of samples**. The simulation below draws two groups from normal distributions whose true means differ by 1 standard deviation, and counts how often a Welch t-test gives p < 0.05.

### 3. What happens when there is no real difference?

If the groups are identical, about 5% of experiments will still be "significant" by chance, which is what p < 0.05 means.

```python
import numpy as np
from scipy import stats

rng = np.random.default_rng(2024)

# ---- 1. A permutation test: how often does shuffling the labels give a difference this big? ----
group_a = np.array([5.1, 4.8, 5.6, 5.0, 5.3, 4.9])
group_b = np.array([5.9, 6.2, 5.7, 6.4, 5.5, 6.0])
observed = group_b.mean() - group_a.mean()
pooled = np.concatenate([group_a, group_b])
count = 0
n_perm = 20000
for _ in range(n_perm):
    rng.shuffle(pooled)
    if abs(pooled[6:].mean() - pooled[:6].mean()) >= abs(observed):
        count += 1
print(f"observed difference = {observed:.2f}, permutation p = {count / n_perm:.4f}, "
      f"Welch t-test p = {stats.ttest_ind(group_b, group_a, equal_var=False)[1]:.4f}")

# ---- 2. Statistical power: how often is a TRUE difference detected? ----
def detection_rate(n, effect, n_sim=4000, alpha=0.05):
    hits, est = 0, []
    for _ in range(n_sim):
        a = rng.normal(0, 1, n)
        b = rng.normal(effect, 1, n)
        p = stats.ttest_ind(b, a, equal_var=False)[1]
        if p < alpha:
            hits += 1
            est.append(b.mean() - a.mean())
    return hits / n_sim, (np.mean(est) if est else float("nan"))

print("true difference = 1.0 standard deviation")
print("n per group   power   average estimated difference among significant results")
for n in (3, 5, 10, 20, 50):
    power, est = detection_rate(n, 1.0)
    print(f"{n:11d}   {power:5.2f}   {est:5.2f}")

# ---- 3. No true difference at all: how often do we still get a 'breakthrough'? ----
power0, _ = detection_rate(5, 0.0, n_sim=20000)
print("no true difference, n = 5 per group: significant in", f"{power0:.1%}", "of experiments")

# Output:
# observed difference = 0.83, permutation p = 0.0044, Welch t-test p = 0.0009
# true difference = 1.0 standard deviation
# n per group   power   average estimated difference among significant results
#           3    0.12    1.91
#           5    0.25    1.66
#          10    0.54    1.30
#          20    0.87    1.09
#          50    1.00    1.00
# no true difference, n = 5 per group: significant in 4.3% of experiments
```

Reading the output:

- **Permutation test.** The two groups differ by 0.83, and shuffling gives a difference that large only 0.44% of the time (p = 0.0044). The Welch t-test gives p = 0.0009, the same conclusion. The two methods ask the question differently, so the numbers do not have to match.
- **Power depends heavily on sample size.** With 3 samples per group, a real difference of 1 standard deviation is found only 12% of the time. With 10 per group it is 54%, with 20 it is 87%, and with 50 it is essentially certain. Small experiments are not just weaker: they mostly miss real effects.
- **The "winner's curse."** Look at the last column. Among the experiments that *did* reach significance with only 3 samples per group, the average estimated difference is 1.91, almost double the true value of 1.0. Small studies that succeed tend to overestimate the effect, which is one reason first results often shrink when repeated.
- **Without any real difference**, n = 5 per group still gave a "significant" result in 4.3% of 20,000 experiments: close to the 5% the threshold is designed to allow (a Welch test is slightly conservative at small sample sizes).

### The lecture's point about responsibility

**Weak conclusions in research or healthcare can waste samples, time, and budget.** The closing line, **responsible bioinformatics requires both computation and evidence-based judgment**, reads the picture as follows: the code can produce a "Breakthrough!" in seconds, and the test is the check that decides whether to believe it.

A few habits follow:

- Decide the comparison and the threshold **before** looking at the results.
- Report the **size** of the effect and a measure of **uncertainty** (such as a confidence interval), not just a p-value.
- Use enough **biological replicates**, and know the power of the design before the experiment.
- Correct for multiple testing when many things are tested (Module 16, §8).

> **Key idea:** statistics separates signal from noise. A convincing-looking difference must be tested, small samples miss real effects and exaggerate the ones they find, and the evidence must be judged, not just computed.

---

## 6. Machine Learning Applications in Bioinformatics

**Machine learning (ML) learns patterns from examples, not just from fixed hand-written rules.** Instead of programming "if gene X is high, call the sample a tumour", you show the algorithm many labelled examples and let it find the rule. This suits biology because **bioinformatics data is noisy, high-dimensional, and hard to interpret manually**: thousands of genes measured on dozens of samples.

**ML supports** the tasks the slide names:

- **disease classification** (tumour or normal, subtype),
- **gene prediction** (finding genes in a genome),
- **biomarker discovery** (which measurements predict an outcome),
- **imaging** (for example, pathology slides),
- **drug discovery**.

The lecture describes it as **a fast research assistant, but poor training leads to poor predictions**, and as a way to **reveal hidden patterns that humans may miss in very large biological datasets.** The goal is to **turn biological data into useful predictions and biological insight.**

### The standard workflow

1. Collect **labelled examples** (for example, samples with a known diagnosis).
2. Split the data: **training** data to learn from, and **testing** data kept aside to measure real performance.
3. Train a model (for example, a logistic regression) on the training data.
4. **Evaluate on data the model has never seen.** Training accuracy is not the test.
5. Interpret and validate the model before relying on it.

### "Poor training leads to poor predictions": a demonstration

The most common mistake in biological ML comes from a problem of shape. With **40 samples and 2,000 genes**, it is easy to find genes that, just by chance, separate the two groups. The code below creates data of **pure noise** and labels with **no relationship** to it, so the only honest accuracy is about 50%. It then makes the same classifier look excellent by one wrong step: choosing the "best" 20 genes using *all* samples, and only then cross-validating. That is called **data leakage**, because information about the test samples leaked into the training step.

```python
import numpy as np
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline

rng = np.random.default_rng(5)
n_samples, n_features = 40, 2000
X = rng.normal(size=(n_samples, n_features))        # 40 patients x 2000 'genes' of pure noise
y = np.array([0] * 20 + [1] * 20)                   # labels with NO relationship to X
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

# WRONG: choose the 20 best genes using ALL samples, then cross-validate
X_leaky = SelectKBest(f_classif, k=20).fit_transform(X, y)
leaky = cross_val_score(LogisticRegression(max_iter=1000), X_leaky, y, cv=cv).mean()

# RIGHT: do the gene selection inside every training fold only
pipe = make_pipeline(SelectKBest(f_classif, k=20), LogisticRegression(max_iter=1000))
honest = cross_val_score(pipe, X, y, cv=cv).mean()
print(f"pure noise, selection before cross-validation : accuracy = {leaky:.2f}")
print(f"pure noise, selection inside cross-validation  : accuracy = {honest:.2f}")

# Now plant a real signal in 10 of the 2000 genes and repeat the honest analysis
X_real = X.copy()
X_real[y == 1, :10] += 1.5
honest_real = cross_val_score(pipe, X_real, y, cv=cv).mean()
print(f"10 informative genes, selection inside cross-validation: accuracy = {honest_real:.2f}")

# a held-out test set gives the same honest answer
perm = rng.permutation(n_samples)
train, test = perm[:28], perm[28:]
model = make_pipeline(SelectKBest(f_classif, k=20), LogisticRegression(max_iter=1000))
model.fit(X_real[train], y[train])
accuracy = model.score(X_real[test], y[test])
print(f"held-out test accuracy on {len(test)} unseen samples = {accuracy:.2f}")

# Output:
# pure noise, selection before cross-validation : accuracy = 0.95
# pure noise, selection inside cross-validation  : accuracy = 0.57
# 10 informative genes, selection inside cross-validation: accuracy = 0.93
# held-out test accuracy on 12 unseen samples = 1.00
```

Reading the output:

| Analysis | Accuracy |
|---|---|
| Pure noise, genes selected **before** cross-validation (leaky) | **0.95** |
| Pure noise, genes selected **inside** each cross-validation fold (honest) | 0.57 |
| 10 truly informative genes, honest pipeline | 0.93 |

The leaky analysis reports 95% accuracy on data that contain **no signal at all.** The honest analysis gives 57%, close to the 50% chance level, which is the truth. When a real signal exists (10 informative genes), the honest pipeline finds it (0.93), so honesty does not stop good models from working. The last line, accuracy on 12 unseen samples, agrees, but with only 12 samples it is a very uncertain number.

### Practical cautions

- **Always evaluate on data not used in any step of training**, including gene selection, scaling, and tuning. Use a pipeline so that every step is fitted on the training data only.
- **Class balance and metrics.** If 95% of samples are healthy, a model that always says "healthy" is 95% accurate and useless. Look at sensitivity, specificity, and similar measures (§8).
- **Batch effects.** A model can learn which machine or hospital a sample came from, instead of the biology.
- **Interpretability.** In biology and medicine, knowing *why* a model predicts something often matters as much as the prediction.

> **Key idea:** ML finds patterns, including false ones. With many more features than samples, evaluate honestly on unseen data and keep every data-dependent step inside the training loop.

---

## 7. AI-Assisted Literature Mining and Hypothesis Generation

**AI literature mining is a practical bioinformatics use of machine learning.** The slide's reason: **biology outputs too much literature for any researcher to track manually.** Software that **scans text to find genes, proteins, diseases, pathways, and links** acts as **a fast research assistant, not a replacement for scientists**, and it **can reveal indirect patterns that inspire new biological hypotheses.**

### Step 1: find the entities

The first task is **named-entity recognition**: spotting which words in a sentence are genes, proteins, diseases, or pathways. Real systems use dictionaries of official names and trained language models, and must handle synonyms and abbreviations.

### Step 2: find the links

The simplest approach to links is **co-occurrence**: if two entities appear in the same sentence, they may be related. More advanced methods extract the relation itself (for example, "inhibits" or "is expressed in").

### Step 3: find the indirect links (hypothesis generation)

The classic idea of **literature-based discovery**, due to Don Swanson in 1986, is the **A–B–C model**. If one body of literature says A is linked to B, and a different body says B is linked to C, but nobody has connected A to C, then **A–C is a candidate hypothesis.** Swanson's famous example connected dietary fish oil with Raynaud's syndrome through shared intermediates (blood viscosity, platelet aggregation, and vascular reactivity), in two bodies of literature that had never cited each other.

The toy example below follows that logic on 12 **invented** sentences. The gene and disease names are made up, so nothing here is a biological finding.

```python
from collections import defaultdict
from itertools import combinations

# Invented sentences. Gene and disease names are made up, so nothing here is a real finding.
corpus = [
    "Patients with neuropathy-Z show raised oxidative stress in nerve tissue.",
    "Impaired iron transport was observed in neuropathy-Z.",
    "Oxidative stress damages nerve cells and slows repair.",
    "PLM3 regulates iron transport in liver cells.",
    "Loss of PLM3 increases oxidative stress in liver cells.",
    "KRN2 controls iron transport in the gut.",
    "KRN2 is expressed in fibrosis-Y tissue.",
    "TRX1 showed no association with fibrosis-Y in this cohort.",
    "TRX1 is a marker of inflammation in fibrosis-Y.",
    "Inflammation is increased in fibrosis-Y.",
    "QOR4 is involved in inflammation of the skin.",
    "Dermatitis-X patients show strong inflammation.",
]
genes = ["PLM3", "KRN2", "TRX1", "QOR4"]
others = ["neuropathy-Z", "fibrosis-Y", "dermatitis-X",
          "iron transport", "oxidative stress", "inflammation"]
terms = genes + others

def mentions(sentence):
    s = sentence.lower()
    return [t for t in terms if t.lower() in s]

cooc = defaultdict(int)
for sentence in corpus:
    for a, b in combinations(sorted(mentions(sentence)), 2):
        cooc[(a, b)] += 1

def together(a, b):
    return cooc.get(tuple(sorted((a, b))), 0)

print("direct co-mentions involving neuropathy-Z:",
      {t: together("neuropathy-Z", t) for t in terms if together("neuropathy-Z", t)})
print("TRX1 and fibrosis-Y co-mentioned:", together("TRX1", "fibrosis-Y"),
      "times (one of them says 'no association')")

# Indirect leads (the 'A-B-C' idea): a gene never mentioned with the disease,
# but linked to the same intermediate terms as the disease.
disease = "neuropathy-Z"
mid_of_disease = {t for t in others if t != disease and together(disease, t)}
print("terms linked to the disease:", sorted(mid_of_disease))
leads = []
for g in genes:
    if together(g, disease) == 0:
        shared = sorted(mid_of_disease & {t for t in others if together(g, t)})
        if shared:
            leads.append((len(shared), g, shared))
for n, g, shared in sorted(leads, reverse=True):
    print(f"lead: {g} ({n} shared intermediates: {', '.join(shared)})")

# Output:
# direct co-mentions involving neuropathy-Z: {'iron transport': 1, 'oxidative stress': 1}
# TRX1 and fibrosis-Y co-mentioned: 2 times (one of them says 'no association')
# terms linked to the disease: ['iron transport', 'oxidative stress']
# lead: PLM3 (2 shared intermediates: iron transport, oxidative stress)
# lead: KRN2 (1 shared intermediates: iron transport)
```

Reading the output:

- The disease `neuropathy-Z` is directly linked to two intermediates, `iron transport` and `oxidative stress`.
- The gene `PLM3` is **never mentioned with the disease**, but it is linked to *both* of the same intermediates, so it becomes the top lead. `KRN2` shares only one intermediate and ranks lower. The genes `TRX1` and `QOR4` share none and are not proposed.
- That is a **hypothesis, not a finding**: it says where a human should look, nothing more.

### Why humans must verify: "AI suggests leads"

The second line of the output shows the problem. `TRX1` and `fibrosis-Y` are counted as co-mentioned twice, but one of those sentences says there was **no association**. A simple co-occurrence counter cannot read negation, and the lecture's closing line applies: **AI suggests leads; humans must verify evidence and test feasibility.** Other pitfalls:

- **Names are ambiguous.** Gene symbols can match ordinary words or other genes, and one gene has many aliases, so names should be mapped to standard identifiers.
- **Frequency is not truth.** A claim repeated in many papers may be repeated from one source, and a published claim may be wrong.
- **Language models can invent.** A general-purpose language model can write fluent summaries and can also produce references or facts that do not exist, so every claim it gives should be traced to the original paper.
- **A lead still needs an experiment.** Even a well-supported lead has to be tested in the lab or in data (feasibility).

> **Key idea:** literature mining reads more than any person can and connects what no single paper connects. It proposes leads, and the scientist checks the evidence and tests the idea.

---

## 8. Clinical and Translational Bioinformatics Applications

**Clinical and translational bioinformatics connects biological data to real healthcare use.** "Translational" means moving findings **from lab discoveries to bedside decisions.** The slide's picture reads left to right: DNA and digital data (`0101`), then **interpretation**, then **better outcomes** for the patient. The field **links genomic, molecular, and clinical information** (for example, sequencing results together with the patient's records), and **translates research complexity into clinical action.**

### Precision medicine and biomarkers

**Precision medicine uses mutations to guide targeted therapy**: the treatment is chosen according to the molecular profile of the patient or the tumour. Well-known examples (from general knowledge, not from the lecture): tumours with HER2 amplification can be treated with HER2-targeted drugs, lung cancers with certain EGFR mutations with EGFR inhibitors, and cancers with BRCA1 or BRCA2 mutations with PARP inhibitors.

The slide's last line: **biomarkers can predict risk, drug response, and prognosis.**

| Biomarker role | Question it answers |
|---|---|
| **Risk** | How likely is this person to develop the disease? |
| **Drug response** | Will this person respond to this treatment? |
| **Prognosis** | How is the disease likely to progress? |

### A good test can still mislead: the role of prevalence

A biomarker test is described by two numbers: **sensitivity** (the fraction of sick people it correctly flags) and **specificity** (the fraction of healthy people it correctly clears). What a patient wants to know, though, is: "the test is positive; how likely am I to be ill?" That is the **positive predictive value (PPV)**, and it depends on how common the disease is.

The code takes a test with 90% sensitivity and 95% specificity (it looks excellent) and applies it to 10,000 people at three different disease prevalences. It then analyses invented trial data that ask whether a biomarker predicts drug response, and finishes with a toy mutation-to-therapy lookup.

```python
from scipy import stats

def confusion(n, prevalence, sensitivity, specificity):
    diseased = round(n * prevalence)
    healthy = n - diseased
    tp = round(diseased * sensitivity); fn = diseased - tp
    tn = round(healthy * specificity);  fp = healthy - tn
    return tp, fn, fp, tn

sens, spec = 0.90, 0.95            # a biomarker test that looks excellent
print("prevalence  TP    FN    FP     TN  | PPV    NPV")
for prevalence in (0.01, 0.10, 0.30):
    tp, fn, fp, tn = confusion(10000, prevalence, sens, spec)
    ppv, npv = tp / (tp + fp), tn / (tn + fn)
    print(f"{prevalence:9.0%}  {tp:4d}  {fn:4d}  {fp:4d}  {tn:5d}  | {ppv:.3f}  {npv:.4f}")

# Does a biomarker predict drug response? (invented trial data)
#                  responded  did not respond
pos = [18, 12]     # biomarker-positive patients
neg = [6, 34]      # biomarker-negative patients
odds_ratio, p = stats.fisher_exact([pos, neg])
print(f"response rate: biomarker+ {pos[0]/sum(pos):.0%}, biomarker- {neg[0]/sum(neg):.0%}")
print(f"odds ratio = {odds_ratio:.1f}, Fisher exact p = {p:.4f}")

# From mutation to therapy: a lookup table is the simplest form of precision-medicine logic
guidelines = {                       # entirely hypothetical gene -> drug-class rules
    ("GENE_A", "V600X"): "therapy class 1",
    ("GENE_B", "amplification"): "therapy class 2",
}
patient_variants = [("GENE_A", "V600X"), ("GENE_C", "R12Q")]
for v in patient_variants:
    print(v, "->", guidelines.get(v, "no matching guideline: needs expert review"))

# Output:
# prevalence  TP    FN    FP     TN  | PPV    NPV
#        1%    90    10   495   9405  | 0.154  0.9989
#       10%   900   100   450   8550  | 0.667  0.9884
#       30%  2700   300   350   6650  | 0.885  0.9568
# response rate: biomarker+ 60%, biomarker- 15%
# odds ratio = 8.5, Fisher exact p = 0.0001
# ('GENE_A', 'V600X') -> therapy class 1
# ('GENE_C', 'R12Q') -> no matching guideline: needs expert review
```

Reading the output:

- **At 1% prevalence**, the 10,000 people include 100 who are ill. The test flags 90 of them (true positives) but also **495 healthy people** (false positives), so only 90 of 585 positive results are real: **PPV = 15.4%.** About 85% of positive results are false alarms, even with an excellent test. A negative result, however, is very reliable (NPV 99.9%).
- At 10% prevalence the PPV is 66.7%, and at 30% it is 88.5%. The same test is much more informative in a population where the condition is common, which is why screening tests are interpreted with the population in mind.
- **Drug response.** In the invented trial, 60% of biomarker-positive patients responded, against 15% of biomarker-negative patients: an odds ratio of 8.5 (Fisher's exact test, p = 0.0001). A difference this large and well tested is what makes a biomarker clinically useful. It would still need confirmation in an independent group of patients (§5, §6).
- **Lookup rules.** The final lines show the simplest form of precision-medicine logic: a known variant maps to a guideline, and an unknown variant gets "needs expert review". Real systems (using resources such as ClinVar, as in Module 16, §6) have far more rules, and the unmatched cases are the hard ones.

### What clinical use demands

Because decisions affect patients, clinical bioinformatics asks more than research does:

- **Validated, versioned pipelines**: this is §1 and §2 with stricter rules, so that every reported result can be traced to the exact data, tools, and parameters.
- **Careful interpretation** with standard evidence classes (Module 16, §6) and expert review of uncertain variants.
- **Statistical and ML validation** on independent data (§5, §6) before a biomarker or model is used.
- **Communication** that a clinician can act on.

> **Key idea:** clinical bioinformatics turns data into decisions for patients, so every earlier piece (reproducible pipelines, sound statistics, honest validation) matters more. Even an excellent test can mislead when the condition is rare.

---

## 9. Key Takeaways

- An end-to-end solution is a chain: reproducible pipeline, workflow tool, suitable computing, code, statistics, machine learning, literature, and clinical application. The chain is only as trustworthy as its weakest link.
- A reproducible pipeline fixes inputs, tools, and parameters in writing (and in files, containers, and version control), so the same recipe gives the same result and reruns easily on new data.
- Workflow tools (Snakemake with rules, Nextflow with processes and channels) run steps in dependency order, skip what is up to date, resume after failures, and make the same pipeline portable across laptop, cluster, and cloud.
- Cloud and HPC provide the capacity for research-scale data. Independent jobs scale almost perfectly, but a single tool is limited by its serial part (Amdahl's law), and failures must be expected.
- Python organises and prepares data (pandas, NumPy), and R explores and communicates (statistics, plots). Write steps as functions with error handling, and record data-cleaning choices, since mean filling shrinks variance and dropping rows loses samples.
- Statistics decides whether a difference is signal: use tests, effect sizes with uncertainty, enough replicates, and multiple-testing correction. Low power misses real effects and exaggerates those it finds (the winner's curse).
- Machine learning finds patterns, including false ones. With far more features than samples, keep all selection and tuning inside the training loop, and evaluate on unseen data. Leakage turned a 50% model into a "95%" model.
- AI literature mining finds entities and links across too much text and suggests hypotheses (the A–B–C idea). It cannot read negation reliably, and humans must verify every lead.
- Clinical and translational bioinformatics links genomic, molecular, and clinical data to decisions. Biomarkers predict risk, response, and prognosis, and a test's predictive value depends on how common the condition is (PPV 15% at 1% prevalence for a test with 90% sensitivity and 95% specificity).

---

## 10. What's Next

The lecture outline continues with data privacy, compliance, and the responsible use of genomic data; building a bioinformatics portfolio and a career roadmap; emerging trends in precision medicine and personalised genomics; and a capstone on designing an end-to-end analysis plan. Those topics were not part of the slides these notes were built from, so they are not covered here and will be added when the material is available.
