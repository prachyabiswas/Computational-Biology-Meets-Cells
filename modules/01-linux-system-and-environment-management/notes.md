# Linux & Conda for Bioinformatics — Notes

> Source: Singleron Free Online Bioinformatics Academy — Module 1 (Linux System & Environment Management)
> Compiled for personal study and GitHub documentation.

---

## Table of Contents

1. [What is Linux?](#1-what-is-linux)
2. [Linux Users & Groups](#2-linux-users--groups)
3. [Directory Structure](#3-directory-structure)
4. [Absolute vs Relative Paths](#4-absolute-vs-relative-paths)
5. [File & Directory Naming Conventions](#5-file--directory-naming-conventions)
6. [Common Commands](#6-common-commands)
7. [Viewing File Contents](#7-viewing-file-contents)
8. [Creating Files](#8-creating-files)
9. [Deleting Files & Directories](#9-deleting-files--directories)
10. [Text Editors — Vim & Nano](#10-text-editors--vim--nano)
11. [Permissions (chmod)](#11-permissions-chmod)
12. [Compression & Archiving (tar, zip)](#12-compression--archiving-tar-zip)
13. [Wildcard Characters](#13-wildcard-characters)
14. [Pipes & Output Redirection](#14-pipes--output-redirection)
15. [Chaining Multiple Commands](#15-chaining-multiple-commands)
16. [Task Submission (nohup, qsub)](#16-task-submission-nohup-qsub)
17. [Conda: Introduction](#17-conda-introduction)
18. [Why Conda Exists](#18-why-conda-exists)
19. [Installing Miniconda](#19-installing-miniconda)
20. [Conda Channels](#20-conda-channels)
21. [Conda Environments](#21-conda-environments)
22. [Common Conda Commands](#22-common-conda-commands)
23. [Common Bioinformatics Software](#23-common-bioinformatics-software)
24. [Hands-on Exercise (with solutions)](#24-hands-on-exercise-with-solutions)

---

## 1. What is Linux?

Linux is a **free and open-source, UNIX-like operating system**. Its kernel was first released by Linus Torvalds in 1991.

Key properties:
- As long as the **GNU General Public License (GPL)** is followed, the entire source code is freely available, modifiable, and redistributable.
- Supports **multiple platforms**.
- Supports **multiple users and multiple tasks** at the same time (unlike a typical single-user desktop OS).

**Linux vs Windows (conceptual difference):**

| Aspect | Windows | Linux |
|---|---|---|
| Navigation | Mouse (click, double-click, drag-drop) | Keyboard-driven commands |
| Storage reference | Disk C, Disk D, named folders | Everything starts from root `/` |
| File discovery | Browsing folders | `cd`, `ls`, path strings |

Why Linux matters in bioinformatics: almost all sequencing/alignment/analysis tools (BWA, STAR, Trinity, CeleScope, etc.) are built to run on Linux servers/clusters, not Windows.

---

## 2. Linux Users & Groups

- Every **process** running on the system is run by a particular **user**.
- Every **file** is owned by a particular user.
- Access to files/directories is **limited** by user identity.
- A **superuser** (root) sits above all regular users and groups — root has unrestricted access.
- Users can be organized into **groups**, and permissions can be granted at the group level.

This matters because on a shared bioinformatics server (e.g., a lab's compute cluster), many people log in as different users, and permissions decide who can read/write/execute what.

---

## 3. Directory Structure

**Core idea: "Everything in Linux is a file."** Even directories and devices are represented as files.

The filesystem is a single tree starting at the **root directory `/`**. Standard top-level directories include:

```
/
├── bin      → essential user binaries (commands)
├── boot     → boot loader files
├── dev      → device files
├── etc      → system configuration files
├── home     → personal directories for each user (e.g., /home/alice)
├── root     → home directory of the root user
├── run      → runtime data
├── sbin     → system binaries (admin commands)
├── tmp      → temporary files
├── usr      → user-installed programs/libraries
└── var      → variable data (logs, etc.)
```

Example from a real shell session:

```bash
$ ls
bin  data  lib    mnt  run   srv  tmp  usr
boot dev   lib64  opt  sbin  sys  var  home  media  root  proc  lost+found
```

---

## 4. Absolute vs Relative Paths

**Absolute path**: starts from the root `/` and lists every directory you pass through. It always works regardless of where you currently are.

```bash
$ pwd
/mnt/sdd/singleron_lecturer/Windrunner/singleron/software
```

**Relative path**: does *not* start with `/`. It is interpreted relative to your **current directory**.

Special symbols:
| Symbol | Meaning |
|---|---|
| `.` | current directory |
| `..` | parent (one level up) directory |
| `~` | the current user's home directory |

```bash
$ cd ..                 # go up one level
$ cd ./singleron_lecturer/Windrunner/singleron/software   # relative path from current location
```

---

## 5. File & Directory Naming Conventions

- All characters are valid **except the forward slash `/`**. If special characters are used, wrap the name in quotes.
- Maximum filename length: **256 characters**.
- Names are **case-sensitive** (`Singleron` ≠ `singleron`).
- Avoid starting a filename with `+`, `-`, or `.` as the first character.
- Unlike Windows, Linux does **not** determine file type from the extension/suffix — it reads the file's **header information** instead. (So renaming `photo.jpg` to `photo.txt` doesn't actually change what kind of file it is.)

Common file types (seen via `ls -l`, first character of the permission string):

| Symbol | Type |
|---|---|
| `-` | Regular file |
| `d` | Directory |
| `l` | Symbolic link |
| `b` | Block device file |
| `c` | Character device file |
| `p` | Pipe |
| `s` | Socket |

---

## 6. Common Commands

A Linux command generally has three parts:

```
command + options + arguments
ls -a /mnt/sdd/singleron_lecturer/Windrunner
```

- Options modify *how* the command behaves (usually prefixed with `-`).
- Arguments are *what* the command acts on (e.g., a file or directory name).

### `ls` — list directory contents
```bash
ls          # multi-column listing
ls -l       # long/detailed listing (permissions, owner, size, date)
ls -a       # show hidden files (those starting with .)
ls -d       # show directory properties instead of its contents
```

### `pwd` — print working directory
Shows the full absolute path of where you currently are.

### `cd` — change directory
```bash
cd ~        # go to home directory
cd -        # go to the previous directory
cd ..       # go to parent directory
cd .        # stay in current directory (no-op)
```

### `cp` — copy
```bash
cp file folder             # copy file into folder
cp file folder/file2       # copy and rename
cp -r folder1 folder2      # copy an entire directory (recursive)
cp -d file1 file2          # copy a symbolic link itself, not its target
```

### `mv` — move / rename
```bash
mv folder1 folder2   # if folder2 exists → moves folder1 INSIDE it
                      # if folder2 doesn't exist → renames folder1 to folder2
mv file folder        # same logic for files
```

### `mkdir` — create a directory
```bash
mkdir class_new
mkdir -p project/1/2/3/4   # -p: create all nested directories at once
mkdir -m 777 software      # -m: set permission while creating
```

### Help commands
| Command | Purpose |
|---|---|
| `whatis` | one-line description of a command |
| `which` | shows the location of the executable |
| `man` | full manual with all options |
| `help` | help for shell built-in commands |
| `info` | similar to `man`, sometimes more detailed |

---

## 7. Viewing File Contents

| Command | Behavior |
|---|---|
| `cat` | prints the whole file starting from line 1 |
| `tac` | prints the file **backwards** (`cat` reversed — literally the name spelled backwards) |
| `nl` | prints with line numbers |
| `more` | shows content page by page (forward only) |
| `less` | like `more`, but you can scroll both forward and backward |
| `head` | shows only the first lines |
| `tail` | shows only the last lines |

Useful flags:
```bash
cat -n file              # display with line numbers
cat file1 file2 file3 > file4     # merge (concatenate) files into file4
cat file1 >> file2        # append file1's content to the end of file2

head -n 5 test_diffgenes.xls   # first 5 lines
tail -n 5 test_diffgenes.xls   # last 5 lines
head -c 100 file               # first 100 bytes
```

Example from a real gene-expression table:
```bash
$ head -n 5 test_diffgenes.xls
gene_id    p_val            avg_logFC         pct.1  pct.2  p_val_adj
FTH1       2.70952885711447e-11  -0.592499992784135   1      1      5.76560645505388e-07
EIF1       1.57964538445031e-10  -0.527195749736074   0.989  0.993  3.36132741357181e-06
```

---

## 8. Creating Files

### `touch`
Two uses:
1. Updates the **timestamp** of an existing file to the current system time (content stays untouched).
2. Creates a new **empty** file if it doesn't already exist.

```bash
touch test1.txt test2.txt        # create two empty files
touch -c test3.txt                # -c: only update timestamp, don't create if missing
touch -t 201808080808 file.txt    # -t: set a specific timestamp (YYYYMMDDhhmm)
```

Check a file's metadata:
```bash
stat test1.txt
```

---

## 9. Deleting Files & Directories

### `rm`
```bash
rm -f file          # force delete, never prompt, ignore missing files
rm -r folder        # recursive delete — needed for non-empty directories. DANGEROUS.
rm -i file          # prompt before every deletion (safer)
rmdir emptyFolder   # only deletes an EMPTY directory
```

> ⚠️ **There is no recycle bin in Linux.** A deleted file is gone — permanently. A common safety habit is adding `alias rm='rm -i'` to your `.bashrc` so every `rm` asks for confirmation.

---

## 10. Text Editors — Vim & Nano

### Vim (`vi`/`vim`)
The most common Linux/Unix text editor. It has **three modes**:

1. **Command mode** — default mode after opening a file; used for navigation and shortcuts.
2. **Insert/Edit mode** — press `a` or `i` from command mode to start typing; press `Esc` to leave.
3. **Last-line mode** — press `:` from command mode to type save/search/replace/exit commands.

```
vim filename     → opens file in Command Mode
   a / i         → enters Insert Mode (Esc to exit)
   :             → enters Last-line Mode
```

Useful shortcuts:
| Key/Command | Action |
|---|---|
| `x` | delete the character under cursor |
| `u` | undo last change; `U` undoes the whole line |
| `Ctrl+R` | redo |
| `A` | append text at the end of the line |
| `dw` | delete from cursor to end of word |
| `:wq` | save and quit |
| `:q!` | discard changes and quit |

Cheat sheet: https://vim.rtorr.com/

### Nano
A simpler, **modeless** editor — you can start typing immediately after opening a file (no mode-switching needed).

| Shortcut | Action |
|---|---|
| `Ctrl+A` / `Ctrl+E` | jump to beginning / end of line |
| `Ctrl+W` | search for text (`Alt+W` to repeat search) |
| `Ctrl+K` / `Ctrl+U` | cut / paste the selected line |
| `Ctrl+O` | save (write out) |
| `Ctrl+X` | exit |
| `Ctrl+G` | open help |

Cheat sheet: https://www.nano-editor.org/dist/latest/cheatsheet.html

---

## 11. Permissions (chmod)

Every file/directory has three permission groups: **owner (u)**, **group (g)**, and **others (o)**. Each group can have three permission types: **read (r)**, **write (w)**, **execute (x)**.

```
drwxr-xr-x
│└┬┘└┬┘└┬┘
│ u  g  o
└── file type (d = directory, - = regular file)
```

### Numeric (octal) notation
Each permission has a value: `r=4`, `w=2`, `x=1`. Sum them per group.

```
rwx r-x r--
7   5   4
```
So `chmod 755 file` means: owner = read+write+execute, group = read+execute, others = read-only.

```bash
chmod 755 test.txt          # set permissions numerically
chmod -R 755 test_folder    # -R: apply recursively to all contents
chmod u=rwx,g=rx,o=r file   # symbolic notation (same result as 754)
```

Related commands:
- `chown` — change file **owner**
- `chgrp` — change file **group**

---

## 12. Compression & Archiving (tar, zip)

### `tar`
```bash
tar -cf archive.tar file1 file2    # create (-c) an archive
tar -tvf archive.tar                # list (-t) contents verbosely
tar -xf archive.tar                 # extract (-x)
tar -xf source.tar -C directory     # extract into a specific directory

tar -czf target.tar.gz file1 file2  # create AND gzip-compress
tar -xzf source.tar.gz              # extract a .tar.gz

tar -cjf target.tar.bz2 file1       # create AND bzip2-compress
tar -xjf source.tar.bz2 -C dir      # extract a .tar.bz2
```
> Note: `-f` must be immediately followed by the archive filename.

### `zip` / `unzip`
```bash
zip -r compressed.zip directory        # compress a directory
unzip -r compressed.zip -d directory   # decompress into a directory
```

---

## 13. Wildcard Characters

| Wildcard | Meaning |
|---|---|
| `*` | matches any string (including empty), 0 or more characters |
| `?` | matches exactly 1 character |
| `[abcd]` | matches any one character among a, b, c, d |
| `[a-z]` | matches any character in the range a–z |
| `{..}` | matches a comma-separated set of generated patterns |
| `[!abcd]` | matches any character **NOT** in the set |

```bash
ls *.sh*            # all files containing ".sh" anywhere in the name
ls *.sh* > sh.txt    # save that list into sh.txt
```

---

## 14. Pipes & Output Redirection

**Pipe (`|`)**: sends the **standard output (stdout)** of one command as the **standard input (stdin)** of the next command.
> Important: a pipe only passes stdout, **not** error output (stderr).

```bash
command1 | command2
command1 | command2 | command3   # can be chained
```

**Redirection (`>`)**: sends a command's output into a file (overwriting it).
```bash
command > file        # overwrite file with command's output
command >> file        # append to file instead of overwriting
```

---

## 15. Chaining Multiple Commands

| Separator | Behavior |
|---|---|
| `;` | run commands one after another regardless of success/failure of the previous one |
| `&&` | run the next command **only if** the previous one succeeded |
| `\|\|` | run the next command **only if** the previous one **failed** (keep trying until one succeeds) |

```bash
cmd1 ; cmd2          # both run no matter what
cmd1 && cmd2         # cmd2 runs only if cmd1 succeeded
cmd1 || cmd2         # cmd2 runs only if cmd1 failed
```

---

## 16. Task Submission (nohup, qsub)

On a shared server, long-running bioinformatics jobs (alignment, assembly) shouldn't block your terminal session or die when you disconnect.

### Direct submission with `nohup`
```bash
nohup [command to run] &
```
`nohup` = "no hang up" → the process keeps running even after you log out. The trailing `&` runs it in the background.

```bash
conda activate celescope2
cd /path/to/shell
nohup sh rna_test.sh &
```

### Cluster submission with `qsub` (SGE scheduler)
```bash
qsub [options] file.sh
```
| Option | Meaning |
|---|---|
| `-V` | pass current environment variables to the job |
| `-cwd` | use current directory as the working directory |
| `-l vf=10g,p=10` | request resources: 10 GB memory, 10 CPUs |
| `-q` | specify a queue |
| `-N` | name the task |

```bash
qsub -cwd -V -l vf=10g,p=10 test.sh
```

---

## 17. Conda: Introduction

**Conda** is an open-source **package manager and environment manager**. It lets you install multiple versions of software packages (and their dependencies) and switch between them easily.

- Originally built for **Python**, but also manages **R, Ruby, Lua, Scala, Java, JavaScript, C/C++, FORTRAN**, and more.
- Works across **Linux, macOS, and Windows**.
- Distributed via **Anaconda** (full distribution) or **Miniconda** (minimal installer).

---

## 18. Why Conda Exists

In early bioinformatics work, installing software was painful:

- `sudo apt` / `sudo yum` → fails if you're not an administrator.
- Downloading a **precompiled binary** → often unavailable for the tool you need.
- **Recompiling from source** → requires understanding the codebase, resolving dependency conflicts, frequent build failures.

**Conda solved this** by letting any user (no admin rights needed) install software and automatically resolve all dependencies into an isolated environment.

---

## 19. Installing Miniconda

Official docs: https://conda.io/projects/conda/en/stable/user-guide/install/index.html

- Miniconda installers exist for Windows, macOS (Intel & Apple Silicon), and Linux (including aarch64/ARM).
- Follow the official instructions for your platform; verify the installer with its published SHA256 hash if desired.

---

## 20. Conda Channels

Think of **channels** like TV channels — conda itself is the "TV," and channels are where it looks to find software to install.

| Step for conda | TV analogy |
|---|---|
| Download conda | Buy a TV |
| Install conda | Mount the TV |
| Configure channels | Search for channels |
| Install packages | Watch TV |

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set show_channel_urls yes
```

**Tips:**
1. Don't add the same channel repeatedly — redundant channels slow down installs. Order matters (later-added channels get **higher** priority).
2. Channel configuration applies to **where the server is located**, not where you (the user) physically are.
3. Your configured channels live in `~/.condarc`.

```bash
conda config --get channels     # view in order added
conda config --show channels    # view current priority order
```

---

## 21. Conda Environments

**Why isolate environments?**
- Conda changes your Python/package setup globally by default — risky.
- **Dependency conflicts**: e.g., Software A needs Python 2.7, Software B needs Python 3.8 — both can't coexist in one environment.
- Helps you manage multiple project types separately: DNA-seq, RNA-seq, ChIP-seq, scRNA-seq, etc., each with its own tool versions.

**Advice:** Never install packages directly into the `base` environment unless you fully understand the consequences — it can break your entire conda installation.

```bash
conda create -n R4            # create a new environment named "R4"
conda activate R4              # enter/activate it
conda deactivate               # exit back to base
conda env list                 # list all environments
conda info --env               # same, with more detail
```

### Deleting / renaming an environment
```bash
conda remove -n R4 --all               # delete environment + all its packages

# Conda has no direct "rename" command — clone, then delete the old one:
conda create -n py2                     # original environment
conda create -n py2_new --clone py2     # clone it under a new name
conda remove -n py2 --all               # remove the old one
```

---

## 22. Common Conda Commands

```bash
conda --help                 # general help
conda config                 # view/edit channel configuration
conda search <package>       # search for a package
conda install <package>      # install a package into the active environment
conda uninstall <package>    # remove a package
conda remove <package>       # same as uninstall
conda update <package>       # update a package to the latest compatible version
conda list                   # list packages in the active environment
```

---

## 23. Common Bioinformatics Software (typically installed via conda)

| Category | Tools |
|---|---|
| **Quality Control (QC)** | fastqc, multiqc, fastp, trimmomatic, cutadapt, trim_galore |
| **Alignment & Quantification** | bwa, hisat2, bowtie, bowtie2, STAR, salmon, subread (featureCounts) |
| **Genome Assembly** | Trinity |
| **Single-cell Transcriptomics** | CeleScope |

---

## 24. Hands-on Exercise (with solutions)

**Task 1 — Download and extract practice data:**
```bash
wget https://singleronbio-opendata.oss-cn-hangzhou.aliyuncs.com/Online_training_Materials/linux/linux_traindata.tar.gz
# if it fails due to certificate issues:
wget --no-check-certificate https://singleronbio-opendata.oss-cn-hangzhou.aliyuncs.com/Online_training_Materials/linux/linux_traindata.tar.gz

tar -xzf linux_traindata.tar.gz
```

**Task 2 — Create a directory named `linux_test`:**
```bash
mkdir linux_test
```

**Task 3 — Copy `catfile.txt` and `dog.txt` into `linux_test` in ONE line:**
```bash
cp catfile.txt dog.txt linux_test/
```

**Task 4 — View the first six lines of `dog.txt`:**
```bash
head -n 6 dog.txt
```

**Task 5 — Merge `catfile.txt` and `dog.txt` into `newfile.txt`:**
```bash
cat catfile.txt dog.txt > newfile.txt
```

**Task 6 — List files containing "file" in the name (inside `linux_test`) and save to `test.txt`:**
```bash
ls linux_test/*file* > test.txt
```

**Task 7 — Create a folder named `test2`:**
```bash
mkdir test2
```

**Task 8 — Copy the `gene` folder into `linux_test`, rename it `GENE`, compress to `.tar.gz`:**
```bash
cp -r gene linux_test/GENE
tar -czf GENE.tar.gz linux_test/GENE
```

**Task 9 — Move the compressed GENE archive into `test2` and extract it:**
```bash
mv GENE.tar.gz test2/
tar -xzf test2/GENE.tar.gz -C test2/
```

**Task 10 — Give `linux_test` read, write, and execute permission for all users:**
```bash
chmod -R 777 linux_test
```

---

## Summary

This module covers the **two foundational skills every bioinformatics pipeline depends on**:

1. **Linux command-line fluency** — navigating, managing, and manipulating files/directories, controlling permissions, and running long jobs reliably on a server.
2. **Conda environment management** — installing scientific software without admin rights and keeping project dependencies isolated and reproducible.

Together, these are prerequisites before moving into actual sequence analysis tools (FastQC, BWA, STAR, Trinity, CeleScope, etc.).
