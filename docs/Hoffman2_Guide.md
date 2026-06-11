# Hoffman2 Cluster — Working Guide (RadCluster)

UCLA IDRE Hoffman2 cluster. Scheduler: **UGE / Grid Engine 8.6.4** (`qsub`, `qrsh`, `qstat`, …).
This guide is grounded in the actual setup used for the RadCluster project (your SSH alias, the
`anaconda3/2023.03` module, the JupyterLab + SSH‑tunnel workflow).

> Quick mental model: **login nodes are for editing/submitting only — never run compute there.**
> All real work runs on **compute nodes** that you request from the scheduler, either *interactively*
> (`qrsh`) or as a *batch job* (`qsub`).

---

## 1. Logging in

### 1.1 Your SSH setup (already configured)

You have a key‑based alias in `~/.ssh/config` (no password needed):

```
Host hoffman2 hoffman2.idre.ucla.edu
    HostName hoffman2.idre.ucla.edu
    User ghoniem
    IdentityFile ~/.ssh/hoffman2_ed25519
    IdentitiesOnly yes
    ServerAliveInterval 60
```

So from a terminal you simply:

```bash
ssh hoffman2
```

You land on one of the **login nodes** (e.g. `login2`). To confirm where you are: `hostname`.

### 1.2 From a machine without the alias

```bash
ssh ghoniem@hoffman2.idre.ucla.edu
```

This uses your UCLA password (and Duo if enabled). Prefer the key‑based alias — it avoids
re‑typing credentials and keeps automation working.

### 1.3 Generating the key (if you ever set up a new machine)

```bash
ssh-keygen -t ed25519 -f ~/.ssh/hoffman2_ed25519        # create a key pair
ssh-copy-id -i ~/.ssh/hoffman2_ed25519.pub ghoniem@hoffman2.idre.ucla.edu
```

Then add the `Host hoffman2` block above to `~/.ssh/config`.

### 1.4 Moving files between your laptop and Hoffman2

> **Run `scp`/`rsync` from your laptop — NOT from inside an SSH session on Hoffman2.** The `hoffman2`
> alias only exists in your laptop's `~/.ssh/config`; the cluster's login nodes can't resolve it
> (`ssh: Could not resolve hostname hoffman2`). Open a fresh local terminal (don't `ssh` in first)
> and run the command there.

```bash
# copy a file up
scp ./input_parameters.xlsx hoffman2:~/RadCluster/RadCluster_1_0/input/

# pull a results directory DOWN (recursive)
scp -r hoffman2:~/RadCluster/RadCluster_1_0/output/<run_dir> ./output/

# rsync is better for large / repeated syncs (only sends changes, resumable)
rsync -avz --progress hoffman2:~/RadCluster/RadCluster_1_0/output/ ./output/
```

For very large transfers use **Globus** (endpoint: *UCLA Hoffman2*).

---

## 2. JupyterLab in your browser (and submitting jobs from it)

JupyterLab must run **on a compute node** (inside a scheduler job), not on a login node. You then
tunnel that node's port to your laptop over SSH. There are two ways to do this.

### 2.1 Recommended: the official `h2jupynb` launcher (one command)

IDRE ships a Python helper that does *everything* — requests a node, starts Jupyter, opens the
tunnel, and launches your browser:

```bash
# download once
curl -O https://raw.githubusercontent.com/rdauria/jupyter-notebook/main/h2jupynb
# or grab the current copy from the IDRE Hoffman2 docs

# launch:  -u user  -t walltime(hr)  -m mem(GB/core)  -c cores  -v python version
python3 h2jupynb -u ghoniem -t 8 -m 8 -c 2 -v 3.10 -l lab
```

Run `python3 h2jupynb -h` to see all flags (GPU, architecture, directory, port, etc.). This is the
cleanest day‑to‑day method. The manual method below is what to fall back on when you want full
control or are reconnecting to an existing job.

### 2.2 Manual method (what we set up this session)

**Step 1 — get an interactive compute node** (here: 2 cores, 8 GB/core, 8‑hour wall clock):

```bash
ssh hoffman2
qrsh -N JUPYNB -l h_rt=8:00:00,h_data=8G -pe shared 2
# you are now ON a compute node, e.g. n1037
```

**Step 2 — start JupyterLab on that node:**

```bash
module load anaconda3/2023.03          # Python 3.10.9 — matches the RadCluster runs
cd $HOME/RadCluster
jupyter lab --no-browser --ip=0.0.0.0 --port=8889 --ServerApp.token=radcluster \
            --ServerApp.root_dir=$HOME/RadCluster
# note the node name (e.g. n1037) and the port/token
```

A reusable launcher script (`~/start_jupyter.sh`) lives in your home directory for exactly this:

```bash
#!/bin/bash
source /u/local/Modules/default/init/bash 2>/dev/null
module load anaconda3/2023.03
cd $HOME/RadCluster
nohup jupyter lab --no-browser --ip=0.0.0.0 --port=8889 \
  --ServerApp.token=radcluster --ServerApp.root_dir=$HOME/RadCluster \
  > $HOME/jupyter_n1037.log 2>&1 &
```

**Step 3 — open the SSH tunnel from your laptop** (replace `n1037` with your actual node):

```bash
ssh -f -N -o ExitOnForwardFailure=yes -L 8889:n1037:8889 hoffman2
```

This forwards `localhost:8889` → login node → `n1037:8889`.

**Step 4 — open in your browser:**

```
http://localhost:8889/lab?token=radcluster
```

To tear the tunnel down later: `ssh -O exit -L 8889:n1037:8889 hoffman2` (or just close the terminal
if you ran it in the foreground).

### 2.3 Submitting NEW jobs from inside JupyterLab

Your Jupyter session is itself just one interactive job. To launch **additional** (e.g. long batch)
jobs, open a **Terminal** tab in JupyterLab (or use `!` in a notebook cell) and use `qsub`:

```bash
qsub my_run.sh                 # submit a batch job script
qstat -u ghoniem               # watch your jobs
qdel <job-id>                  # cancel one
```

A minimal batch submit script `my_run.sh`:

```bash
#!/bin/bash
#$ -N radcluster_run            # job name
#$ -l h_rt=24:00:00,h_data=8G   # 24 h wall clock, 8 GB per core
#$ -pe shared 4                 # 4 cores on one node
#$ -cwd                         # run from current directory
#$ -o run.$JOB_ID.log           # stdout
#$ -j y                         # merge stderr into stdout
#$ -M ghoniem@ucla.edu          # email
#$ -m bea                       # mail at begin/end/abort

source /u/local/Modules/default/init/bash
module load anaconda3/2023.03
cd $HOME/RadCluster/RadCluster_1_0
export OMP_NUM_THREADS=$NSLOTS
python codes/Python_Testing/run_simulation.py
```

> **Why batch over interactive for production runs:** the partial run we looked at stopped at the
> solver's internal `timeout_s`, but interactive `qrsh` sessions *also* die the moment your laptop
> sleeps or the SSH session drops. A `qsub` batch job keeps running on the cluster regardless of your
> connection — that's the right tool for long simulations.

---

## 3. Navigating & using the cluster

### 3.1 Node types

| Node class | What it is | How you get it |
|---|---|---|
| **Login nodes** | Editing, compiling, submitting, light file ops. **No compute.** | Where `ssh hoffman2` lands you |
| **Shared compute** (`pod_smp.q`, etc.) | General community pool. Time‑limited. | `qrsh`/`qsub` with `-l h_rt=…` (≤ 24 h) |
| **Group / `highp` nodes** | Nodes your PI's group owns. Longer walltime, priority. | add `-l highp` (only if your group owns nodes) |
| **High‑memory** | Large‑RAM nodes for big problems | request large `h_data` and/or `-l highmem` |
| **GPU nodes** | CUDA GPUs for accelerated codes | `-l gpu,…` (see §3.4) |

Your current Jupyter session runs on a shared‑pool node (`pod_smp.q@n1037`).

### 3.2 Requesting resources (the key `qsub`/`qrsh` flags)

| Flag | Meaning | Example |
|---|---|---|
| `-l h_rt=HH:MM:SS` | **Wall‑clock limit**. Shared pool max ≈ **24 h**; group `highp` up to **14 days**. | `-l h_rt=24:00:00` |
| `-l h_data=NG` | **Memory per core/slot** (not total). Total = `h_data × slots`. | `-l h_data=8G` |
| `-pe shared N` | **N cores on a single node** (shared‑memory / OpenMP). | `-pe shared 4` |
| `-pe dc* N` | N cores spread across nodes (distributed / MPI). | `-pe dc* 32` |
| `-l highp` | Use your group's owned nodes (longer walltime). | `-l highp` |
| `-l arch=intel-gold*` | Pin a CPU architecture. | `-l arch=intel-gold*` |
| `-l exclusive` | Whole node to yourself (no sharing). | `-l exclusive` |

**Memory rule of thumb:** total memory = `h_data` × number of slots. Need 64 GB on 4 cores →
`-pe shared 4 -l h_data=16G`.

**Your actual account limits** (verified 2026‑06‑04 via `mygroup`):

- **Max cores in a single non‑`highp` job:** `slots = 600`
- **Max simultaneous jobs:** `500`
- **Max tasks per array job:** `75000`
- Your group does **not** own dedicated `highp` nodes — you run in the **shared community pool**, so
  plan around the **~24 h** shared‑pool walltime cap (checkpoint long runs, or split into resumable
  chunks). Re‑check with `mygroup` if your group later buys nodes.

### 3.3 Interactive vs. batch

```bash
# Interactive shell on a compute node (good for debugging, Jupyter):
qrsh -l h_rt=2:00:00,h_data=4G -pe shared 2

# Batch submission (good for long production runs — survives disconnects):
qsub my_run.sh
```

### 3.4 CPU vs. GPU

CPU is the default. To request a **GPU**:

```bash
# interactive, 1 GPU + 1 CPU core, 4 h:
qrsh -l h_rt=4:00:00,h_data=8G,gpu,cuda=1 -pe shared 1
```

In a batch script the GPU request goes on a `#$ -l` line, e.g.:

```bash
#$ -l h_rt=8:00:00,h_data=16G,gpu,cuda=1
```

`cuda=N` is the number of GPUs per node. To target a **specific architecture**, add a compute‑
capability constraint, e.g. `-l gpu,cuda=1,cuda_cc=8.0` (Ampere or newer). Then load the matching
toolkit (`module load cuda`) and verify inside the job with `nvidia-smi`.

**GPU models actually present on Hoffman2** (verified 2026‑06‑04 via `qhost -F`):

| GPU model | GPU memory | Compute cap. | Notes |
|---|---|---|---|
| **NVIDIA H200 NVL / H100 NVL** | 141 GB / 94 GB | 9.0 | Newest, scarce — top‑tier training/HPC |
| **NVIDIA A100** (SXM4‑40 GB & PCIe‑80 GB) | 40 / 80 GB | 8.0 | High‑end, in demand |
| **NVIDIA L40S** | 45 GB | 8.9 | **Most plentiful** (~13 nodes) — best availability |
| **NVIDIA RTX A6000** | 48 GB | 8.6 | Large‑memory workstation GPU |
| **Tesla V100‑PCIE** | 32 GB | 7.0 | Older HPC card |
| **GeForce RTX 2080 Ti** | 11 GB | 7.5 | Consumer, small memory |
| **Tesla P4** | 7.5 GB | 6.1 | Inference card, low memory |
| **GeForce GTX 1080 Ti** | 11 GB | 6.1 | Old consumer |
| **Tesla K40m** | 11 GB | 3.5 | Legacy — too old for recent CUDA/PyTorch |

Practical picks: for **best availability** request an **L40S** (`cuda_cc=8.9`); for **max GPU memory
or newest features** target **A100/H100/H200** (`cuda_cc>=8.0`, expect longer queue waits); avoid the
K40m (compute cap 3.5 is below most modern frameworks' minimum).

> RadCluster's current solvers are **CPU‑based (SUNDIALS/CVODE)** — you only need GPU nodes if/when a
> code path is GPU‑accelerated (e.g. a future ML‑surrogate or CUDA linear solver).

### 3.5 Software modules

Hoffman2 uses environment **modules** — nothing scientific is on your `PATH` until you load it.

```bash
module avail                 # list everything available
module avail anaconda        # filter
module load anaconda3/2023.03   # what RadCluster uses (Python 3.10.9)
module load cmake gcc cuda      # e.g. for the C++ SUNDIALS solver / GPU
module list                  # what's loaded now
module purge                 # unload everything
```

Put your usual `module load` lines at the top of every batch script (login‑shell state does **not**
carry into a job).

### 3.6 Storage

| Location | Quota | Currently used | Notes |
|---|---|---|---|
| `$HOME` (`/u/home/g/ghoniem`) | **60 GB** (5,000,000 files) | 8.07 GB / 22,179 files | Code, configs, inputs — back up; don't fill it |
| `$SCRATCH` (`/u/scratch/g/ghoniem`) | **2,000 GB** | ~0 GB | Large temporary I/O — **auto‑purged** (files deleted after a set age). Not for keeping results |
| Group/project space | — | — | Only if your group buys an allocation |

*(Quotas verified 2026‑06‑04 via `myquota`.)* Check anytime:

```bash
myquota          # your home/scratch usage vs. limits
```

> For RadCluster: keep the repo in `$HOME`, but if a run produces huge `results_y.npy` / `.pkl`
> files, write them to `$SCRATCH` during the run and `rsync` the keepers back to `$HOME` or your
> laptop afterward (scratch is purged).

### 3.7 Everyday command cheat‑sheet

```bash
# jobs
qsub script.sh                 # submit batch job
qrsh -l h_rt=2:00:00 -pe shared 2   # interactive node
qstat -u ghoniem               # your jobs (r=running, qw=queued/waiting, Eqw=error)
qstat -j <job-id>              # full detail: resources, usage, why it's pending
qdel <job-id>                  # cancel a job
myjobs                         # friendly summary of your jobs

# nodes / cluster
qhost                          # all nodes: cores, load, memory
qhost -F gpu                   # nodes with GPU resources
mygroup                        # your group + any owned (highp) nodes
myquota                        # disk usage

# inside a job
echo $NSLOTS                   # cores granted (use for OMP_NUM_THREADS / -j)
echo $JOB_ID  $SGE_O_WORKDIR   # job id / submit directory
nvidia-smi                     # GPU status (GPU jobs only)
```

### 3.8 Reading job state

- `r` — running
- `qw` — queued, waiting for resources (be patient, or request fewer cores / less time)
- `Eqw` — error state; inspect with `qstat -j <id>` and clear with `qmod -cj <id>` or `qdel`
- `dr` / `dt` — being deleted

---

## 4. RadCluster‑specific notes

- **Repo on Hoffman2:** `~/RadCluster` (mirrors this local repo). Active module: `RadCluster_1_0/`.
- **Python env:** `module load anaconda3/2023.03` (Python 3.10.9 — matches run provenance).
- **Run output:** `RadCluster_1_0/output/<timestamp>_<config>/` — `provenance.md`, `summary.csv`,
  `diagnostics.txt`, `results_t.npy`, `results_y.npy`, `plots/`. A `_PARTIAL` suffix + `run_status:
  interrupted` means the solver hit its internal `timeout_s` before finishing.
- **Long production runs:** prefer `qsub` batch jobs over interactive Jupyter so they survive
  disconnects; and raise the solver `timeout_s` if you want to reach high dose.
- **C++ SUNDIALS solver:** needs `module load cmake gcc` to build under `cpp_utils/` → `build/`.

---

*Last updated: 2026‑06‑04. Some site limits (quotas, max walltimes, available GPU models) change over
time — confirm against the current IDRE Hoffman2 documentation and `qhost` / `myquota` output.*
