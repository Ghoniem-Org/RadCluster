#!/bin/bash
# hoffman2_array.sh -- run the digital-twin campaign as a UGE array job.
#
#   qsub hoffman2_array.sh
#
# ONE ARRAY TASK == ONE SUBTASK OF MACHINE INDEX 3.  Hoffman2 is a single
# participant in a four-machine campaign (0 MacBook Pro, 1 Matrix-PC,
# 2 Nasr Workstation, 3 Hoffman2); its 381-row share is split a second time
# across the 16 array tasks by --subtask, and each writes its own
# results/<design>_machine3_t<K>.jsonl.  Tasks never touch the same bytes, so
# they cannot merge-conflict and a task killed at h_rt costs only its own rows.
#
# WHY AN ARRAY OF SMALL TASKS RATHER THAN ONE BIG JOB
#   A single 64-core job queues behind everything on a shared cluster and, when
#   it hits h_rt, dies whole.  Sixty-four 1-core tasks backfill into gaps, start
#   almost immediately, and die one at a time.  The work is embarrassingly
#   parallel with OMP_NUM_THREADS=1 (plan S11(m): OpenMP does not engage under
#   bin_moment anyway), so there is nothing to gain from a shared-memory job.
#
# RESUMPTION IS THE DESIGN, NOT A FALLBACK
#   Rows append as they complete and a restart skips row_ids already present, so
#   resubmitting this identical script is always safe and always makes progress.
#   Expect to submit it several times; see the resubmission loop at the bottom.
#
# ---------------------------------------------------------------------------
#  NOTHING TO SET.  Every value below is either detected (the machine index, via
#  machines.json) or frozen by plan S12(q)/(r)/(t).  This block used to carry
#  five <<<SET>>> placeholders; each was an opportunity to give two machines the
#  same index, which is the one misconfiguration the merge cannot distinguish
#  from a machine that never reported.
# ---------------------------------------------------------------------------

#$ -cwd
#$ -j y
#$ -o joblogs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -N radcluster_t2
#  h_rt: 24:00:00 is the CEILING on the shared pool -- MEASURED 2026-08-03, not
#  assumed: qsub -w v accepts 24h but refuses 336h with "no permission for
#  cluster queue *_pod.q" on every owned queue. ghoniem is in no group queue,
#  so 24 h is the hard limit and --stop-after-s below must sit under it.
#$ -l h_rt=24:00:00,h_data=4G
#  'shared' PE confirmed present (228 PEs on this cluster). If you raise this,
#  raise --workers to match AND re-run campaign_layout.py -- the weights depend
#  on it.
#$ -pe shared 4
#  Array range = number of tasks. MUST equal --of-subtasks below (16).
#$ -t 1-16

set -euo pipefail

# ── settings from campaign_layout.py ────────────────────────────────────────
REPO=${REPO:-$HOME/RadCluster}
DESIGN=design/T3_rev6.csv
# NOTHING TO SET.  Hoffman2 is machine index 3 of 4 and run_ensemble detects
# that from machines.json (env SGE_ROOT + node login<N>/n<N>), so the index,
# --of and --weights are no longer typed here.  That removes the failure this
# block used to invite: a mistyped index means two machines compute the same
# rows and some rows are computed by nobody.
#
# The 16 array tasks share machine index 3 and split its 381 rows a second time
# by --subtask, 23-24 rows each, each writing its own results file.

# ── grid: the rev-6 BASELINE, fixed by the author 2026-08-05 (plan S11(q)) ──
# Not "from gate S11(i)-3" any more -- truncation was withdrawn as an
# admissibility criterion, so there is no convergence gate left to set I and V.
# They are declared constants, sized off the measured EUROFER97 band:
#   I=30000 -> d_ceiling(<100>)  = 39.7 nm = 2.27x the measured <100> mean
#   V= 5000 -> d_ceiling(cavity) =  4.8 nm = 2.16x the measured cavity mean
# active_window (not full_system): every timing on record -- the Mac's 5580.6 s
# reference and the 26059.6 s Hoffman2 comparison -- was measured in it, so
# switching modes would invalidate both STOP_AFTER_S and TIMEOUT_S below.
GRID="--equations bin_moment --I 30000 --V 5000 \
      --i-discrete 50 --v-discrete 5 --i-bin 25 --v-bin 25 \
      --i-mobile-default 50 --v-mobile-default 5 \
      --dose 1.0 --rtol 1e-6 --solver-mode active_window"

# ── walltime arithmetic, from the MEASURED Hoffman2 speed ───────────────────
# 2026-08-04: job 14223957 ran the reference row here in 26059.6 s against the
# Mac's 5580.6 s -> Hoffman2 is 4.67x SLOWER PER CORE (--speed hoffman2:0.214),
# uncontended (99.4 % cpu/wall) on a Xeon Gold 6240.  Observables matched the
# Mac to printed precision, so this is hardware, not a divergent trajectory.
#
# At the rev-6 baseline grid that is ~3130 s x 4.67 = ~14600 s/row HERE.
# The previous STOP_AFTER_S=81000 assumed ~1 h/row and reserved 5400 s; that
# would have SIGKILLed a row in flight on nearly every task.  Reserve one whole
# row instead: a row starting just under 68000 finishes by ~82600 < 86400.
# STOP_AFTER_S + TIMEOUT_S must stay under h_rt (86400), or the scheduler
# SIGKILLs a row that the timeout would otherwise have retired gracefully --
# and a SIGKILL loses the partial trajectory that plan S12(s) exists to keep.
#   56000 + 30000 = 86000 < 86400.
# At ~17650 s/row here (3780 Mac-s x 4.67) each of the 4 workers starts a row at
# roughly 0 / 17650 / 35300 / 52950 and the last finishes by ~70600.
STOP_AFTER_S=56000
# 30000 s = 6424 Mac-equivalent seconds, above the 5276 s slowest of the three
# calibrated rows, so most rows finish fully.  A row that does not is KEPT: the
# solver is asked to finalize, flushes its trajectory, and the row contributes
# to every rung of the dose ladder it reached (plan S12(s)).
TIMEOUT_S=30000

# ── environment ─────────────────────────────────────────────────────────────
# ALL VERSIONS PINNED, and all three corrected from the 2026-08-03 probe:
#   'module load gcc'    would give 7.5.0 (default), and the SYSTEM gcc is
#                        4.8.5 -- neither compiles the C++17 in cpp_utils
#   'module load cmake'  would give 3.19.5 (default); CMakeLists wants >= 3.10
#                        but the system cmake is 2.8.12 and fails outright
#   'module load python/3.11'  DOES NOT EXIST. Stock python modules top out at
#                        3.9.6; anaconda3/2023.03 gives Python 3.10.9 with
#                        numpy 1.23.5 / scipy 1.10.0 / pandas 1.5.3.
# The other machines record Python 3.11.7, so this is a KNOWN mismatch.
# merge_and_sobol does not gate on it (python is not in the split check) and the
# numerics live in the C++ solver -- the agreement gate below is what actually
# decides comparability. Recorded per row either way.
. /u/local/Modules/default/init/bash
module purge
module load gcc/11.3.0 cmake/3.30.0 anaconda3/2023.03
# SUNDIALS is vendored by hoffman2_setup.sh: the probe found no sundials module
# and nothing under /u/local/apps, so find_package(SUNDIALS REQUIRED) needs this.
export CMAKE_PREFIX_PATH="$HOME/opt/sundials-7.1.1:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$HOME/opt/sundials-7.1.1/lib64:$HOME/opt/sundials-7.1.1/lib:${LD_LIBRARY_PATH:-}"

cd "$REPO/RadCluster_2_1"

# Serial workers, one row each.  Not a tuning knob: it also removes OpenMP
# reduction-order nondeterminism, which is what makes results from different
# builds poolable at all (plan S11(l)).
export OMP_NUM_THREADS=1

echo "=== task $SGE_TASK_ID -> subtask $(( SGE_TASK_ID - 1 ))/16 of machine 3 on $(hostname) ==="

# ── solver: each machine builds its own; git does not carry build/ ──────────
# Serialise the build across tasks landing on the same node, or 64 tasks race
# on the same build directory and corrupt it.
BUILD="$REPO/RadCluster_2_1/build"
if [ ! -x "$BUILD/solver" ]; then
  echo "*** solver not built. Run hoffman2_setup.sh ONCE on a login node first --"
  echo "*** it vendors SUNDIALS 7.1.1 (absent on this cluster) and builds the"
  echo "*** solver. Building it from inside an array task would have 16 tasks"
  echo "*** racing on one build directory and would waste the queue wait."
  exit 1
fi

# ── agreement gate: do NOT contribute rows from an unverified build ─────────
# Hoffman2's compiler and SUNDIALS differ from every other participant, so this
# is a THIRD binary. plan S11(l) measured Mac vs MATRIX-PC2 at 2.95e-08 and
# judged them poolable, but that was measured, not assumed -- and it must be
# measured here too. Only task 1 runs it; the rest trust the exit status file.
cd "$REPO/RadCluster_2_1/digital_twin"
mkdir -p joblogs
if [ "$SGE_TASK_ID" -eq 1 ]; then
  python check_machine.py > joblogs/check_machine.hoffman2.log 2>&1 && echo ok > .agree.ok || echo FAIL > .agree.ok
fi
for _ in $(seq 1 60); do [ -f .agree.ok ] && break; sleep 10; done
if [ "$(cat .agree.ok 2>/dev/null)" != "ok" ]; then
  echo "*** check_machine FAILED or timed out -- refusing to contribute rows."
  echo "*** See joblogs/check_machine.hoffman2.log. A build that does not agree"
  echo "*** cannot be pooled with the other machines (plan S11(l))."
  exit 1
fi

# ── run ─────────────────────────────────────────────────────────────────────
python run_ensemble.py \
    --design "$DESIGN" \
    --machine auto \
    --subtask $(( SGE_TASK_ID - 1 )) --of-subtasks 16 \
    $GRID \
    --timeout-s "$TIMEOUT_S" \
    --stop-after-s "$STOP_AFTER_S"

echo "=== task $SGE_TASK_ID finished, rc=$? ==="

# ---------------------------------------------------------------------------
# AFTERWARDS, FROM THE LOGIN NODE (compute nodes usually have no outbound net):
#
#   cd $REPO && git pull --rebase
#   git add RadCluster_2_1/digital_twin/results/*_machine*.jsonl \
#           RadCluster_2_1/digital_twin/results/*_machine*.manifest.json
#   git commit -m "T2 results, hoffman2 tasks 1-64"
#   git push
#
# Push the MANIFESTS too, not just the rows: campaign_layout.py sizes the next
# round from session_wall_s + completed, and MATRIX-PC2's manifest never
# arrived, which is exactly why its weight had to be declared instead of
# measured (plan S11(n)).
#
# TO RESUME after tasks hit h_rt -- safe to run repeatedly, it only ever adds:
#   until python merge_and_sobol.py --design $DESIGN --results results/ \
#            | grep -q "0 MISSING"; do qsub hoffman2_array.sh; sleep 3600; done
# ---------------------------------------------------------------------------
