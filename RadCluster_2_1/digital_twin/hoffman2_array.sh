#!/bin/bash
# hoffman2_array.sh -- run the digital-twin campaign as a UGE array job.
#
#   qsub hoffman2_array.sh
#
# ONE ARRAY TASK == ONE MACHINE INDEX.  Each task writes its own
# results/<design>_machine<K>.jsonl, so tasks never touch the same bytes and can
# never merge-conflict -- the same property that lets four workstations share a
# campaign, applied to N cluster tasks.
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
#  VALUES YOU MUST SET -- all marked <<<SET>>>.  Do not guess the physics ones:
#  take them from campaign_layout.py's generated command and from gate S11(i)-3.
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
#  Array range = number of tasks. MUST equal the K you passed to
#  campaign_layout.py --split hoffman2:K
#$ -t 1-16

set -euo pipefail

# ── settings from campaign_layout.py ────────────────────────────────────────
REPO=<<<SET: /u/home/.../RadCluster>>>
DESIGN=design/T3_rev6.csv
OF=<<<SET: total machine indices, from campaign_layout.py>>>
WEIGHTS=<<<SET: the --weights string, byte-identical everywhere>>>
BASE=<<<SET: this group's first machine index, from campaign_layout.py>>>

# ── grid, from gate S11(i)-3 -- NOT the S11(f) placeholder ──────────────────
GRID="--equations bin_moment --I <<<SET>>> --V 10000 \
      --i-discrete 50 --v-discrete 5 --i-bin 25 --v-bin 25 \
      --i-mobile-default 50 --v-mobile-default 5 \
      --dose 1.0 --rtol 1e-6 --solver-mode full_system"

# Park 1.5 rows before h_rt so the scheduler never SIGKILLs a row in flight.
# h_rt 24 h = 86400 s; at ~1 h/row that is 86400 - 5400.
STOP_AFTER_S=81000
TIMEOUT_S=<<<SET: p99 row cost measured ON HOFFMAN2, not on the Mac>>>

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

MACHINE=$(( BASE + SGE_TASK_ID - 1 ))
echo "=== task $SGE_TASK_ID -> machine index $MACHINE of $OF on $(hostname) ==="

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
    --machine "$MACHINE" --of "$OF" --weights "$WEIGHTS" \
    --workers 4 \
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
