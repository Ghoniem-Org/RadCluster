#!/bin/bash
# hoffman2_setup.sh -- ONE-TIME setup of the RadCluster toolchain on Hoffman2.
#
# Run once on a LOGIN NODE (it is a build, not compute -- allowed, and it needs
# the outbound network that compute nodes do not have).  Takes ~15-25 min,
# almost all of it SUNDIALS.
#
#   ssh ghoniem@hoffman2.idre.ucla.edu
#   bash hoffman2_setup.sh
#
# WHY SUNDIALS IS BUILT FROM SOURCE
#   The 2026-08-03 probe found NO sundials module and nothing under
#   /u/local/apps/sundials.  cpp_utils/CMakeLists.txt does
#   find_package(SUNDIALS REQUIRED) for 7.1.1 with CVODE + nvecserial, so it
#   must be vendored.  This is the single biggest difference between bringing
#   Hoffman2 online and bringing up another Mac.
#
# ENVIRONMENT, AS MEASURED (not assumed) -- probe of 2026-08-03, login4:
#   UGE 8.6.4 (qsub/qstat/qhost; no Slurm)   CentOS 7, Linux 3.10
#   default gcc 4.8.5 / cmake 2.8.12 -- BOTH TOO OLD, C++17 is required
#   gcc/11.3.0 and cmake/3.30.0 available as modules
#   NO python/3.11 module (max stock 3.9.6); anaconda3/2023.03 = Python 3.10.9
#   $SCRATCH=/u/scratch/g/ghoniem, 2 TB quota, empty; $HOME 60 GB, 4.85 used
#   login node reaches github (HTTP 200)
set -euo pipefail

SUNDIALS_VER=7.1.1
PREFIX=${PREFIX:-$HOME/opt/sundials-$SUNDIALS_VER}
REPO=${REPO:-$HOME/RadCluster}
BUILD_TMP=${BUILD_TMP:-$SCRATCH/sundials_build}

. /u/local/Modules/default/init/bash
# Versions are PINNED.  Unversioned 'module load gcc' gives 7.5.0 and
# 'cmake' gives 3.19.5; the defaults drift and a silent toolchain change is
# exactly the class of error the provenance block exists to catch.
module purge
module load gcc/11.3.0 cmake/3.30.0 anaconda3/2023.03
echo "== toolchain =="; gcc --version | head -1; cmake --version | head -1; python3 -V

# ── 1. SUNDIALS ─────────────────────────────────────────────────────────────
if [ -f "$PREFIX/include/cvode/cvode.h" ]; then
  echo "== SUNDIALS already at $PREFIX -- skipping"
else
  echo "== building SUNDIALS $SUNDIALS_VER -> $PREFIX"
  mkdir -p "$BUILD_TMP" && cd "$BUILD_TMP"
  [ -f sundials-$SUNDIALS_VER.tar.gz ] || \
    curl -fL -o sundials-$SUNDIALS_VER.tar.gz \
      "https://github.com/LLNL/sundials/releases/download/v$SUNDIALS_VER/sundials-$SUNDIALS_VER.tar.gz"
  rm -rf sundials-$SUNDIALS_VER && tar xzf sundials-$SUNDIALS_VER.tar.gz
  cmake -S sundials-$SUNDIALS_VER -B build-sun \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
        -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF \
        -DENABLE_LAPACK=ON        `# Woodbury preconditioner path needs it` \
        -DENABLE_OPENMP=OFF       `# campaign runs OMP_NUM_THREADS=1 (plan S11(m))`
  cmake --build build-sun -j 8
  cmake --install build-sun
fi

# ── 2. repo ─────────────────────────────────────────────────────────────────
if [ -d "$REPO/.git" ]; then
  echo "== repo present, pulling"; git -C "$REPO" pull --ff-only || \
    echo "   (pull skipped -- resolve by hand if it diverged)"
else
  echo "== cloning repo -> $REPO"
  git clone https://github.com/<<<SET: your org>>>/RadCluster.git "$REPO"
fi

# ── 3. solver ───────────────────────────────────────────────────────────────
echo "== building the RadCluster solver"
export CMAKE_PREFIX_PATH="$PREFIX:${CMAKE_PREFIX_PATH:-}"
cmake -S "$REPO/RadCluster_2_1/cpp_utils" -B "$REPO/RadCluster_2_1/build" \
      -DCMAKE_BUILD_TYPE=Release
cmake --build "$REPO/RadCluster_2_1/build" -j 8
ls -la "$REPO/RadCluster_2_1/build/solver"

# ── 4. python deps ──────────────────────────────────────────────────────────
# anaconda3/2023.03 already carries numpy 1.23.5 / scipy 1.10.0 / pandas 1.5.3.
# Install into a --user prefix rather than the shared module.
echo "== python deps"
python3 -m pip install --user -q -r "$REPO/requirements.txt" 2>&1 | tail -3 || \
  echo "   (some deps may already be satisfied by the anaconda module)"

# ── 5. THE GATE: does this build agree with the other machines? ─────────────
# Hoffman2's gcc 11.3.0 on CentOS 7 is a THIRD binary, after the Mac's clang and
# MATRIX-PC2's. Plan S11(l) measured Mac vs MATRIX-PC2 at 2.95e-08 and judged
# them poolable -- but that was MEASURED. Do not contribute rows until this
# passes here too.
echo
echo "== check_machine (agreement gate) =="
cd "$REPO/RadCluster_2_1/digital_twin"
python3 check_machine.py; rc=$?
if [ $rc -eq 0 ]; then echo ok > .agree.ok; else echo FAIL > .agree.ok; fi
echo "   check_machine rc=$rc"
echo
echo "NOTE: check_machine compares 12 quantities at rtol=1e-9 with NO absolute"
echo "floor. Two of them (delta_FP, pile) are legitimately ~1e-8/1e-9, where a"
echo "relative comparison is meaningless -- plan S11(l) saw a 54% 'disagreement'"
echo "between two values that were both ~1e-9 against a 0.05 bar. If this failed"
echo "ONLY on those, it is the known defect (gate S11(i)-7), not a bad build."
echo
echo "Now measure the row cost HERE before setting weights or --timeout-s:"
echo "  cores are Xeon Gold 6338N @ 2.20GHz, likely slower per core than the Mac,"
echo "  so campaign_layout.py --speed hoffman2:<factor> must be measured, not guessed."
