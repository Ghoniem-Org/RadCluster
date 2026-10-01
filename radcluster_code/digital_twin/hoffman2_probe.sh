#!/bin/bash
# hoffman2_probe.sh -- READ-ONLY reconnaissance of a cluster login node.
#
# Writes nothing, submits nothing, installs nothing.  Every command here either
# prints an environment fact or reports that it could not.  Run it before
# filling in the <<<SET>>> values in hoffman2_array.sh -- those are guesses
# until this says otherwise, and a wrong h_rt or module name costs a whole
# array job's queue wait to discover.
#
#   ssh -S ~/.ssh/cm-hoffman2 host 'bash -s' < hoffman2_probe.sh
#
set -u   # NOT -e: a missing tool is a RESULT here, not a failure.

hr(){ printf '\n== %s %s\n' "$1" "$(printf '%.0s-' $(seq 1 $((60-${#1}))))"; }
have(){ command -v "$1" >/dev/null 2>&1; }
try(){ if have "$1"; then "$@" 2>&1 | head -"${HEAD:-6}"; else echo "  ($1 not found)"; fi; }

hr "identity and host"
echo "  whoami   : $(whoami)"
echo "  hostname : $(hostname -f 2>/dev/null || hostname)"
echo "  home     : $HOME"
echo "  shell    : $SHELL"
echo "  uname    : $(uname -srm)"
[ -r /etc/os-release ] && . /etc/os-release && echo "  os       : ${PRETTY_NAME:-?}"

hr "scheduler"
# Hoffman2 has historically been UGE/SGE (qsub), but clusters migrate to Slurm.
# The job script is written for UGE, so this decides whether it needs rewriting.
for s in qsub qstat qhost sbatch squeue sinfo; do
  if have $s; then echo "  FOUND $s -> $(command -v $s)"; else echo "  absent $s"; fi
done
if have qstat; then echo "  --- qstat -help (first lines) ---"; HEAD=4 try qstat -help; fi
if have sbatch; then echo "  --- sinfo ---"; HEAD=8 try sinfo; fi

hr "queue / resource limits  (what h_rt and core count may I ask for?)"
if have qconf; then
  echo "  queues:"; HEAD=20 try qconf -sql
  for q in $(qconf -sql 2>/dev/null | head -6); do
    echo "  --- $q h_rt/h_data ---"
    qconf -sq "$q" 2>/dev/null | grep -E "h_rt|h_data|slots" | sed 's/^/     /'
  done
else
  echo "  (qconf absent -- cannot read queue config)"
fi
if have myquota; then echo "  --- myquota ---"; HEAD=12 try myquota; fi
if have quota;   then echo "  --- quota ---";   HEAD=8  try quota -s; fi

hr "modules"
for init in /u/local/Modules/default/init/bash \
            /u/local/Modules/default/init/modules.sh \
            /etc/profile.d/modules.sh; do
  [ -r "$init" ] && { . "$init" 2>/dev/null && echo "  sourced $init"; break; }
done
if have module || type module >/dev/null 2>&1; then
  for m in gcc cmake python intel sundials cuda; do
    echo "  --- module avail $m ---"
    (module avail "$m") 2>&1 | grep -vi "^--*$" | head -8 | sed 's/^/     /'
  done
else
  echo "  (no module command)"
fi

hr "toolchain already on PATH"
HEAD=2 try gcc --version
HEAD=2 try g++ --version
HEAD=1 try cmake --version
HEAD=1 try make --version
HEAD=1 try git --version
for p in python3 python3.11 python3.12 python; do
  have $p && echo "  $p -> $($p -V 2>&1) at $(command -v $p)"
done

hr "SUNDIALS (the solver needs CVODE; is it provided or must we vendor it?)"
found=0
for d in /u/local/apps/sundials /u/local/apps/SUNDIALS /usr/local/include/sundials \
         /usr/include/sundials; do
  [ -e "$d" ] && { echo "  FOUND $d"; found=1; }
done
have pkg-config && pkg-config --list-all 2>/dev/null | grep -i sundials | sed 's/^/  pkg-config: /'
[ $found -eq 0 ] && echo "  none in the usual places -- expect to build SUNDIALS too"

hr "filesystems and scratch"
echo "  \$SCRATCH  : ${SCRATCH:-(unset)}"
echo "  \$TMPDIR   : ${TMPDIR:-(unset)}"
for d in "$HOME" "${SCRATCH:-}" /u/scratch /u/home; do
  [ -n "$d" ] && [ -d "$d" ] && df -h "$d" 2>/dev/null | tail -1 | sed "s|^|  $d -> |"
done

hr "outbound network FROM THIS LOGIN NODE"
# Matters for transport: if the login node can reach github, results can be
# pushed from here after an array job drains.  Compute nodes usually cannot,
# which is why the job script does NOT try to push.
for host in github.com; do
  if have curl; then
    code=$(curl -s -o /dev/null -w '%{http_code}' -m 8 "https://$host" 2>/dev/null)
    echo "  https://$host -> HTTP ${code:-FAIL}"
  elif have wget; then
    wget -q --spider -T 8 "https://$host" && echo "  https://$host -> OK" || echo "  https://$host -> FAIL"
  else
    echo "  (no curl/wget)"
  fi
done
have git && { GIT_TERMINAL_PROMPT=0 timeout 12 git ls-remote https://github.com/git/git HEAD >/dev/null 2>&1 \
  && echo "  git https clone: OK" || echo "  git https clone: FAILED/blocked"; }

hr "login-node CPU (NOT what jobs get -- for reference only)"
if [ -r /proc/cpuinfo ]; then
  echo "  model : $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
  echo "  cores : $(grep -c ^processor /proc/cpuinfo)"
fi
have qhost && { echo "  --- qhost (compute node inventory, first rows) ---"; HEAD=12 try qhost; }

hr "existing RadCluster checkout, if any"
for d in "$HOME/RadCluster" "$HOME/Documents/RadCluster" "${SCRATCH:-/nonexistent}/RadCluster"; do
  [ -d "$d" ] && echo "  FOUND $d (git: $(git -C "$d" rev-parse --short HEAD 2>/dev/null || echo 'not a repo'))"
done
echo "  (none listed above means nothing is checked out yet)"

hr "done"
echo "  Probe was read-only: no files written, no jobs submitted."
