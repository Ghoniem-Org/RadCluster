"""
gitsync.py — multi-machine coordination for the verification campaign, over git.

The machines share nothing but the GitHub remote, so the remote is the lock.
Two operations matter:

  claim(run_id)   -- atomically take ownership of a run, or discover someone
                     else already has it
  publish(...)    -- record a finished run and push it, so every other machine
                     (and the user, from anywhere) sees the result immediately

WHY A SEPARATE WORKTREE.  Every git operation here runs in a dedicated worktree
under `.sync/`, never in the checkout the solver is running from.  This is not
tidiness.  `run_ensemble.py` carries a 20-line comment about the day a
`git pull --rebase --autostash` replaced a results file underneath a running
worker: a file handle names an inode, not a path, so the worker kept happily
appending to an unlinked inode and 5 rows / 6.2 core-hours were written into
nothing for three hours before anyone noticed.  A campaign that pulls and
pushes after every run would re-arm that trap on a much shorter fuse.  Confining
git to a worktree the solver never reads means a pull cannot move anything the
solver is holding open.

WHY ONE FILE PER RUN.  A claim is `claims/<run_id>.json`, written by exactly one
machine and never edited by another.  Two machines therefore cannot produce a
CONTENT conflict -- only a push race, where the loser fetches, replays its one
commit on top and pushes again.  That is a retry, not a merge.  The alternative
(one shared ledger file appended by everyone) conflicts on every concurrent
finish and needs a real merge strategy.

The branch is separate from `main` so a day of automated commits does not land
in the history of the paper.
"""
from __future__ import annotations

import errno
import fcntl
import json
import os
import platform
import random
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_BRANCH = "campaign-verification"
SYNC_DIR = HERE / ".sync"          # gitignored worktree
PUSH_RETRIES = 6


class SyncError(RuntimeError):
    pass


def _pid_alive(pid) -> bool:
    """Is a process with this pid running?  Used to tell a live sibling worker
    from a claim left behind by one that died."""
    try:
        os.kill(int(pid), 0)
    except (OSError, TypeError, ValueError) as exc:
        return getattr(exc, "errno", None) == errno.EPERM   # alive, not ours
    return True


def _git(args, cwd, check=True, timeout=180):
    p = subprocess.run(["git", *args], cwd=str(cwd), capture_output=True,
                       text=True, timeout=timeout)
    if check and p.returncode != 0:
        raise SyncError(f"git {' '.join(args)} failed in {cwd}:\n"
                        f"{p.stdout}\n{p.stderr}")
    return p


MACHINE_ID_FILE = Path.home() / ".radcluster_machine_id"


def machine_id() -> str:
    """Stable, UNIQUE identity for this host.

    It was platform.node().split(".")[0], which is neither.  This machine's
    hostname changed from "MacBook-Pro" to "Mac.san.rr.com" mid-campaign (a
    router handing out a domain), so it silently started claiming as "Mac" --
    the id the OTHER machine already used.  Two hosts sharing an id is not
    cosmetic: claim() trusts pid liveness to tell a live sibling worker from a
    dead one, and a pid is only meaningful on the machine that issued it.  A
    stale-looking claim could then be taken over while its real owner was still
    computing, which is precisely the duplicate work the board exists to stop.

    So the id is persisted on first use and never derived from a name that can
    change underneath us.  RADCLUSTER_MACHINE still overrides, for naming a
    machine deliberately.
    """
    env = os.environ.get("RADCLUSTER_MACHINE")
    if env:
        return env
    try:
        cached = MACHINE_ID_FILE.read_text(encoding="utf-8").strip()
        if cached:
            return cached
    except OSError:
        pass
    import uuid
    ident = f"{platform.node().split('.')[0]}-{uuid.uuid4().hex[:4]}"
    try:
        MACHINE_ID_FILE.write_text(ident + "\n", encoding="utf-8")
    except OSError:
        pass                 # unwritable home: fall back to per-process id
    return ident


class CampaignSync:
    """Git-backed claim board for one campaign branch."""

    def __init__(self, repo_root: Path, branch: str = DEFAULT_BRANCH,
                 sync_dir: Path = SYNC_DIR, remote: str = "origin",
                 offline: bool = False):
        self.repo = Path(repo_root).resolve()
        self.branch = branch
        self.dir = Path(sync_dir).resolve()
        self.remote = remote
        self.offline = offline
        self.machine = machine_id()

    # WORKTREE MUTEX.  Several workers on ONE machine share this worktree, and
    # concurrent `git add`/`commit`/`push` in a single worktree race on
    # .git/index.lock -- one of them fails, and which one is luck.  Every
    # sequence of git commands therefore runs under an flock on a sidecar file.
    # Cross-machine exclusion is still the remote's job; this only serialises
    # the local processes that share these files.
    @contextmanager
    def _lock(self, timeout_s=300):
        self.dir.parent.mkdir(parents=True, exist_ok=True)
        lf = self.dir.parent / ".sync.lock"
        deadline = time.time() + timeout_s
        fh = open(lf, "a+")
        try:
            while True:
                try:
                    fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError:
                    if time.time() > deadline:
                        raise SyncError(f"worktree lock busy for {timeout_s}s")
                    time.sleep(0.2 + random.random() * 0.3)
            yield
        finally:
            try:
                fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
            finally:
                fh.close()

    # ── setup ────────────────────────────────────────────────────────────────
    def ensure(self):
        """Create the worktree and branch if absent; make it current."""
        with self._lock():
            if not (self.dir / ".git").exists():
                self._create_worktree()
            self.refresh()
        (self.dir / "claims").mkdir(exist_ok=True)
        (self.dir / "results").mkdir(exist_ok=True)

    def _create_worktree(self):
        self.dir.parent.mkdir(parents=True, exist_ok=True)
        if not self.offline:
            _git(["fetch", self.remote], self.repo, check=False)
        remote_ref = f"{self.remote}/{self.branch}"
        has_remote = _git(["rev-parse", "--verify", "-q", remote_ref],
                          self.repo, check=False).returncode == 0
        has_local = _git(["rev-parse", "--verify", "-q", self.branch],
                         self.repo, check=False).returncode == 0
        if has_remote:
            _git(["worktree", "add", "--checkout", str(self.dir),
                  "-B", self.branch, remote_ref], self.repo)
        elif has_local:
            _git(["worktree", "add", str(self.dir), self.branch], self.repo)
        else:
            # First machine to start the campaign creates the branch, ORPHANED.
            # Rooting it at HEAD instead would check the entire repository out a
            # second time -- measured at 210 MB here, duplicated on every
            # machine, to hold a directory of small JSON files.  The link back
            # to the code is kept where it belongs: each claim records the
            # `code_sha` it ran at.
            #
            # --no-checkout means those 210 MB are never written at all: the
            # index is populated but no file lands, and `rm -rf --cached` then
            # empties the index before the first commit.
            _git(["worktree", "add", "--detach", "--no-checkout",
                  str(self.dir), "HEAD"], self.repo)
            _git(["checkout", "--orphan", self.branch], self.dir)
            _git(["rm", "-rf", "--cached", "."], self.dir, check=False)
            # `checkout --orphan` carries the index into the working tree, so
            # the files land anyway; emptying the index does not remove them.
            # This runs BEFORE claims/ and results/ are created, so there is
            # nothing of ours for it to take -- and the worktree is a dedicated,
            # gitignored directory, so `clean -fdx` here cannot reach the user's
            # checkout.
            _git(["clean", "-fdxq"], self.dir, check=False)
        for sub in ("claims", "results"):
            (self.dir / sub).mkdir(exist_ok=True)
            keep = self.dir / sub / ".gitkeep"
            if not keep.exists():
                keep.write_text("")
        self._commit("campaign: initialise claim board", paths=["claims", "results"])
        self._push()

    def refresh(self):
        """Bring the worktree to the remote's state, discarding nothing of ours.

        Our own commits are preserved by rebasing them on top; the worktree
        holds only campaign bookkeeping, so there is never local work here that
        did not come from this module.
        """
        if self.offline:
            return
        if _git(["fetch", self.remote, self.branch], self.dir,
                check=False).returncode != 0:
            return                       # branch not on the remote yet
        r = _git(["rebase", f"{self.remote}/{self.branch}"], self.dir, check=False)
        if r.returncode != 0:
            _git(["rebase", "--abort"], self.dir, check=False)
            _git(["reset", "--hard", f"{self.remote}/{self.branch}"], self.dir)

    # ── claim board ──────────────────────────────────────────────────────────
    def claims(self) -> dict:
        """Every claim currently on the board, by run_id."""
        out = {}
        d = self.dir / "claims"
        if not d.exists():
            return out
        for f in sorted(d.glob("*.json")):
            try:
                out[f.stem] = json.loads(f.read_text(encoding="utf-8"))
            except Exception:
                out[f.stem] = {"run_id": f.stem, "status": "unreadable"}
        return out

    def claim(self, run_id: str, meta: dict = None) -> bool:
        """Try to take `run_id`.  True if we own it, False if someone else does.

        The remote decides.  We write the claim, commit and push; if the push is
        rejected we refresh and look again -- if the run is now claimed by
        another machine we lose it cleanly, and if it is still free we retry.
        """
        for attempt in range(PUSH_RETRIES):
          with self._lock():
            self.refresh()
            existing = self.claims().get(run_id)
            if existing:
                st = existing.get("status")
                owner, opid = existing.get("machine"), existing.get("pid")
                if st == "done":
                    return False              # already computed
                if owner != self.machine:
                    return False              # another machine owns it
                # SAME MACHINE.  This used to return True on the reasoning that
                # a claim by this machine must be ours resuming after a restart.
                # With several workers on one box that is false, and it handed
                # TWO of them the same run: both took T1_B1 and started
                # identical multi-hour solves.  The owning PROCESS is what
                # decides -- a live sibling keeps it, a dead one leaves a stale
                # claim we may take over.
                if opid and int(opid) == os.getpid():
                    return True               # genuinely ours
                if opid and _pid_alive(opid):
                    return False              # a live worker on this box owns it
                # else: stale claim from a process that died -> take it over
            rec = {
                "run_id": run_id, "status": "running",
                "machine": self.machine, "host": platform.node(),
                "pid": os.getpid(),
                "claimed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "code_sha": self._head_sha(),
            }
            rec.update(meta or {})
            self._write_claim(rec)
            self._commit(f"campaign: claim {run_id} [{self.machine}]",
                         paths=[f"claims/{run_id}.json"])
            if self._push():
                return True
          # Lost the race (or a transient failure): back off and re-examine.
          # Outside the lock, so a sibling can make progress meanwhile.
          time.sleep(1.0 + random.random() * 2.0 * (attempt + 1))
        raise SyncError(f"could not settle a claim for {run_id} after "
                        f"{PUSH_RETRIES} attempts")

    def next_unclaimed(self, run_ids: list[str]) -> str | None:
        self.refresh()
        held = self.claims()
        for rid in run_ids:
            c = held.get(rid)
            if c is None:
                return rid
            # A run whose owner died mid-flight stays 'running' forever; that is
            # deliberate.  Reassigning it automatically would risk two machines
            # computing it at once.  `--release` is the explicit way back.
            if c.get("status") == "failed" and c.get("machine") == self.machine:
                continue
        return None

    # ── publishing results ───────────────────────────────────────────────────
    def publish(self, run_id: str, summary: dict, artifacts: dict = None) -> bool:
        """Record a finished run and push it.  Called after EVERY run."""
        with self._lock():
            return self._publish_locked(run_id, summary, artifacts)

    def _publish_locked(self, run_id, summary, artifacts):
        self.refresh()
        rec = self.claims().get(run_id, {"run_id": run_id})
        rec.update(summary)
        rec["machine"] = self.machine
        rec["published_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        self._write_claim(rec)
        paths = [f"claims/{run_id}.json"]
        rdir = self.dir / "results" / run_id
        rdir.mkdir(parents=True, exist_ok=True)
        for name, content in (artifacts or {}).items():
            safe = os.path.basename(str(name))
            tgt = rdir / safe
            if isinstance(content, (bytes, bytearray)):
                tgt.write_bytes(content)
            else:
                tgt.write_text(str(content), encoding="utf-8")
            paths.append(f"results/{run_id}/{safe}")
        status = rec.get("status", "?")
        self._commit(f"campaign: {run_id} {status} [{self.machine}]", paths=paths)
        for attempt in range(PUSH_RETRIES):
            if self._push():
                return True
            time.sleep(1.0 + random.random() * 2.0 * (attempt + 1))
            self.refresh()
        return False

    def release(self, run_id: str) -> bool:
        """Drop a claim so another machine can take the run."""
        with self._lock():
            return self._release_locked(run_id)

    def _release_locked(self, run_id):
        self.refresh()
        f = self.dir / "claims" / f"{run_id}.json"
        if not f.exists():
            return False
        f.unlink()
        self._commit(f"campaign: release {run_id} [{self.machine}]",
                     paths=[f"claims/{run_id}.json"])
        return self._push()

    # ── plumbing ─────────────────────────────────────────────────────────────
    def _write_claim(self, rec):
        d = self.dir / "claims"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{rec['run_id']}.json").write_text(
            json.dumps(rec, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8")

    def _head_sha(self):
        r = _git(["rev-parse", "--short", "HEAD"], self.repo, check=False)
        return r.stdout.strip() if r.returncode == 0 else "unknown"

    def _commit(self, msg, paths):
        _git(["add", "-A", *paths], self.dir, check=False)
        if not _git(["diff", "--cached", "--quiet"], self.dir,
                    check=False).returncode:
            return False                 # nothing staged
        _git(["commit", "-m", msg], self.dir)
        return True

    def _push(self) -> bool:
        if self.offline:
            return True
        return _git(["push", self.remote, f"HEAD:{self.branch}"], self.dir,
                    check=False).returncode == 0
