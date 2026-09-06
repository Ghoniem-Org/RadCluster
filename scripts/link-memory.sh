#!/usr/bin/env bash
#
# link-memory.sh — point this machine's Claude Code memory at the copy in git.
#
# WHY.  Claude Code keeps per-project memory in
# ~/.claude/projects/<slug>/memory/, where <slug> is the project's ABSOLUTE PATH
# with the separators replaced.  That makes it machine-local twice over: it
# depends on your home directory and on where you cloned the repo, and it sits
# in the same directory as the session transcripts (tens of MB of .jsonl) that
# must never be synced.  So the memory is committed at .claude/memory/ and each
# machine symlinks the harness path to it.  Run this once per machine, from
# anywhere inside the checkout.
#
# Reverting is just: rm the symlink.  The files stay in git either way.
#
set -euo pipefail

repo="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
src="$repo/.claude/memory"

[ -d "$src" ] || { echo "error: no $src -- wrong checkout?" >&2; exit 1; }

# The slug transform: absolute path, separators to hyphens (leading / included,
# so the slug starts with '-').  Verified against this repo on macOS 2026-09-06.
slug="$(printf '%s' "$repo" | sed 's|/|-|g')"
dest="$HOME/.claude/projects/$slug/memory"

# If Claude Code has already run here under a DIFFERENT slug rule (paths
# containing dots or other characters it may sanitise further), say so rather
# than silently creating a second, unread directory.
if [ -d "$HOME/.claude/projects" ]; then
  # NOTE the `--` and -F on both greps: every slug begins with '-' (the leading
  # path separator), which grep otherwise parses as an option bundle.
  existing="$(ls -1 "$HOME/.claude/projects" 2>/dev/null | grep -iF -- "$(basename "$repo")" || true)"
  if [ -n "$existing" ] && ! printf '%s\n' "$existing" | grep -qxF -- "$slug"; then
    echo "warning: computed slug '$slug' is not among the existing project dirs:"
    printf '  %s\n' $existing
    echo "         Claude Code may sanitise this path differently.  Link the one"
    echo "         it actually uses instead, then tell the repo about it."
  fi
fi

if [ -e "$dest" ] && [ ! -L "$dest" ]; then
  # A real directory is already there: preserve whatever it holds instead of
  # clobbering memories written before this machine was linked up.
  backup="$dest.local-$(date +%Y%m%d%H%M%S)"
  mv "$dest" "$backup"
  echo "note: moved pre-existing local memory to $backup"
  echo "      merge anything worth keeping into $src and commit it."
fi

mkdir -p "$(dirname "$dest")"
ln -sfn "$src" "$dest"
echo "linked $dest"
echo "    -> $src"
