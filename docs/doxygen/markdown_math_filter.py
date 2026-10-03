"""Input filter for Doxygen: rewrite GitHub-style Markdown math as Doxygen formulas.

Doxygen runs this on radcluster_code/CLAUDE.md (see FILTER_PATTERNS in the Doxyfile) and reads the result
from standard output. GitHub renders ``$...$`` and ``$$...$$``; Doxygen needs ``\\f$...\\f$`` and
``\\f[...\\f]``. Nothing else in the file is changed, and the file on disk is never written.

Fenced code blocks and inline code spans are left alone, so a dollar sign in a shell command stays a dollar
sign. The number of lines is preserved, so Doxygen's line numbers still point at the source file.

    python3 docs/doxygen/markdown_math_filter.py radcluster_code/CLAUDE.md
"""
import re
import sys

_DISPLAY = re.compile(r"\$\$(.+?)\$\$", re.DOTALL)
_INLINE = re.compile(r"(?<![\\$])\$(?!\$)([^$\n]+?)(?<![\\$])\$(?!\$)")
_CODE_SPAN = re.compile(r"(`+)(.+?)\1")
_FENCE = re.compile(r"^\s*(```|~~~)")


def _convert_inline(text):
    """Convert ``$...$`` outside inline code spans."""
    pieces = []
    last = 0
    for match in _CODE_SPAN.finditer(text):
        pieces.append(_INLINE.sub(lambda m: r"\f$" + m.group(1) + r"\f$", text[last:match.start()]))
        pieces.append(match.group(0))
        last = match.end()
    pieces.append(_INLINE.sub(lambda m: r"\f$" + m.group(1) + r"\f$", text[last:]))
    return "".join(pieces)


def convert(markdown):
    """Return ``markdown`` with its math rewritten as Doxygen formulas."""
    out = []
    prose = []

    def flush():
        if prose:
            block = "".join(prose)
            block = _DISPLAY.sub(lambda m: r"\f[" + m.group(1) + r"\f]", block)
            out.append("".join(_convert_inline(line) for line in block.splitlines(keepends=True)))
            prose.clear()

    in_fence = False
    for line in markdown.splitlines(keepends=True):
        if _FENCE.match(line):
            flush()
            in_fence = not in_fence
            out.append(line)
        elif in_fence:
            out.append(line)
        else:
            prose.append(line)
    flush()
    return "".join(out)


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as source:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stdout.write(convert(source.read()))
