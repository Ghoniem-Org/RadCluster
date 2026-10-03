# Code reference (Doxygen)

Doxygen builds a browsable reference of RadCluster's code from its docstrings and header comments.

## Build it

Install [Doxygen](https://www.doxygen.nl/download.html) (1.9 or later; 1.9.8 was used to check this configuration) and Python 3, then from the repository root:

```bash
doxygen docs/doxygen/Doxyfile
```

Open `docs/doxygen/build/html/index.html`. The build takes a few seconds and writes about 18 MB (840 files). `docs/doxygen/build/` is not tracked by git.

## What it covers

| Covered | Source |
|---|---|
| the Python package: abstract core, EUROFER-97 declaration, inputs, rates, orchestration, post-processing | `radcluster_code/py_utils/` |
| the C++ solver: CVODE driver, preconditioners, sparse Jacobian, EUROFER-97 kernels | `radcluster_code/cpp_utils/` |
| the physics and numerics checks, and the figure generators | `radcluster_code/codes/` |
| the calibration and verification campaign | `radcluster_code/digital_twin/` |
| the physics and solver reference, as a page with its equations | `radcluster_code/CLAUDE.md` |
| the main page | `docs/doxygen/mainpage.md` |

Not covered: the notebooks, the input workbook, the experimental database, the derivations and papers in `docs/Formulation/`, and the campaign's design and result files. Build directories, caches, run outputs and `archive/` are excluded.

The code is not written with Doxygen commands. Python docstrings are shown as written, and every class, function and file is listed whether or not it has a comment (`EXTRACT_ALL`). The C++ headers already use `/** ... */` comments, which Doxygen reads directly.

## Files

| File | Contents |
|---|---|
| `Doxyfile` | the configuration: only the settings that differ from Doxygen's defaults |
| `mainpage.md` | the main page of the reference |
| `markdown_math_filter.py` | rewrites the `$...$` and `$$...$$` equations of `radcluster_code/CLAUDE.md` as Doxygen formulas while Doxygen reads it; the file on disk is not changed |

When the version changes, update `PROJECT_NUMBER` in `Doxyfile`.

## Equations

The equations are typeset in the browser by MathJax, which the pages load from a public server (jsDelivr). Without a network connection the rest of the reference works and the equations appear as LaTeX source.

## Publishing

`.github/workflows/documentation.yml` builds the reference on every push to `main` and keeps the HTML as a workflow artifact. It publishes to GitHub Pages only when the workflow is run by hand from the Actions tab, and only after Pages has been enabled in the repository settings with "GitHub Actions" as its source. Nothing is published until then.

## The wiki

The pages of the [wiki](https://github.com/Ghoniem/RadCluster/wiki) are kept in `docs/wiki/`, one Markdown file per page. Edit them there, and copy them to the wiki's own repository to publish:

```bash
git clone https://github.com/Ghoniem/RadCluster.wiki.git
cp docs/wiki/*.md RadCluster.wiki/
cd RadCluster.wiki && git add -A && git commit -m "Update from docs/wiki" && git push
```

## Known warnings

A build gives no warnings. They are written to `docs/doxygen/build/warnings.log`.
