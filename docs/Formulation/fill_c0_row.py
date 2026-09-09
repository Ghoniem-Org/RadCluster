#!/usr/bin/env python3
"""Fill Table 4's finer-than-C1 row in campaign_results.tex from the board.

The rung is whichever of C0 (i_d 1600) / C0b (800) / C0c (600) has completed --
they are alternatives, not a ladder, added in decreasing cost as each proved
slower than expected.  If none has, the row says so rather than being dropped:
a missing rung is a fact about the campaign, not a gap to hide.
"""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
SYNC = HERE.parent.parent / "RadCluster_2_1/digital_twin/verification/.sync/claims"
D = {"N111": 12.568, "N100": 1.065, "d100": 5.741, "Nv": 2.948, "dcav": 2.544}

def pct(v, k):
    return f"{v:.3f} (${(v-D[k])/D[k]*100:+.1f}\\%$)"

row = None
for rid, lbl, res in (("T4_C0_nc", "C0", "1/2.5"), ("T4_C0b_nc", "C0b", "1/5"),
                      ("T4_C0c_nc", "C0c", "1/6.7")):
    f = SYNC / f"{rid}.json"
    if not f.exists():
        continue
    j = json.loads(f.read_text())
    v = (j.get("scored") or {}).get("2") or {}
    if j.get("status") != "done" or not v or v.get("missing"):
        continue
    row = (f"{lbl} (${res}$)   & {j.get('N_eq','---')} "
           f"& {pct(v['N_loops_111']/1e21,'N111')} & {pct(v['N_loops_100']/1e21,'N100')} "
           f"& {pct(v['d_100_nm'],'d100')} & {pct(v['N_voids']/1e20,'Nv')} "
           f"& {pct(v['d_cavity_nm'],'dcav')} \\\\")
    break

if row is None:
    row = ("C0 (finer than C1) & 1226--3214 & \\multicolumn{5}{l}{\\dnr{} --- the "
           "no-coarsening variant is markedly stiffer the more of the frozen "
           "spectrum is resolved discretely; C1 ran in \\SI{12}{\\second}, "
           "$i_d=600$--1600 did not complete} \\\\")

tex = HERE / "campaign_results.tex"
s = tex.read_text()
s = re.sub(r"^C0 \(\$1/2\$\).*\\\\$", lambda _m: row, s, count=1, flags=re.M)
tex.write_text(s)
print("filled:", row[:110])
