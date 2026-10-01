# Reference Audit: *A Graph-Based Cluster Dynamics Framework for Irradiated Materials* (Ghoniem, 2026)

**Scope:** This audit covers all 156 references: the main text has refs 1–135 and the Supplementary Material adds refs 136–156. Two checks were run:

- **Metadata:** every entry was compared with its authoritative record (Crossref, doi.org, the publisher or OSTI). The fields compared were authors, title, journal, volume, year, pages and DOI.
- **Citation use:** every in-text citation (about 305 citing sentences, extracted from the PDF) was checked against what the cited paper actually contains.

**Audit date:** 16 Sep 2026. Page numbers are PDF pages.

## Summary

| Result | Count | References |
|---|---|---|
| Metadata OK | 111 | — |
| Minor errors (pages, issue, missing DOI, typos, truncated title) | 23 | 1, 2, 4, 5, 12, 24, 35, 41, 42, 45, 46, 48, 49, 77, 99, 105, 107, 108, 118, 125, 131, 133, 136 |
| **Major errors** (wrong authors/title/DOI, or two papers conflated) | **16** | 10, 21, 28, 39, 40, 43, 88, 89, 90, 91, 96, 97, 113, 119, 123, 124 |
| **Not found: probably fabricated** | **5** | 30, 32, 33, 34, 73 |
| Not verifiable online (internal tech report) | 1 | 11 |
| **Citations that misattribute content** | **27 citing sentences** | 8, 10, 22, 23, 24, 28, 30, 32–34, 37, 38, 39, 48, 77, 78, 79, 80, 89, 90, 95, 96, 118, 144, 153 |

Citation-use verdicts overall: about 149 SUPPORTED, 92 PLAUSIBLE, 27 MISMATCH. PLAUSIBLE means the topic matches but a specific number could not be confirmed from the abstract or record, and applies mostly to parameter tables.

---

## 1. Likely fabricated references (remove or replace)

| Ref | As printed | Finding |
|---|---|---|
| [30] | Ferreira & Wirth, *Multispecies cluster dynamics of H, He, V and SIAs in tungsten…*, NME 39 (2024) 101631 | The DOI resolves to M. Islam, a DFT study of SiC polymorphs, in NME **38** (2024) 101631. No paper with this title or these authors exists. |
| [32] | Chen, Caro, Capolungo, Hyde, Marian, *Algorithmic and computational implementations of cluster dynamics…*, CMS 209 (2022) 111431 | The DOI resolves to Martínez-Agustín et al., on Cahn-Hilliard/Swift-Hohenberg pattern formation in CMS 210. No such paper exists. |
| [33] | Ke, Wirth, Devanathan, *Modeling tools for radiation effects…*, Annu. Rev. Mater. Res. 53 (2023) 59–90 | The DOI returns 404, and no such article is in the vol. 53 table of contents. |
| [34] | Zheng, Maurer, Niu, Chookajorn, *Cluster dynamics modeling of microstructure evolution in nuclear materials: a review*, Prog. Nucl. Energy 158 (2023) 104615 | The DOI returns 404 and no record exists. |
| [73] | Singh, Golubov, Trinkaus, Risø-R-1644(EN) (2008) | No record of this report number or title in DTU Orbit, OSTI/ETDEWEB or web search. |

All five are cited only in the p.2 "broad materials portfolio" and review sentences, except [73], which is used on p.17 and supp pp.80–85. Suggested real substitutes, which you should verify before use:

- For the W multispecies case: Blondel et al., FST 71 (2017) 84 [63], which you already cite.
- For CD reviews: Kohnert, Wirth & Capolungo [20] and Xiong et al. [35], both already cited. Also "Microstructure modeling of nuclear structural materials: Recent progress and future directions", Comput. Mater. Sci. (2023), OSTI 2208824.
- For [73]: Trinkaus, Heinisch, Barashev, Golubov, Singh [72] alone, or Heinisch, Singh & Golubov, JNM 283–287 (2000) 737 [147].

[11] (Ghoniem & Cho, UCLA-ENG-7845, 1978) has no online record. Confirm it from your own files.

## 2. Major metadata errors (real paper, wrong entry)

| Ref | Error | Correct entry |
|---|---|---|
| [10] | Hayns 1976: the title does not exist and the DOI 10.1016/0022-3115(76)90163-5 returns 404 | M. R. Hayns, *On the group method for the approximate solution of a hierarchy of rate equations describing nucleation and growth kinetics*, JNM 59 (1976) 175–182, doi:10.1016/0022-3115(76)90132-X. The "transient stages" sentence on p.1 fits Hayns, JNM 56 (1975) 267–274 better. |
| [21] | Wrong title ("tritium release from irradiated vanadium") | Ghoniem, Alhajji, Kaletta, *The effect of helium clustering on its transport to grain boundaries*, JNM 136 (1985) 192–206 |
| [28] | Wrong authors and pages | J. Boisse, C. Domain, C. S. Becquart, JNM 455 (2014) **10–15** |
| [39] | Wrong title | Stewart, Osetskiy, Stoller, *Atomistic studies of formation and diffusion of helium clusters and bubbles in BCC iron*, JNM 417 (2011) 1110–1114 |
| [40] | Wrong year and source type | Wagner & Kampmann, *Homogeneous second phase precipitation*, in P. Haasen (ed.), *Phase Transformations in Materials*, Materials Science and Technology Vol. 5, VCH, Weinheim, **1991**, pp. 213–303 |
| [43] | Wrong authors (title, journal and DOI are correct) | **K. Xu, B. G. Thomas, Y. Wu, H. Wang, H. Kong, Z. Wu**, Metals 10 (2020) 1685 |
| [88] | The DOI points to an Al-amorphous-alloy paper, the authors are wrong, and a bib note leaked into the entry | Z. Jiao, S. Taller, K. Field, G. Yeli, M. P. Moody, G. S. Was, *Microstructure evolution of T91 irradiated in the BOR60 fast reactor*, **JNM 504 (2018) 122–134**, doi:10.1016/j.jnucmat.2018.03.024 (verified) |
| [89] | The DOI points to a SiC corrosion paper, the authors are wrong, and a bib note leaked into the entry | C. Zheng, E. R. Reese, K. G. Field, E. Marquis, S. A. Maloy, D. Kaoumi, *…HT9 after neutron irradiation: effect of dose*, **JNM 523 (2019) 421–433**, doi:10.1016/j.jnucmat.2019.06.019 (verified via OSTI) |
| [90] | Three records conflated: the DOI is Kim et al. 2021, vol./art. no. 551/152971 is an ODS paper, and the title is Zhong & Tan 2024 | The bib-note content (T91/NF616, ATR, 8.2 dpa) matches W. Zhong, T. A. Saleh, L. Tan, *Neutron irradiation induced defects and clustering in NF616 and T91*, JNM 552 (2021) 153001, doi:10.1016/j.jnucmat.2021.153001 |
| [91] | Title and data are from Weiß et al. (= [99]); the DOI points to a U–Mo FIB paper; no such Klimenkov paper exists | Delete it and cite [99] (Weiß, Gaganidze, Aktaa, JNM 426 (2012) 52–58, doi:10.1016/j.jnucmat.2012.03.027, verified). It is cited on p.29 and in Table S21 (supp p.120). |
| [96] | Wrong authors: four are invented and M. Hou is missing | D. A. Terentyev, L. Malerba, M. Hou, PRB 75 (2007) 104108 (verified) |
| [97] | Wrong authors, volume and pages | D. Terentyev, N. Anento, A. Serra, V. Jansson, H. Khater, G. Bonny, JNM **408** (3) (2011) **272–284** |
| [113] | Two papers conflated; the DOI does not resolve | Either De Backer, Sand, Nordlund, Lunéville, Simeone, Dudarev, *Subcascade formation and defect cluster size scaling…*, EPL **115** (2016) 26001, or De Backer, Domain, Becquart, Lunéville, Simeone, Sand, Nordlund, *A model of defect cluster creation in fragmented cascades in metals based on morphological analysis*, J. Phys.: Condens. Matter 30 (2018) 405701 |
| [119] | Wrong authors | C. S. Becquart, R. Ngayam-Happy, P. Olsson, C. Domain |
| [123] | Wrong co-author | M. M. Rahman, **F. El-Mellouhi**, N. Mousseau |
| [124] | Wrong authors | **C. Barouh, T. Schuler, C.-C. Fu, T. Jourdan**, PRB 92 (2015) 104102 |

Also, bib `note` text has leaked into the printed entries for [88]–[91]. Remove those `note` fields.

## 3. Citations that don't support the sentence

| Page | Sentence (abridged) | Ref(s) | Problem and suggestion |
|---|---|---|---|
| p.1 | "void nucleation of quenched-in vacancies [7, 8]" | [8] | Katz & Wiedersich treat nucleation under simultaneous vacancy and interstitial supersaturation, which is irradiation, not quenching. Split the citation. |
| p.1 | "transient stages … mid-1970s [10]" | [10] | The cited title does not exist (see §2). |
| p.2 | "void swelling and He effects on cavity stability in austenitic steels [12, 21, 22]" | [22] | Ortiz & Caturla 2007 is a generic study of cascade defect evolution, not austenitic or He cavities. [21] is He transport, and [12] is not austenitic-specific either. |
| p.2 | "SIA-cluster and Cr-precipitate kinetics in F/M steels [23–27]" | [23], [24] | [23] is an MD primary-damage study. [24] is an FeNiCr interatomic potential for *austenitic* plasticity. Neither is a CD study of F/M steels. |
| p.2 | "He and H-isotope retention in W [28–30]" | [28], [30] | [28] is a DFT/MD self-trapping study, not a retention CD study. [30] is fabricated. |
| p.2 | CD reviews "[32–35]" | [32]–[34] | Fabricated (see §1). |
| p.2 | "tracking HeℓVm … prohibitive [18, 36, 37]" | [37] | Schäublin & Chiu is MD of He-bubble hardening and says nothing about CD state-space cost. |
| p.2 | "reduced mean-He descriptions used heuristically [38, 39]" | [38], [39] | [38] is Fokker–Planck theory of planar interstitial loops with no helium. [39] is atomistic He-cluster diffusion, not a reduced CD. Cite CD papers that actually use a mean-He closure. |
| p.3 | "critical size about 23 interstitials at 330 °C [48]" | [48] | Dudarev et al. give no such number, only temperature regimes. State that it is your own evaluation of their energy expression. |
| p.23 | "glissile SIA clusters [74–79]" | [78], [79] | Wolfer 2007 covers dislocation bias factors and Kohnert & Capolungo 2019 covers the sink strength/bias of DD networks. Neither addresses glissile clusters. Keep [74]–[77]. |
| p.25 | Grouping with a zeroth moment and slope per group removes mass leakage "[77]" | [77] | This is from Golubov, Ovcharenko, Barashev, Singh, PMA 81 (2001) 643, i.e. **[16]**, not the 2000 production-bias paper. |
| p.25 | "He migrates with Em ≈ 0.06 eV [80]" | [80] | Fu et al. 2005 (Nat. Mater.) covers vacancy and SIA clusters. The He value is from **Fu & Willaime, PRB 72 (2005) 064117 = [116]**. Also re-check the "10^6 times faster than vacancies at 400 °C" figure against your Em values. |
| p.29 | "Burgers vector switches a/2⟨111⟩→a⟨100⟩ with dose [88–90]" | [89], [90] | The HT9 study reports a⟨100⟩ dominance at *both* doses, with no switch. [90] is a conflated entry. Rephrase. |
| p.31 | "C and N trap clusters; solutes limit flights to about 22 nm [95–97]" | [95], [96] | Both are pure-metal MD studies of SIA-cluster motion, with no C/N or solute trapping and no 22 nm value. Only [97] (C–cluster interaction) is relevant. Find the source of the 22 nm figure. |
| supp p.67 | Table S3 EUROFER rate parameters "[28, 80, 115–126]" | [28], [118] | Both are **tungsten** papers. Remove them or replace with Fe sources. |
| supp p.80 | Toroidal cross-section σ = π² r_c r_n "[144]" | [144] | Not in Jansson et al. Mark it as your own derivation or cite the correct source. |
| supp p.90 | C/N/solute trapping shortens 1D flights "[95, 144, 153]" | [95], [144], [153] | None of these treats solute trapping. [153] covers intrinsic self-trapping in pure Fe. |

Softer points, which are defensible but worth tightening:

- [1] and [3] are cited for the "capacitance of arbitrary sink shape" generalization on p.17. They give only the spherical result.
- On p.34, the 375–455 °C loop-character crossover is attributed to [100]. It more likely comes from the companion paper [101].
- The Table S20 value for 300 °C/15 dpa attributed to [87] (4.70×10²¹) differs from the primary source [92] (4.9×10²¹).
- The p.3 citation of [44] reads as backing a claim about existing CD codes. Move it next to the preconditioner statement.
- A false "[111]" citation was detected on supp p.91. It is only the Burgers vector "[111] + ½[11̄1̄] → [100] [49]", but check that the LaTeX source has no stray `\cite` there.

## 4. Minor fixes

| Ref | Fix |
|---|---|
| [1] | Add doi:10.1515/zpch-1918-9209 and capitalize the German nouns |
| [2] | pp. 1–89; add doi:10.1103/RevModPhys.15.1 |
| [4] | Issue (2); pp. 471–478 |
| [5] | pp. 335–351 |
| [12] | Second author "D. D. Cho", consistent with [11] |
| [24] | "S. Poncé" (currently garbled as "Póné") |
| [35] | pp. 5785–5802 |
| [41] | Add issue (3–4) |
| [46] | Add doi:10.1214/aoms/1177729893 |
| [48], [49] | The printed DOIs contain a line-break space; check the bib |
| [99] | Add doi:10.1016/j.jnucmat.2012.03.027; "EUROFER97"; "Journal of Nuclear Materials" |
| [107] | Full title ends "…in α-Fe: A molecular dynamics study of 50 keV cascades" |
| [108] | Full title ends "…: A molecular dynamics study" |
| [125] | Full title ends "…: Harmonic and anharmonic contributions" |
| [131] | "Lumping Analysis in…" (no leading "A") |
| [133] | Prentice Hall, Englewood Cliffs, NJ |
| [136] | "GMRES" |

Minor capitalization or missing issue numbers only: [42], [45], [77], [105], [118].

---

## Appendix: per-reference verification detail

The detail below comes from the automated checks: Crossref/DOI records plus abstracts. Verdicts marked PLAUSIBLE mean the topic matches but the specific numbers were not confirmed against the full text.

### [1] Smoluchowski 1917
- **Metadata:** MINOR — Title, author, journal, vol. 92, pp. 129–168 match. De Gruyter/Crossref list the issue as vol. "92U", issue 1, published 1 Nov 1918 (1917 is the traditional citation year, since the paper was submitted in 1916/17). No DOI given; the record's DOI is 10.1515/zpch-1918-9209. Capitalization: "Theorie", "Koagulationskinetik", "Lösungen" are nouns and should be capitalized in German.
- **Source checked:** https://api.crossref.org/works/10.1515/zpch-1918-9209 ; https://www.degruyterbrill.com/document/doi/10.1515/zpch-1918-9209/html
- **Citations:**
  - p.1: SUPPORTED — the paper introduces the coagulation population balance with the 4πDR diffusion-limited kernel.
  - p.17: PLAUSIBLE — the diffusion-limited kernel for spheres comes from Smoluchowski. The generalization to "electrostatic capacitance of the absorbing body" for non-spherical sinks is later work (Ham [5], Waite), not Smoluchowski.
  - supp p.68: SUPPORTED — the Smoluchowski population balance is an infinite coagulation network.
  - supp p.76: SUPPORTED — the 4πrD rate constant for a sphere is Smoluchowski's result.
  - supp p.85: SUPPORTED — Smoluchowski's pair coagulation uses summed diffusivities (D1+D2) and summed radii.
- **Recommended fix:** Add doi:10.1515/zpch-1918-9209, capitalize the German nouns, and optionally note the vol. 92U issue (1917/1918).

### [2] Chandrasekhar 1943
- **Metadata:** MINOR — Author, title, journal, vol. 15, issue 1 and year are correct. The pages are 1–89, but the entry gives only "1". The DOI 10.1103/RevModPhys.15.1 is missing. Capitalization: "Reviews of Modern Physics".
- **Source checked:** https://link.aps.org/doi/10.1103/RevModPhys.15.1 ; https://ui.adsabs.harvard.edu/abs/1943RvMP...15....1C/abstract (via WebSearch)
- **Citations:**
  - p.1: SUPPORTED — Chapter III of the review recasts Smoluchowski's coagulation and diffusion-limited capture theory in stochastic terms.
- **Recommended fix:** "Reviews of Modern Physics 15 (1) (1943) 1–89. doi:10.1103/RevModPhys.15.1".

### [3] Waite 1957a
- **Metadata:** OK — Title, author, Phys. Rev. 107 (2) 463–470 (1957) and DOI all match. The line-break inside the DOI is only an artifact of PDF extraction.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRev.107.463
- **Citations:**
  - p.1: SUPPORTED — Waite's theory treats non-random (correlated pair) initial distributions in diffusion-limited reactions.
  - p.17: PLAUSIBLE — Waite gives the kernel 4πr0D for point-like or spherical reactants. "Electrostatic capacitance" for arbitrary sink shapes is not the focus of this paper.
  - supp p.78: SUPPORTED — Waite uses the reaction radius r0 as the sum of radii, with relative diffusivity D = DA + DB.
  - supp p.85: SUPPORTED — the paper uses summed diffusivities for two mobile species.
- **Recommended fix:** none

### [4] Waite 1957b
- **Metadata:** MINOR — Title, author, volume and year are correct. The issue (2) is missing, and the pages are 471–478 rather than "471" alone. The DOI is correct.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRev.107.471
- **Citations:**
  - p.1: SUPPORTED — the paper applies the diffusion-limited theory to annealing of radiation damage in germanium.
- **Recommended fix:** "Physical Review 107 (2) (1957) 471–478."

### [5] Ham 1958
- **Metadata:** MINOR — Author (Frank S. Ham), title, J. Phys. Chem. Solids 6 (4) and 1958 are correct. The pages are 335–351, but only "335" is given. The DOI is correct.
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3697(58)90053-2
- **Citations:**
  - p.1: SUPPORTED — Ham solves diffusion-limited growth of precipitates, including spherical and ellipsoidal particles in a cell model. That work is the basis of capture by extended, non-point sinks.
- **Recommended fix:** Pages "335–351".

### [6] Dienes & Damask 1958
- **Metadata:** OK — Title, authors, J. Appl. Phys. 29 (12) 1713–1721 (1958) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1063/1.1723032
- **Citations:**
  - p.1: SUPPORTED — the classic radiation-enhanced diffusion rate-equation model, which combines defect production, vacancy–interstitial recombination and sink loss. "First" is a historical judgement that is broadly accepted.
- **Recommended fix:** none

### [7] Kiritani 1973
- **Metadata:** OK — Title, author, J. Phys. Soc. Jpn. 35 (1) 95–107 (1973) and DOI all match. The issue number (1) is missing but optional.
- **Source checked:** https://api.crossref.org/works/10.1143/JPSJ.35.95
- **Citations:**
  - p.1: SUPPORTED — the paper analyses clustering of supersaturated (quenched-in) vacancies.
  - p.25 (lineage): SUPPORTED — Kiritani introduced the size-grouping method in this paper.
  - p.25 (geometric grouping broadens the SDF and drifts the mean size): PLAUSIBLE — Kiritani did introduce geometric grouping. The criticism of its artifacts comes from later work (e.g. Golubov et al. [16]), not from this paper.
- **Recommended fix:** none (optionally add issue 1)

### [8] Katz & Wiedersich 1971
- **Metadata:** OK — Title, authors, J. Chem. Phys. 55 (3) 1414–1425 (1971) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1063/1.1676236
- **Citations:**
  - p.1: MISMATCH (partial) — the sentence groups [7, 8] under "void nucleation of quenched-in vacancies". Katz & Wiedersich treat nucleation in materials supersaturated with both vacancies and interstitials, i.e. under irradiation, not quenching. Kiritani [7] fits the quench case; [8] is the irradiation (two-species) nucleation theory.
- **Recommended fix:** Reword, e.g. "Void nucleation from quenched-in vacancies [7] and under simultaneous vacancy–interstitial supersaturation [8] followed".

### [9] Gruber 1967
- **Metadata:** OK — AIP lists "Calculated Size Distributions for Gas Bubble Migration and Coalescence in Solids", E. E. Gruber, J. Appl. Phys. 38 (1) 243, which matches the entry. The Crossref API was rate-limited (HTTP 429), so the end page 250 and the DOI 10.1063/1.1708962 could not be confirmed from Crossref. They are consistent with the AIP URL (vol. 38, issue 1, p. 243).
- **Source checked:** WebSearch → https://pubs.aip.org/aip/jap/article-abstract/38/1/243/3952/ (Crossref api.crossref.org/works/10.1063/1.1708962 attempted, 429)
- **Citations:**
  - p.1: SUPPORTED — the paper computes bubble size distributions from migration and coalescence.
- **Recommended fix:** none (optionally add issue 1)

### [10] Hayns 1976
- **Metadata:** MAJOR — No paper titled "The transient stages of damage accumulation in irradiated materials with point defect clusters" by M. R. Hayns appears in Crossref. I listed all of Hayns's J. Nucl. Mater. papers from 1975–1977 and searched the web. The given DOI 10.1016/0022-3115(76)90163-5 does not resolve (Crossref 404). The volume and start page (J. Nucl. Mater. 59, p. 175, 1976) belong to a different Hayns paper: "On the group method for the approximate solution of a hierarchy of rate equations describing nucleation and growth kinetics", J. Nucl. Mater. 59 (1976) 175–182, doi:10.1016/0022-3115(76)90132-X. The end page 187 is also wrong. The title as printed looks fabricated or conflated. A related Hayns paper on transient loop evolution is "The nucleation and early growth of interstitial dislocation loops in irradiated materials", J. Nucl. Mater. 56 (1975) 267–274, doi:10.1016/0022-3115(75)90042-2.
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3115(76)90163-5 (404); https://api.crossref.org/works?query.author=Hayns&filter=issn:0022-3115,from-pub-date:1975,until-pub-date:1977 ; WebSearch (no hit for the title)
- **Citations:**
  - p.1 ("transient stages of damage accumulation were modeled in the mid-1970s"): MISMATCH as cited, because the cited title does not exist. If the entry is corrected to Hayns 1976 (group method), the sentence should credit the grouping and nucleation kinetics rather than "transient stages". Hayns 1975 (J. Nucl. Mater. 56) fits "early/transient loop evolution" better.
  - p.25 (grouping lineage): SUPPORTED only if the entry is corrected to the 1976 group-method paper, which is exactly about grouping.
- **Recommended fix:** Replace with "M. R. Hayns, On the group method for the approximate solution of a hierarchy of rate equations describing nucleation and growth kinetics, J. Nucl. Mater. 59 (1976) 175–182. doi:10.1016/0022-3115(76)90132-X". Optionally add Hayns, J. Nucl. Mater. 56 (1975) 267–274 for the p.1 "transient stages" statement, and reword that sentence.

### [11] Ghoniem & Cho 1978 (UCLA-ENG-7845)
- **Metadata:** NOT FOUND (unverifiable) — A UCLA technical report has no DOI. Web searches found no record of UCLA-ENG-7845 or this title, and the UCLA Matrix lab 1970–79 publication page lists only the 1979 phys. stat. sol. (a) paper by these authors. The report may well exist (the author is Ghoniem himself), but it cannot be confirmed online. The co-author appears as "D. D. Cho" in [12].
- **Source checked:** WebSearch; https://www.seas.ucla.edu/matrix/html/pub_70-79.html
- **Citations:**
  - p.1: PLAUSIBLE — the topic is consistent with the follow-up journal paper [12]. Priority ("first fully coupled treatment") cannot be confirmed from an unavailable report.
- **Recommended fix:** The author should confirm the report number and date from their own records. If it is not archived, consider citing only [12]. Use a consistent co-author name ("D. D. Cho").

### [12] Ghoniem & Cho 1979
- **Metadata:** MINOR — Title, phys. stat. sol. (a) 54, 171–178 (1979) and DOI are correct, and issue 1 is missing. Crossref and UCLA list the second author as "D. D. Cho"; the entry has "D. Cho", which is also inconsistent with [11].
- **Source checked:** https://api.crossref.org/works?query.bibliographic=Ghoniem+Cho+simultaneous+clustering...; https://www.seas.ucla.edu/matrix/html/pub_70-79.html
- **Citations:**
  - p.1: SUPPORTED — the paper covers simultaneous vacancy and interstitial clustering under irradiation with rate equations. The priority claim is the author's own judgement ("plausibly").
  - p.2 ("void swelling and helium effects on cavity stability in austenitic stainless steels and Fe–Cr–Ni model alloys [12, 21, 22]"): PLAUSIBLE/weak — the 1979 paper treats clustering of point defects (vacancy and interstitial clusters). It is not a helium-effects or cavity-stability study, so it supports only the "void/vacancy clustering" part of the sentence.
- **Recommended fix:** Change the author to "N. M. Ghoniem, D. D. Cho" and add issue (1). Limit the p.2 citation of [12] to the void-clustering part of the sentence.

### [13] Ghoniem & Sharafat 1980
- **Metadata:** OK — Title, authors, J. Nucl. Mater. 92 (1) 121–135 (1980) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3115(80)90148-8
- **Citations:**
  - p.1: SUPPORTED — this is a Fokker–Planck description of interstitial loop size evolution during irradiation.
- **Recommended fix:** none

### [14] Bacon, Gao & Osetsky 2000
- **Metadata:** OK — Title, authors (Yu. N. Osetsky), J. Nucl. Mater. 276 (1–3) 1–12 (2000) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/S0022-3115(99)00165-8
- **Citations:**
  - p.1: SUPPORTED — an MD review of primary damage on picosecond time scales, which supports the contrast with CD time scales.
  - supp p.65 (fission vs fusion spectrum defect-production table in α-Fe): PLAUSIBLE — the paper gives MD cascade defect production vs PKA energy for bcc Fe and other metals. It does not treat fission or fusion neutron spectra directly, so any spectrum-averaged values attributed to it are unverified.
- **Recommended fix:** none

### [15] Marian et al. 2017
- **Metadata:** OK — Title, all 12 authors in the correct order, Nucl. Fusion 57 (9) 092008 (2017) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1088/1741-4326/aa5e8d
- **Citations:**
  - p.1: PLAUSIBLE — the multiscale-modelling review discusses MD time-scale limits and the role of mean-field/CD methods. However, it is tungsten-specific, which is an odd choice for a general or EUROFER97 time-scale statement.
- **Recommended fix:** none (optionally add a general multiscale review)

### [16] Golubov et al. 2001
- **Metadata:** OK — Title, authors in order, Philos. Mag. A 81 (3) 643–658 (2001) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1080/01418610108212164
- **Citations:**
  - p.2 (cluster-size distribution feeds property models): PLAUSIBLE — this is a grouping-method paper; the general statement is compatible but not its focus.
  - p.2 (grouping closures introduce mass-conservation errors unless designed): SUPPORTED — the paper builds a grouping scheme that conserves defect number and mass, contrasting it with earlier schemes.
  - p.25 (grouping lineage): SUPPORTED.
- **Recommended fix:** none

### [17] Surh, Sturgeon & Wolfer 2004
- **Metadata:** OK — Title, authors, J. Nucl. Mater. 325 (1) 44–52 (2004) and DOI match. An erratum exists: J. Nucl. Mater. (2005), ScienceDirect PII S0022311505000954.
- **Source checked:** https://www.sciencedirect.com/science/article/abs/pii/S0022311503004847
- **Citations:**
  - p.2: PLAUSIBLE — this is a hybrid master-equation/Fokker–Planck void swelling model compared against experiments. It fits "evolution of the size distribution" and "swelling property prediction" in general.
- **Recommended fix:** none (optionally cite the erratum)

### [18] Caturla, Ortiz & Fu 2008
- **Metadata:** OK — Title, authors (M. J. Caturla, C. J. Ortiz, C.-C. Fu), C. R. Physique 9 (3–4) 401–408 (2008) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.crhy.2007.09.004
- **Citations:**
  - p.2 (meeting point of atomistic parameterization and engineering prediction): PLAUSIBLE — a review of kinetic (rate-theory and KMC) modelling of He and point-defect accumulation parameterized by ab initio and MD.
  - p.2 (explicit He–V grid is prohibitive): PLAUSIBLE — the paper discusses He–vacancy cluster kinetics. The explicit statement about the computational cost of the (ℓ, m) grid cannot be confirmed from the abstract.
- **Recommended fix:** none

### [19] Marian & Bulatov 2011
- **Metadata:** OK — Title, authors, J. Nucl. Mater. 415 (1) 84–95 (2011) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2011.05.045
- **Citations:**
  - p.2: SUPPORTED — stochastic CD bridges atomistic inputs and mesoscale evolution.
  - p.13: SUPPORTED — SCD samples the master equation with a Gillespie-type SSA.
  - supp p.68: SUPPORTED.
- **Recommended fix:** none

### [20] Kohnert, Wirth & Capolungo 2018
- **Metadata:** OK — Title, authors, Comput. Mater. Sci. 149, 442–459 (2018) and DOI match. Minor capitalization: "A review".
- **Source checked:** https://api.crossref.org/works/10.1016/j.commatsci.2018.02.049
- **Citations:**
  - p.2 (meeting point of scales): SUPPORTED — a CD methods review.
  - p.2 (shaping reviews, fragmented algorithmic framework): SUPPORTED/PLAUSIBLE — it is a review. The "common diagnosis" is the authors' synthesis.
  - p.2 (grouping/binning closures and mass conservation): SUPPORTED — the review covers grouping methods and their conservation properties.
- **Recommended fix:** none
### [21] Ghoniem 1985
- **Metadata:** MAJOR — Wrong title. The manuscript gives "The effect of helium clustering on tritium release from irradiated vanadium". Crossref gives "The effect of helium clustering on its transport to grain boundaries". Authors (N.M. Ghoniem, J.N. Alhajji, D. Kaletta), JNM 136(2–3) (1985) 192–206 and the DOI all match. The vanadium/tritium title is invented.
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3115(85)90007-8 ; https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/0022-3115(85)90007-8
- **Citations:**
  - p.2: PLAUSIBLE — The real paper is an early rate-theory/cluster-dynamics treatment of helium clustering and transport to grain boundaries, so it fits "helium effects". I could not get the abstract to confirm that it covers void swelling or cavity stability in austenitic or Fe–Cr–Ni alloys.
- **Recommended fix:** Change the title to "The effect of helium clustering on its transport to grain boundaries". Keep the rest (add issue 2–3).

### [22] Ortiz 2007
- **Metadata:** OK — Title, authors (C.J. Ortiz, M.J. Caturla), PRB 75(18) 184101 (2007) and DOI all match. Only a lowercase "role" and a line-break space inside the DOI.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.75.184101 ; Semantic Scholar API (abstract not available)
- **Citations:**
  - p.2: MISMATCH (probable) — This is a generic rate-theory/kinetic Monte Carlo study of cascade defect evolution (intracascade clustering, correlated recombination). Its subject is not void swelling or helium effects on cavity stability in austenitic steels or Fe–Cr–Ni alloys. No abstract was available, and the title shows no helium or austenitic focus.
- **Recommended fix:** Move it to a context about cascade damage and rate theory. For austenitic/helium cavity stability, cite suitable papers instead, e.g. Ghoniem, Sharafat, Williams & Mansur, JNM (1983), or Stoller & Odette.

### [23] Stoller 2000
- **Metadata:** OK — Title, author (R.E. Stoller), JNM 276(1–3) (2000) 22–32 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/S0022-3115(99)00204-4
- **Citations:**
  - p.2: MISMATCH — This is a molecular-dynamics (MD) study of primary damage in displacement cascades in Fe. It is not a cluster-dynamics study of SIA-cluster or Cr-precipitate kinetics in ferritic/F-M steels.
  - supp p.65: SUPPORTED — The paper gives MD cascade defect-production characteristics versus energy and temperature in alpha-Fe, which fits a table comparing defect production.
- **Recommended fix:** Remove it from the p.2 list of cluster-dynamics applications and keep it in the supplementary table.

### [24] Bonny 2011
- **Metadata:** MINOR — The 4th author's name is garbled: manuscript has "S. Póné", Crossref has "S. Poncé". Title, MSMSE 19(8) 085008 (2011) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1088/0965-0393/19/8/085008
- **Citations:**
  - p.2: MISMATCH — This paper builds an interatomic potential for the austenitic FeNiCr alloy for plasticity studies. It is neither a cluster-dynamics paper nor about ferritic or F-M steels.
- **Recommended fix:** Fix the author name to "S. Poncé" and remove the reference from this sentence. For Cr-precipitate/SIA kinetics in ferritic steels, cite a real CD study of Fe–Cr (e.g. Terentyev/Malerba/Bonny Fe–Cr alpha-prime kinetics, or Jourdan et al.).

### [25] Malerba 2021
- **Metadata:** OK — Title, all 11 authors in order, NME 29 (2021) 101069 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.nme.2021.101069
- **Citations:**
  - p.2: PLAUSIBLE — This review compiles physical mechanisms and parameters for microstructure-evolution models (OKMC/CD) in pure Fe. It is relevant to SIA-cluster kinetics, but it covers pure Fe, not Cr precipitates or F-M steels.
  - supp p.65: PLAUSIBLE — The review covers primary damage (cascade defect production) in Fe. I could not confirm fission-vs-fusion spectrum comparisons from the metadata.
- **Recommended fix:** none (optionally add Part II on Fe–Cr alloys for the Cr-precipitate claim)

### [26] Gao 2021a
- **Metadata:** OK — Title, authors (Jie Gao, E. Gaganidze, B. Kaiser, J. Aktaa), JNM 547 (2021) 152822 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2021.152822 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311521000453 (search listing)
- **Citations:**
  - p.2 (list): PLAUSIBLE — A He–V cluster-dynamics study of Eurofer97, which is a ferritic-martensitic steel. It is about He–V clustering, not SIA-cluster or Cr-precipitate kinetics.
  - p.2 (Gao et al. Eurofer97 dual-beam): SUPPORTED — The title states He+/Fe3+ dual-beam irradiation of Eurofer97.
- **Recommended fix:** none

### [27] Gao 2021b
- **Metadata:** OK — Title, authors, JNM 557 (2021) 153212 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2021.153212 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311521004359 (abstract)
- **Citations:**
  - p.2 (list): SUPPORTED — A cluster-dynamics study of He bubbles, C15 clusters and dislocation loops in Eurofer97, which fits SIA-cluster kinetics in F-M steels.
  - p.2 (dual-beam): SUPPORTED — The abstract states Fe3+/He+ dual-ion-beam irradiation of Eurofer97.
  - p.34: SUPPORTED — The abstract confirms coupled He bubbles, C15 clusters and loops under dual-beam Fe3+/He+ irradiation (603–773 K).
- **Recommended fix:** none

### [28] Becquart 2014 (actually Boisse 2014)
- **Metadata:** MAJOR — Author list and pages are wrong. Manuscript: "C.S. Becquart, C. Domain", pp. 311–315. Crossref: J. Boisse, C. Domain, C.S. Becquart, JNM 455(1–3) (2014) pp. 10–15. The title matches, apart from capitalization.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2014.02.031
- **Citations:**
  - p.2: PLAUSIBLE — A DFT/MD study of He self-trapping and trap mutation in W. It is relevant to He retention in W but is not a cluster-dynamics study.
  - supp p.67 (Table S3, reaction-rate parameters for EUROFER CD): MISMATCH — This is a tungsten paper and cannot supply rate parameters for EUROFER/Fe.
- **Recommended fix:** Correct to "J. Boisse, C. Domain, C.S. Becquart, J. Nucl. Mater. 455 (2014) 10–15". Remove it from Table S3 and cite Fe-specific sources for those parameters.

### [29] Faney 2014
- **Metadata:** OK — Title, authors (T. Faney, B.D. Wirth), MSMSE 22(6) 065010 (2014) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1088/0965-0393/22/6/065010
- **Citations:**
  - p.2: SUPPORTED — A spatially dependent cluster-dynamics model of He retention and evolution in He-irradiated W.
- **Recommended fix:** none

### [30] Ferreira 2024
- **Metadata:** NOT FOUND — The DOI 10.1016/j.nme.2024.101631 points to a different paper: M. Islam, "A comprehensive investigation on the physical properties of SiC polymorphs ... DFT study", Nucl. Mater. Energy 38 (2024) 101631. A Crossref bibliographic search and web searches found no "Multispecies cluster dynamics of hydrogen, helium, vacancies and self-interstitials in tungsten" by Ferreira & Wirth. The entry is likely fabricated.
- **Source checked:** https://api.crossref.org/works/10.1016/j.nme.2024.101631 ; Crossref query.bibliographic search (no match) ; WebSearch
- **Citations:**
  - p.2: MISMATCH — The cited paper does not exist, and the DOI target is a SiC DFT paper.
- **Recommended fix:** Replace it with a verified multispecies W cluster-dynamics paper, e.g. Blondel et al., "Benchmarks and tests of a multidimensional cluster dynamics model of helium implantation in tungsten", Fusion Sci. Technol. 71 (2017) 84, or a Xolotl H/He paper (Blondel, Bernholdt, Hammond, Wirth). Verify the chosen reference before citing.

### [31] Veshchunov 2008
- **Metadata:** OK — Title, author (M.S. Veshchunov), JNM 374(1–2) (2008) 44–53 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2007.06.021
- **Citations:**
  - p.2: PLAUSIBLE — The paper models coalescence of grain-face gas bubbles in UO2, which is related to bubble evolution in fuels. It is about coalescence, not nucleation, and is not strictly a cluster-dynamics paper.
- **Recommended fix:** none (optionally add a fuel cluster-dynamics reference, such as fission-gas CD by Veshchunov or Skorek)

### [32] Chen 2022
- **Metadata:** NOT FOUND — The DOI 10.1016/j.commatsci.2022.111431 points to Martínez-Agustín et al., "3D pattern formation from coupled Cahn-Hilliard and Swift-Hohenberg equations...", Comput. Mater. Sci. 210 (2022) 111431. A Crossref search for Chen/Marian in Comput. Mater. Sci. 2022 and an exact-title web search found no such paper. The entry is likely fabricated.
- **Source checked:** https://api.crossref.org/works/10.1016/j.commatsci.2022.111431 ; Crossref query.bibliographic ; WebSearch (exact title)
- **Citations:**
  - p.2: MISMATCH — The cited paper does not exist.
- **Recommended fix:** Replace it with the real cluster-dynamics review A.A. Kohnert, B.D. Wirth, L. Capolungo, "Modeling microstructural evolution in irradiated materials with cluster dynamics methods: A review", Comput. Mater. Sci. 149 (2018) 442–459, doi:10.1016/j.commatsci.2018.02.049 (verified on OSTI 1477651).

### [33] Ke 2023
- **Metadata:** NOT FOUND — The DOI 10.1146/annurev-matsci-080921-072619 returns 404 at doi.org, and the Crossref lookup was rate-limited. The table of contents of Annu. Rev. Mater. Res. vol. 53 (2023) has no article by Ke, Wirth or Devanathan and no radiation-effects modeling review. Pages 53–79 and 81–104 belong to other articles. The entry is likely fabricated.
- **Source checked:** https://doi.org/10.1146/annurev-matsci-080921-072619 (404) ; https://www.annualreviews.org/content/journals/matsci/53/1 ; WebSearch
- **Citations:**
  - p.2: MISMATCH — The cited paper does not exist.
- **Recommended fix:** Delete it, or replace it with a verified multiscale radiation-modeling review.

### [34] Zheng 2023
- **Metadata:** NOT FOUND — Crossref returns 404 for DOI 10.1016/j.pnucene.2023.104615. Web searches for the title and author set (Zheng, Maurer, Niu, Chookajorn) in Prog. Nucl. Energy 158 found nothing. The author set looks implausible for this topic. The entry is likely fabricated.
- **Source checked:** https://api.crossref.org/works/10.1016/j.pnucene.2023.104615 (404) ; WebSearch
- **Citations:**
  - p.2: MISMATCH — The cited paper does not exist.
- **Recommended fix:** Delete it. Real CD reviews already exist, e.g. Kohnert et al. 2018 (see [32]) and Xiong et al. 2024 [35].

### [35] Xiong 2024
- **Metadata:** MINOR — Pages are off by one: manuscript has 5786–5803, Crossref has 5785–5802. Title, all 9 authors, JOM 76(10) (2024) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1007/s11837-024-06717-w
- **Citations:**
  - p.2: SUPPORTED — A review of cluster dynamics for radiation damage whose abstract discusses limitations of CD models and possible improvements, which fits "shaped the field / common diagnosis".
- **Recommended fix:** Change the pages to 5785–5802.

### [36] Trinkaus 2003
- **Metadata:** OK — Title, authors (H. Trinkaus, B.N. Singh), JNM 323(2–3) (2003) 229–242 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2003.09.001
- **Citations:**
  - p.2: PLAUSIBLE — A major review of He accumulation, He–V clustering and bubble nucleation in metals, relevant to the complexity of the He–V state space. I could not confirm an explicit statement that tracking the full (l, m) grid is prohibitive.
- **Recommended fix:** none

### [37] Schäublin 2007
- **Metadata:** OK — Title, authors (R. Schäublin, Y.L. Chiu), JNM 362(2–3) (2007) 152–160 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2007.01.187
- **Citations:**
  - p.2: MISMATCH — This is an MD simulation of how He bubbles affect dislocation pinning and hardening in Fe. It says nothing about the computational cost of tracking He_l V_m clusters in cluster dynamics.
- **Recommended fix:** Remove it. Cite instead He–V CD papers that discuss state-space size, e.g. Marian & Bulatov (2011), Kohnert & Wirth, or Dunn et al. (2015 grouping method).

### [38] Ghoniem 1989
- **Metadata:** OK — Title, author (N.M. Ghoniem), PRB 39(16) 11810–11819 (1989) and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.39.11810
- **Citations:**
  - p.2 (reduced mean-helium descriptions): MISMATCH — The paper is a stochastic (Fokker–Planck) theory of planar interstitial clustering and dislocation loops. It contains no helium or mean-helium-per-cavity closure.
  - p.2 (Fokker–Planck representations): SUPPORTED — It is a stochastic Fokker–Planck treatment of cluster size evolution.
- **Recommended fix:** Remove it from the mean-helium sentence and keep it for the Fokker–Planck sentence.

### [39] Stewart 2011
- **Metadata:** MAJOR — Wrong title. The manuscript gives "Atomistic studies of helium defect properties in bcc iron: comparison of He–V cluster models". Crossref gives "Atomistic studies of formation and diffusion of helium clusters and bubbles in BCC iron". Crossref spells the second author "Yuri Osetskiy". JNM 417(1–3) (2011) 1110–1114 and the DOI match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2010.12.217 ; ORNL/ResearchGate listings via WebSearch
- **Citations:**
  - p.2: MISMATCH — This is an atomistic MD/statics study of how He clusters and bubbles form and diffuse in Fe. It is not a reduced cluster-dynamics description tracking a mean He load per cavity.
- **Recommended fix:** Correct the title. For "mean helium per cavity" reduced descriptions, cite actual CD works using that approximation, e.g. Ghoniem et al. (1983) or Stoller & Odette, or Golubov et al. rate-theory He bubble models.

### [40] Wagner 1981 (actually 1991)
- **Metadata:** MAJOR — The year and source type are wrong. The chapter is R. Wagner, R. Kampmann, "Homogeneous second phase precipitation", in P. Haasen (ed.), Phase Transformations in Materials, Materials Science and Technology: A Comprehensive Treatment, Vol. 5, VCH, Weinheim, 1991, pp. 213–303. It is not a journal, and it is not from 1981; "Materials Science and Technology" here is a book series, not the journal of that name. A later Wiley-VCH online reissue (with P.W. Voorhees added, pp. 213–304) is dated 2006/2013.
- **Source checked:** https://help.thermocalc.com/2025b/content/precipitation/theory-overview.htm ; Northwestern Scholars and Stanford SearchWorks listings via WebSearch (the Crossref query was rate-limited)
- **Citations:**
  - p.2: PLAUSIBLE — The chapter is the standard reference on classical nucleation and growth, including Zeldovich/Fokker–Planck continuum treatments and the KWN method. It fits "Fokker–Planck representations" in a general way.
- **Recommended fix:** Change to "R. Wagner, R. Kampmann, Homogeneous second phase precipitation, in: P. Haasen (Ed.), Phase Transformations in Materials, Materials Science and Technology Vol. 5, VCH, Weinheim, 1991, pp. 213–303."
### [41] Russell 1984
- **Metadata:** MINOR — all core fields match (K.C. Russell, Prog. Mater. Sci. 28, 1984, 229–434, DOI correct). Issue number (3–4) is missing.
- **Source checked:** https://api.crossref.org/works/10.1016/0079-6425(84)90001-X ; https://www.osti.gov/etdeweb/biblio/5720477
- **Citations:**
  - p.2: PLAUSIBLE — Russell's review covers nucleation and growth theory under irradiation, including continuum (Zeldovich/Fokker–Planck-type) treatments of the size distribution. It is a broad phase-stability review, not a paper on Fokker–Planck closures, so it works as background but is not a primary Fokker–Planck CD reference.
- **Recommended fix:** Add issue "(3–4)". For the Fokker–Planck claim, consider also citing a dedicated FP-closure paper (e.g., Ghoniem & Sharafat 1980, or Surh/Sturgeon/Wolfer 2004).

### [42] Ovcharenko 2003
- **Metadata:** MINOR — authors, title, journal, volume 152(2), pages 208–226, 2003 and DOI all match. Only capitalization differs: the record has "Master equations"; the DOI prefix is lowercase "s0010" (harmless).
- **Source checked:** https://api.crossref.org/works/10.1016/s0010-4655(02)00816-0
- **Citations:**
  - p.2: SUPPORTED — GMIC++ is a grouping method for large master-equation sets. The mass-conservation remark is the manuscript's own critique and is reasonable for grouping schemes.
  - p.25: SUPPORTED — a standard reference in the grouping/coarse-graining lineage.
- **Recommended fix:** none (optional: "Master equations" capitalization per record).

### [43] Cui 2020 (as printed)
- **Metadata:** MAJOR — the title, journal, volume, issue, article number, year and DOI (Metals 10(12) 1685, 2020, 10.3390/met10121685) all match a real paper, but the authors are wrong. The manuscript lists "S. Cui, M. Mamivand, D. Morgan". Crossref and the MDPI page give Kun Xu, Brian G. Thomas, Yueyue Wu, Haichuan Wang, Hui Kong, Zhaoyang Wu (Anhui Univ. of Technology / Colorado School of Mines). No Cui–Mamivand–Morgan paper with this title was found.
- **Source checked:** https://api.crossref.org/works/10.3390/met10121685 ; https://www.mdpi.com/2075-4701/10/12/1685 ; WebSearch
- **Citations:**
  - p.2: SUPPORTED — the paper compares size-grouping methods for cluster dynamics (Al3Sc precipitation) against ungrouped solutions, which fits "grouping or binning closures". It also proposes a log-linear within-group distribution to improve accuracy.
- **Recommended fix:** Change authors to "K. Xu, B. G. Thomas, Y. Wu, H. Wang, H. Kong, Z. Wu". Keep title, journal and DOI.

### [44] Saad 2003
- **Metadata:** OK — Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., SIAM, 2003, DOI 10.1137/1.9780898718003 all match. The hyphenation "So- ciety" is only a PDF line break.
- **Source checked:** https://api.crossref.org/works/10.1137/1.9780898718003
- **Citations:**
  - p.3: PLAUSIBLE — Saad supports general statements about preconditioning and iterative solvers for structured sparse systems. It says nothing about CD codes or conservation diagnostics, and the sentence mostly makes a claim about the CD literature. The citation fits only the "standard preconditioners" part.
- **Recommended fix:** Move [44] right after "standard preconditioners perform poorly on this structure", so it is not read as backing the claim about existing CD codes.

### [45] Woodbury 1950
- **Metadata:** MINOR — the report exists: M. A. Woodbury, Inverting modified matrices, Memorandum Report 42, Statistical Research Group, Princeton University, 1950 (4 pp.). The author's first name is Max A. Woodbury (initials fine). "Tech. Rep." is the entry-type label only. No DOI exists.
- **Source checked:** https://catalog.princeton.edu/catalog/4461270 ; https://en.wikipedia.org/wiki/Woodbury_matrix_identity (WebSearch)
- **Citations:**
  - p.3: SUPPORTED — the original source of the Woodbury low-rank-update inverse identity used in an SMW preconditioner.
- **Recommended fix:** none (optionally "Memorandum Report 42" without "Tech. Rep." duplication).

### [46] Sherman & Morrison 1950
- **Metadata:** MINOR — authors, title, Ann. Math. Stat. 21(1), 124–127, 1950 all match. The DOI is missing: 10.1214/aoms/1177729893.
- **Source checked:** https://api.crossref.org/works?query.bibliographic=Sherman%20Morrison%20Adjustment%20of%20an%20inverse%20matrix...
- **Citations:**
  - p.3: SUPPORTED — the original rank-one inverse-update formula.
- **Recommended fix:** Add doi:10.1214/aoms/1177729893.

### [47] Golub & Van Loan 2013
- **Metadata:** OK — Matrix Computations, 4th ed., Johns Hopkins University Press, Baltimore, 2013 confirmed. Optional additions: ISBN 978-1-4214-0794-4, or SIAM-hosted DOI 10.1137/1.9781421407944.
- **Source checked:** https://www.press.jhu.edu/books/title/10678/matrix-computations ; https://epubs.siam.org/doi/book/10.1137/1.9781421407944
- **Citations:**
  - p.3: SUPPORTED — the standard reference for SMW updates and banded/bordered matrix algorithms.
- **Recommended fix:** none

### [48] Dudarev 2008
- **Metadata:** MINOR — authors, title, PRL 100(13) 135503, 2008 and DOI all match. The only problem is the DOI line break "PhysRevLett. 100.135503", which leaves a space inside the DOI.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevLett.100.135503 ; full text https://scientific-publications.ukaea.uk/wp-content/uploads/Published/PRLVOL100p135503.pdf
- **Citations:**
  - p.3: MISMATCH (value unverified) — the paper gives no explicit critical size of "about 23 interstitials at 330 °C". It describes temperature regions: a/2<111> favoured below about 350 °C, <100> favoured above about 550 °C, and an intermediate region. A critical size can be derived from its energy expression, but that number comes from the authors' own evaluation, not from [48].
  - p.30: SUPPORTED — the paper gives the prismatic-loop free energy as perimeter × (prelogarithmic ln term plus core-traction and anharmonic core terms).
  - supp p.88: SUPPORTED — the functional form matches. Check the signs: the full text writes Ftot = P[F̂ ln(4R/eρ) − FΔ − Fc], with the core terms subtracted, while Eq. S61 appears to add them (the extracted text may be garbled).
  - supp p.94 (δ, Tc): SUPPORTED — the paper gives core cutoff ρ = 0.4 nm (4.0 Å) and Tc ≈ 912 °C ≈ 1185 K.
  - supp p.94 (core terms): SUPPORTED — the paper tabulates core-traction energies 0.345 (<111>) and 0.387 (<100>) eV/Å and anharmonic core terms 0.46 (<111>) and 0.33 (<100>). The units of Fc should be checked against the paper, which lists them as best-fit values.
- **Recommended fix:** Fix the DOI space. On p.3, either state that the ~23-SIA / 330 °C critical size is derived from the [48] energetics ("evaluated from the formation energy of [48]"), or remove the number. Check the sign convention of the core terms in Eq. S61 against [48].

### [49] Marian 2002
- **Metadata:** MINOR — authors, title, PRL 88(25) 255507, 2002 and DOI all match. The only problem is the DOI line break "PhysRevLett.88. 255507", which leaves a space inside the DOI.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevLett.88.255507 ; WebSearch (APS abstract page and OSTI record listed but blocked from fetch)
- **Citations:**
  - p.4: PLAUSIBLE — [49] is the source of the junction reaction between comparable 1/2<111> loops. "Thermally activated rotation of single loops" is not from [49]; that is the Eyre–Bullough / later rotation literature. As placed, [49] seems to cover both mechanisms.
  - p.31: SUPPORTED — the MD mechanism in [49] is two similar-size 1/2<111> loops reacting through a 1/2<110> intermediate to form a <100> loop (from my knowledge of the paper; the abstract fetch was blocked).
  - supp p.91: SUPPORTED — the Burgers-vector reaction 1/2[111] + 1/2[11̄1̄] → [100] is the reaction in [49].
- **Recommended fix:** Fix the DOI space. On p.4, add a separate citation for loop rotation (e.g., Eyre & Bullough, Philos. Mag. 12 (1965) 31) or rephrase so that [49] covers only the junction reaction.

### [50] Ghoniem 2026 (RadCluster software)
- **Metadata:** OK — the GitHub repository Ghoniem/RadCluster exists and is public. Its description, "Cluster dynamics modeling of Eurofer microstructure", matches the entry. There is no DOI; consider archiving a release on Zenodo for a citable DOI and version.
- **Source checked:** https://github.com/Ghoniem/RadCluster
- **Citations:**
  - p.4: SUPPORTED — this is the code the paper describes.
  - p.6: SUPPORTED — self-reference to the code-variable mapping.
  - p.25: SUPPORTED — self-reference to the code's numerical scheme (the scheme itself was not audited).
- **Recommended fix:** Add a version/commit or Zenodo DOI and access date.

### [51] Horn & Jackson 1972
- **Metadata:** OK — F. Horn, R. Jackson, General mass action kinetics, Arch. Ration. Mech. Anal. 47(2) 81–116, 1972, DOI 10.1007/BF00251225 all match.
- **Source checked:** https://api.crossref.org/works/10.1007/BF00251225
- **Citations:**
  - p.6: SUPPORTED — the founding paper of the mass-action CRN formalism and complex balancing.
  - supp p.66: SUPPORTED — the stoichiometric form ẋ = S J(x) under mass action.
  - supp p.68: SUPPORTED — the complex-balancing theorems of [51] assume finite networks with constant rate constants.
- **Recommended fix:** none

### [52] Feinberg 1987
- **Metadata:** OK — the title, Chem. Eng. Sci. 42(10) 2229–2268, 1987 and DOI 10.1016/0009-2509(87)80099-4 all match.
- **Source checked:** https://api.crossref.org/works/10.1016/0009-2509(87)80099-4
- **Citations:**
  - p.6: SUPPORTED — the deficiency-zero and deficiency-one theorems.
  - supp p.66: SUPPORTED.
  - supp p.68: SUPPORTED — deficiency theorems for finite mass-action networks.
- **Recommended fix:** none

### [53] Feinberg 2019
- **Metadata:** OK — Foundations of Chemical Reaction Network Theory, Applied Mathematical Sciences, Springer, Cham, 2019, DOI 10.1007/978-3-030-03858-8 all match. Crossref did not return the volume number; vol. 202 matches the Springer series listing as I know it.
- **Source checked:** https://api.crossref.org/works/10.1007/978-3-030-03858-8
- **Citations:**
  - p.6: SUPPORTED — a comprehensive monograph on deficiency theory.
  - supp p.66: SUPPORTED.
- **Recommended fix:** Remove the space in "doi: 10.1007/..." (formatting).

### [54] Temkin 1996
- **Metadata:** OK — O. N. Temkin, A. V. Zeigarnik, D. Bonchev, Chemical Reaction Networks: A Graph-Theoretical Approach, CRC Press, 1996 (ISBN 0849328675) confirmed. There is no DOI; a 2020 Routledge/CRC reissue also exists.
- **Source checked:** https://www.routledge.com/Chemical-Reaction-Networks-A-Graph-Theoretical-Approach/Temkin-Zeigarnik-Bonchev/p/book/9780367448479 ; Amazon listing (WebSearch)
- **Citations:**
  - p.7: SUPPORTED — the book reviews graph-theoretical treatments of reaction mechanisms.
- **Recommended fix:** none (optional ISBN).

### [55] Klamt 2009
- **Metadata:** OK — S. Klamt, U.-U. Haus, F. Theis, Hypergraphs and cellular networks, PLoS Comput. Biol. 5(5) e1000385, 2009 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1371/journal.pcbi.1000385
- **Citations:**
  - p.7: SUPPORTED — the paper argues that reactions with multiple reactants and products need directed hypergraphs, not simple graphs.
  - supp p.66: SUPPORTED — the same concept: directed hyperedges with two tails.
- **Recommended fix:** none

### [56] Goodwin et al. Cantera 3.2.0 (2025)
- **Metadata:** OK — the Zenodo record 10.5281/zenodo.17620923 is Cantera 3.2.0 (published 17 Nov 2025). The authors Goodwin, Moffat, Schoegl, Speth, Weber are in the correct order and the title matches.
- **Source checked:** https://zenodo.org/records/17620923 (the Crossref API returned 429, so Zenodo was used instead)
- **Citations:**
  - p.7: SUPPORTED — Cantera separates mechanism declaration (YAML input) from the solvers.
- **Recommended fix:** none

### [57] Damian 2002
- **Metadata:** OK — authors in order, title, Comput. Chem. Eng. 26(11) 1567–1579, 2002 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1016/S0098-1354(02)00128-X
- **Citations:**
  - p.7: SUPPORTED — KPP generates solver code from a declarative mechanism description.
- **Recommended fix:** none

### [58] Hucka 2003
- **Metadata:** OK — title, Bioinformatics 19(4) 524–531, 2003 and DOI match. The first eight authors are in the correct order and "et al." covers the remaining 36 of 44 authors, which is acceptable.
- **Source checked:** https://api.crossref.org/works/10.1093/bioinformatics/btg015
- **Citations:**
  - p.7: SUPPORTED — SBML is a declarative, simulator-independent format for network models.
- **Recommended fix:** none

### [59] Loman 2023
- **Metadata:** OK — the eight authors in order, title, PLOS Comput. Biol. 19(10) e1011530, 2023 and DOI all match.
- **Source checked:** https://api.crossref.org/works/10.1371/journal.pcbi.1011530
- **Citations:**
  - p.7: SUPPORTED — Catalyst.jl is a symbolic reaction-network DSL that is kept separate from the numerical solvers.
- **Recommended fix:** none

### [60] Gao 2016
- **Metadata:** OK — C. W. Gao, J. W. Allen, W. H. Green, R. H. West, the title, Comput. Phys. Commun. 203, 212–225, 2016 and DOI 10.1016/j.cpc.2016.02.013 all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.cpc.2016.02.013
- **Citations:**
  - p.7: SUPPORTED — RMG builds mechanisms automatically from reaction-family templates and rate rules.
  - supp p.68: SUPPORTED — RMG organizes mechanisms into reaction families with rate rules. It was built mainly for combustion and pyrolysis but is used more widely; "combustion mechanisms" is a fair description.
- **Recommended fix:** none. The text also cites [61] alongside RMG, which is outside this batch.
### [61] Liu 2021
- **Metadata:** OK. Authors (14, same order), title, J. Chem. Inf. Model. 61(6) 2686–2696 (2021) and DOI all match. "gen- eration" is only a PDF line-break hyphen.
- **Source checked:** https://api.crossref.org/works/10.1021/acs.jcim.0c01480
- **Citations:**
  - p.7: SUPPORTED. RMG v3.0 builds kinetic networks automatically from reaction-family templates.
- **Recommended fix:** none

### [62] Martín-Bragado 2013
- **Metadata:** OK. Authors (Martin-Bragado, Rivera, Valles, Gomez-Selles, Caturla), title, Comput. Phys. Commun. 184(12) 2703–2710 (2013) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.cpc.2013.07.011
- **Citations:**
  - p.7: SUPPORTED. MMonCa is an object kinetic Monte Carlo code for irradiation damage, where reaction events come from general object and interaction rules. That fits "rule-based generation of the reaction set" in a loose sense.
- **Recommended fix:** none

### [63] Blondel 2017
- **Metadata:** OK. Authors (Blondel, Bernholdt, Hammond, Hu, Maroudas, Wirth), title, Fusion Sci. Technol. 71(1) 84–92 (2017) and DOI match.
- **Source checked:** https://api.crossref.org/works/10.13182/FST16-109
- **Citations:**
  - p.7: SUPPORTED. This paper benchmarks the Xolotl spatially resolved (multidimensional) cluster dynamics code for He in W.
- **Recommended fix:** none

### [64] van Kampen 2007
- **Metadata:** OK. N.G. van Kampen, Stochastic Processes in Physics and Chemistry, 3rd ed., North-Holland (Elsevier), 2007. The DOI 10.1016/B978-0-444-52965-7.X5000-4 is confirmed on ScienceDirect. The space inside the DOI ("doi:10.1016/ B978…") is a line-break artefact; check the .bib file has none.
- **Source checked:** https://www.sciencedirect.com/book/9780444529657/stochastic-processes-in-physics-and-chemistry (the Crossref fetch was rate-limited); WebSearch (publisher and retailer listings)
- **Citations:**
  - p.13: SUPPORTED. The standard reference for the chemical master equation, with propensities and stoichiometric state changes.
- **Recommended fix:** none

### [65] Gillespie 1976
- **Metadata:** OK. Daniel T. Gillespie, J. Comput. Phys. 22(4) 403–434 (1976), DOI 10.1016/0021-9991(76)90041-3; title matches. "D. T. Gillespie" is correct.
- **Source checked:** https://www.sciencedirect.com/science/article/abs/pii/0021999176900413 ; https://www.mindat.org/reference.php?id=15876812 (the Crossref fetch was rate-limited)
- **Citations:**
  - p.13: SUPPORTED. This is the original stochastic simulation algorithm, an exact method that samples the chemical master equation.
  - supp p.68: SUPPORTED. Same point: the SSA is the stochastic counterpart used by stochastic cluster dynamics.
- **Recommended fix:** none

### [66] Berg 1977
- **Metadata:** OK. H.C. Berg, E.M. Purcell, "Physics of chemoreception", Biophys. J. 20(2) 193–219 (1977); DOI matches.
- **Source checked:** https://api.crossref.org/works/10.1016/S0006-3495(77)85544-6
- **Citations:**
  - p.17: SUPPORTED. Berg and Purcell write the diffusive current to an absorber of any shape as 4πDCc∞, where C is the electrostatic capacitance.
  - supp p.76: SUPPORTED. Same capacitance analogy for the far field and the current.
  - supp p.77: SUPPORTED. The paper uses the flat-disc capacitance 2a/π (for a disc-shaped absorbing patch). Note that applying it to a loop is the manuscript's own analogy.
- **Recommended fix:** none

### [67] Seeger 1977
- **Metadata:** OK. A. Seeger, U. Gösele, "Steady-state diffusion of point defects to dislocation loops", Phys. Lett. A 61(6) 423–425 (1977); DOI matches.
- **Source checked:** https://api.crossref.org/works/10.1016/0375-9601(77)90355-3
- **Citations:**
  - p.17: PLAUSIBLE. The toroidal-geometry diffusion solution for loops supports "uses the loop radius". The part about elastic interaction and a capture efficiency above unity only for SIAs rests mainly on [68, 69]. From the title and abstract level I cannot confirm that this 3-page letter treats elastic bias.
  - supp p.77: SUPPORTED, but the exact expression is unverified. The torus result C ≈ πr/ln(8r/rc) is the classical toroidal absorber form this paper is known for. The prefactor and the (0.5–0.7) range could not be checked against the full text.
- **Recommended fix:** none. Optionally, cite [67] only for geometry and [68, 69] for the bias.

### [68] Brailsford 1981
- **Metadata:** OK. A.D. Brailsford, R. Bullough, "The theory of sink strengths", Phil. Trans. R. Soc. Lond. A 302(1465) 87–137 (1981); DOI matches.
- **Source checked:** https://api.crossref.org/works/10.1098/rsta.1981.0158
- **Citations:**
  - p.17: SUPPORTED. The canonical treatment of sink strengths, including dislocation and loop bias (capture efficiencies).
  - supp p.77: SUPPORTED. It covers sink efficiency with geometric and elastic-bias parts.
- **Recommended fix:** none

### [69] Jourdan 2015
- **Metadata:** OK. T. Jourdan (sole author), J. Nucl. Mater. 467 (2015) 286–301; title and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2015.09.046 ; https://api.openalex.org/works/doi:10.1016/j.jnucmat.2015.09.046
- **Citations:**
  - p.17: SUPPORTED. The paper is about dislocation and loop bias factors (size-dependent capture efficiencies) in rate-equation cluster dynamics.
  - supp p.77: SUPPORTED. Same topic: loop capture efficiency including elastic bias.
- **Recommended fix:** none

### [70] Borodin 1998
- **Metadata:** OK. V.A. Borodin, "Rate theory for one-dimensional diffusion", Physica A 260(3–4) 467–478 (1998); DOI matches.
- **Source checked:** https://api.crossref.org/works/10.1016/S0378-4371(98)00338-0 ; https://api.openalex.org/works/doi:10.1016/S0378-4371(98)00338-0 (no abstract in either record)
- **Citations:**
  - p.17: SUPPORTED. Borodin's 1D rate theory gives a sink strength for purely 1D walkers that is quadratic in sink density. This is the known result of the paper; the abstract was not available to confirm.
  - supp p.80: PLAUSIBLE. The pure-1D limit k² = 6κ_m κ with κ = σN matches the Barashev/Trinkaus notation more closely than Borodin's. The exact form could not be confirmed in [70].
  - supp p.83: SUPPORTED. An analytical treatment of 1D reaction kinetics.
  - supp p.84: PLAUSIBLE. The "marked interception" result and the cubic loss frequency derived from it could not be confirmed from the available metadata.
- **Recommended fix:** none

### [71] Barashev 2001
- **Metadata:** OK. A.V. Barashev, S.I. Golubov, H. Trinkaus, Phil. Mag. A 81(10) 2515–2532 (2001); title and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1080/01418610108217161
- **Citations:**
  - p.17: SUPPORTED. The paper derives reaction kinetics of 1D-gliding SIA clusters with voids and dislocations. For the pure 1D case the sink strength depends quadratically on sink density.
  - supp p.80: SUPPORTED. The κ = σN form of the 1D sink strength comes from this line of work.
  - supp p.83: SUPPORTED. An analytical treatment.
  - supp p.84: PLAUSIBLE. The specific "marked interception" argument could not be confirmed.
- **Recommended fix:** none

### [72] Trinkaus 2002
- **Metadata:** OK. H. Trinkaus, H.L. Heinisch, A.V. Barashev, S.I. Golubov, B.N. Singh, "1D to 3D diffusion-reaction kinetics of defects in crystals", Phys. Rev. B 66(6) 060105 (2002); DOI matches. "(R)" (Rapid Communication) is correct.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.66.060105
- **Citations:**
  - p.17: SUPPORTED. This is the 1D-to-3D master-curve interpolation by exactly these five authors.
  - supp p.81: SUPPORTED. Same master curve.
  - supp p.83: SUPPORTED. An analytical treatment of mixed 1D/3D kinetics.
  - supp p.85: SUPPORTED. The master curve was checked against kMC of gliders with occasional direction changes (Heinisch et al.).
- **Recommended fix:** none

### [73] Singh 2008 (Risø-R-1644)
- **Metadata:** NOT FOUND. No record of a Risø report "Risø-R-1644(EN)" or of this title by Singh, Golubov and Trinkaus. Searched WebSearch (exact report number and title phrases), OSTI/ETDEWEB results, DTU Orbit (the person page returned 404) and WorldWideScience. The only related items found are journal papers: Heinisch, Singh, Golubov, "The effects of one-dimensional glide on the reaction kinetics of interstitial clusters" (J. Nucl. Mater. 283–287, 2000), and Trinkaus, Singh, Golubov, "Progress in modelling the microstructural evolution in metals under cascade damage conditions" (J. Nucl. Mater. 283–287, 2000). The entry may be fabricated or garbled.
- **Source checked:** WebSearch ("Risø-R-1644", title phrases); https://www.osti.gov/etdeweb/biblio/20356908 ; https://www.osti.gov/etdeweb/biblio/20356781 ; https://worldwidescience.org/topicpages/d/diffusion+reaction+kinetics.html
- **Citations:**
  - p.17: PLAUSIBLE, source not verified. The master curve is already covered by [72].
  - supp p.80: PLAUSIBLE, source not verified. Pure-1D limit.
  - supp p.81: PLAUSIBLE, source not verified.
  - supp p.85: PLAUSIBLE, source not verified.
- **Recommended fix:** Remove the entry unless the author has a physical copy with a confirmed report number. Otherwise cite [72] alone, or add a verified source such as H. Trinkaus, B.N. Singh, S.I. Golubov, J. Nucl. Mater. 283–287 (2000) 89–98 (DOI 10.1016/S0022-3115(00)00225-7; verify volume, pages and DOI before use) or Heinisch, Singh, Golubov, J. Nucl. Mater. 283–287 (2000) 737.

### [74] Woo 1990
- **Metadata:** OK. C.H. Woo, B.N. Singh, Phys. Status Solidi (b) 159(2) 609–616 (1990); title and DOI match.
- **Source checked:** https://api.openalex.org/works/doi:10.1002/pssb.2221590210 (Crossref and Wiley were rate-limited or blocked); WebSearch hit on the Wiley DOI page
- **Citations:**
  - p.23: SUPPORTED. This paper introduced production bias, where cascade-produced SIA clusters supply extra interstitials.
- **Recommended fix:** none

### [75] Woo 1992
- **Metadata:** OK. C.H. Woo, B.N. Singh, Phil. Mag. A 65(4) 889–912 (1992); title and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1080/01418619208205596
- **Citations:**
  - p.23: SUPPORTED. Production bias from clustering of cascade point defects, including glissile SIA clusters.
- **Recommended fix:** none

### [76] Singh 1992
- **Metadata:** OK. B.N. Singh, A.J.E. Foreman, Phil. Mag. A 66(6) 975–990 (1992); title and DOI match.
- **Source checked:** https://api.crossref.org/works/10.1080/01418619208248002
- **Citations:**
  - p.23: SUPPORTED. Production bias and transient void swelling driven by the interstitial-cluster supply.
- **Recommended fix:** none

### [77] Golubov 2000
- **Metadata:** MINOR. Authors, title, J. Nucl. Mater. 276, 78–89 (2000) and DOI match. The issue (1–3) is missing, and the capital "T" in "Towards" differs slightly from the record (style only).
- **Source checked:** https://api.crossref.org/works/10.1016/S0022-3115(99)00171-3
- **Citations:**
  - p.23: SUPPORTED. The generalized production bias model with 1D-gliding SIA clusters in fcc and bcc metals.
  - p.25: MISMATCH. The grouping method that carries a zeroth moment and a slope (first moment) per group, removing mass leakage, is from S.I. Golubov, A.M. Ovcharenko, A.V. Barashev, B.N. Singh, "Grouping method for the approximate solution of a kinetic equation describing the evolution of point-defect clusters", Phil. Mag. A 81 (2001) 643–658. It is not from the 2000 production-bias paper.
- **Recommended fix:** Add issue 1–3. For the p.25 grouping-method sentence, replace [77] with Golubov et al., Phil. Mag. A 81 (2001) 643 (DOI 10.1080/01418610108214306; verify), and/or Jourdan, Bocquet, Soisson, Comput. Mater. Sci. 49 (2010) 2 if that is the intended source.

### [78] Wolfer 2007
- **Metadata:** OK. W.G. Wolfer, "The dislocation bias", J. Comput.-Aided Mater. Des. 14(3) 403–417 (2007); DOI matches.
- **Source checked:** https://link.springer.com/article/10.1007/s10820-007-9051-3 (Crossref was rate-limited)
- **Citations:**
  - p.23: MISMATCH. The paper derives dislocation bias factors for the preferential absorption of SIAs over vacancies. It says nothing about glissile SIA clusters as an extra interstitial supply. It could be cited for why dislocation bias alone is not enough, but not for the glissile-cluster claim.
- **Recommended fix:** Move [78] to a sentence about dislocation bias, or reword to "beyond the dislocation bias [78, 79], most plausibly glissile SIA clusters [74–77]".

### [79] Kohnert 2019
- **Metadata:** OK. A.A. Kohnert, L. Capolungo, Phys. Rev. Materials 3(5) 053608 (2019); title and DOI match.
- **Source checked:** https://api.openalex.org/works/doi:10.1103/PhysRevMaterials.3.053608 ; WebSearch hit on link.aps.org (the APS page returned 403)
- **Citations:**
  - p.23: MISMATCH. The paper computes sink strength and dislocation bias of 3D dislocation-dynamics networks. It does not concern glissile SIA clusters or production bias.
- **Recommended fix:** Same as [78]: cite it for the dislocation bias, not for the glissile-cluster supply.

### [80] Fu 2005
- **Metadata:** OK. C.-C. Fu, J. Dalla Torre, F. Willaime, J.-L. Bocquet, A. Barbu, Nature Mater. 4 (2005) 68–74; title and DOI match. The issue (1) is missing (optional).
- **Source checked:** https://www.nature.com/articles/nmat1286 ; WebSearch (ADS record 2005NatMa...4...68F) (Crossref was rate-limited)
- **Citations:**
  - p.25: MISMATCH. The paper (ab initio plus kMC of vacancy and SIA clusters in α-Fe) gives no helium or hydrogen migration energy. Em ≈ 0.06 eV for interstitial He in α-Fe is from C.-C. Fu, F. Willaime, "Ab initio study of helium in α-Fe: dissolution, migration, and clustering with vacancies", Phys. Rev. B 72 (2005) 064117. If "h" means hydrogen, cite a DFT H-in-Fe study such as Jiang & Carter, Phys. Rev. B 70 (2004) 064102. Also check "10^6 times faster than vacancies at 400 °C": with E_v^m ≈ 0.67 eV (Fu et al.) the ratio of Boltzmann factors at 673 K is about exp(0.61/0.058) ≈ 4×10^4, not 10^6.
  - supp p.67 (Table S3): PLAUSIBLE, values unverified. Fu et al. 2005 is the standard source for vacancy and SIA (cluster) migration and binding energies in Fe (e.g. E_v^m ≈ 0.67 eV, SIA ≈ 0.34 eV). The table values attributed to it were not checked one by one.
- **Recommended fix:** For the p.25 Em ≈ 0.06 eV statement, replace [80] with Fu & Willaime, Phys. Rev. B 72 (2005) 064117 (DOI 10.1103/PhysRevB.72.064117; verify), and recheck the "10^6" factor. Keep [80] for the Table S3 point-defect energies.
### [81] Koiwa 1974
- **Metadata:** OK. Crossref matches: M. Koiwa, J. Phys. Soc. Jpn. 37(6) (1974) 1532–1536.
- **Source checked:** https://api.crossref.org/works/10.1143/JPSJ.37.1532
- **Citations:**
  - p.25: SUPPORTED. The paper is a critique of the grouping method for vacancy-clustering size distributions, so it belongs in the coarse-graining lineage.
- **Recommended fix:** none

### [82] Ghoniem 1988
- **Metadata:** OK. Crossref matches: N.M. Ghoniem, J. Nucl. Mater. 155–157 (1988) 1123–1127. Title uses "Fokker-Planck".
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3115(88)90480-1
- **Citations:**
  - p.25: PLAUSIBLE. A moments solution of the Fokker–Planck equation reduces the size distribution to a few moments, which is a form of coarse-graining. Its main subject, though, is determining the bias factor, not grouping.
- **Recommended fix:** none

### [83] Huang & Ghoniem 1995
- **Metadata:** OK. Crossref matches: H. Huang, N.M. Ghoniem, Phys. Rev. E 51(6) (1995) 5251–5260.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevE.51.5251
- **Citations:**
  - p.25: SUPPORTED. The moment method for multidimensional Fokker–Planck equations is a reduced description of the size distribution.
- **Recommended fix:** none

### [84] Jourdan 2014
- **Metadata:** OK. Crossref matches: T. Jourdan, G. Bencteux, G. Adjanor, J. Nucl. Mater. 444(1–3) (2014) 298–313. There is a stray space in the DOI string ("jnucmat. 2013"), which is a line-break artifact.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2013.10.009
- **Citations:**
  - p.25: SUPPORTED. Known content: the paper couples discrete rate equations for small clusters to a Fokker–Planck description for large sizes, and adds grouping. Crossref gave no abstract.
- **Recommended fix:** none

### [85] Kohnert & Wirth 2017
- **Metadata:** OK. Crossref matches: A.A. Kohnert, B.D. Wirth, Modelling Simul. Mater. Sci. Eng. 25(1) 015008. Crossref lists 2016 as the online date; the issue is dated 2017, which is acceptable.
- **Source checked:** https://api.crossref.org/works/10.1088/1361-651X/25/1/015008
- **Citations:**
  - p.25: SUPPORTED. The abstract says it "advances a framework for grouping arbitrary cluster dynamics problems".
- **Recommended fix:** none

### [86] Terrier 2017
- **Metadata:** OK. Crossref matches: P. Terrier, M. Athènes, T. Jourdan, G. Adjanor, G. Stoltz, J. Comput. Phys. 350 (2017) 280–295.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jcp.2017.08.015
- **Citations:**
  - p.25: SUPPORTED. The paper is a hybrid deterministic/stochastic coupling approach, with the stochastic part covering the large-size tail.
- **Recommended fix:** none

### [87] Dethloff 2018
- **Metadata:** OK. Crossref matches: C. Dethloff, E. Gaganidze, J. Aktaa, Nucl. Mater. Energy 15 (2018) 23–26.
- **Source checked:** https://api.crossref.org/works/10.1016/j.nme.2018.05.015 ; https://doaj.org/article/1b80e08a57304f2f9396c1b8ede68e52
- **Citations:**
  - p.28 (analysis cutoff and diffraction conditions): SUPPORTED. The review "identifies significant inconsistencies in earlier investigations" of loop size distributions and recommends standardizing TEM procedures.
  - p.33 (most density below the visibility limit): PLAUSIBLE. This fits a critical review of loop size distributions, but the abstract does not state it.
  - supp p.111 (ranges of loop densities and diameters): PLAUSIBLE. These are compiled values and cannot be checked from the abstract.
  - supp p.118 (⟨100⟩ fraction of 27–78 %): PLAUSIBLE. The spread across analyses is the review's theme, but the numbers are unverified.
  - supp p.118–119 Table S20 rows (300/330 °C, 15/32 dpa; f⟨100⟩ = 77, 72, 45, 27, 31 %; densities and diameters): PLAUSIBLE, numbers unverified. The 300 °C/15 dpa row gives 1.40/4.70 ×10^21 m^-3, while the primary paper [92] gives 1.4/4.9 ×10^21. Check which analysis each row comes from.
- **Recommended fix:** none to the metadata. Check the table values against the review's own tables.

### [88] "Gao 2018" (actually Jiao 2018)
- **Metadata:** MAJOR. The DOI 10.1016/j.actamat.2018.01.016 resolves to a different paper: Bi et al., "Novel deformation-induced polymorphic crystallization and softening of Al-based amorphous alloys", Acta Mater. 147 (2018) 90–99. The title and the bib note (T91, BOR-60, 376–524 °C, 15.4/35.1 dpa, bimodal cavities) belong to Z. Jiao, S. Taller, K. Field, G. Yeli, M.P. Moody, G.S. Was, J. Nucl. Mater. 504 (2018) 122–134, doi:10.1016/j.jnucmat.2018.03.024. The author list (Gao, Liu, Schreiber, Toloczko, Lear, Cole, Was), the journal (Acta Materialia) and the volume/pages (142, 240–254) are all wrong. Bib-note text has also leaked into the entry.
- **Source checked:** https://api.crossref.org/works/10.1016/j.actamat.2018.01.016 ; https://api.crossref.org/works?query.bibliographic=Microstructure+evolution+of+T91+irradiated+in+the+BOR60+fast+reactor ; https://api.crossref.org/works/10.1016/j.jnucmat.2018.03.024 ; https://www.nomelab.com/publication/2018-01-01_jiao2018microstructure/
- **Citations:**
  - p.29 (dominant Burgers vector switches from a/2⟨111⟩ to a⟨100⟩ as dose increases): PLAUSIBLE, but weak. The abstract reports a⟨100⟩ loops at 376–415 °C and network dislocations at 460–524 °C, so the change it describes is with temperature. It does not report a dose-driven switch from a/2⟨111⟩, and it is a neutron study, not an ion study.
- **Recommended fix:** Replace with: Z. Jiao, S. Taller, K. Field, G. Yeli, M.P. Moody, G.S. Was, Microstructure evolution of T91 irradiated in the BOR60 fast reactor, J. Nucl. Mater. 504 (2018) 122–134, doi:10.1016/j.jnucmat.2018.03.024. Delete the bib note. Soften the sentence so it does not claim a within-study dose switch.

### [89] "Gao 2019" (actually Zheng 2019)
- **Metadata:** MAJOR. The DOI 10.1016/j.jnucmat.2019.03.026 resolves to a different paper: Shin et al., "Factors affecting the hydrothermal corrosion behavior of chemically vapor deposited silicon carbides", J. Nucl. Mater. 518 (2019) 350–356. The title and bib note are from C. Zheng, E.R. Reese, K.G. Field, E. Marquis, S.A. Maloy, D. Kaoumi, J. Nucl. Mater. 523 (2019) 421–433, doi:10.1016/j.jnucmat.2019.06.019. The authors (Gao, Liu, Schreiber, Lear, Was), volume (516) and pages (312–322) are wrong. The note says "377 °C"; the record says 650 ± 23 K (about 377 °C), which is fine. Bib-note text has leaked into the entry.
- **Source checked:** https://doi.org/10.1016/j.jnucmat.2019.03.026 (via sciencedirect pii S0022311518315587) ; https://www.osti.gov/pages/biblio/1529011 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311519302156
- **Citations:**
  - p.29 (Burgers vector switch from a/2⟨111⟩ to a⟨100⟩ with dose): MISMATCH. The abstract says loops were "predominantly in the a⟨100⟩ orientation at both doses" (17.1 and 35.1 dpa), so no switch with dose is reported.
- **Recommended fix:** Replace with: C. Zheng, E.R. Reese, K.G. Field, E. Marquis, S.A. Maloy, D. Kaoumi, Microstructure response of ferritic/martensitic steel HT9 after neutron irradiation: effect of dose, J. Nucl. Mater. 523 (2019) 421–433, doi:10.1016/j.jnucmat.2019.06.019. Delete the bib note. Either cite it only for a⟨100⟩ dominance at higher dose or rephrase the claim.

### [90] Zhong & Tan 2021 (conflated entry)
- **Metadata:** MAJOR. Three different records are mixed together:
  - The DOI 10.1016/j.jnucmat.2020.152634 resolves to B.K. Kim et al., "Effects of helium on irradiation response of reduced-activation ferritic-martensitic steels...", J. Nucl. Mater. 545 (2021) 152634.
  - The cited volume and article number (551, 152971) belong to Cao & Zhou, an ODS steel paper.
  - The title "The microstructure effects on irradiation response of ferritic–martensitic steels" is W. Zhong, L. Tan, J. Nucl. Mater. 593 (2024) 154990, doi:10.1016/j.jnucmat.2024.154990. That paper covers 9Cr-NbMo and 9Cr-Ta in HFIR at 400/490 °C and 14.7 dpa, not T91 in ATR.
  - The bib note (T91/NF616, ATR, up to 8.2 dpa, 292–431 °C, loops in all samples, Ni-rich clusters) matches W. Zhong, T.A. Saleh, L. Tan, "Neutron irradiation induced defects and clustering in NF616 and T91", J. Nucl. Mater. 552 (2021) 153001, doi:10.1016/j.jnucmat.2021.153001.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2020.152634 ; https://doi.org/10.1016/j.jnucmat.2021.152971 (sciencedirect pii S002231152100194X) ; https://api.crossref.org/works/10.1016/j.jnucmat.2024.154990 ; https://www.sciencedirect.com/science/article/abs/pii/S002231152400093X ; https://www.sciencedirect.com/science/article/abs/pii/S0022311521002245
- **Citations:**
  - p.29 (Burgers vector switch to a⟨100⟩ with dose): MISMATCH.
    - The intended NF616/T91 ATR paper reports 1/2⟨111⟩ loops as dominant up to 8.2 dpa. That fits only the low-dose side of the claim; no switch is shown.
    - The 2024 Zhong & Tan title is about microstructure optimization, softening and cavities, not loop character.
- **Recommended fix:** Replace with: W. Zhong, T.A. Saleh, L. Tan, Neutron irradiation induced defects and clustering in NF616 and T91, J. Nucl. Mater. 552 (2021) 153001, doi:10.1016/j.jnucmat.2021.153001. Delete the bib note. Rephrase the claim, since [88–90] together show a⟨100⟩ at higher dose and temperature and 1/2⟨111⟩ at low dose only across different studies.

### [91] "Klimenkov 2012" (conflated with Weiß 2012 = [99])
- **Metadata:** MAJOR.
  - The DOI 10.1016/j.jnucmat.2012.01.022 resolves to B.D. Miller et al., "Advantages and disadvantages of using a focused ion beam to prepare TEM samples from irradiated U–10Mo monolithic nuclear fuel", J. Nucl. Mater. 424 (2012) 38–42.
  - The title "Quantitative characterization of microstructural defects in up to 32 dpa neutron irradiated EUROFER97" belongs to O.J. Weiß, E. Gaganidze, J. Aktaa, J. Nucl. Mater. 426 (2012) 52–58, doi:10.1016/j.jnucmat.2012.03.027. That paper is already ref [99].
  - The bib-note data (WBDF TEM, 3.4 and 4.8 nm at 15/32 dpa, 330–335 °C, voids) are also from Weiß et al.
  - No Klimenkov/Materna-Morris/Möslang paper with this title exists. Their nearest paper, "Characterization of radiation induced defects in EUROFER 97 after neutron irradiation", J. Nucl. Mater. 417 (2011) 124–126, doi:10.1016/j.jnucmat.2010.12.261, covers 16.3 dpa at 250–450 °C in HFR, not BOR-60 or 32 dpa.
  - "J. Nucl. Mater. 428 (2012) 48–54" does not match either paper. The entry is effectively a fabricated duplicate of [99].
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2012.01.022 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311512001481 ; https://api.crossref.org/works/10.1016/j.jnucmat.2010.12.261 ; https://inis.iaea.org/records/t3443-3gt31
- **Citations:**
  - p.29 (swelling below 0.01–0.1 % at 32 dpa, 330 °C): PLAUSIBLE only if the source is Weiß et al. The Weiß abstract reports void densities at least an order of magnitude below loop densities, but no swelling number. As attributed to Klimenkov, it is a MISMATCH.
  - supp p.111 (compiled ranges): PLAUSIBLE, unverified, and depends on the conflated source.
  - supp p.118 (1.4→1.7 ×10^22 m^-3, 3.4→4.8 nm): SUPPORTED by Weiß et al. [99]. It is not a separate Klimenkov source, so this is a duplicate citation.
  - supp p.118–119 Table S20 (330 °C 15/32 dpa "all" rows): SUPPORTED via Weiß et al. The duplicate attribution should go.
  - supp p.120 (cavity diameter 1.6–2.6 nm, density 1.1–1.4 ×10^21 m^-3, swelling <0.01 %, BOR-60): PLAUSIBLE, numbers unverified. They are probably from Weiß et al. The abstract confirms void analysis but not these values.
  - supp p.120 Table S21 rows (330/335 °C BOR-60): PLAUSIBLE, same caveat.
- **Recommended fix:** Delete [91] and cite [99] (Weiß et al. 2012, doi:10.1016/j.jnucmat.2012.03.027) everywhere [91] appears. If a Klimenkov source was really intended, cite Klimenkov et al., J. Nucl. Mater. 417 (2011) 124–126, doi:10.1016/j.jnucmat.2010.12.261, and only for the 16.3 dpa, 250–450 °C HFR data.

### [92] Dethloff 2016
- **Metadata:** OK. Crossref matches: C. Dethloff, E. Gaganidze, J. Aktaa, Nucl. Mater. Energy 9 (2016) 471–475. "EURO- FER" is a hyphenation artifact.
- **Source checked:** https://api.crossref.org/works/10.1016/j.nme.2016.05.009 ; https://www.sciencedirect.com/science/article/pii/S2352179115300508
- **Citations:**
  - p.29 (swelling below 0.01–0.1 % at 32 dpa): PLAUSIBLE. The paper covers HFR at 15 dpa/300 °C and BOR-60 at 15 and 32 dpa/330 °C and reports void density and size. A swelling value was not confirmed.
  - supp p.111 (compiled ranges): PLAUSIBLE.
  - supp p.118 (⟨100⟩ fraction range): SUPPORTED for 78 %, since 4.9/(4.9+1.4) = 78 % at HFR 15 dpa. The paper says ⟨100⟩ loops are more frequent.
  - supp p.118 Table S20 row (300 °C, 15 dpa): PLAUSIBLE with a likely error. The paper gives ⟨100⟩ at 4.9 ×10^21 m^-3 with mean diameter 4.1 nm, and ½⟨111⟩ at 1.4 ×10^21 m^-3 with 2.8 nm. The table gives ⟨111⟩ 4.2 nm and ⟨100⟩ 2.8 nm, so the diameters look swapped and 4.2 should be 4.1.
  - supp p.120 (higher HFR cavity density attributed to helium): SUPPORTED. The paper attributes enhanced void formation in HFR to about 10 appm He, compared with negligible He in BOR-60.
  - supp p.120 Table S21 row (HFR, "~2–3", "<0.1 %"): PLAUSIBLE. The paper reports a mean void size of about 2.3 nm and a much higher density than BOR-60. The swelling value is unverified.
- **Recommended fix:** none to the metadata. Check and correct the Table S20 diameters (⟨100⟩ 4.1 nm, ½⟨111⟩ 2.8 nm).

### [93] Coppola & Klimenkov 2019
- **Metadata:** OK. Crossref matches: R. Coppola, M. Klimenkov, Metals 9(5) (2019) 552.
- **Source checked:** https://api.crossref.org/works/10.3390/met9050552
- **Citations:**
  - p.29 (SANS at "2.5 and 8 dpa"): SUPPORTED in substance, with a minor inaccuracy. The doses are 2.7 and 8.4 dpa (the supplement uses the correct values).
  - supp p.111 (0.9–3 nm, TEM and SANS): PLAUSIBLE.
  - supp p.120 (mean diameter 0.9→1.3 nm at 2.7–8.4 dpa and 300 °C; 2.6 nm at 16.3 dpa and 250 °C; volume fractions 0.001–0.006): PLAUSIBLE. The abstract confirms the doses, temperatures, volume fractions of 0.001–0.006 and "average radii increasing with the dose". It reports radii, so check that 0.88/1.32/2.58 nm are diameters and not radii.
  - supp p.120 Table S21 rows: PLAUSIBLE, with the same radius-versus-diameter check.
- **Recommended fix:** Change "2.5 and 8 dpa" to "2.7 and 8.4 dpa" on p.29. Check the radius/diameter convention in the SANS values.

### [94] Arakawa 2007
- **Metadata:** OK. Crossref matches: K. Arakawa, K. Ono, M. Isshiki, K. Mimura, M. Uchikoshi, H. Mori, Science 318(5852) (2007) 956–959.
- **Source checked:** https://api.crossref.org/works/10.1126/science.1145386
- **Citations:**
  - p.31 (long-range glide in high-purity iron): SUPPORTED. The paper shows in-situ one-dimensional diffusion of ½⟨111⟩ nanoloops in high-purity α-Fe. Note the paper shows 1D diffusion, not necessarily "long-range" in a quantified sense.
  - supp p.83 (1D transport by in-situ microscopy): SUPPORTED.
  - supp p.90 (glide distances far exceed obstacle spacing; genuinely 1D): PLAUSIBLE. 1D motion is supported; the comparison with obstacle spacing is the manuscript's inference.
- **Recommended fix:** none

### [95] Osetsky 2003
- **Metadata:** OK. Crossref matches: Yu.N. Osetsky, D.J. Bacon, A. Serra, B.N. Singh, S.I. Golubov, Philos. Mag. 83(1) (2003) 61–91. The issue number is missing.
- **Source checked:** https://api.crossref.org/works/10.1080/0141861021000016793
- **Citations:**
  - p.31 (C and N trap clusters strongly, solutes are weaker obstacles, flights limited to about 22 nm): MISMATCH. The paper is MD of self-interstitial cluster transport in pure Fe and Cu and does not treat C/N trapping or solute obstacles. The 22 nm value is not from this paper.
  - supp p.83 (1D transport of SIA clusters established atomistically): SUPPORTED.
  - supp p.90 (C/N trapping and substitutional solutes shorten flights): MISMATCH for the solute claim. At most it supports the pure-metal 1D motion and Burgers vector changes.
- **Recommended fix:** Keep [95] for 1D transport only. Cite a solute or impurity study (for example [97] with the correct authors, or Terentyev et al. on Cr/C interactions) for the trapping and flight-length claim, and give a traceable source for "22 nm".

### [96] Terentyev 2007 (author list wrong)
- **Metadata:** MAJOR. Crossref lists the authors as D.A. Terentyev, L. Malerba, M. Hou (three authors). The manuscript lists "D. Terentyev, L. Malerba, G. Bonny, A. T. Al Motasem, M.-C. Marinica, F. Willaime", with Hou missing and four authors added. Title, journal, volume, year and article number (Phys. Rev. B 75 (2007) 104108) are correct.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.75.104108 (the APS abstract page returned 403)
- **Citations:**
  - p.31 (C/N trapping, solute obstacles, about 22 nm flights): MISMATCH, based on known content. The paper is MD on the dimensionality (1D vs 3D, Burgers vector rotation) of SIA cluster motion in pure bcc Fe, with no C/N or solute trapping. The 22 nm value is unverified.
  - supp p.90 (Burgers vector reoriented after a shorter flight in alloys): PLAUSIBLE only for the reorientation, change-of-direction part, which the paper covers in pure Fe. It does not support the alloy or solute attribution.
- **Recommended fix:** Correct the authors to: D.A. Terentyev, L. Malerba, M. Hou, Dimensionality of interstitial cluster motion in bcc-Fe, Phys. Rev. B 75 (2007) 104108. Cite it only for dimensionality and reorientation.

### [97] "Terentyev, Malerba, Hou 2011" (actually Terentyev, Anento, Serra, Jansson, Khater, Bonny 2011)
- **Metadata:** MAJOR. The DOI 10.1016/j.jnucmat.2010.11.053 resolves to the correct title. The real record is D. Terentyev, N. Anento, A. Serra, V. Jansson, H. Khater, G. Bonny, J. Nucl. Mater. 408(3) (2011) 272–284. The manuscript has the wrong authors (Terentyev, Malerba, Hou; the author sets of [96] and [97] appear swapped or mixed) and the wrong volume and pages (419, 336–347).
- **Source checked:** https://doi.org/10.1016/j.jnucmat.2010.11.053 (redirected to https://www.sciencedirect.com/science/article/abs/pii/S0022311510007610). The Crossref API was rate-limited.
- **Citations:**
  - p.31 (carbon traps clusters strongly): SUPPORTED. The abstract says "carbon–vacancy complexes act as strong traps for ½⟨111⟩ loops". The nitrogen part and the about 22 nm flight length are not from this paper.
  - supp p.90 (C and N trap gliding clusters strongly): SUPPORTED for carbon; nitrogen is not covered.
- **Recommended fix:** D. Terentyev, N. Anento, A. Serra, V. Jansson, H. Khater, G. Bonny, Interaction of carbon with vacancy and self-interstitial atom clusters in α-iron studied using metallic–covalent interatomic potential, J. Nucl. Mater. 408 (3) (2011) 272–284, doi:10.1016/j.jnucmat.2010.11.053.

### [98] Arakawa 2006
- **Metadata:** OK. Crossref matches: K. Arakawa, M. Hatanaka, E. Kuramoto, K. Ono, H. Mori, Phys. Rev. Lett. 96(12) (2006) 125506. Capitalize "Burgers".
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevLett.96.125506 (the abstract pages at APS and PubMed were not retrievable)
- **Citations:**
  - p.31 (interacting loops adopt the Burgers vector of the larger partner): PLAUSIBLE. The title and known content are in-situ TEM of ½⟨111⟩ loops in Fe changing Burgers vector through elastic interaction with neighbouring loops, not external dislocations. The "larger partner" detail could not be confirmed from a fetched abstract.
  - supp p.92 (same claim): PLAUSIBLE, same caveat. The paper concerns ½⟨111⟩→½⟨111⟩ changes, not ⟨100⟩ formation, so its relevance to templated ⟨100⟩ absorption is an analogy.
- **Recommended fix:** Capitalize "Burgers" in the title. Check the "larger partner" statement against the paper's text.

### [99] Weiß 2012
- **Metadata:** MINOR. Title, authors, journal, volume, issue, pages and year all match. The DOI is missing (10.1016/j.jnucmat.2012.03.027). "eurofer97" should be "EUROFER97" and "Journal of nuclear materials" should be "Journal of Nuclear Materials". This entry duplicates the title wrongly given in [91].
- **Source checked:** https://www.sciencedirect.com/science/article/abs/pii/S0022311512001481 ; https://www.osti.gov/etdeweb/biblio/21607125 (search listing)
- **Citations:**
  - p.33 (most of the number density below the TEM visibility limit): PLAUSIBLE. The abstract says the observable defects cannot fully explain the measured hardening, which implies invisible defects, but it does not state the distribution directly.
  - supp p.118 (1.4→1.7 ×10^22 m^-3, 3.4→4.8 nm at 15→32 dpa, 330–340 °C): SUPPORTED. The values match the abstract exactly.
  - supp p.118–119 Table S20 "all" rows: SUPPORTED.
- **Recommended fix:** Add doi:10.1016/j.jnucmat.2012.03.027 and fix the capitalization. Merge the [91] citations into [99].

### [100] Gao 2022
- **Metadata:** OK. Crossref matches: J. Gao, E. Gaganidze, J. Aktaa, J. Nucl. Mater. 559 (2022) 153409.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2021.153409 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311521006292
- **Citations:**
  - p.33 (formation free energy to 1100 K from anisotropic elasticity with Varshni-fitted elastic constants): SUPPORTED in substance. The abstract confirms anisotropic elastic self-energies with temperature-dependent stiffness tensors, valid to 1100 K. The "Varshni" functional form is not in the abstract and is unverified.
  - p.34 (Gao et al. crossover at 375–455 °C for n = 100): PLAUSIBLE, unverified. The abstract gives only the parameterization. Size-dependent ⟨100⟩/⟨111⟩ crossover temperatures are more likely in the companion cluster-dynamics paper [101], so check whether this row should cite [101].
  - supp p.88 (form ε0(T) n ln n + ε1(T) n): PLAUSIBLE. It matches a "simple analytical formula" built from elastic (n ln n) and core terms; the exact functional form was not confirmed.
- **Recommended fix:** none to the metadata. Check that the 375–455 °C crossover values come from [100] and not [101].
### [101] Gao 2022
- **Metadata:** OK. Authors (Jie Gao, Ermile Gaganidze, Jarir Aktaa), title, Acta Materialia 233 (2022) 117983 all match.
- **Source checked:** https://api.crossref.org/works/10.1016/j.actamat.2022.117983 ; WebSearch (the same group's companion paper "Parameterization on formation free energy of dislocation loops up to 1100 K in bcc iron", J. Nucl. Mater. 2021, S0022311521006292)
- **Citations:**
  - p.33: PLAUSIBLE. This is a cluster dynamics study of the relative ½⟨111⟩/⟨100⟩ loop populations in α-Fe, which matches. The loop free energy up to 1100 K comes from the group's companion 2021 JNM paper, cited here as [100]. That fits the sentence's structure, but I could not check from the fetched record that [101] uses the anisotropic-elasticity/Varshni parameterization (no abstract was available).
  - p.34: PLAUSIBLE. The title says the paper looks at how C15 cluster stability affects the ⟨100⟩ population, which matches the topic. The specific claim that C15 clusters "on collapse supply most ⟨100⟩ nuclei" is consistent with the paper but not confirmed from an abstract.
- **Recommended fix:** none

### [102] Marinica 2012
- **Metadata:** OK. Authors, title, PRL 108(2) 025501 (2012) match. The only problem is a line-break artifact in the DOI ("10.1103/ PhysRevLett...").
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevLett.108.025501
- **Citations:**
  - p.34: SUPPORTED. This is the original paper proposing C15 Laves-phase interstitial clusters in bcc Fe.
- **Recommended fix:** none (remove the stray space in the DOI)

### [103] Calder 1993
- **Metadata:** OK. A.F. Calder, D.J. Bacon, JNM 207 (1993) 25–45, title matches.
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3115(93)90245-T
- **Citations:**
  - supp p.65: SUPPORTED. This is a classic MD study of displacement cascades in α-Fe, a suitable source for defect-production characteristics. The individual table values are not verified.
- **Recommended fix:** none

### [104] Stoller 1997
- **Metadata:** OK. R.E. Stoller, G.R. Odette, B.D. Wirth, "Primary damage formation in bcc iron", JNM 251 (1997) 49–60.
- **Source checked:** https://api.openalex.org/works/doi:10.1016/S0022-3115(97)00256-0 (Crossref returned 429)
- **Citations:**
  - supp p.65: SUPPORTED. MD primary damage in bcc Fe (defect survival and clustered fractions against PKA energy).
- **Recommended fix:** none

### [105] Stoller 1998
- **Metadata:** MINOR. Crossref gives the title as "Point Defect Cluster Formation in Iron Displacement Cascades Up to 50 keV", MRS Proceedings vol. 540, 1998, article 679. The manuscript's title, volume and year match. The page range 679–684 cannot be confirmed from Crossref, which gives only the start page. The symposium volume is usually dated 1999 in print.
- **Source checked:** https://api.crossref.org/works/10.1557/PROC-540-679
- **Citations:**
  - supp p.65: SUPPORTED. The abstract describes MD cascades in Fe up to 50 keV at 100–900 K, including clustered fractions and cluster size distributions.
- **Recommended fix:** none required (optionally cite as MRS Symp. Proc. 540 (1999) 679)

### [106] Gao 2000
- **Metadata:** OK. F. Gao, D.J. Bacon, Yu.N. Osetsky, P.E.J. Flewitt, T.A. Lewis, JNM 276 (2000) 213–220, title matches.
- **Source checked:** https://api.crossref.org/works/10.1016/S0022-3115(99)00180-4
- **Citations:**
  - supp p.65: SUPPORTED. MD on sessile interstitial clusters formed in α-Fe cascades, which is relevant to defect-production characteristics.
- **Recommended fix:** none

### [107] Soneda 2001
- **Metadata:** MINOR. The title is truncated. The full title is "Vacancy loop formation by 'cascade collapse' in α-Fe: A molecular dynamics study of 50 keV cascades". The manuscript writes "α-iron" and drops "of 50 keV cascades". Authors, Phil. Mag. Lett. 81(9) (2001) 649–659 and DOI are correct.
- **Source checked:** https://api.crossref.org/works/10.1080/09500830110062799
- **Citations:**
  - supp p.65: SUPPORTED. MD of 50 keV cascades in α-Fe with vacancy loop formation, which is relevant to defect production.
- **Recommended fix:** Restore the full title: "...in α-Fe: a molecular dynamics study of 50 keV cascades".

### [108] Terentyev 2006
- **Metadata:** MINOR. The subtitle is missing. The full title is "Effect of the interatomic potential on the features of displacement cascades in α-Fe: A molecular dynamics study". Authors (7, in order), JNM 351 (2006) 65–77 and DOI are correct.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2006.02.020
- **Citations:**
  - supp p.65: SUPPORTED. MD cascades in α-Fe comparing interatomic potentials (defect production and clustering).
- **Recommended fix:** Add ": a molecular dynamics study" to the title.

### [109] Calder 2010
- **Metadata:** OK. A.F. Calder, D.J. Bacon, A.V. Barashev, Yu.N. Osetsky, Phil. Mag. 90(7–8) (2010) 863–884, title matches.
- **Source checked:** https://api.crossref.org/works/10.1080/14786430903117141
- **Citations:**
  - supp p.65: SUPPORTED. MD of high-energy cascades in Fe and the origin of large SIA clusters.
- **Recommended fix:** none

### [110] Souidi 2011
- **Metadata:** OK. A. Souidi, M. Hou, C.S. Becquart, L. Malerba, C. Domain, R.E. Stoller, JNM 419 (2011) 122–133, title matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1016/j.jnucmat.2011.08.049 (Crossref returned 429)
- **Citations:**
  - supp p.65: PLAUSIBLE. The paper couples MD/BCA primary damage in Fe with OKMC long-term evolution, so it is related to defect production. It is more about how primary damage affects microstructure evolution than a direct source of production numbers.
- **Recommended fix:** none

### [111] Bjorkas 2012
- **Metadata:** OK. C. Björkas, K. Nordlund, M.J. Caturla, PRB 85(2) 024105 (2012), title matches.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.85.024105
- **Citations:**
  - supp p.65: PLAUSIBLE. MD cascade defect distributions in α-Fe used as OKMC input for damage accumulation. The topic is relevant, but the paper is not mainly a production-statistics study.
  - supp p.91: MISMATCH (extraction artifact, not a real citation). The "[111]" in "[111] + ½[11̄1̄] → [100] [49]" is a Burgers-vector Miller index, and the sentence actually cites [49]. Check that the LaTeX source does not contain a stray \cite there.
- **Recommended fix:** none for the entry. At supp p.91, confirm the Burgers vector is typeset as math and not as a citation.

### [112] Zarkadoula 2013
- **Metadata:** OK. E. Zarkadoula, S.L. Daraszewicz, D.M. Duffy, M.A. Seaton, I.T. Todorov, K. Nordlund, M.T. Dove, K. Trachenko, J. Phys.: Condens. Matter 25(12) 125402 (2013), title matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1088/0953-8984/25/12/125402
- **Citations:**
  - supp p.65: SUPPORTED. Large-scale MD of high-energy (up to 0.5 MeV) cascades in Fe, which is relevant to defects produced by fusion-spectrum PKAs.
- **Recommended fix:** none

### [113] De Backer 2018
- **Metadata:** MAJOR. The entry merges two different papers:
  - The title "A model of defect cluster creation in fragmented cascades in metals [based on morphological analysis]" belongs to A. De Backer, C. Domain, C.S. Becquart, L. Lunéville, D. Simeone, A.E. Sand, K. Nordlund, J. Phys.: Condens. Matter 30(40) (2018) 405701, doi:10.1088/1361-648X/aadb4e.
  - The author list (De Backer, Sand, Nordlund, Lunéville, Simeone, Dudarev) and the article number 26001 belong to "Subcascade formation and defect cluster size scaling in high-energy collision events in metals", EPL 115(2) (2016) 26001, doi:10.1209/0295-5075/115/26001.
  - The DOI given, 10.1209/0295-5075/124/26001, did not resolve in OpenAlex (404). "EPL 124 (2018)" is wrong for both papers.
- **Source checked:** https://api.openalex.org/works/doi:10.1209/0295-5075/124/26001 (404); https://api.openalex.org/works/doi:10.1088/1361-648X/aadb4e ; https://api.openalex.org/works/doi:10.1209/0295-5075/115/26001 ; WebSearch (iopscience, PubMed 30124201)
- **Citations:**
  - supp p.65: PLAUSIBLE. Both candidate papers deal with subcascade fragmentation and cluster size scaling in high-energy cascades (Fe and W), which suits a fission-versus-fusion comparison. Which paper is meant is ambiguous.
- **Recommended fix:** Pick one paper. For cluster size scaling with PKA energy, use De Backer, Sand, Nordlund, Lunéville, Simeone, Dudarev, EPL 115 (2016) 26001, doi:10.1209/0295-5075/115/26001. For the morphological model, use De Backer, Domain, Becquart, Lunéville, Simeone, Sand, Nordlund, J. Phys.: Condens. Matter 30 (2018) 405701, doi:10.1088/1361-648X/aadb4e.

### [114] Nordlund 2018
- **Metadata:** OK. The 13 authors in order, title, JNM 512 (2018) 450–479 and DOI all match.
- **Source checked:** https://api.openalex.org/works/doi:10.1016/j.jnucmat.2018.10.027
- **Citations:**
  - supp p.65: SUPPORTED. A comprehensive review of primary radiation damage, including NRT/arc-dpa and cascade cluster statistics.
- **Recommended fix:** none

### [115] Fu 2004
- **Metadata:** OK. C.-C. Fu, F. Willaime, P. Ordejón, PRL 92(17) 175503 (2004), title matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1103/PhysRevLett.92.175503
- **Citations:**
  - supp p.67: PLAUSIBLE. This is a standard DFT source for SIA and di-SIA formation, binding and migration energies in α-Fe (for example the SIA migration energy of about 0.34 eV). The specific Table S3 values are not verified.
- **Recommended fix:** none

### [116] Fu 2005
- **Metadata:** OK. C.-C. Fu, F. Willaime, PRB 72(6) 064117 (2005), title matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1103/PhysRevB.72.064117
- **Citations:**
  - supp p.67: PLAUSIBLE. This is the standard DFT source for He solution, migration and He–vacancy binding energies in Fe, which fits He-related rate parameters. The specific table values are not verified.
- **Recommended fix:** none

### [117] Becquart 2005
- **Metadata:** OK. C.S. Becquart, C. Domain, J. Foct, "Ab initio calculations of some atomic and point defect interactions involving C and N in Fe", Philosophical Magazine 85(4–7) (2005) 533–540. The DOI resolves to this paper. There is only a spacing artifact in the DOI line.
- **Source checked:** https://api.openalex.org/works/doi:10.1080/02678370412331320152
- **Citations:**
  - supp p.67: PLAUSIBLE. DFT C/N–vacancy and C–SIA interactions in Fe could supply carbon-trapping parameters for EUROFER. Whether Table S3 actually uses C/N interaction energies is not verified. If it has no C/N parameters, the citation is unneeded.
- **Recommended fix:** none (check that Table S3 has a C/N-related entry)

### [118] Becquart 2007
- **Metadata:** MINOR. C.S. Becquart, C. Domain, "Ab initio calculations about intrinsic point defects and He in W", NIMB 255(1) 23–26. OpenAlex gives 2006 (online date), while the print issue is Feb 2007, so the manuscript's year is acceptable. Title, authors, volume, pages and DOI match.
- **Source checked:** https://api.openalex.org/works/doi:10.1016/j.nimb.2006.11.006
- **Citations:**
  - supp p.67: MISMATCH. The paper covers point defects and He in tungsten, not iron. It cannot justify reaction-rate parameters for EUROFER/α-Fe unless Table S3 explicitly uses W data, for example for the W solute in EUROFER, and that seems unlikely.
- **Recommended fix:** Remove it from the Table S3 citation list, or replace it with an Fe source for He/point defects, e.g. Fu & Willaime PRB 72 (2005) 064117 [116] (already cited) or Becquart & Domain on He in Fe.

### [119] Becquart 2018
- **Metadata:** MAJOR. The author list is wrong. The manuscript gives "C. S. Becquart, F. Soisson". The record gives C.S. Becquart, R. Ngayam-Happy, P. Olsson, C. Domain. F. Soisson is not an author. Title, JNM 500 (2018) 92–109 and DOI match. OpenAlex lists the year as 2017 (online date), but volume 500 is dated March 2018.
- **Source checked:** https://api.openalex.org/works/doi:10.1016/j.jnucmat.2017.12.022
- **Citations:**
  - supp p.67: PLAUSIBLE. DFT binding of SIAs and small SIA clusters to solutes (Cr, Mn, Ni, Cu, Si, P, etc.) in Fe, which is relevant to solute-trapping parameters. The specific values are not verified.
- **Recommended fix:** Change the authors to "C. S. Becquart, R. Ngayam-Happy, P. Olsson, C. Domain".

### [120] Borodin 2007
- **Metadata:** OK. V.A. Borodin, P.V. Vladimirov, "Diffusion coefficients and thermal stability of small helium–vacancy clusters in iron", JNM 362(2–3) (2007) 161–166.
- **Source checked:** https://api.openalex.org/works/doi:10.1016/j.jnucmat.2007.01.019
- **Citations:**
  - supp p.67: PLAUSIBLE. Gives diffusivities and dissociation energetics of small HenVm clusters in Fe, which suits He-V cluster rate parameters. The specific values are not verified.
- **Recommended fix:** none
### [121] Seletskaia 2005
- **Metadata:** OK — Seletskaia, Osetsky, Stoller, Stocks; PRL 94(4) 046403 (2005); title, DOI match.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevLett.94.046403 ; https://api.openalex.org/works/doi:10.1103/PhysRevLett.94.046403 (abstract)
- **Citations:**
  - supp p.67 (Table S3, group cite [28,80,115–126]): PLAUSIBLE — DFT study of He interstitial/defect energetics in Fe (tetrahedral site preferred, magnetic effects); appropriate source for He formation/binding energies, but the specific Table S3 values cannot be matched from the abstract; the group citation does not say which parameter comes from this paper.
- **Recommended fix:** none for metadata; say in Table S3 which parameter(s) come from [121], for example the He interstitial formation energy.

### [122] Olsson 2007
- **Metadata:** OK — Olsson, Domain, Wallenius; PRB 75(1) 014110 (2007); title, DOI match.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.75.014110
- **Citations:**
  - supp p.67 (Table S3 group cite): PLAUSIBLE — ab initio Cr–vacancy/SIA interaction energies in bcc Fe fit a EUROFER (Fe–9Cr) parameter table; which values come from it cannot be checked.
- **Recommended fix:** none for metadata; say in the table which parameters come from it.

### [123] Rahman 2023
- **Metadata:** MAJOR — wrong co-author: the manuscript lists "M. M. Rahman, L. K. Béland, N. Mousseau"; the record lists **Md Mijanur Rahman, Fedwa El-Mellouhi, Normand Mousseau**. Title, journal, vol. 7(9), 093602 (2023) and DOI are correct.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevMaterials.7.093602 ; https://api.openalex.org/works/doi:10.1103/PhysRevMaterials.7.093602
- **Citations:**
  - supp p.67 (Table S3 group cite): PLAUSIBLE — k-ART study of 1–8 vacancy clusters in α-Fe gives formation energies and activation/migration barriers, which suits vacancy-cluster binding/migration inputs. The specific numbers cannot be checked.
- **Recommended fix:** Change the authors to "M. M. Rahman, F. El-Mellouhi, N. Mousseau".

### [124] Schuler 2015 (actually Barouh et al. 2015)
- **Metadata:** MAJOR — the author list belongs to a different paper. The manuscript lists "T. Schuler, M. Nastar, F. Soisson"; the record for this title/DOI (PRB 92(10) 104102, 2015) lists **Caroline Barouh, Thomas Schuler, Chu-Chun Fu, Thomas Jourdan**. Title, journal, volume, article number and DOI match the Barouh et al. paper.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.92.104102
- **Citations:**
  - supp p.67 (Table S3 group cite): PLAUSIBLE — covers vacancy-mediated diffusion of interstitial solutes (C, N, O) and vacancy–solute binding in α-Fe; relevant only if Table S3 has solute–vacancy or C/N parameters. The specific values cannot be checked.
- **Recommended fix:** Change to "C. Barouh, T. Schuler, C.-C. Fu, T. Jourdan, Predicting vacancy-mediated diffusion of interstitial solutes in α-Fe, Phys. Rev. B 92 (2015) 104102." If a Schuler–Nastar–Soisson paper was meant instead, use that paper's own metadata and DOI.

### [125] Chiesa 2009
- **Metadata:** MINOR — the title is truncated. The full title is "Free energy of a ⟨110⟩ dumbbell interstitial defect in bcc Fe: Harmonic and anharmonic contributions". Authors, PRB 79(21) 214109 (2009) and DOI are correct.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.79.214109
- **Citations:**
  - supp p.67 (Table S3 group cite): PLAUSIBLE — gives SIA formation free energy/entropy in Fe, which suits SIA formation parameters. Specific values cannot be checked.
- **Recommended fix:** Add the missing subtitle ": Harmonic and anharmonic contributions".

### [126] Nguyen-Manh 2006
- **Metadata:** OK — Nguyen-Manh, Horsfield, Dudarev; PRB 73(2) 020101(R) (2006); title, DOI match.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevB.73.020101
- **Citations:**
  - supp p.67 (Table S3 group cite): PLAUSIBLE — DFT SIA configurations and formation energies in bcc metals, including the Fe ⟨110⟩ dumbbell; suits SIA parameters. Specific values cannot be checked.
- **Recommended fix:** none

### [127] Érdi & Tóth 1989
- **Metadata:** OK — the book exists: P. Érdi, J. Tóth, *Mathematical Models of Chemical Reactions: Theory and Applications of Deterministic and Stochastic Models*, Manchester University Press, 1989 (ISBN 0719022088; Princeton University Press co-edition, ISBN 0691085323). There is no DOI.
- **Source checked:** WebSearch → https://www.abebooks.com/9780719022081/ , https://www.amazon.com/dp/0691085323 , Semantic Scholar record
- **Citations:**
  - supp p.66: SUPPORTED — standard monograph on deterministic mass-action kinetics, with the stoichiometric formulation ẋ = S·rate(x).
- **Recommended fix:** none (optionally add "(Nonlinear Science series)").

### [128] Koch 2010
- **Metadata:** OK — I. Koch, "Petri Nets – A Mathematical Formalism to Analyze Chemical Reaction Networks", Molecular Informatics 29(12) 838–843 (2010); DOI matches.
- **Source checked:** https://api.crossref.org/works/10.1002/minf.201000086
- **Citations:**
  - supp p.66: PLAUSIBLE — a review of Petri-net (bipartite species/reaction) representations of CRNs, a reasonable source for the graph view. The specific "directed hyperedges with two tails" wording is hypergraph language rather than Petri-net language, but the two are equivalent.
- **Recommended fix:** none

### [129] Faeder 2009
- **Metadata:** OK — Faeder, Blinov, Hlavacek, "Rule-Based Modeling of Biochemical Systems with BioNetGen", Methods Mol. Biol. vol. 500 (Systems Biology), Humana Press, 2009, pp. 113–167; DOI matches. The book editor (I. V. Maly) could be added.
- **Source checked:** https://api.crossref.org/works/10.1007/978-1-59745-525-1_5 ; https://api.openalex.org/works/doi:10.1007/978-1-59745-525-1_5
- **Citations:**
  - supp p.68: SUPPORTED — BioNetGen generates reaction networks from rules over structured species.
- **Recommended fix:** none (optionally add editor I. V. Maly).

### [130] Danos & Laneve 2004
- **Metadata:** OK — V. Danos, C. Laneve, "Formal molecular biology", Theor. Comput. Sci. 325(1) 69–110 (2004); DOI matches.
- **Source checked:** https://api.crossref.org/works?query.bibliographic=Formal+molecular+biology+Danos+Laneve (the result lists DOI 10.1016/j.tcs.2004.03.065); WebSearch
- **Citations:**
  - supp p.68: SUPPORTED — this is the paper that introduced the κ (Kappa) calculus. Later Kappa-language papers (e.g. Danos et al. 2007, CONCUR) could also be cited.
- **Recommended fix:** none

### [131] Wei & Kuo 1969
- **Metadata:** MINOR — the title is slightly off. The manuscript has "A lumping analysis in monomolecular reaction systems. …"; the record has "Lumping Analysis in Monomolecular Reaction Systems. Analysis of the Exactly Lumpable System" (no leading "A"). Authors (J. Wei, J. C. W. Kuo), Ind. Eng. Chem. Fundam. 8(1) 114–123 (1969) and DOI are correct.
- **Source checked:** https://api.openalex.org/works/doi:10.1021/i160029a019 ; https://pubs.acs.org/doi/10.1021/i160029a019 (via WebSearch listing). Crossref fetch was rate-limited (HTTP 429).
- **Citations:**
  - supp p.68: SUPPORTED — this is the classic paper on exact lumping of reaction systems, which the manuscript's graph coarsening is compared to.
- **Recommended fix:** Remove the leading "A" from the title.

### [132] Hagberg 2008
- **Metadata:** OK — Hagberg, Schult, Swart, "Exploring Network Structure, Dynamics, and Function using NetworkX", Proc. 7th Python in Science Conf. (SciPy2008), Pasadena, pp. 11–15 (eds. Varoquaux, Vaught, Millman).
- **Source checked:** WebSearch → https://proceedings.scipy.org/articles/TCWV9851 , https://www.osti.gov/biblio/960616
- **Citations:**
  - supp p.69: SUPPORTED — this is the standard citation for NetworkX (MultiDiGraph).
- **Recommended fix:** none (optionally add editors).

### [133] Meyer 1988
- **Metadata:** MINOR — the book exists: B. Meyer, *Object-Oriented Software Construction*, 1st ed., Prentice Hall, 1988 (ISBN 0136290493). Prentice Hall is usually located at Englewood Cliffs, NJ (or Prentice Hall International, Hemel Hempstead), not "New York".
- **Source checked:** WebSearch → https://en.wikipedia.org/wiki/Object-Oriented_Software_Construction , https://www.amazon.com/dp/0136290493 , https://openlibrary.org/books/OL2033837M
- **Citations:**
  - supp p.70: SUPPORTED — Meyer's 1988 book is the source of the open–closed principle.
- **Recommended fix:** Change the publisher location to "Englewood Cliffs, NJ" (or drop the location).

### [134] Hindmarsh 2005
- **Metadata:** OK — all seven authors in the correct order; ACM TOMS 31(3) 363–396 (2005); DOI matches.
- **Source checked:** https://researchr.org/publication/HindmarshBGLSSW05 ; https://computing.llnl.gov/projects/sundials/publications (Crossref rate-limited)
- **Citations:**
  - supp p.71: SUPPORTED — the SUNDIALS/CVODE paper describing the variable-order BDF integrator. The post-step nonnegativity floor is the manuscript's own implementation choice.
- **Recommended fix:** none

### [135] Gardner 2022
- **Metadata:** OK — Gardner, Reynolds, Woodward, Balos; ACM TOMS 48(3) (2022) 1–24; DOI 10.1145/3539801 matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1145/3539801 ; https://computing.llnl.gov/projects/sundials/publications
- **Citations:**
  - supp p.71: SUPPORTED — modern SUNDIALS paper (vector, matrix and linear-solver abstractions).
- **Recommended fix:** none (optionally give the ACM article number instead of pages 1–24).

### [136] Saad & Schultz 1986
- **Metadata:** MINOR — capitalization only: "Gmres" should be "GMRES". Authors, SIAM J. Sci. Stat. Comput. 7(3) 856–869 (1986) and DOI are correct.
- **Source checked:** https://api.openalex.org/works/doi:10.1137/0907058 ; https://epubs.siam.org/doi/10.1137/0907058 (WebSearch)
- **Citations:**
  - supp p.71: SUPPORTED — SPGMR is SUNDIALS' scaled, preconditioned GMRES, and this is the original GMRES paper.
- **Recommended fix:** Protect the acronym in BibTeX ({GMRES}).

### [137] Davis 2010
- **Metadata:** OK — T. A. Davis, E. Palamadai Natarajan, "Algorithm 907: KLU, A Direct Sparse Solver for Circuit Simulation Problems", ACM TOMS 37(3) article 36, pp. 1–17 (2010); DOI matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1145/1824801.1824814 ; https://dl.acm.org/doi/10.1145/1824801.1824814 (WebSearch listing)
- **Citations:**
  - supp p.71: SUPPORTED — the original KLU paper.
- **Recommended fix:** none

### [138] Curtis 1974
- **Metadata:** OK — Curtis, Powell, Reid, "On the Estimation of Sparse Jacobian Matrices", J. Inst. Math. Appl. (now IMA J. Appl. Math.) 13(1) 117–119 (1974); DOI matches. Only minor capitalization differences.
- **Source checked:** https://api.openalex.org/works/doi:10.1093/imamat/13.1.117
- **Citations:**
  - supp p.71: SUPPORTED — this is the CPR method: grouping structurally orthogonal columns so one finite difference covers several columns. The graph-coloring view was formalized later by Coleman & Moré (1983, SIAM J. Numer. Anal. 20:187); consider citing that paper too.
- **Recommended fix:** none (optionally also cite Coleman & Moré 1983).

### [139] Hager 1989
- **Metadata:** OK — W. W. Hager, "Updating the Inverse of a Matrix", SIAM Review 31(2) 221–239 (1989); DOI matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1137/1031049
- **Citations:**
  - supp p.72: SUPPORTED — a standard review of the Sherman–Morrison–Woodbury identity and its use for low-rank updates.
- **Recommended fix:** none

### [140] Gösele & Seeger 1976
- **Metadata:** OK — U. Gösele, A. Seeger, "Theory of bimolecular reaction rates limited by anisotropic diffusion", Phil. Mag. 34(2) 177–193 (1976); DOI matches.
- **Source checked:** https://api.openalex.org/works/doi:10.1080/14786437608221934 (reconstructed abstract)
- **Citations:**
  - supp p.78: PLAUSIBLE — the paper does treat diffusion-limited reaction rates with anisotropic diffusion and non-spherical reaction volumes, reducing to the isotropic 3D case as a limit. Coordinate rescaling that turns the reaction surface into an ellipsoid is the standard method of that paper. The specific leading-order result using the geometric-mean D̄ = (DxDyDz)^(1/3) cannot be confirmed from the abstract; check it against the paper's text.
  - supp p.84: SUPPORTED — the paper supports the claim that anisotropic or unequal diffusivities need a diffusion-tensor transformation (mapping to an ellipsoidal/elliptical reaction surface). The abstract also discusses reactions between 1D migrating defects on non-parallel paths.
- **Recommended fix:** none for metadata; check the geometric-mean statement against the paper's equations.
### [141] Adjanor 2022a
- **Metadata:** OK. Crossref gives G. Adjanor, "Complete characterization of sink-strengths for mutually 1D-mobile defect clusters: Extension to diffusion anisotropy analog cases", J. Nucl. Mater. 572 (2022) 153970. Everything matches. The space inside the DOI ("10.1016/ j.jnucmat...") comes from PDF line wrapping. Check that the .bib has no space.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2022.153970 ; https://arxiv.org/abs/1808.10362 (preprint abstract)
- **Citations:**
  - supp p.78: SUPPORTED. The paper derives 1D-1D sink strengths that depend on the ratio of the two diffusivities through a power law.
  - supp p.79: PLAUSIBLE. A loss that goes as c_A^2 c_B + c_B^2 c_A fits a 1D-1D sink strength that scales with the square of concentration. The abstract does not state the exact prefactor 6*pi^2*R^4.
  - supp p.81: PLAUSIBLE. The paper covers random rotations (the "1DR" cases). The fraction f_v = (v-1)/v is not in the abstract, so it is unverified.
  - supp p.83: SUPPORTED. This matches the title and scope of [141, 142].
  - supp p.84 (cubic loss frequency, Eq. S16): PLAUSIBLE. Same reasoning as p.79; the exact expression is unverified.
  - supp p.84 (D_A != D_B maps to diffusion toward an ellipse, diffusion-tensor transformation): SUPPORTED. This is the "diffusion anisotropy analog" in the title.
  - supp p.84 (conditional probability p_x for intersecting vs. skew pairs): PLAUSIBLE. The preprint abstract does not mention it. Check the full text.
- **Recommended fix:** none. Remove the line-wrap space in the DOI if it is in the .bib.

### [142] Adjanor 2022b
- **Metadata:** OK. Crossref gives G. Adjanor, "Complete characterization of sink-strengths for 1D to 3D mobilities of defect clusters: Bridging between limiting cases with effective sink-strengths calculations", J. Nucl. Mater. 572 (2022) 154010. Everything matches. "ef- fective" is a hyphenation artifact from the PDF.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2022.154010
- **Citations:**
  - supp p.78: SUPPORTED. The paper extends sink strengths across the 1D-to-3D mobility range.
  - supp p.83: PLAUSIBLE. The title covers 1D-to-3D mobility in general, while the sentence says "mutually 1D-mobile". The mutually 1D case is mainly in [141].
- **Recommended fix:** none

### [143] Jourdan 2025
- **Metadata:** OK. Crossref gives T. Jourdan, G. Adjanor, "Computation of sink strengths in complex microstructures", J. Nucl. Mater. 616 (2025) 156021. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2025.156021 ; https://www.sciencedirect.com/science/article/pii/S0022311525004155
- **Citations:**
  - supp p.78: PLAUSIBLE. The abstract describes a lifetime-based OKMC method for microstructures with several sink types, some of them mobile, including the 1D-to-3D SIA transition and the effect of vacancy mobility. It relates to "arbitrary ratios D_A/D_B" only loosely. The sentence is cut off in the extraction, so the rest of the claim cannot be checked.
- **Recommended fix:** none

### [144] Jansson 2013
- **Metadata:** OK. Crossref gives V. Jansson, L. Malerba, A. De Backer, C.S. Becquart, C. Domain, J. Nucl. Mater. 442(1-3) (2013) 218-226. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2013.08.052 ; https://arxiv.org/pdf/1407.7218 (full text)
- **Citations:**
  - supp p.80: MISMATCH. The full text has no orientation-averaged toroidal cross-section sigma = pi^2 r_c r_n. It only reports that random vs. parallel loop orientation changes the 1D sink strength by about 11%. The formula needs another source, or should be marked as the authors' own derivation.
  - supp p.81: SUPPORTED, with a transcription issue. The paper defines d_j = (sqrt(3)/2) a_0 as the first-nearest-neighbour distance, and l_ch = d_j sqrt(n_ch) with n_ch the number of jumps before a change of direction. The manuscript's "d_j = 3a/2" is either a PDF garbling or a real error. It must read sqrt(3)a/2, and the "glide jump distance" wording should match.
  - supp p.90: MISMATCH. The paper does not discuss carbon, nitrogen or substitutional-solute trapping. It treats n_ch as a free parameter in pure-Fe OKMC. It fits at most as the source of the l_ch parameterization, not of the trapping claim.
- **Recommended fix:** Check the p.81 equation reads d_j = sqrt(3) a/2. Drop [144] as support for the toroidal cross-section on p.80 (cite a correct source or derive it) and for solute trapping on p.90.

### [145] Trinkaus 1992
- **Metadata:** OK. Crossref gives H. Trinkaus, B.N. Singh, A.J.E. Foreman, "Glide of interstitial loops produced under cascade damage conditions: Possible effects on void formation", J. Nucl. Mater. 199(1) (1992) 1-5. Everything matches. The DOI suffix is lowercase "l" in the manuscript and "L" in Crossref; DOIs are case-insensitive.
- **Source checked:** https://api.crossref.org/works/10.1016/0022-3115(92)90433-l
- **Citations:**
  - supp p.83: PLAUSIBLE. This is the founding paper on 1D glide of cascade-produced SIA loops and its effect on void formation. Rafting and dislocation decoration are covered mainly in the 1997 paper [146], so the pairing is acceptable.
- **Recommended fix:** none

### [146] Trinkaus 1997
- **Metadata:** OK. Crossref gives H. Trinkaus, B.N. Singh, A.J.E. Foreman, "Mechanisms for decoration of dislocations by small dislocation loops under cascade damage conditions", J. Nucl. Mater. 249(2-3) (1997) 91-102. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/S0022-3115(97)00230-4
- **Citations:**
  - supp p.83: SUPPORTED. The paper covers dislocation decoration, loop-loop and loop-dislocation interactions, and rafting under cascade damage.
- **Recommended fix:** none

### [147] Heinisch 2000
- **Metadata:** OK. Crossref gives H.L. Heinisch, B.N. Singh, S.I. Golubov, "The effects of one-dimensional glide on the reaction kinetics of interstitial clusters", J. Nucl. Mater. 283-287 (2000) 737-740. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/S0022-3115(00)00258-0
- **Citations:**
  - supp p.83: SUPPORTED. It is a KMC study of the reaction kinetics of 1D-gliding SIA clusters.
  - supp p.85: SUPPORTED. It tests analytic mixed 1D/3D sink-strength expressions against KMC with occasional direction changes. The link to the specific closure Eq. (S28) is unverified.
  - supp p.90: PLAUSIBLE. It is relevant to correlated 1D kernels, but the l^-1 decay of the segment kernel cannot be confirmed from this paper.
- **Recommended fix:** none

### [148] Heinisch 2007
- **Metadata:** OK. Crossref gives H.L. Heinisch, H. Trinkaus, B.N. Singh, "Kinetic Monte Carlo studies of the reaction kinetics of crystal defects that diffuse one-dimensionally with occasional transverse migration", J. Nucl. Mater. 367-370 (2007) 332-337. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2007.03.034
- **Citations:**
  - supp p.83: SUPPORTED. It is a KMC study of the reaction kinetics of 1D-diffusing defects.
  - supp p.85: SUPPORTED. It checks analytic expressions against KMC for 1D migration with occasional transverse or direction changes.
- **Recommended fix:** none

### [149] Huang 2019
- **Metadata:** OK. Crossref gives Shu Huang, Jaime Marian, "Rates of diffusion controlled reactions for one-dimensionally-moving species in 3D space", Philos. Mag. 99(20) (2019) 2562-2583. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1080/14786435.2019.1635721
- **Citations:**
  - supp p.83: SUPPORTED. It gives analytical reaction rates for 1D-moving species.
  - supp p.90: PLAUSIBLE. It is relevant to 1D reaction kernels; the specific l^-1 statement is unverified.
- **Recommended fix:** none

### [150] Torney 1983
- **Metadata:** OK. Crossref gives D.C. Torney, H.M. McConnell, "Diffusion-limited reaction rate theory for two-dimensional systems", Proc. R. Soc. Lond. A 387 (issue 1792) (1983) 147-170. The issue number is missing from the manuscript, which is optional.
- **Source checked:** https://api.crossref.org/works/10.1098/rspa.1983.0055
- **Citations:**
  - supp p.84: SUPPORTED. The abstract says the 2D rate function goes to zero asymptotically as (ln t)^-1, which matches the logarithmically slow long-time 2D absorbing-disc result.
- **Recommended fix:** none (optionally add issue 1792)

### [151] Redner 2001
- **Metadata:** OK. The DOI 10.1017/CBO9780511606014 redirects to Cambridge Core "A Guide to First-Passage Processes" (S. Redner, Cambridge University Press, 2001). Crossref returned HTTP 429, so this was confirmed through the doi.org redirect and a WebSearch of the Cambridge Core and BU pages.
- **Source checked:** https://doi.org/10.1017/CBO9780511606014 (redirects to cambridge.org/core/books/guide-to-firstpassage-processes/59066FD9754B42D22B028E33726D1F07); WebSearch results (cambridge.org, physics.bu.edu/~redner)
- **Citations:**
  - supp p.84: SUPPORTED. The book covers 2D first passage and absorption by a disc, including the logarithmic long-time behaviour.
- **Recommended fix:** none

### [152] Kohnert 2015
- **Metadata:** OK. Crossref gives Aaron A. Kohnert, Brian D. Wirth, "Cluster dynamics models of irradiation damage accumulation in ferritic iron. II. Effects of reaction dimensionality", J. Appl. Phys. 117(15) (2015) 154306. Everything matches; only the capitalization of "Effects" differs.
- **Source checked:** https://api.crossref.org/works/10.1063/1.4918316
- **Citations:**
  - supp p.90: SUPPORTED/PLAUSIBLE. The abstract describes 1D reaction kinetics in cluster dynamics that depend on hop length. The specific l^-1 decay is unverified.
- **Recommended fix:** none (capitalize "Effects")

### [153] Terentyev 2008
- **Metadata:** OK. Crossref gives D.A. Terentyev, T.P.C. Klaver, P. Olsson, M.-C. Marinica, F. Willaime, C. Domain, L. Malerba, "Self-trapped interstitial-type defects in iron", Phys. Rev. Lett. 100(14) (2008) 145503. Everything matches; the manuscript has "D. Terentyev" without the middle initial, which is fine.
- **Source checked:** https://api.crossref.org/works/10.1103/PhysRevLett.100.145503
- **Citations:**
  - supp p.90: MISMATCH. The paper is about intrinsic self-trapping: small SIA clusters in pure Fe take non-parallel, sessile configurations that are immobile. It is not about trapping of gliding clusters by carbon, nitrogen or substitutional solutes, or about Burgers-vector reorientation by solutes. This verdict comes from the title and the known content of the paper; no abstract could be fetched (OSTI blocked, PubMed page empty).
- **Recommended fix:** Remove [153] from the solute-trapping sentence, or reword it to cite [153] only for intrinsic self-trapping or immobilization of SIA clusters in Fe.

### [154] Klimenkov 2011
- **Metadata:** OK. Crossref gives M. Klimenkov, E. Materna-Morris, A. Möslang (Crossref itself misspells the name as "Mölsang"), "Characterization of radiation induced defects in EUROFER 97 after neutron irradiation", J. Nucl. Mater. 417(1-3) (2011) 124-126. The manuscript entry is correct.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2010.12.261 ; https://inis.iaea.org/records/t3443-3gt31 (abstract)
- **Citations:**
  - supp p.111: PLAUSIBLE. The data are EUROFER97 neutron irradiation at 16.3 dpa and 250-450 C with loops and He bubbles. The numeric ranges quoted are aggregated over several references and cannot be tied to this paper.
  - supp p.118/119 (Table S20 rows at 16.3 dpa, 250-450 C): PLAUSIBLE. The dose (16.3 dpa), temperature range, 1/2<111> Burgers vector (f<100> = 0) and a loop-density maximum near 300 C all agree with the abstract. The individual sizes and densities (7 nm / 2e21 m^-3, 14 nm / 4e21, 35 nm / 3e20, and the 400-450 C values) are unverified. Also, a table row with a density of 1e19 m^-3 at 450 C is labelled "all" loops, while the abstract mentions He bubbles; check whether the 400-450 C rows are loops or cavities.
- **Recommended fix:** none for the metadata. Verify the Table S20 values against the paper's table.

### [155] Klimenkov 2020
- **Metadata:** OK. Crossref gives M. Klimenkov, U. Jäntsch, M. Rieth, A. Möslang, "Correlation of microstructural and mechanical properties of neutron irradiated EUROFER97 steel", J. Nucl. Mater. 538 (2020) 152231. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2020.152231 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311519316460 (abstract)
- **Citations:**
  - supp p.111: PLAUSIBLE. This is EUROFER97 neutron loop and void data. The aggregated numeric ranges are unverified.
  - supp p.118 (ii) (density falls about 3 orders of magnitude from 300 to 415 C; sizes from 5-8 to 100-180 nm; f<100> from 10-27% to 79-87%): PLAUSIBLE. The abstract gives 250-415 C, 1/2<111> loops at 250-300 C and a preference for <100> loops at higher temperature, which agrees qualitatively. The specific numbers are unverified. The abstract describes an "average dose of 16.3 dpa", while the manuscript gives sample-specific doses of 13.4-18.1 dpa. The abstract also says 1/2<111> loops were found at 250 and 300 C, yet the table gives 10-27% <100> at those temperatures. Both points should be checked against the full text.
  - supp p.118 (iii): PLAUSIBLE. It is one of several sources for the scatter in the <100> fraction.
  - supp p.118/119 (Table S20 rows): PLAUSIBLE (unverified values).
  - supp p.119 (Fig. S6 caption): PLAUSIBLE. It depends on the Table S20 values above.
- **Recommended fix:** none for the metadata. Check the per-temperature doses (13.4-18.1 dpa vs. the 16.3 dpa average) and the <100> fractions at 250-300 C against the paper.

### [156] Chauhan 2021
- **Metadata:** OK. Crossref gives Ankur Chauhan, Qian Yuan, Christian Dethloff, Ermile Gaganidze, Jarir Aktaa, "Post-irradiation annealing of neutron-irradiated EUROFER97", J. Nucl. Mater. 548 (2021) 152863. Everything matches.
- **Source checked:** https://api.crossref.org/works/10.1016/j.jnucmat.2021.152863 ; https://www.sciencedirect.com/science/article/abs/pii/S0022311521000866 (abstract)
- **Citations:**
  - supp p.111: PLAUSIBLE. The paper covers EUROFER97 irradiated at 330 C to 15 dpa, with 1/2<111> and <100> loops and cavities studied by TEM. The aggregated numeric ranges are unverified.
  - supp p.118 (iii): PLAUSIBLE for the 15 dpa <100> fraction. The abstract reports only 330 C / 15 dpa (plus annealing at 550 C for 3 h). Chauhan 2021 should not be cited for the "27-31% at 32 dpa" part; in Table S20 those rows come from [87].
  - supp p.118/119 (Table S20 rows at 330 C, 15 dpa: 6.2 nm, 7 nm, 8 nm, 13 nm): PLAUSIBLE. The condition matches. The values are unverified, and the 8 nm and 13 nm "all" rows at 2e21 and 4.1e21 m^-3 are probably post-annealed states or cavities. Label them if they are annealed (550 C/3 h) and not as-irradiated.
- **Recommended fix:** none for the metadata. In Table S20, say which [156] rows are as-irradiated and which are after annealing. Cite [156] only for the 15 dpa part of statement (iii).
