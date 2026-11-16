# Integration-site analysis — handoff

Written for an agent or analyst picking this up without the conversation that produced it. It covers
why the analysis exists, what was tried and rejected, what is assumed, and where the traps are. The
code explains *what* it does; this explains *why*, and what already went wrong.

Status as of 2026-11-16: pipeline runs end to end, 9 sites called, report renders. Nothing is
committed to git — the project has a standing rule that commits happen only when asked.

---

## 1. The question

Mouse bone marrow transduced with a lentiviral vector (`LeGO_SFFV_mCherry_DRE`), sorted on mCherry
protein, sequenced with 10x 3′ v3 GEX + hashtags. Eight tubes:

| tube | compartment | sort | notes |
|---|---|---|---|
| Rx_1, Rx_2 | stroma | — | excluded from this analysis by request |
| Rx_3 | HSC | mCherry− | |
| Rx_4 | HSC | mCherry+ | carries most of the signal |
| Rx_5 | HSC | wild-type | **negative control — never saw the vector** |
| Rx_6 | MNC | mCherry+ | |
| Rx_7 | MNC | mCherry− | |
| Rx_8 | MNC | wild-type | **negative control** |

The transgene assay answers *whether* a cell carries the vector. This analysis asks *where it
integrated*, which turns a genomic locus into a **clone label** — cells sharing an integration site
descend from one transduced ancestor. That is the route to clonality without a dedicated barcode
library.

## 2. The idea, and why it works

Transcription of the cassette is supposed to terminate at the polyA site in the 3′ LTR's R region. It
leaks. A minority of transcripts read through to the end of the provirus and on into flanking host
DNA. Such a read aligns to the reporter contig as far as its last base and soft-clips the rest — and
**that clipped tail is host genomic sequence naming the integration locus**. The read still carries
its 10x cell barcode, so the locus arrives attached to a cell.

Two pieces of biology make this tractable, and both took effort to establish:

**(a) The junction is at exactly one position.** Integrase performs 3′ processing on the blunt
reverse-transcription product: it cuts off the terminal `GT`, leaving a recessed `CA` that attacks
host DNA. Every provirus is joined to the genome at that `CA`. The reporter contig ends
`…TCTCTAG**CA**` — that exact processed terminus. So **contig position 3236 is the integration
junction itself**, and a clip anchored anywhere else cannot be one. This is enzymology, not a
threshold, which is why it is the first filter and costs no sensitivity.

**(b) 3′ capture concentrates reads where we need them.** The cell barcode and UMI attach at the polyA
tail, and only the barcode-bearing fragment is amplified, so reads fall within ~300–600 bases of the
tail. The LTR is at the 3′ end of the transcript, so it carries 10–20× the reads of the coding
sequence. The same bias that makes mCherry nearly undetectable makes this analysis possible.

## 3. The reference contig — design decisions and what was rejected

`data_hub/10x_resources/reference_additions/LeGO_SFFV_mCherry_DRE_transcript.fa`, 3,236 bp = plasmid
2395–5630, built into the Cell Ranger reference by `scripts/slurm/cellranger/mkref.sh`.

| decision | reason |
|---|---|
| start at plasmid 2395 (mCherry's ATG) | the 5′ UTR gets no reads in a 3′ library, and this makes mCherry land at contig 1–708 so the annotation checks itself |
| end at plasmid 5630 (last base of 3′ LTR) | nothing past the LTR integrates. The polyA cleavage site is earlier, at contig **3152**; the 84 bases beyond are the read-through substrate this analysis depends on |
| **only one LTR copy** | the plasmid's two LTRs are byte-identical over all 181 bases. With both present every read landing in them multi-maps and is discarded — removing exactly the 3′ end where the reads are. The GTF cannot fix this: the GTF governs counting, not placement. |
| plasmid kept separately, not as reference | ~2.6 kb is bacterial backbone that never enters a cell; used only to recognise vector-derived clips |

**Rejected approaches, so you don't re-propose them:**

- *Masking the 5′ LTR in the full plasmid* — doesn't simplify anything and doesn't generalise to other
  plasmid architectures.
- *Separate contigs per reporter* — loses reads spanning the P2A and loxP junctions.
- *A GTF-generator script* — written, then deleted as over-engineered. Provenance now lives in the
  GTF's own comment header so it travels with the file.

The genes annotated on the contig are `mCherry` (1–708), `DRE` (781–1830), `WPRE-LTR` (1939–3236).
They are three windows onto **one mRNA**, not three genes. Do not sum or average them.

## 4. The method as it stands

Two channels, both off the same BAMs. The 3' channel below is the primary one and defines the site
list; the 5' channel (§4a) is a smaller, independent confirmation.

```
BAM → extract clips → bwa → cluster → sites
```

Filter chain, in order (each clip counted once, under the first reason it meets):

| filter | value | rationale |
|---|---|---|
| anchor ends exactly at junction | 3236 | the integrase `CA`. Removes ~660 reads/tube that clip the transcript's own polyA at the cleavage site — 90% of those are ≥80% adenine, against 1% at the junction |
| not vector | clip starts `GTGGCGCCCGAACAGG` | see §6, this one is load-bearing |
| not polyA | <80% A | |
| has corrected CB and UB | | |
| discovery tier | clip ≥25 bp | long enough to place alone genome-wide; **defines** the site list |
| assignment tier | clip ≥8 bp | placed against the mini-reference only; **adds cells**, see §4b |
| MAPQ | ≥30 | on the realigned clip, discovery only |
| cluster | single-linkage, 20 bp | junctions scatter a few bases because clips start wherever the read began |
| site floor | ≥3 cells and ≥3 UMIs | the one genuinely arbitrary parameter |

Discovery and assignment are different questions. Discovery searches 2.7 Gb; assignment searches a
few kb of already-called flanks and additionally requires an exact position. See §4b.

## 4a. The 5' channel and the TSD check

The provirus's other end. Both LTRs are identical and the contig keeps one copy, so a read crossing
the **host → 5' LTR** junction also lands on that copy, clipping host to its *left*. Anchor: contig
**2736**, the processed 5' terminus (`...GGGGGG AC TGGAAGG...`; integrase removes the `AC`, leaving
the `TG` that partners the `CA` at 3236).

**These anchors smear and the repair is mandatory.** The contig continues upstream of 2736, so host
bases matching it are absorbed by greedy extension and the reported start slides left; a `GGGGGG` run
at 2728–2733 makes this common. Anchoring on 2736 exactly returns **1 clip**; a ten-base window with
the absorbed bases added back returns **150**. The 3' channel is immune only because the contig stops
at 3236. Skip the repair and junctions land tens of bases out.

`read_bwa_hits()` takes a `channel` argument because the convention inverts: the junction-adjacent
base is the clip's *first* at the 3' end and its *last* at the 5' end, so the two strand cases swap.

**Target-site duplication.** Integrase cuts with a 5 bp stagger and repair duplicates those bases
onto both flanks, so a provirus's two junction coordinates differ by **4**. IS01, IS03, IS04 and IS08
all pass at exactly 4. Junctions >50 bp from any site are reported as a distinct locus, not a
failure. Written to `tsd_check.csv`.

This is the strongest check available and it has already paid for itself — see the minus-strand entry
in §8, which no 3'-channel filter could have caught.

## 4b. The mini-reference (cell assignment)

Discovery asks "is there an integration anywhere in 2.7 Gb?"; assignment asks "does this clip match
one of the flanks already called?". Treating those as one question, with one clip-length threshold,
throws away most of the data.

`write_site_flanks()` cuts ±200 bp around each **discovery-tier** site — seeding it from anything the
assignment step produced would make the method confirm itself — and `read_mini_hits()` places clips
against it. The discriminator is **position, not score**: a clip belongs to a site only if its
junction-adjacent base lands on that site's junction base, at the offset carried in the sequence
name. One accepted landing site out of ~400 per flank is what makes an 8 bp clip decisive.

Result: cells labelled rise from 339 to **552**. Both wild-type tubes stay at **zero**.

**Two traps, both hit:**

- *Low-complexity flanks attract clips that belong nowhere.* JH584304.1:35,417 has a 23-base
  homopolymer against 5–6 for every genuine site, and the mini-reference inflated it 3 → 14 cells —
  all in mCherry-negative tubes, none in the two high-signal ones, 4.7× the gain real sites saw.
  `max_homopolymer = 12` now keeps such a site out of the mini-reference; it is still reported from
  discovery, with zero assigned cells.
- *The two output tables can silently disagree.* The site list and its coordinates come from
  discovery, but the cell counts must come from assignment. `update_site_counts()` rewrites them in
  the assign stage; without it the site table reported only the cells that happened to carry a 25 bp
  clip.

## 5. Assumptions — be explicit about these

1. **The contig's last base is the processed viral terminus.** 3236 is hardcoded in `IS_DEFAULTS`. If
   the cassette is ever re-cut with a different 3′ boundary, this must move. Everything depends on it.
2. **Only intact proviruses are visible.** Integrants with deleted or rearranged LTR ends present no
   clean `CA` and are invisible. Site counts are a floor, not a census.
3. **Detection is read-through-limited, so absence is weak evidence.** A cell is labelled only if one
   of its rare read-through molecules happened to be sequenced. Most members of a real clone carry no
   label. **Clone sizes here are samples, not counts.**
4. **Per-cell confidence is near-uniform and low.** 96% of assignments rest on a single UMI. Confidence
   lives at the *site* level (202 independent cells for IS01), not the cell level. Do not use read
   depth to pick "high-confidence cells" — the variation is mostly PCR.
5. **Repetitive flanks are silently lost.** A junction whose host sequence is repetitive yields MAPQ-0
   clips. Sites are enriched for integration into unique sequence, which partly confounds the
   intragenic enrichment used as a sanity check in §7.
6. **HSC and MNC are not directly comparable.** Labelled barcodes survive QC at 82% in HSC but 58% in
   MNC, against ~69% overall in both. The HSC direction is explicable (detection favours deeply
   sequenced cells, which also pass QC); MNC is not explained.
7. **Hashtags identify the same animals across compartments.** The data supports this — IS01 is
   B0303 in both HSC and MNC — but it is an assumption inherited from the sample sheet.

## 6. The ψ artefact — the trap that no quality filter catches

About 10% of junction clips are **not** host. They start `GTGGCGCCCGAACAGGGACTTGAAAGCGAAAGG`, which
is HIV-1 sequence immediately 3′ of U5, at plasmid 213 — i.e. the **5′** LTR → ψ junction.

Why they land at the 3′ end: the two LTRs are byte-identical, and the contig deliberately keeps only
one copy. A molecule spanning the 5′ LTR's junction with ψ therefore matches the contig's 3′ terminus
perfectly, is placed there **with MAPQ 255**, and clips into vector sequence. Structurally it is
indistinguishable from a host junction. No MAPQ, coverage or cell-count filter sees it, because the
anchor genuinely *is* a unique alignment. **Only the clip's sequence gives it away.** The
de-duplication that makes the contig usable is the same thing that makes this artefact invisible.

This is a worked example of a general trap: if the vector architecture changes, look for *any* other
sequence that can present the LTR terminus as a left arm.

**Origin: resolved 2026-11-16 — packaged vector genomic RNA.** Three hypotheses, all tested:

| hypothesis | verdict |
|---|---|
| internal oligo-dT priming at plasmid 444 | refuted — no coverage cliff there (below) |
| 2-LTR circles | refuted — 0 of 1,618 junction clips carry the U5→U3 circle junction, 0 reciprocal in the 5' clips |
| plasmid carry-over | refuted — coverage stops dead at the bacterial backbone |

The backbone scan is the decisive one. Probing unplaced reads across the plasmid in Rx_4: ψ 490,
leader 80, RRE/cPPT 570, SFFV promoter 8, and then `backbone` 0, `backbone` 2, `ori` 0, `AmpR` 0,
`AmpR` 0. Plasmid DNA would cover the backbone — it is the same molecule. Wild-type Rx_5 is ~0
everywhere. What *is* covered is exactly the packaged genome's footprint: R–U5–ψ–RRE–cPPT through to
the 3' LTR, no backbone. The low SFFV signal fits too, since a promoter sits upstream of the mRNA.

The original priming scan, kept because it is what refuted hypothesis 1:

```
plasmid   Rx_4    Rx_5 (WT)   Rx_6
220        306        0        233
260        490        0        543    ← peak
444                            (A-rich tract, hypothesised priming site)
460        136        1        177    ← no cliff
700         80        0         80
                                       scanned: 141M / 156M / 138M unplaced reads
```

Coverage runs straight through the A-rich tract, so priming there is not the mechanism. The
wild-type tube is clean (3 reads in 156M), so the material is vector-derived rather than
contamination. Either way the ψ clips are removed by sequence, so this never affected the results —
it was a question about what the library contains, now answered.

## 7. Results, and the validations

9 sites called, 8 usable, **552 distinct cells labelled**.

| id | locus | gene | cells | note |
|---|---|---|---|---|
| IS01 | chr13:83,681,311 | *Mef2c* | 202 | |
| IS02 | chr9:53,079,916 | *Ddx10* | 145 | |
| IS03 | chr19:38,933,742 | *Hells* | 80 | |
| IS04 | chr7:83,792,638 | *Abhd17c* | 74 | |
| IS05 | chr9:64,711,166 | *Gm36459* | 44 | **Rx_4 (HSC) only** |
| IS06 | chr1:67,556,980 | intergenic | 7 | |
| IS08 | chr18:17,071,116 | intergenic | 6 | |
| IS09 | chr2:31,009,148 | *Fnbp1* | 6 | |
| IS07 | JH584304.1:35,417 | intergenic | **0** | **artefact, gated out.** Unplaced scaffold, 23-base homopolymer in its flank, zero cells in both high-signal tubes, no 5' junction. Excluded from the mini-reference by `max_homopolymer`; reported from discovery only. |

**Site IDs are stable but not rank-ordered.** They are assigned from discovery-tier counts with ties
broken by locus, so IS07 sits between IS06 and IS08 despite having no assigned cells. Do not infer
importance from the number.

Three checks the pipeline never used as input, so they are genuinely independent:

1. **Integration-site preference — underpowered, do not cite as a validation.** 6 of 8 usable sites
   are intragenic and *Mef2c*, *Ddx10*, *Hells* are actively transcribed in HSC, which is the
   expected lentiviral preference. But the honest null is not "the genic fraction of the genome": a
   site can only be called where a clip places uniquely, and unique sequence is itself enriched in
   genes. Against 3,000 random positions filtered through the *same* mappability requirement, 49.2%
   are intragenic, against 75% for the called sites — one-sided binomial **p = 0.135**. Consistent
   with the literature, not evidence for it at n = 8. The other three validations are unaffected and
   carry the result.
2. **FACS agreement.** 220 of 223 labelled cells are in the mCherry⁺ fraction. The sort knows nothing
   about soft-clipped reads. The 3 in mCherry⁻ are expected — that gate is on protein, and a
   provirus-carrying cell below the gate is a true positive the sort got wrong.
3. **Replicate confinement — the strongest.** `dm2_class` (hashtag demultiplex) assigns each cell to
   one of four mice. An integration event happens once, in one animal, so a real clone *must* sit in
   a single replicate. Every site is **98–100% confined** to one:

   | replicate | clones |
   |---|---|
   | B0301 | IS04, IS08 |
   | B0302 | none |
   | B0303 | IS01, IS02, IS03, IS06, IS07, IS09 |
   | B0304 | IS05 |

   And the mapping agrees across compartments. Clone burden is strikingly uneven between animals.

4. **Target-site duplication (§4a).** The two junctions of a provirus must differ by 4 bp. IS01,
   IS03, IS04 and IS08 all do — two independent read populations from opposite ends agreeing to the
   base. It also surfaced two candidates the 3' channel missed, chrX:100,691,117 and
   chr14:62,948,717, and left chr19:38,932,385 unexplained: 25 reads, 1.4 kb from IS03, more than
   IS03's own 5' junction.

**Both wild-type tubes contribute zero cells to every site** — no false positives across 28,010 cells
that never saw the vector. Rx_1/2/3/7 contribute 1–2 cells scattered across the large sites; that is
the index-hopping floor, ~1–2% of real signal.

## 8. Troubleshooting log — mistakes already made

Read this before debugging anything. Several of these cost hours.

**Analysis errors**

- **Exact string matching missed real sites.** The first pass searched one representative 40-mer per
  clip group against the genome FASTA with exact matching. A single sequencing error in the chosen
  representative makes the search fail silently and report "no match in GRCm39". This hid
  chr9:64,711,166 entirely and mis-set two others. *Always align; never string-match.* The
  representative 40-mer quoted at the time (`…CCAAACTGTT…`) was itself a 1-of-70 misread; consensus
  is `…CCAACCTGTT…`.
- **The polyA-cleavage clips were nearly fed to the aligner.** ~660 reads/tube clip at 3144–3152. 90%
  are polyA, but the other 10% is real-looking sequence with good MAPQ, real barcodes and real UMIs
  that would have become spurious sites. Nothing downstream catches them. The anchor-position filter
  is what removes them.
- **Cleavage site was recorded as 3153, actually 3152.** The old justification used read *start*
  positions in a 20-base bin — indirect, 16 reads, and incapable of single-base precision. The polyA
  clips measure the 3′ end directly with ~660 reads. Corrected in the README and GTF.
- **"Ambient is ruled out by the wild-type tubes" — wrong.** A WT tube has no vector anywhere, so it
  bounds only cross-tube and mismapping error. A sorted-negative tube containing transduced cells has
  its own internal ambient pool. Different controls bound different things.
- **Minus-strand junctions were reported ~48 bp off, and lost cells.** `read_bwa_hits()` took bwa's
  `POS` as the junction. A 3′ clip's junction-adjacent base is its *first*, which is `POS` only when
  the clip aligned forward; reverse-complemented it sits at the alignment's far end. Worse than a
  wrong coordinate: `POS` = junction − reflen + 1, and reflen varies with clip length, so clips from
  one minus-strand junction scattered over ~37 bp and the 20 bp clustering window split them,
  dropping fragments below the cell floor. Fixing it moved IS01, IS04 and IS06 and recovered 20
  assignments (323 → 343). Invisible by inspection, because each fragment still clustered tightly.
  Found only by the 5′ TSD check — plus-strand sites came out at the expected 4 bp while minus-strand
  sat at 45. **This is the argument for building the second channel: it is the only thing that caught
  it.**

- **Site IDs silently repointed between runs.** `label_sites()` ordered by cell count and three sites
  tied at 3, so the tie broke on grouping order: IS07–IS09 named different loci after a rerun, and
  every document citing them became wrong with nothing to flag it. Ties now break by chromosome and
  position.
- **The site table and the cell table disagreed about the same quantity.** `summarise_sites()` was
  fed discovery hits while `summarise_cells()` was fed assignment hits, so the site table under-
  reported cells. Both now come from the same hits via `update_site_counts()`.
- **Lowering a threshold raises the noise floor faster than the signal.** Dropping the assignment
  floor to 8 bp gained 67% more cells, but 14 of them were one low-complexity artefact — the false
  site gained 4.7× while real sites gained 1.6×. Always check the per-tube breakdown and the
  wild-type control before banking a yield improvement, not just the headline count.

**Infrastructure traps**

- **SLURM partition: use `short`, never `short_el9`/`long_el9`.** renv keys its library path on the OS
  release (`renv/library/linux-rocky-8.10/...`). On a Rocky 9 node it finds no library, bootstraps an
  empty one, and the first `library()` call fails. Cost one failed job and left a stray
  `renv/library/linux-rocky-9.7/` that had to be deleted. `short` and `long` are the *same* 115 nodes
  and differ only in time limit, so `short` is the whole pool for a short job.
- **`cellranger multi` refuses an `--output-dir` that exists but isn't its own pipestance.** Don't
  `mkdir` the output dir to stash a config; put configs beside the results.
- **Quarto must be invoked via `tools/quarto.sh`.** Calling `env/bin/quarto` directly fails because
  conda's Quarto needs its `activate.d/quarto.sh`. This was once misdiagnosed as a broken install.

**ggplot traps (all cost a silently wrong figure)**

- **`colour = NA` drops points.** Whether mapped through a scale or passed as a fixed parameter,
  ggplot treats an NA colour as a missing value and removes the row *at draw time*. This silently
  deleted 142 of 149 points. Use **`stroke = 0`** to suppress a `shape = 21` border.
- **`"transparent"` + `alpha =` makes an opaque white ring.** `"transparent"` is `#FFFFFF00` — white
  with alpha 0 — and `alpha=` overwrites the alpha channel, giving `#FFFFFFD9`.
- **`ggplot_build()` does not reveal dropped rows.** The removal happens later, in the draw step. A
  test using `ggplot_build()` reported "safe" when it wasn't.
- **Character facet/colour variables drop empty levels.** The wild-type column vanished from a facet
  grid because no labelled cell is wild-type — the panel that most needed showing. Use factors with
  explicit levels plus `drop = FALSE`. Same cause made clone colours shift between panels.
- **`scale_size` rescales per panel.** Without shared `limits`, a dot of a given size means different
  things in each facet.
- **Lab distance metric was transposed.** `dist()` on a 3×n matrix compares the three *channels*, not
  the colours. This produced a palette ranking that was exactly backwards. Lab matrices must be n×3.

## 9. Open questions

1. **MNC labelled barcodes survive QC at 58% vs 68% overall**, opposite to HSC. Unexplained.
2. **IS05 is HSC-only with 44 cells, and its cells look spatially clustered** in the HSC UMAP. Not
   tested for significance; 20 points on a 50,000-cell map can cluster by chance.
3. **The `≥3 cells and ≥3 UMIs` floor is arbitrary.** A TSD-verified site arguably needs no floor at
   all: IS09 passes on one 5' read, and chrX:100,691,117 passes while sitting below the floor in
   both channels. The waterfall is printed so it can be revisited.
4. **Clonality proper** — site counts per compartment as a diversity measure is the intended
   destination and has not been attempted. Assumption 3 is the obstacle: clone sizes are samples.
5. **chr19:38,932,385** — 25 reads, 1.4 kb from IS03, not explained by anchor absorption. Test
   whether those clips share barcodes with IS03's (same clone, same provirus) and whether the host
   side ends at a splice donor (a spliced vector-host fusion).
6. **chrX:100,691,117 and chr14:62,948,717** — 5' junctions with no called 3' site. Real sites the
   3' channel missed, or artefacts?

### Closed

- **ψ artefact origin** — packaged vector genomic RNA (§6). Priming, 2-LTR circles and plasmid all
  refuted; the backbone scan was decisive.
- **IS07 / JH584304.1** — artefact, now gated out by flank complexity (§7).
- **Multi-copy founders** — no. Every informative pair of sites co-occurs *below* chance (IS01+IS02:
  6 observed vs 26 expected in Rx_4; 2 vs 24 in Rx_6), which is the signature of disjoint clones
  rather than independent integrations. The ~2% of cells at more than one site read as doublets.

## 10. Files

| | |
|---|---|
| `src/R/integration-sites.R` | all logic, as functions. Attaches no packages — `source()` it and step through. All thresholds in `IS_DEFAULTS`. |
| `scripts/R/call-integration-sites.R` | driver, two entry points (`extract`, `call`) because bwa runs between them |
| `scripts/slurm/analysis/integration-sites.sh` | the job. ~1 min, 4.7 GB peak. |
| `scripts/slurm/setup/bwa-index.sh` | one-off index build, ~50 min |
| `notebooks/integration-sites/integration-sites.qmd` | the report |
| `output/objects/integration-sites/*.csv` | 7 output tables, including `tsd_check.csv` |
| `data_hub/10x_resources/reference_additions/README.md` | contig derivation and verification — **outside git, not version controlled** |

```bash
sbatch scripts/slurm/analysis/integration-sites.sh
tools/quarto.sh render notebooks/integration-sites/integration-sites.qmd
```

Not done: the report is **not registered in `_targets.R`** (the report block there is the user's), and
`tools/catalogue-data.sh` has not been run on the new Cell Ranger inputs.
