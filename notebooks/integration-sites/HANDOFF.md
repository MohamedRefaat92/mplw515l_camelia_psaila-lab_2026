# Integration-site analysis — handoff

Written for an agent or analyst picking this up without the conversation that produced it. It covers
why the analysis exists, what was tried and rejected, what is assumed, and where the traps are. The
code explains *what* it does; this explains *why*, and what already went wrong.

Status as of 2026-11-18: pipeline runs end to end on the cluster (job 35742986, 1 h 49 m), **15
sites** called, 1,213 cell-site assignments over 1,186 cells. Nothing is committed to git — the
project has a standing rule that commits happen only when asked.

What changed on 2026-11-18, since it changes most of the numbers below: the genome-mapped channel
was already calling sites, but its reads were being discarded afterwards, so clone sizes were built
from the 3′ channel alone. They now carry both channels. The site count went 9 → 15 and the labelled
cell count 552 → 1,186. Site labels are now pinned to loci by `metadata/integration_site_registry.csv`
rather than assigned by rank — see §7.

**Three results changed direction on the same day, and the pattern is worth more than any of them.**
Each came from using a denominator that already had the answer built into it:

1. *Co-occurrence.* Comparing observed pairs against `|A|·|B| / labelled cells` assumes complete
   detection, so "below chance" was guaranteed. The fix was a control, not a correction: sites in
   different mice cannot share a cell, so cross-mouse pairs measure the doublet floor directly.
2. *Leakage.* "B0304 leaks at 14% against B0303's 0.2%" divided labelled cells by labelled cells.
   B0304 has 92% of its cells in the mCherry⁻ tubes. Per cell the animal effect is 2.4-fold, and the
   ranking inverts.
3. *The 100%-mCherry⁻ sites.* Called "selection" on the assumption that a site found in those tubes
   must read 100%. It need not — discovery and assignment both pool every tube — so the asymmetry is
   real, just thin at 2–6 cells.

**Never divide labelled cells by labelled cells.** All three had the shape of a strong result.

**13 sites, not 15, for any repertoire statement.** IS14 and IS15 clear the cell floor only on the
two channels' reads summed, and have no 5′ junction either. `integration_sites.csv` carries
`called_by` and `pooled_only` so this does not have to be remembered.

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
| assignment tier | clip ≥10 bp | matched against the mini-reference only, by exact position, no aligner; **adds cells**, see §4b. The floor was 8 until the decoy flanks measured where they stop returning hits |
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
assignment step produced would make the method confirm itself — and `assign_clips_exact()` places
clips against it, with no aligner. The discriminator is **position, not score**: a clip belongs to a
site only if its junction-adjacent base lands on that site's junction base, at the offset carried in
the sequence name. One accepted landing site out of ~400 per flank is what makes a 10 bp clip
decisive. (bwa used to do this placement, via a `read_mini_hits()` that no longer exists: its `-T 12`
silently enforced a 12 bp floor while `min_assign` advertised 8, so 169 reads were written and
dropped with no error. The floor is now whatever `min_assign` says, and `min_assign` is 10 because
that is where the decoy flanks stop returning hits.)

Result: cells labelled rise from 339 to **552** on the 3′ channel alone, and to **1,186** once the
genome-mapped channel's cells are folded in as well. The wild-type tubes stayed at zero through all
of that except for two cells — one in Rx_1 and one in Rx_5, both at IS04 — which is the measured
cross-tube rate, 0.16% of labelled cells.

**A third trap, hit on 2026-11-18: the cache that only breaks on restart.** The genomic stage writes
one CSV a tube so the expensive scan can be resumed. Read back with guessed column types, that cache
is a landmine — the two wild-type tubes legitimately return **zero rows**, `readr` types every column
in a header-only file as character, and binding them to a tube with real rows dies with *"Can't
combine `pos` <character> and `pos` <double>"*. The failure lands at the start of the rerun, after
the first run has already been paid for, and it is caused by the two tubes that are working
correctly. `read_genomic_cache()` declares the schema instead. **Any cache whose contents can
legitimately be empty needs its types declared, not inferred.**

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

15 sites called — **13 defensibly**, see `pooled_only` below — 12 in the mini-reference, **1,186 distinct cells labelled**. `mouse` is the
dominant hashtag class, which is 96–100% of each site's cells.

| id | locus | gene | cells | mouse | note |
|---|---|---|---|---|---|
| IS02 | chr9:53,079,916 | *Ddx10* | 342 | B0303 | |
| IS01 | chr13:83,681,311 | *Mef2c* | 318 | B0303 | |
| IS04 | chr7:83,792,638 | *Abhd17c* | 178 | B0301 | |
| IS03 | chr19:38,933,742 | *Hells* | 130 | B0303 | also the spliced-exon channel, §7.3 |
| IS10 | chr13:5,676,485 | *Gm35330* | 99 | B0304 | **invisible to the 3′ channel** — its clips are MAPQ 0 |
| IS05 | chr9:64,711,166 | *Gm36459* | 57 | B0304 | **Rx_4 (HSC) only**; 2.6 UMIs a cell, the highest of any site |
| IS11 | chr1:16,580,387 | *Stau2* | 25 | B0304 | new from the genomic channel |
| IS09 | chr2:31,009,148 | *Fnbp1* | 19 | B0301 | |
| IS12 | chr11:111,656,449 | intergenic | 12 | B0304 | new from the genomic channel |
| IS06 | chr1:67,556,980 | intergenic | 10 | B0303 | |
| IS08 | chr18:17,071,116 | intergenic | 10 | B0303 | |
| IS13 | chr2:126,704,897 | *Trpm7* | 6 | B0303 | all six cells in mCherry⁻ tubes |
| IS07 | JH584304.1:35,417 | intergenic | 3 | B0303 | **artefact.** Unplaced scaffold, 23-base homopolymer in its flank, no 5′ junction, and all three cells in mCherry⁻ tubes. Held out of the mini-reference by `max_homopolymer`; its cells come from the genomic channel, which does not consult the flank |
| IS14 | chr10:19,804,209 | intergenic | 2 | B0304 | **`pooled_only`** — clears the cell floor only on both channels' reads summed, no 5′ junction, held out of the mini-reference. Exclude from repertoire counts |
| IS15 | chr9:82,806,975 | *Phip* | 2 | B0304 | as IS14 |

`called_by` in `integration_sites.csv` carries this per site: 9 sites called by both channels, 4 by the genome-mapped channel alone, 2 (IS14, IS15) only on the pooled reads. **The defensible repertoire count is 13.**

**Site IDs are pinned to loci, not to ranks.** `metadata/integration_site_registry.csv` holds the
mapping; a locus in it keeps its number for good and a new locus takes the lowest number never used.
This replaced rank ordering after the genomic channel changed the ranking and silently repointed IS01
through IS04 at different loci. Do not infer importance from the number, and do not renumber the
registry.

**What bounds the error in each tier.** The decoy flanks (100 of them, zero hits at 10 bp and above)
bound the *assignment* tier. The genome-mapped channel has its own decoy now — `genomic_decoy_rate()`,
run as the `genomic-decoy` stage — which repeats the scan with the terminus permuted to the same base
composition. Genome-wide across all eight tubes (job 35768797, 2 h): **15 decoy reads against 1,422
real, 1.1%**, scattered over 7 positions where the real reads concentrate at 47. No decoy position is
dense enough to be called a site, so the rate is **1.1% at the read level and zero at the site
level**. Rx_8 makes the point best: a wild-type tube with zero real reads and 5 decoy reads, so in
material with no provirus the real 13-mer finds nothing and a random one finds five.

Composition-match the decoy *k*-mer, not the whole terminus. The real 13-mer `AAAATCTCTAGCA` has six
adenines; permuting all 18 bases leaves the tested 13 with three, and A-rich *k*-mers are commoner in
the genome, so that version understates the rate in the direction that flatters the channel. The
first attempt made exactly that mistake.

**Open items closed on 2026-11-18, with what each turned out to be:**

| item | outcome |
|---|---|
| host-gene UMI confound | **real**, and strand-dependent. Junction reads carry `GN:Z:<host gene>` with `xf` bit 8 set, so they count as host-gene UMIs — but only when the provirus is intragenic *and* in the gene's orientation. `cis_confound` flags the four: IS03, IS04, IS10, IS11. 20.9% of carriers' *Abhd17c* UMIs at IS04, 23.5% of *Gm35330* at IS10, 12.1% of *Stau2* at IS11, 2.3% of *Hells* at IS03. The earlier "177 spurious *Ddx10* UMIs" was wrong — IS02 is antisense |
| 5′ anchoring | the absorbed-base repair is sound. 2 of 76 placed hits start with a soft clip, 40 are clean end-to-end, NM = 0 throughout |
| 5′ floor | `min_assign` is advertised for this channel and inert — only the discovery tier is written to `five_prime.fa`, so no clip under 25 bp ever reaches an aligner. Not a silent discard; they are simply never sent |
| `cluster_window = 20` | hides nothing. All 1,396 genomic reads inside a called site sit at offset exactly 0 |
| TSD on IS10–IS15 | unreachable, not failing. The whole 5′ channel is 76 reads and none lands within 2 kb of those six loci |
| 1,186 → figures | 5 cells in the stroma tubes (in neither embedded panel) and 33 in no demultiplexer output, leaving 1,148 before QC. The 33 are missing rather than filtered and are worth a look |
| ambient (#43) | holds, and the new channel is a check on it: junction-bearing barcodes went 568 → 669 while the non-cell count stayed at **15**. Signal added without background |
| metadata queries 1–2 | computed from `demux_res_*`; see INFORMATION-GAPS. Query 3, the FACS mCherry⁺ fraction, is still outside the data |

**The clip table (step 1 of the statistical-assignment plan), added 2026-11-18.**
`clip_table.csv.gz` is one row per molecule — tube + barcode + UMI — with no floor applied, and
`clip_matches.csv.gz` is every clip against every flank at every mismatch count. Three things about
it matter more than the files:

- **One matcher, several views.** `clip_site_matches()` is the only code that decides whether a clip
  matches a flank; `assign_clips_exact()` is a two-line filter over it. It records `nm` and `budget`
  per row rather than a 0/1 verdict, which is what lets a view ask a different question — exact
  matching is `nm == 0`, a lower floor is one argument — and lets the gate prove the views agree.
- **The gate is an equality that stops the run.** Filtered at the current floor, the match table must
  return 15 sites / 1,213 assignments / 1,186 cells. It failed on first run and caught two real bugs:
  a union across a molecule's reads that let a 5-base read's matches through a floor applied to its
  30-base sibling (52 spurious assignments), and a dplyr trap where `clip_len = max(clip_len)`
  followed by `clip = first(clip[clip_len == max(clip_len)])` in one `summarise()` selects an
  arbitrary row, because by the second line `clip_len` is already the scalar.
- **The gate belongs on reads, not molecules.** The assignment tier is a per-read filter, so only a
  per-read view can reproduce it to the row. The molecule table is a convenience built on top.

With it, the floor question is measured rather than argued: under exact matching the 100 decoy flanks
are clean from **6 bp** up and dirty below (1 of 11 molecules at 5 bp, 10 of 17 at 4 bp, all of them
at 1–2 bp). Lowering 10 → 6 would add 73 assignments over 47 new cells, 4.0%, with no wild-type
additions. **The floor stays at 10**; the measurement is recorded so the choice is costed.

**Writes are verified now, because they were silently truncating.** On 2026-11-18 the `/well/ludwig`
fileset hit its 391 TB limit with ~23 GB left and writes into `output/` began returning short files —
no header, no error, no warning. A deliberate 5,000-row write produced 88 lines. Four outputs were
corrupted that way including `cell_site_assignments.csv`, which the report reads. Every table now
goes through `write_csv_checked()`, which counts the lines back and errors on a mismatch; it caught
three further truncations and halted instead of leaving plausible-looking corruption. The two large
artefacts are gzipped (18× smaller) for the same reason. **Do not remove either guard** — the
filesystem is shared and the condition will recur.

**The detection model (step 3), added 2026-11-18 — this is where the factor of three is.**
`detection_model.csv`: per-site λ from a zero-truncated Poisson on per-cell molecule counts, the
detection-corrected clone size with a bootstrap interval, and an observed-against-fitted check.

| site | detected | molecules/cell | detection | clone size (95% CI) |
|---|---|---|---|---|
| IS05 | 57 | 2.63 | **91%** | 63 (60–67) |
| IS02 | 342 | 1.41 | 52% | 660 (593–753) |
| IS01 | 318 | 1.17 | 28% | 1,154 (942–1,480) |
| IS11 / IS06 / IS08 | 25 / 10 / 10 | 1.16 / 1.80 / 1.40 | 26 / 73 / 51% | 95 / 14 / 20 |
| nine others | 451 | ≤ 1.09 | — | — |

Three things not to get wrong with it. **The molecule count must come from the clip table**, not from
summing `umis_three + umis_genomic` over the cell table — that double-counts the 83 molecules in both
channels *and* drops every sub-floor molecule, and the two errors cancel badly enough that the wrong
answer looks reasonable (mean 1.271 against the correct 1.236). **`NA` is the estimator declining**,
not a gap: at a truncated mean of 1.00 the fit is degenerate and the clone size goes to infinity.
**Read the fit from obs against exp, never from `var_mean_ratio`** — a truncated Poisson with small λ
is underdispersed by construction (0.41, 0.17), so comparing it to 1 invents a problem. The Poisson
fits: 240 IS02 cells at one molecule against 232 expected, 76 at two against 85. The overdispersion
the plan expected came from the shared-UMI artefact and disappeared with it.

**Repertoire completeness, answered without the FACS record.** The denominator needed is the mCherry⁺
cells *sequenced*, which is in hand: 43,332 across Rx_4 and Rx_6, of which 1,131 are labelled — a
**2.6% assumption-free floor** — and ~2,900 after detection correction, so the fifteen sites account
for roughly **7% of the sorted mCherry⁺ population**. The missing 93% leave a trace that is
accounted for rather than mysterious: of 1,422 junction-bearing cells, 236 carry a molecule matching
no called site, and of the 76 long enough to place, 67 land at MAPQ < 30 (repetitive flank), 15 do
not place, and 29 sit at seven uncalled positions with 1–2 cells each — all below `min_cells = 3`.
So the missing repertoire visible here is seven micro-loci plus the blind spot, not a hidden
population. Step 4's posteriors are not needed for this.

**Step 4 (posteriors over K+1 hypotheses) is not built, and here is the number that decided it.**
Across the whole dataset **25 molecules match more than one called site exactly**, all with clips of
1–2 bases; above 2 bp every matching molecule matches exactly one site. The mouse-restricted prior —
the plan's "largest free gain", cutting K from 15 to the 2–7 sites in the cell's own mouse — resolves
**8** of those 25. So the posterior framework would operate on eight molecules at lengths the plan
itself says not to chase. Its remaining value is supplying weights and a detection covariate to
Phase 6, which does not exist yet; build it when that consumer does, and keep the plan's EM warning
verbatim when you do (priors from clone sizes and clone sizes from posteriors is a loop whose fixed
point pulls ambiguous mass onto the largest clone — the inflation already visible in all seven
cross-mouse landings sitting on IS01 and IS02).

**Three changes on 2026-11-18, in response to "every number in the report must be reproducible".**

1. **The floor is 6, not 10, and it comes with `exact_below = 15`.** These are one decision, not two:
   under the `L %/% 10` mismatch budget a 6-base clip is allowed one mismatch, and 21 of the 35
   molecules at 6 bp then hit a decoy flank; under exact matching, none do. The per-base decoy curve
   is clean from 6 bp up and dirty below (1 of 9 molecules at 5 bp, 10 of 17 at 4 bp). The change
   gains 57 three-prime assignments at 6–9 bp and **loses 5** at 10–12 bp that had rested on a
   mismatch — the losses are the point, not a regret. `assign_keep()` holds the rule in one place so
   the equality gate cannot pass against a stale version of it.
2. **Cell evidence is restricted to Cell Ranger's called barcodes.** 33 of what the cell table called
   cells were not cells — the BAM carries corrected `CB` tags for droplets Cell Ranger did not call,
   and this pipeline read them without checking. They are the ambient population. The filter is
   applied where cells are counted and **deliberately not** to the clip and genomic tables, because
   `ambient_rate()` is computed from exactly the rows it removes. 34 removed on the current run.
3. **Five measurements that were prose are now outputs.** The floor curve, the completeness figures,
   the unassigned-molecule decomposition and the ambient rate were computed in scratch scripts and
   typed into the notebook; they would have gone stale silently the first time anything upstream
   moved, while twelve other chunks recomputed around them. They are now
   `clip_floor_curve.csv`, `repertoire_completeness.csv`, `unassigned_molecules.csv` and
   `ambient_rate.csv`, read by render chunks. The mono-replicate control needs the hashtag calls, so
   it is a notebook chunk rather than a pipeline output — and it had to be **moved** in the document,
   because it was written into the floor section where `site_replicate` and `embedding` do not yet
   exist and would have failed at render.

Current counts: **15 sites, 1,222 assignments over 1,195 cells**, gate passing.

**The assignment floor is ONE OBJECT, not two constants (2026-11-18).** `IS_DEFAULTS$assign_floor` is
`list(min_len = 6L, exact_below = 15L)` and there is a load-time assertion that refuses an incoherent
pair. This is a guard, not tidiness: under the `L %/% 10` mismatch budget **21 of the 35 molecules at
6 bp match a decoy flank — a 60% error rate** — so `min_len = 6` is defensible only with exact
matching beneath it. Asking a reader not to separate the two numbers is not a guard; making them one
value is. The assertion was tested by feeding it an incoherent policy, and it rejects at `source()`.

**What licenses the floor, and what does not.** 0.31 expected spurious assignments across the 931
accepted molecules, against 1,222 real — **0.02 per site**, which changes no clone size, interval or
figure. Report that absolute count, not a percentage: a rate on the 25 assignments the 6 bp band adds
reads as a tolerance being brushed when it is 0.26 events. *General rule, since it recurs: a marginal
error rate is the right frame when the increment will be interpreted on its own, an absolute count
when it will be pooled.* Two things that look like supporting evidence and are not: the
mono-replicate control's zero cross-mouse errors in that band (only ~0.16 were expected, so zero is
consistent with anything), and the wild-type tubes (as always).

**`n_decoy` is 1,000, and decoy comparisons are capped at `exact_below`.** At 100 flanks the curve
showed zero hits at 6 bp, which read as clean and was an absence of evidence — chance predicts 1.67,
so *p* = 0.19. At 1,000 the measurement exists: **15 observed against 13.5 predicted, obs/exp = 1.11**,
which validates the matcher against theory and is the form the standing check should read. Capping
decoy comparisons at 15 bp cost nothing — chance matching at 30 bases is 4⁻³⁰ — and took `assign`
from ~13 minutes back to **3 m 52 s** with the decoy count unchanged. `ratio_interpretable` marks the
lengths where obs/exp can be read at all; below an expectation of ~5 events a single chance hit sends
the ratio to 30, and a standing check that alarms on Poisson noise is the monotonicity test's mistake
again.

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
