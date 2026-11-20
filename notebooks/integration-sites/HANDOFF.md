# Integration-site analysis — handoff

Written for an agent or analyst picking this up without the conversation that produced it. It covers
why the analysis exists, what was tried and rejected, what is assumed, and where the traps are. The
code explains *what* it does; this explains *why*, and what already went wrong.

Status as of 2026-11-18: pipeline runs end to end on the cluster (job 35742986, 1 h 49 m), **15
sites** called, **1,226 cell-site assignments over 1,199 cells**, clip-table gate PASS. Committed
through the floor change, the called-cell filter and the review responses; the project has a standing
rule that commits happen only when asked.

**Every number in this file comes from a CSV, and the CSVs move.** These figures were rewritten
twice in one day — when the assignment floor went 10 → 6, and when assignments were restricted to
called cells — and both times prose quoting them had to be chased through the report and through
here. Where a section below quotes a table, regenerate it from the named CSV rather than editing the
digits. The provenance banner in the report now verifies each of those CSVs against the manifest the
pipeline wrote beside them (§11), so a figure here that disagrees with the report means this file is
stale, not that the data moved unnoticed.

What changed on 2026-11-18, since it changes most of the numbers below: the genome-mapped channel
was already calling sites, but its reads were being discarded afterwards, so clone sizes were built
from the 3′ channel alone. They now carry both channels. The site count went 9 → 15 and the labelled
cell count 552 → 1,195. Site labels are now pinned to loci by `metadata/integration_site_registry.csv`
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

**EVERY CLONE IS PERFECTLY CONFOUNDED WITH ITS ANIMAL, so any test comparing a site's cells against
other cells must stratify by replicate.** Each site is mono-replicate by construction - that is what
makes the hashtag a clonality control in the first place - and the same fact means a site's cells and
its mouse's cells are the same set. Tubes do not help: they are hashtag-multiplexed, so Rx_4 holds
four mice. This has already produced one wrong result. The spatial-clustering test was first run with
a tube-stratified null and reported IS05, IS10, IS04, IS03 and IS01 as transcriptionally confined;
re-run with tube x animal x depth, only IS05 and IS10 survive, and the other three were animal
effects. IS05 and IS10 are both B0304, which is 1,614 of ~20,500 called cells in Rx_4 - so the
original null was comparing a tenth-sized, separately-processed population against a background that
was 92% other animals. Multiplicity correction would not have caught it (IS04's Holm value is 0.042
and survives at 0.05); only the stratum did.

**12 sites, not 17, for any repertoire statement.** Four are excluded now, not three: IS16 joins
them, admitted by rule C′ as an integration but untested as a clone. Three of the four fail the
flank-complexity gate (`flank_gated`), so their flanks cannot support quantification: IS07 at
`JH584304.1:35,417`, IS14 at `chr10:19,804,209`, IS15 at `chr9:82,806,975`. IS14 and IS15
additionally clear the cell floor only on the two channels' reads summed (`pooled_only`) and have no
5′ junction. The fourth, IS16, is excluded for a different reason: rule C′ admits it on a
strand-matched target-site duplication, which establishes that an integration happened there and
says nothing about whether a clone did — its 3 cells are all B0301 at *p* = 0.10.

`site_evidence.csv` carries `defensible` and, beside it, `evidence_tier`. Prefer the tier. A bare
count has been quoted in this project as 15, 13, 12 and 16, each correct under a caveat recorded
somewhere else, and the drift is not carelessness but the predictable result of a one-dimensional
summary of multi-dimensional evidence. **Quote "17 across four tiers, 12 defensible"**, where `defensible` now means
flank-clean AND mono-replicate beyond chance, not merely above the floor. Derive the *reason* per site from the table too: the first attempt rendered
"IS07: NA; IS14: NA; IS15: NA" by subsetting twice.

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
dropped with no error. The floor is now whatever `assign_floor` says, and it is 6 with exact
matching required beneath 15 bases - one policy, not two constants.)

Result: cells labelled rise from 339 to **552** on the 3′ channel alone, and to **1,195** once the
genome-mapped channel's cells are folded in as well. The wild-type tubes stayed at zero through all
of that except for one called cell, the Rx_1 event at IS04: the Rx_5 cell that used to make it two
was an index-hopping barcode that Cell Ranger never called, so restricting to called cells removed
it. The measured cross-tube rate is **1 in 1,226 assignments**.

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

17 sites called — **12 defensibly**, across four `evidence_tier` values (A: 4, B: 6, C: 6, D: 1); IS07's flank cannot support quantification; IS14, IS15, IS16 and IS17 fail the mono-replicate clonality test — **1,199 distinct cells labelled**. `mouse` is the
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

`called_by` in `integration_sites.csv` carries this per site: 9 sites called by both channels, 4 by the genome-mapped channel alone, 2 (IS14, IS15) only on the pooled reads, and 1 (IS16) by rule C′. **The defensible repertoire count is 12**, and `site_evidence.csv` carries `evidence_tier` beside `defensible` so the count travels with its basis. This line read 13 for a week after the third exclusion was added at the top of this file — one place corrected, one place not, which is what an unversioned number in prose does.

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
| 1,186 → figures | 5 cells in the stroma tubes (in neither embedded panel) and 33 in no demultiplexer output, leaving 1,148 before QC. **Closed, and the 33 are gone:** see the demultiplexer row below |
| the demultiplexer's cell list | **closed.** Cell Ranger calls 161,168 cells, the demultiplexer emits 161,091 — the gap that showed as 20,502 against 20,514 in Rx_4. It is one-directional (no barcode is in the demux output that Cell Ranger did not call) and the 77 missing are **exactly** the 77 cells with zero HTO counts across all four tags: `identical(sort(zero_HTO_called), sort(dropped))` is TRUE. deMULTIplex2 cannot place a cell with no tag reads, so it drops it before classification — which is why they carry no call rather than a "negative" call. **None of the labelled cells is among them** (0 observed against 0.57 expected by chance), so nothing in the report is affected. The "33 labelled cells in no demultiplexer output" predated the called-cell filter and is now 0 of 1,199 |
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
  reproduce the pipeline's own output row for row — currently 17 sites / 1,226 assignments / 1,199
  cells, but the gate compares against what the run produced rather than against a number written
  here, which is why it survived the floor change. It failed on first run and caught two real bugs:
  a union across a molecule's reads that let a 5-base read's matches through a floor applied to its
  30-base sibling (52 spurious assignments), and a dplyr trap where `clip_len = max(clip_len)`
  followed by `clip = first(clip[clip_len == max(clip_len)])` in one `summarise()` selects an
  arbitrary row, because by the second line `clip_len` is already the scalar.
- **The gate belongs on reads, not molecules.** The assignment tier is a per-read filter, so only a
  per-read view can reproduce it to the row. The molecule table is a convenience built on top.

With it, the floor question is measured rather than argued, and the measurement moved the floor.
**It is now 6 bp with exact matching required below 15 bp**, held in one object,
`IS_DEFAULTS$assign_floor = list(min_len = 6L, exact_below = 15L)`, with a load-time assertion that
refuses an incoherent pair — tested by feeding it one. The two numbers are a single decision: under
the `L %/% 10` mismatch budget, 21 of the 35 molecules at 6 bp match a decoy, a 60% error rate, so a
6-bp floor is defensible only with exact matching beneath it. Do not treat `min_len` as adjustable
on its own.

Three stages of being wrong about this are worth knowing, because each looked like a result:

- **8 bp**, from the wrong instrument (decoys that could not have matched anyway).
- **10 bp**, from the right instrument read at too coarse a resolution.
- **"clean from 6 bp up" at 100 decoy flanks**, which was an *absence of evidence*: chance alone
  predicts only 1.67 hits across the 35 molecules at that length, so zero carried *p* = 0.19.

At `n_decoy = 1000` the measurement exists: **15 observed against 13.5 predicted at 6 bp**,
obs/exp = 1.11, which validates the matcher against the analytic 4⁻ᴸ rate rather than merely failing
to contradict it. The floor's licence is reported as an absolute count — **0.31 expected spurious
across 931 accepted molecules** against 1,226 real, 0.02 per site — because both percentage framings
invite overreading. Decoy comparisons are capped at `exact_below`, since chance matching at 30 bases
is 4⁻³⁰ and comparing 1,000 flanks against the 748 longer molecules was paying for nothing: that
capping is what kept `assign` at 3 m 52 s while the decoy count went up ninefold.

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
| IS05 | 53 | 2.75 | **92%** | 58 (56–61) |
| IS02 | 342 | 1.41 | 52% | 660 (591–757) |
| IS01 | 320 | 1.17 | 27% | 1,167 (952–1,530) |
| IS11 / IS06 / IS08 | 25 / 11 / 10 | 1.16 / 1.73 / 1.40 | 26 / 70 / 51% | 95 / 16 / 20 |
| nine others | 461 | ≤ 1.09 | — | — |

(Regenerate from `detection_model.csv` rather than editing: these moved once already when the floor
changed to 6 and the called-cell filter landed, and the prose quoting them had to be chased.)

Three things not to get wrong with it. **The molecule count must come from the clip table**, not from
summing `umis_three + umis_genomic` over the cell table — that double-counts the 83 molecules in both
channels *and* drops every sub-floor molecule, and the two errors cancel badly enough that the wrong
answer looks reasonable (mean 1.271 against the correct 1.236). **`NA` is the estimator declining**,
not a gap: at a truncated mean of 1.00 the fit is degenerate and the clone size goes to infinity.
**Read the fit from obs against exp, never from `var_mean_ratio`** — a truncated Poisson with small λ
is underdispersed by construction (0.41, 0.17), so comparing it to 1 invents a problem. The Poisson
fits: 240 IS02 cells at one molecule against 232 expected, 76 at two against 85. The overdispersion
the plan expected came from the shared-UMI artefact and disappeared with it.

**Repertoire completeness, answered without the FACS record — and the denominator was wrong once.**
The first answer used the mCherry⁺ cells *sequenced*, 43,332 across Rx_4 and Rx_6, giving a 2.6%
floor and ~7% after detection correction. That denominator assumes the sort gate is pure, and it is
not: a free, junction-independent discriminator — vector transcript from the three reporter genes in
the reference — says **65.2% of Rx_4 and 43.2% of Rx_6 called cells carry any, 53.6% pooled, against
4 of 26,822 in the wild-type control**. So the gate is about half pure.

`repertoire_completeness()` therefore reports **two** denominators: sorted cells as the upper bound,
vector-positive cells as the defensible one — a **5.3% floor and NO usable ceiling**; see §9. The
12.5% ceiling was an artefact of imputing one median detection rate ten times without an error bar,
and is withdrawn.

**What the unaccounted remainder is, and what cannot be said about it.** Of 1,422 junction-bearing
cells, 236 carry a molecule matching no called site; of the 76 long enough to place, 67 land at
MAPQ < 30 (repetitive flank), 15 do not place, and 29 sit at seven uncalled positions with 1–2 cells
each, all below `min_cells = 3`. That census is worth having. It does **not** license the conclusion
this file used to draw from it — "the missing repertoire is seven micro-loci plus the blind spot, not
a hidden population" — because a molecule census has no power over a population defined by yielding
no junction molecule. The claim spent the detection correction twice and is withdrawn. Step 4's
posteriors would not rescue it either.

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

Current counts: **17 sites, 1,226 assignments over 1,199 cells**, gate passing.

**The assignment floor is ONE OBJECT, not two constants (2026-11-18).** `IS_DEFAULTS$assign_floor` is
`list(min_len = 6L, exact_below = 15L)` and there is a load-time assertion that refuses an incoherent
pair. This is a guard, not tidiness: under the `L %/% 10` mismatch budget **21 of the 35 molecules at
6 bp match a decoy flank — a 60% error rate** — so `min_len = 6` is defensible only with exact
matching beneath it. Asking a reader not to separate the two numbers is not a guard; making them one
value is. The assertion was tested by feeding it an incoherent policy, and it rejects at `source()`.

**What licenses the floor, and what does not.** 0.31 expected spurious assignments across the 931
accepted molecules, against 1,226 real — **0.02 per site**, which changes no clone size, interval or
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

Figures here are read from the CSV named beside them. The four orphan 5' junctions were open
questions until 2026-11-19; `orphan_evidence.csv` now carries the evidence and a verdict per
junction, so they are listed under *Closed* with what each turned out to be.

1. **`Rnf17` in IS05's cells, which nothing explains.** The strongest hit in IS05's differential
   expression by a wide margin: **+4.5 log2, detected in 24% of its cells and 0 of 240 matched
   controls**, *p* = 4e-14, Holm 4e-10 across 9,479 genes tested. It is a piRNA-pathway gene with no
   known haematopoietic role, it is on a different chromosome from IS05 so it is not *cis*, and its
   *p* is unchanged by matching on transcriptional cluster so it is not a state effect. A gene absent
   from 240 matched cells and present in a quarter of a clone's is either real clone-specific
   expression or an artefact not yet identified. **Flagged, not theorised.**
2. **Direction of causation for IS05.** Its cells occupy a Vwf⁺ HSC state (see *Closed*), but a
   provirus that pushed a cell into that state and a clone that arose in such a cell predict the same
   picture. The *cis* scan rules out only the simplest version of the first. Nothing in this data
   separates them.
3. **The chr19 fusion's cells are not IS03's cells — 1 of 14 overlap.** The hashtag leg says these
   are B0303's cells and IS03 is B0303's site, but at the cell level the fusion junction and IS03's
   own channels label almost disjoint sets, and 13 of the 14 fusion cells carry no site label at all.
   Either the junction detects cells IS03's channels missed (consistent with its detection being too
   low for the ZTP) or they are B0303 cells carrying something else at that locus. **This is the one
   item flagged to come back**; it does not change the verdict that the junction is not an
   integration site.
4. **Whether IS10's, IS11's and IS03's host-gene elevation is real or residual artefact — now
   leaning strongly to artefact.** It survives subtracting every extractable junction UMI (0.8–2.4
   log2, *p* < 1e-06), and the clip table holds only molecules whose soft clip was long enough to
   extract, so the correction is a lower bound. What is new is that **the residual occurs only at
   SAME-SENSE sites** (see *Closed*), which is precisely where the artefact can operate and nowhere
   else. A genuine cis effect has no reason to track orientation that perfectly. Do not quote these
   three as cis effects. Settling it properly still needs an orthogonal measurement — host-gene
   expression in unlabelled cells of the same clone — which needs a clone membership this assay
   cannot supply.
5. **Is the antisense decrease real? One small effect and one null.** On the valid basis
   (unadjusted, since subtraction is invalid for antisense) IS02's *Ddx10* is −0.27 at *p* = 0.008
   and IS01's *Mef2c* is −0.18 at *p* = 0.9 — so one small decrease and one that is not significant,
   not two clean effects. These are the only two measurements here the artefact cannot inflate, and
   antisense insertion reducing a gene's expression is an ordinary mechanism, but two sites with one
   null is not a result and both are B0303's largest clones, so a shared state effect is not
   excluded either.
6. **Clonality proper** — site counts per compartment as a diversity measure. Still not attempted,
   and the obstacle is now better characterised rather than removed: see the completeness entry
   below. Per-compartment comparison is additionally confounded with tube (IS05 is Rx_4 only,
   IS13 Rx_3 only), so it is not a matter of running a statistic on the current table.

### Closed

- **ORIENTATION PREDICTS THE SIGN OF THE HOST-GENE EFFECT, 6 of 6 — which says what the effect IS.**
  Four host genes are inflated in their site's own cells and two are *deflated*, which was logged as
  unexplained. Put the insertions' orientations beside them and the split is perfect:

  | site | host gene | orientation | `cis_confound` | log2 FC, unadjusted | corrected on the valid basis |
  |---|---|---|---|---|---|
  | IS11 | *Stau2* | same sense | TRUE | +2.40 | **+2.45** (subtracted) |
  | IS10 | *Gm35330* | same sense | TRUE | +2.20 | **+1.69** (subtracted) |
  | IS04 | *Abhd17c* | same sense | TRUE | +1.27 | +0.65, *p* = 0.7 (subtracted) |
  | IS03 | *Hells* | same sense | TRUE | +1.00 | **+0.92** (subtracted) |
  | IS02 | *Ddx10* | **antisense** | FALSE | **−0.27**, *p* = 0.008 | −0.27 (unadjusted; see below) |
  | IS01 | *Mef2c* | **antisense** | FALSE | **−0.18**, *p* = 0.9 | −0.18 (unadjusted; see below) |

  **THE SUBTRACTION IS ONLY VALID AT SAME-SENSE SITES, and applying it everywhere manufactured a
  result.** The first version of this subtracted junction UMIs from the host gene at every site,
  including antisense ones where those reads were never counted to that gene - so it deleted real
  expression from exactly the cells under test. The damage is visible as a monotone drift: at IS02
  the "correction" drove *Ddx10* from −0.27 to −0.74 to −0.92 as more was taken away, which is the
  signature of over-correction and not of a stronger effect. It was about to be written up as the
  headline. `strip_junctions()` now subtracts only where `cis_confound` is TRUE, and
  `orientation_corrected()` picks the valid basis per site so nothing quotes one column across all
  six.

  With 4 same-sense and 2 antisense, the chance of the negatives landing exactly on the antisense
  pair is 1/choose(6,2) = **0.067**. *Suggestive, not significant*, and on six sites it cannot be
  more — but the direction is what a mechanism predicts rather than what noise predicts.

  **Same-sense:** a read-through transcript continues in the host gene's own direction, so Cell
  Ranger counts it as that gene's expression. The elevation IS the artefact, by construction — and
  `cis_confound`, derived from annotation alone with no knowledge of any of this, flags exactly those
  four. **Antisense:** the read-through runs against the gene and cannot be counted to it, so no
  inflation is possible and what remains is uncontaminated.

  *The consequence for the other open question:* the residual that survived full subtraction appears
  ONLY at same-sense sites. A real cis effect would not track orientation that way, since a provirus
  can perturb a neighbour whichever way it points. Incomplete correction is now the better-supported
  reading of IS10, IS11 and IS03.

- **THE FLOOR WAS COUNTING DROPLETS, NOT CELLS, and that was a bigger problem than its value.** A
  junction read carries whatever barcode its droplet had; the cell table, the clone sizes, the
  detection model and every completeness figure restrict to Cell Ranger's CALLED barcodes. So the
  threshold was applied to one quantity and read as another. Relaxing to two cells exposed it: three
  of the four positions admitted arrived in the site table with **one** cell and one with **none**,
  because each was one real cell plus one empty droplet. `call_sites()` and `cluster_support()` now
  take the called-barcode table and count called cells; `cluster_support.csv` reports `barcodes`
  beside `cells` so the gap is a column rather than a surprise.
- **The floor is now `>=2 cells and >=2 UMIs`, chosen from the trade-off rather than argued.** The
  old value rested on one argument about one cluster (chr11:69,817,658 has 8 of 10 UMIs in a single
  barcode) which is a statement about that cluster, not about a threshold. What settles a threshold is
  gain against loss at every value, and the loss needs a quality axis. The axis is **vector
  transcript**: a cell carrying none cannot emit a read-through junction, so the share of a cluster's
  cells carrying any measures whether it could be what it claims, measured per cell and independent of
  the junction data. On called cells, at `min_umis = 2`:

  | min_cells | sites | cells | mean vector share | single-tube | both channels |
  |---|---|---|---|---|---|
  | 1 | 19 | 996 | 0.89 | 3 | 15 |
  | **2** | **17** | **994** | **0.87** | **1** | **14** |
  | 3 | 14 | 988 | 0.85 | 1 | 12 |
  | 4-6 | 13 | 985 | 0.86 | 1 | 11 |

  No cliff, and relaxing *improves* mean vector share (0.85 -> 0.87) while adding no single-tube site;
  both admissions are seen by both channels. So the relaxed value is the defensible one and the
  evidence tier carries the difference. Tightening to 4 buys nothing.

  `min_umis = 2` is **implied, not inert**: two called cells cannot share one molecule. The only
  non-trivial setting would be 3, and the data rejects it sharply - at 3 the floor admitted the WORSE
  of the two best marginal clusters and rejected the better. chr11:69,817,658 passes a 3-UMI floor on
  10 UMIs in one cell of one tube (vector share 0.5); chr1:156,241,543 fails it on 2 UMIs across two
  tubes and both channels with 2 of 2 cells vector-positive. Selecting molecules-in-one-cell over
  cells-in-two-tubes is backwards.

  **The report's prose no longer depends on the value.** The section reads the floor from
  `IS_DEFAULTS` and the figures from `site_floor_curve.csv`, which now carries `mean_vector_share` and
  `min_vector_share` per grid point, so changing the floor changes the text.
- **`defensible` now requires clonality to be demonstrated, because the floor had been standing in
  for it.** Once chrX cleared the relaxed floor directly, every flag that had excluded IS16 switched
  off - not `pooled_only`, not `tsd_rescued`, not `flank_gated` - and it became defensible because the
  THRESHOLD moved. Its three cells are all B0301 against that replicate's 46.8% share, *p* = 0.10,
  which establishes nothing. `site_evidence()` now runs `replicate_confinement()` per site and
  `defensible` requires *p* < 0.05. It costs nothing where it matters (every other defensible site
  passes; IS13 at *p* = 2.8e-05 on 7 of 7) and removes exactly the one that should not have been in.
  **The defensible count is 12 before and after the floor change** - which is the test of the fix.
- **A floor change also revealed a lost denominator.** With chrX clearing the floor, rule C' promotes
  nothing, and `tsd_rescued_clusters()` took an early return that skipped setting the candidate count
  - so the driver printed "  orphan 5' junction(s) tested, 0 promoted" with the numbers blank. 0 of N
  is only a statement with N beside it, and N went missing exactly in the case where it matters most.
  Fixed; `tsd_chance_rate()` also now coerces an empty denominator to 0 rather than emitting
  zero-length columns.
- **Rule C' now promotes nothing, and that is the right state for it.** chrX:100,691,113 clears the
  relaxed floor on its own, so the rule fires on no cluster in this run. It is kept, with its warrant
  measured (2.05e-07 genome-wide, 5.4e-07 against the transcribed footprint), against the next run.
- **IS05's confinement is explained, and it is not *cis*.** Matched differential expression - animal,
  tube and depth matched, NOT cluster, since cluster is the mediator - gives **8 genes at Holm < 0.05
  over 9,479 tested**. The up-regulated set is *Vwf*, *Clu*, *Clec11a*, *Selenom*, *Neo1*: the Vwf+
  platelet-biased long-term HSC programme, which is coherent for a clone in a myeloproliferative model
  rather than a grab-bag. **Nothing in +/-1 Mb survives** (21 of 45 window genes testable, best
  *Dennd4a* at Holm-within-window 0.11), agreeing with the `cis_confound = FALSE` flag reached from
  annotation alone.
  *The split is the answer:* **4 of the 8 survive matching on cluster** (*Rnf17*, *Neo1*, *Selenom*,
  *Vwf*) and 4 do not (*Clec11a*, *Clu*, *Plek*, *Egr3*). So the confinement is part cell-state - 42%
  of IS05's cells sit in one cluster against a 24% background, and that cluster is Vwf-positive - and
  part within-state. What it does NOT establish is direction; see the open list.
  *A design note worth keeping:* `matched_de()` draws the control set ONCE and tests every gene
  against it. Redrawing per gene adds sampling noise across genes and makes the multiplicity
  correction wrong, because each test would be against a different background.

- **ONE site is transcriptionally confined: IS05. Four other results were withdrawn getting there.**
  The question was open since IS05's cells looked clustered on the HSC UMAP. Three successive
  versions of the test gave three different answers, and the two wrong ones are recorded because the
  errors are reusable.

  *Version 1, tube-stratified null:* IS05 13.8, IS10 11.4, IS04 3.4, IS03 2.5, IS01 2.2 — five
  "confined" sites. Wrong, because tubes are hashtag-multiplexed and every site here is
  mono-replicate, so **every clone is perfectly confounded with its animal** and no tube stratum can
  separate them. IS05 and IS10 are both B0304, 1,614 of ~20,500 called cells in Rx_4.

  *Version 2, tube × animal × depth-decile:* IS05 6.5, IS10 5.1 (Holm 0.009 both); IS04 1.6, IS03
  2.0, IS01 1.4 — **three results were animal effects.** Multiplicity would not have caught them:
  IS04's Holm value under the tube null is 0.042 and clears 0.05. The stratum did.

  *Version 3, excluding the shared B0304 cluster:* both survivors over-occupied one transcriptional
  cluster, so the test was re-run with that cluster removed from member set and null pool alike.
  **IS05 survives at z = 4.4 on 28 of its 48 cells (Holm 0.048); IS10 collapses from 5.1 to 1.2
  (Holm 0.63).** IS10's enrichment WAS the shared cluster. It is demoted, and the
  two-survivors-one-mouse worry resolves by one of them not being a finding.

  **So: IS05 only.** Confined beyond its animal, its tube, its depth and the shared B0304 state; 28
  cells after the exclusion, which is the denominator to quote. `cis_confound = FALSE` and its host
  gene is absent from the PCA's variable features, so it is also the one site immune to the junction-
  read artefact below.

  *Two side notes.* IS02 — the largest clone, 137 HSC cells — sits at z = −0.3, indistinguishable
  from its own animal at matched depth. State that as a fact about IS02; it is n = 1 for any claim
  relating clone size to confinement. And B0304 as a whole IS enriched for itself (z = 47.7 with no
  clone label) — which does **not** bear on the result, because the animal-stratified null draws from
  B0304, so the animal effect is absorbed into the null by construction. That is why IS05's null rose
  0.0018 → 0.0048 when the stratum changed. Measuring the animal effect outside the null confirms the
  version-1 error and touches nothing else.

  *And the ranking statistic was wrong in a smaller way.* It was `(e - n)/d + n`, the right statistic
  plus a site-specific offset, and unreconstructable from its own printed inputs — two plausible
  guesses give different winners. Now `excess_per_detection = (enrichment - null)/detection`, stated
  in the table. Subtract first, then divide: the signal scales with d, the null does not. **It is NA
  for IS10 (no fitted λ), so the survivors were never rankable against each other anyway.** Holm is
  pre-registered and printed alone; Bonferroni was removed because Holm dominates it and showing both
  at the one site where they differ hands the reader the decision after the fact.
- **The junction-read *cis* artefact is large, cannot reach the clustering statistic, and is only
  partly removable.** A junction read overlaps the host gene, so at an intragenic same-sense site the
  junction molecules are counted as host-gene UMIs in exactly that site's cells. Against
  animal/tube/depth-matched controls (NOT cluster-matched — cluster is a mediator):

  | site | host gene | no adjustment | genome-channel UMIs removed | all junction UMIs removed | in the PCA? |
  |---|---|---|---|---|---|
  | IS10 | *Gm35330* | +2.20 | +2.08 | **+1.69** (5e-07) | no |
  | IS11 | *Stau2* | +2.40 | +2.25 | **+2.45** (2e-06) | yes |
  | IS04 | *Abhd17c* | +1.27 | +0.97 | **+0.65** (*p* = 0.7) | no |
  | IS03 | *Hells* | +1.00 | +0.98 | **+0.92** (1e-08) | yes |

  **It cannot explain any clustering**, and the reason is one column: the PCA is built on 3,000
  variable features and the relevant host genes are not among them. Check membership before
  rebuilding a counts matrix — it is cheaper and more direct.

  **IS04's effect is entirely artefactual** (p = 0.9 after full subtraction). The other three survive
  at 0.8–2.4 log2, and that is where this data stops: the clip table holds only molecules whose soft
  clip was long enough to extract, so a read-through below that threshold is counted by Cell Ranger
  and not subtracted. Genuine cis effect and incomplete correction are **not separable here**, because
  the artefact and the label are the same event — every cell of a site is by definition a cell that
  yielded a junction molecule there. One route is closed for good: 270 UMIs against a median
  nCount_RNA of 13,394 is under 0.01% of library size, so no *other* gene is affected.
- **The *Klf6* hypothesis is unsupported, and the PC-loading screen that produced it is deleted.**
  IS10 sits 239 kb from *Klf6* — enhancer range, haematopoietic TF, MPN-to-AML model. Matched DE
  without cluster matching gives **+0.26 log2, raw *p* = 0.039, Holm 0.11 across the 6 genes in ±1 Mb
  that are detectable** — below the 0.3 log2 bar and not significant corrected. The off-chromosome
  control is +0.01 at *p* = 0.81, so **the design resolves nothing below ~0.2–0.3 log2**; quote that,
  not the point estimate. The window scan's only Holm-significant gene is *Gm35330*, IS10's own host
  gene, i.e. the artefact. So nothing in reach is a clean candidate — and nothing needs to be, since
  IS10's confinement did not survive version 3.
  *An earlier run of this test matched on cluster and returned +0.09, p = 0.36, with the control at
  −0.33, p = 0.06.* Cluster is the mediator; matching on it bounds only the effect net of state.
  Dropping it raised the estimate and cleaned the control, which is the signature of having
  conditioned on a mediator. `cis_loading_check()` is removed: PCs are global, and it selected the
  gene it then tested.
- **The splice generalisation, reported as the statistic that has power rather than the one that
  does not.** Zero of 28 discovery clusters sit on an exon end with a GT donor — and that is **not
  evidence**, because 0.28% of the genome lies within ±5 bp of an exon end, so 0.08 is expected and a
  zero is what chance predicts. This family cannot be powered at 28 clusters (a bare GT scan expects
  ~1.75 chance hits). The powered measurement runs the other way: **1 of the 4 five-prime orphans is
  a splice artefact**, *p* ≈ 0.01 for ≥1 against the same background. So the 5′ channel's exclusion
  rests on mechanism plus a **measured contamination rate of roughly 1 in 4** — not on a measured
  absence of phantoms in the trusted channels, which is unmeasurable at this scale. Run no further
  variants. Nothing here has power over an unannotated donor.
- **The chr19 fusion: three legs, not four, and they are not independent.** It was written up as
  "four predictions, all of which held". Priced: the exon-boundary leg is load-bearing (*p* ≈ 0.003
  for one position); **the "exactly one `AG|ATTC` in the vector contig" leg is worth nothing** — the
  contig is 3,236 bp, so a 6-mer is expected 1.58 times counting both strands and P(≥1) = **0.79**;
  the GT donor is a consistency check, since ~99% of introns begin GT given an exon end. It cannot be
  strengthened either: all 32 reads carry exactly four bases beyond the junction, and four bases
  cannot pin an acceptor in 3.2 kb. **And the two surviving sequence legs are not independent of each
  other** — exon boundary and IS03-in-the-same-intron are both statements about position within
  *Hells*. The independent leg is hashtag identity: 13 of 13 non-multiplet cells are B0303 against
  that replicate's 0.382 share, *p* = 3.6e-06, using no annotation, alignment coordinate or sequence.

  **The cell-level check came out the wrong way and is reported, not explained.** If these were
  transcripts of the IS03 provirus, IS03's own channels should have labelled most of the same cells.
  **Overlap is 1 of 14**, and only that one fusion cell carries any site label at all. The reading
  that fits is that the fusion junction detects cells IS03's channels missed — consistent with IS03's
  detection being too low for the ZTP — but that is an interpretation of a 13-cell discrepancy, and
  "these are B0303 cells carrying something else at that locus" is not excluded. What it does settle
  is narrower and sufficient: admitting the junction would have created a sixteenth clone out of a
  transcript whose cells are almost all unlabelled, rather than double-counting IS03. **Bring this
  one back.**
- **The sub-clone fraction was quoted off a stale stratum.** §*The discriminator that uses detection
  depth* said 0.105 against 0.493 at k = 2 and "a quarter to a third carry both". The depth variable
  had been corrected to DISTINCT molecules and the prose had not: 0.105 is 12/114, the fingerprint of
  the superseded `umis_three + umis_genomic` denominator, which pushed single-molecule cells into
  k = 2. The table itself was already right. Corrected and now computed inline: **12 of 101, observed
  0.119 against 0.424 expected, ratio 0.28 at k = 2; 0.44 against 0.64, ratio 0.69 at k = 3** — and
  **at k = 4 it falls back to 0.11 on 9 cells**, which a dose-response cannot do and noise at that
  count easily can. The prose had stopped at k = 3. Fisher at k = 2 is unchanged at *p* = 0.00056.
  The honest range for the sub-clone fraction is a quarter to two-thirds on strata of 101, 16 and 9
  cells, which leaves its size essentially unconstrained.
- **`chrX:100,691,113` is IS16, admitted by a rule rather than by hand.** The floor is applied
  WITHIN a discovery channel, so chrX failed at 2 cells in the 3′ channel while IS14 and IS15, which
  reach 3 on the two discovery channels summed, passed. Pooling all three channels is not the fix:
  the largest 5′ orphan here is the chr19 *Hells* fusion at 14 cells, and a pooled-across-everything
  floor would admit it as a site larger than six real ones. **Rule C′** admits a sub-floor cluster
  only when a 5′ junction supplies a 3′ partner at signed Δ = ±4 *whose sign agrees with the 3′
  alignment's strand* (`strand_ok`, the same test `tsd_check()` applies) and the pooled count then
  reaches the floor. The warrant is not that this rejects the *Hells* fusion — that is a consequence
  of the rule, not its basis — but that a TSD is a signature of integrase 3′ processing and of
  nothing else: splicing, mismapping and chimeric reads can each put two features 4 bp apart, none
  can make the sign agree with the strand. And the specificity is computed, not asserted:
  `tsd_chance_rate.csv` gives 2.05 × 10⁻⁷ expected chance admissions, from the genome read out of the
  reference's own `.fai` (2.728 Gb), 28 discovery clusters, a 5 bp accepting window and 4 orphans
  tested — and against the transcribed footprint rather than the whole genome, which is the honest
  denominator since lentivirus prefers gene bodies, **5.4e-07** (the union of protein-coding gene
  bodies is 1.04 Gb, 38% of this reference, so x2.6; it had been guessed at "5-10x, so ~1e-06", and
  gene bodies are too large a share of the mouse genome for that).

  **The sequence half is NOT independent confirmation and that claim is withdrawn.** `tsd_obs = 1`
  says the duplicated span appears in every 5′ read — but those reads *align* to that span, so their
  bases match the reference there by construction of a clean alignment. `tsd_three_way()` now
  compares the 5′ clip's host tail, the 3′ clip's host head and the reference, which is the form
  originally specified; it is still not independent evidence, because both clips align to the same
  span and their agreement is the Δ = 4 relationship restated in bases. It is a check on the
  pipeline — and a worthwhile one, since this analysis has had three coordinate bugs — not evidence
  about biology. **IS16 rests on the coordinate relationship and its chance rate alone**, which still
  clears at 5e-07, so the decision is unchanged. **IS16 is tier D and NOT defensible.** Its 3 cells are all
  B0301, which is 46.8% of called cells in those tubes, so *p* = 0.10. C′ supports that an
  integration happened, not that a clone exists.
- **The 12.5% completeness ceiling is withdrawn.** It came from filling in detection at the ten
  sites where the ZTP is degenerate with the MEDIAN of the six that fit — a model (detection is
  exchangeable across sites), unstated, with no uncertainty, fitted on a sample selected for high
  λ. Replaced by a regression of log λ on log mean vector expression, which is measured per cell
  independently of the junction data. The slope is **1.04 ± 0.72** — proportionality, exactly what
  the mechanism predicts, since a junction read is a read-through of the vector transcript — but
  *p* = 0.22 and R² = 0.34 on six points, so the covariate is not distinguishable from no covariate.
  The prediction intervals span two orders of magnitude and summed across sixteen sites put the
  corrected count ABOVE the denominator (`ceiling_uninformative = TRUE`). **This data bounds
  completeness from below only. Quote the 5.3% floor and no ceiling.**
  *And Chao1 is not the alternative:* `min_cells = 3` discards every singleton and doubleton before
  the site table exists, so f₁ = f₂ = 0 BY CONSTRUCTION and Chao1 and ACE are undefined rather than
  biased. What this assay has that ecology does not is a measured detection process, so completeness
  is reported stratified by clone size instead: P(discovered) = 1 − pbinom(`min_cells` − 1, N, d),
  in `discovery_power.csv`. At d = 0.39 a clone of 3 is found 6% of the time, of 5 30%, of 10 82%,
  of 20 99.5%. **So the sixteen sites are close to the complete set of large clones, and the 84%
  unaccounted share needs no extra hypothesis** — a repertoire of thousands of small clones is
  invisible here by arithmetic. Silent proviruses remain an equally consistent explanation; the FACS
  mCherry⁺ fraction per animal still separates them.
- **`chr19:38,932,385` is not an integration site — it is a spliced host-vector fusion.** Four
  predictions, all of which held. The junction sits on the LAST BASE of a *Hells* exon
  (`exon_end_dist = 0`), not near it. The intron after that base begins **GT**. The bases beyond the
  junction are the same four in every read that has any — **ATTC** — and they are not the next
  *Hells* exon, which begins ATGCTT. The vector contig contains exactly one `AG|ATTC`, an acceptor
  at **position 875**. And IS03, at chr19:38,933,742, lies **inside that very intron**
  (38,932,386–38,933,834), 92 bp before its acceptor. The independent instrument agrees: 13 of 13
  non-multiplet cells are replicate **B0303**, IS03's own mouse, *p* = 3.6 × 10⁻⁶. So these reads
  are a transcript from a provirus already called, not a second integration.
  *A side result worth more than the question:* only **1** of those 14 cells is assigned to IS03, so
  the splice junction detects cells IS03's own channels missed — independent evidence that IS03's
  per-cell detection is low, consistent with its clone size being inestimable.
  `chr19:38,932,221` (1 read, 1 cell) carries the same ATTC and is the same fusion, sub-threshold.
- *(superseded by the IS16 entry above — kept because the reasoning that was WRONG here is worth
  not repeating.)* This entry said chrX was "the same situation as IS14 and IS15, which are called".
  It is not. IS14 and IS15 reach 3 cells on the two DISCOVERY channels summed; chrX reaches 3 only
  by counting the 5′ channel, which is never allowed to call. The axis on which chrX looked stronger
  (it has a TSD, they do not) is not the axis the floor uses, and arguing from the wrong axis is how
  a hand-promotion nearly got made. Rule C′ admits it on the TSD explicitly, which is a different
  and defensible thing.
- **`chr14:62,948,717` is unsupported.** Four reads, intergenic, no soft clip at all (32M and 26M,
  both ending exactly at the junction), no 3' partner, and its 3 non-multiplet cells split across
  mice (2 of 3 in B0304, *p* = 0.49) — which a single clone cannot do.
- **The `≥3 cells and ≥3 UMIs` floor is really a one-part floor, and the second part is inert.**
  `cluster_support.csv` is every discovery cluster with its support and the floor NOT applied;
  `site_floor_curve.csv` is the site count over a 6×6 grid. Sites admitted, rows `min_cells`,
  columns `min_umis`:

  | | 1 | 2 | 3 | 4 | 5 | 6 |
  |---|---|---|---|---|---|---|
  | **1** | 28 | 21 | 16 | 14 | 14 | 14 |
  | **2** | 20 | 20 | 16 | 14 | 14 | 14 |
  | **3** | 15 | 15 | **15** | 13 | 13 | 13 |
  | **4** | 13 | 13 | 13 | 13 | 13 | 13 |

  Three things fall straight out. **`min_umis = 3` rejects nothing** that `min_cells = 3` has not
  already rejected — 1, 2 and 3 all give 15 — so the advertised two-part floor is one part plus
  decoration. **`min_cells` is the load-bearing axis**, running 28 / 20 / 15 / 13 / 13 / 13 as it
  goes 1 to 6, and the count is still falling at 3 and flat from 4. And **the floor costs almost no
  cells**: coverage goes 1,040 → 1,022 from (1,1) to (3,3), so the threshold is not about how many
  cells get labelled but about how many thinly-supported positions get a name.

  That reduces the whole question to `min_cells` 3 against 2, which is worth exactly one cluster:
  **chr11:69,817,658** — 11 reads, **10 UMIs**, 2 cells, 1 tube, both channels at MAPQ 255. More
  molecular evidence than most called sites, and rejected. But 8 of those 10 UMIs sit in a single
  cell, a stroma-panel barcode, so it is one integration observed repeatedly rather than a clone -
  and UMIs within one cell are not independent evidence that an integration is present in the
  population, which is the job `min_cells` exists to do. **Leave the floor at 3.** The case for
  simplifying it to `min_cells` alone is cosmetic and would change no result.
- **The MNC/HSC QC asymmetry is a depth bias meeting two different QC windows.** Labelled barcodes
  survive at **81.3% in HSC and 58.8% in MNC**, while overall survival is near-identical (69.2% and
  67.6%) — enriched among survivors in one compartment, depleted in the other, which cannot be a
  property of the clones. It is a property of what it takes to *detect* one: a junction read needs
  depth, so labelled cells are systematically deeper (HSC median 4,725 genes against 3,728; MNC
  3,494 against 2,104). Per-criterion failure rates, labelled against unlabelled, say the rest:

  | criterion | HSC labelled / all | MNC labelled / all |
  |---|---|---|
  | below `nFeature_min` | **3.4% / 23.5%** | 1.0% / 10.2% |
  | above `nFeature_max` | 10.8% / 3.0% | 9.1% / 3.0% |
  | outside the ribosomal window | 1.3% / 7.6% | **18.3% / 10.6%** |
  | multiplet | 7.7% / 6.5% | 17.9% / 14.5% |

  HSC sets `nFeature_min` at 2,500, so nearly a quarter of all cells fail it and labelled cells
  almost never do: that one term is worth ~20 points and swamps everything else. MNC sets it at
  1,200, where the same advantage is worth only ~9, while labelled cells pay ~8 on the ribosomal
  window (1–20%), ~6 on the upper gene bound and ~3 on multiplets. The sign of the net flips. So
  **nothing is wrong with the clones or the labelling** — the two compartments' QC windows simply
  interact differently with a depth-biased subset. Whether MNC's ribosomal window should be revisited
  is a QC question, not an integration-site one.
- **ψ artefact origin** — packaged vector genomic RNA (§6). Priming, 2-LTR circles and plasmid all
  refuted; the backbone scan was decisive.
- **IS07 / JH584304.1** — called, but excluded from quantification: its flank fails the complexity
  gate (`flank_gated`), so it is one of the three sites outside the defensible 12 (§7).
- **Multi-copy founders** — no. Every informative pair of sites co-occurs *below* chance (IS01+IS02:
  6 observed vs 26 expected in Rx_4; 2 vs 24 in Rx_6), which is the signature of disjoint clones
  rather than independent integrations. The ~2% of cells at more than one site read as doublets.

## 10. Files

| | |
|---|---|
| `src/R/integration-sites.R` | all logic, as functions. Attaches no packages — `source()` it and step through. All thresholds in `IS_DEFAULTS`. |
| `scripts/R/call-integration-sites.R` | driver, five stages (`extract`, `genomic`, `genomic-decoy`, `sites`, `assign`), separate because bwa runs between them and because a control must not sit on the critical path of a run someone is waiting for |
| `scripts/R/export-hashtag-calls.R` | exports the demultiplexer's calls to `work/hashtag_calls.csv.gz`, which is what lets the report resolve replicates without loading a 2 GB Seurat object |
| `scripts/R/record-output-manifest.R` | records existing outputs in `manifest.tsv`. Not part of a run — the pipeline writes it as it goes. For outputs that predate it, or a lost manifest. |
| `scripts/slurm/analysis/integration-sites.sh` | the job. 1 h 49 m for the full five stages (job 35742986); the genome-mapped scan dominates. |
| `scripts/slurm/setup/bwa-index.sh` | one-off index build, ~50 min |
| `notebooks/integration-sites/integration-sites.qmd` | the report, registered in `_targets.R` |
| `output/objects/integration-sites/*.csv` | 25 output tables plus `manifest.tsv`; `work/` holds the per-tube genomic cache, the called-barcode list and the hashtag calls. Three are new: `tsd_chance_rate.csv` (the warrant for rule C′), `detection_regression.csv` (λ predicted from vector expression, replacing the median imputation) and `discovery_power.csv` (completeness stratified by clone size) |
| `data_hub/10x_resources/reference_additions/README.md` | contig derivation and verification — **outside git, not version controlled** |

```bash
sbatch scripts/slurm/analysis/integration-sites.sh
tools/quarto.sh render notebooks/integration-sites/integration-sites.qmd   # renders AND archives
```

Render from `tools/quarto.sh`, not `env/bin/quarto`: conda's quarto needs the variables its
activation script sets, and `./env` is never activated here. A bare `quarto render` fails three times
in a row on three different missing tools, with errors that read like a broken install rather than a
relocated one. The wrapper also archives the result, which matters because
`output/reports/latest/` is overwritten on every render and the dated copy under
`output/reports/archive/` is the only permanent one.

Not done: `tools/catalogue-data.sh` has not been run on the new Cell Ranger inputs.

## 11. Provenance — what is verified, and what that is worth

Every number in the report has to be reproducible from data and code that is available, which means
something has to check that the files behind those numbers are the pipeline's own output. Three
layers do, and they catch different failures:

1. **`write_csv_checked()` — the write happened in full.** Counts the lines back and errors on a
   mismatch. This exists because the fileset filled and writes began *truncating* silently (§8); do
   not remove it.
2. **`manifest.tsv` — the file is still the one the pipeline wrote.** `record_in_manifest()` records
   size, md5, timestamp and the full commit SHA in the same call that line-counted the write. The
   report's provenance banner compares each declared input against it: `matches pipeline manifest`
   is a pass, `DIFFERS from pipeline manifest` is a red banner, because every number downstream of
   an unaccounted file is then computed from something nobody can account for. A re-run stage
   replaces its own rows and leaves the others.
3. **The declaration in `_targets.R` — the list of files being checked is complete.** The report's
   target lists *every* CSV the report reads, not a representative few. A file left off that list is
   a file whose numbers reach the report unverified, and the report also stops re-rendering when it
   moves — a stale report that looks current.

Two cautions. The manifest records the **full 40-character SHA** deliberately: a short hash is often
a valid number, and `8e55518` parses as 8 × 10⁵⁵⁵¹⁸, i.e. `Inf`, which is how the first printout of
it came out. Read that file with every column as character. And
`scripts/R/record-output-manifest.R` attests **less** than the pipeline does — it records whatever is
on disk now, so it can bless a file that was already wrong. Re-run the pipeline rather than that
script if there is any doubt.

What is *not* verified: the banner's "Pipeline" check reads amber because the report has never been
built by `tar_make()`. Building it would rebuild twenty upstream targets, the whole Seurat chain
among them, so the inputs are read from the target's declaration instead
(`declared_inputs()`). That is an honest amber, not a missing feature — but if the Seurat chain is
ever rebuilt for other reasons, build the report target in the same run and the check turns green.
