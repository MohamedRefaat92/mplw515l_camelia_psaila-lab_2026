# Integration sites — open information gaps

Facts about how the experiment was done that the data cannot supply, each of which changes a specific
decision in the analysis. Recorded here rather than in a conversation so the dependency is visible to
whoever picks this up. See `HANDOFF.md` for the analysis itself.

Status key: **open** — needed, not yet answered. **answered** — with the date and what changed.

---

## 1. Time between transduction and harvest — **open**

**Why it matters.** It decides whether unintegrated vector forms are still plausible. Linear vector
DNA, 1-LTR and 2-LTR circles persist for days after transduction and then dilute out with division.
If harvest was weeks later they are effectively gone; if days, they are live hypotheses for the ψ
coverage in `HANDOFF.md` §6, which is currently "refuted priming hypothesis, mechanism unresolved".

**What changes.** Whether the Phase 5 ψ tests are worth running at all, and how to read a positive.

**How it could be answered without asking:** the unprocessed-terminus discriminator. A clip anchored
at contig **2734** retains the `AC` that integrase removes, so it cannot come from an integrated
provirus — it is unintegrated linear vector or plasmid. The ratio of 2734 to 2736 anchors is a direct
estimate of how much unintegrated material is left. Cheap, and worth running before asking.

## 2. Intended MOI, or measured vector copy number — **answered from the data; useful only as confirmation**

**Answered from the data, 2026-11-18 — and in the opposite direction to the first attempt.** This
entry said on 2026-11-16 that multi-copy founders are *absent*, on the grounds that every informative
pair of sites co-occurs below chance. That was wrong, and the reason it was wrong is instructive: the
expectation it compared against assumed complete detection, so "below chance" was guaranteed and
carried no information.

The hashtag demultiplex supplies the control that does carry information. Every site is
mono-replicate, so **two sites in different mice cannot share a cell** and any co-occurrence between
them is a doublet or an ambient pick-up. Across Rx_4 and Rx_6: cross-mouse pairs give 1 observed
co-occurrence against 263.5 expected under independence; same-mouse pairs give 26 against 202.1. A
doublet process cannot tell the two classes apart, so applying the same-mouse rate to the cross-mouse
pairs predicts 33.9 and one was seen (*p* ≈ 7 × 10⁻¹⁴). **Doublets and ambient are excluded as the
dominant mechanism, and multi-copy founders in B0303 are now the favoured explanation.** All 26
co-occurrences are in that one mouse, among IS01, IS02, IS06 and IS08.

**The objection that would have killed this, tested.** A droplet holding cells from two mice carries
two hashtags, which is what a demultiplexer calls a doublet — so if such droplets had been dropped
before the test, cross-mouse pairs would have been stripped along the exact axis of the comparison.
They were not: the test runs on every labelled barcode, multiplets included (11.6% of them), and the
one cross-mouse event is itself a called multiplet. Within the multiplet-enriched subset, cross-mouse
co-occurrence is still 36-fold below independence.

**Detection depth then separates the two readings that co-occurrence counts cannot.** Conditioning on
how many junction molecules a cell yielded, the fraction carrying both IS01 and IS02 goes 0 → 10.5%
→ 31.8% with depth, while the matched cross-mouse control stays at exactly 0 at every depth (12 of
114 against 0 of 91 at two molecules, Fisher *p* = 0.0014). Doublets produce no such gradient. But
the gradient reaches only about a quarter of what full sharing predicts, so the reading is a
**double-insertion sub-clone inside B0303** rather than a single founder — which also explains the
inconsistent implied sizes (2,312 from IS01+IS02 against 365 from IS01+IS06) that co-occurrence alone
could only record as a caveat.

MOI is therefore no longer the prior the co-occurrence test needs — the test stands on its own
control. It would still be useful as outside confirmation, to say whether this was by design, and to
bound what fraction of a clone can carry two proviruses.

**Why it mattered.** It set the prior for the site co-occurrence test. At low MOI almost every
transduced cell carries one provirus and co-occurrence should be rare; at high MOI multi-copy
founders are expected, and two sites in one cell is then a clone rather than a doublet. The data now
answers that question directly, so the prior is no longer load-bearing.

**What changes.** How to read the `ambiguous` flag in `cell_site_assignments.csv`. The pipeline
deliberately refuses to decide — it flags cells at more than one site and leaves the call downstream,
because a multi-copy clone and a doublet are indistinguishable from the data alone. MOI is the
outside information that breaks the tie.

Also bears on replicate B0303 carrying seven of the fifteen sites and 819 of the 1,186 labelled
cells: that is either a more heavily transduced animal or one with more clonal expansion, and MOI
distinguishes them. The co-occurrence result above says the two are not mutually exclusive — the
sharing between IS01, IS02, IS06 and IS08 is in that same animal.

## 3. Is genomic DNA available from these animals? — **open**

**Why it matters.** A bulk integration-site assay (LAM-PCR, LM-PCR, or similar) on the same gDNA
would give an independently derived site list.

**What changes.** It would convert the hardest remaining problem into a solved one. The
`≥3 cells and ≥3 UMIs` floor exists only because nothing independent confirms a weak site; with a
trusted whitelist the floor could be dropped entirely and the mini-reference seeded from it rather
than from this assay's own discovery tier, removing the last circularity. It would also settle IS09,
chrX:100,691,117 and chr14:62,948,717 outright.

Lower priority than it was: the 5′ TSD check now provides *some* independent confirmation. But TSD
only reaches sites with 5′ reads, which is four of nine.

## 4. Has Cre been applied to any sample? — **open**

**Why it matters.** Both loxP sites flank the entire transcription unit. Excision removes the vector
transcript, and with it the read-through that this whole analysis depends on — while the integration
site itself persists in the genomic DNA.

**What changes.** If any sample saw Cre, absence of a clone label in that sample means nothing, and
cross-sample comparisons of clone representation are invalid for it. This would be a hard
confounder, not a caveat: it breaks the link between "carries the provirus" and "detectable here".

Nothing in the current results suggests it has been applied — the mCherry⁺ fractions label normally —
but it has not been ruled out, and it would be silent if it had.

---

## Answered

### Per-animal, per-tube cell counts — **computed 2026-11-18, no longer a query**

Two of the outstanding metadata requests were for the hashtag composition of each tube. It is in the
targets store already, in `demux_res_hsc`, `demux_res_mnc` and `demux_res_stroma` as
`final_assign`, keyed on the Seurat cell name. The table, over every barcode the demultiplexer
called:

| tube | B0301 | B0302 | B0303 | B0304 | multiplet | negative |
|---|---|---|---|---|---|---|
| Rx_1 Stroma WT | 1,409 | 1,620 | 1,103 | 1,425 | 77 | 749 |
| Rx_2 Stroma Het | 3,077 | 4,714 | 1,420 | 5,244 | 659 | 3,369 |
| Rx_3 HSC Het mCh⁻ | 2,258 | 630 | 4,989 | **15,291** | 1,386 | 511 |
| Rx_4 HSC Het mCh⁺ | 11,401 | 277 | 5,847 | 1,614 | 1,342 | 21 |
| Rx_5 HSC WT | 5,089 | 5,885 | 6,599 | 5,552 | 1,974 | 1,712 |
| Rx_6 MNC mCh⁺ | 6,721 | 657 | 11,710 | **201** | 3,270 | 252 |
| Rx_7 MNC mCh⁻ | 6,959 | 2,229 | 4,735 | 5,526 | 3,543 | 262 |
| Rx_8 MNC WT | 3,097 | 4,076 | 3,486 | 4,128 | 2,478 | 517 |

**This table overturned a result.** B0304 contributes 20,817 cells to the two mCherry⁻ tubes and
1,815 to the two mCherry⁺ ones — 92% sorted negative, against ~34% for B0301 and ~36% for B0303. The
report had read B0304's labelled-cell composition as a leakage finding; normalised by these
denominators the animal effect in the mCherry⁻ fraction is a factor of 2.4, not 70, and relative to
each animal's own mCherry⁺ rate B0304 leaks the least of the three. Print the denominators before
dividing.

What is still outstanding from the metadata requests is the third one: **the mCherry⁺ fraction at
sort, per animal**, which is a number from the FACS record rather than from the data. It is now the
input to a detection-rate cross-check rather than to the clone-size estimate, since the
zero-truncated Poisson supplies that from the molecule counts directly.


### Read 2 length — **answered 2026-11-16**: 150 bp, uniform.

Measured from the BAMs. Fixes the 5′ channel's clip-length ceiling at 125 bp (150 minus a 25 bp
minimum anchor) and makes clip length a proxy for insert size over the detectable band.

### Are the BAMs paired-end? — **answered 2026-11-16**: no.

0 reads flagged paired, 0 with non-zero TLEN. 10x 3′ v3 puts barcode and UMI on R1 and the insert on
R2, and Cell Ranger aligns R2 only. This killed a proposed gate that would have measured the insert
size distribution from the BAMs; it is not recorded there.

---

## How to use this file

When one is answered, move it to **Answered** with the date and say what it changed. If an answer
makes a recommendation in `HANDOFF.md` or the task list wrong, fix that too — a stale recommendation
is worse than an open question, because it looks settled.
