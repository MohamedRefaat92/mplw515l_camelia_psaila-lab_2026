# Two parameter decisions needing guidance

Both arose from implementing task T1 of the 2026-11-18 report review. Every figure here is read from
`clip_floor_curve.csv` or a timed run, not from notes. Written against commit `4279063` plus the
uncommitted T1 work.

Neither is blocking: the pipeline runs, the gate passes, 15 sites / 1,222 assignments / 1,195 cells.
Both are cases where I can state the trade-off precisely and have no principled basis for choosing.

---

## Decision 7 — `n_decoy = 1000` puts a control on the critical path

### What changed and why

T1 said the floor decision rested on about one event, and offered either cutting 1,000–10,000 decoys
or taking 7 bp and stopping. `n_decoy` went 100 → 1,000. It was the right call and it worked: at 100
flanks the curve showed **zero** decoy hits at 6 bp, which read as clean and was not — chance alone
predicts only 1.67 hits across the 35 molecules at that length, so a zero carried *p* = 0.19.

At 1,000 flanks the measurement exists: **15 observed against 13.5 predicted at 6 bp**, so chance
matching at the floor is measured and agrees with the analytic 4⁻ᴸ rate to 11%. That turned the
floor's justification from an underpowered absence into a validated rate.

### What it costs

`sites` + `assign` went from **~2 minutes to 14 m 13 s**. The matcher compares every clip against
every flank in both orientations, so cost is roughly linear in flank count: 112 → 1,012 flanks is
~9×. It now runs on **every** invocation of the pipeline, and `assign` is the second-slowest stage
after the genome-mapped scan.

### The power/cost curve, at 6 bp (35 molecules)

| `n_decoy` | per-molecule chance | expected hits | *p*(observe 0) |
|---|---|---|---|
| 100 | 4.8% | 1.7 | 0.19 — no power |
| **300** | 13.6% | 4.8 | **0.008** |
| **1,000 (current)** | 38.6% | 13.5 | 1.3 × 10⁻⁶ |
| 3,000 | 76.9% | 26.9 | 2.1 × 10⁻¹² |

300 decoys already give a decisive zero at roughly a third of the added cost.

### Options

1. **Keep 1,000 in the main run.** Simplest; costs 12 extra minutes per run forever.
2. **Drop to 300.** *p*(0) = 0.008 is ample for this purpose; cost maybe 4–5 min rather than 12.
3. **Move the decoy curve to an opt-in stage**, as was already done for the genome-mapped channel's
   decoy (`genomic-decoy`). The main run keeps a small decoy set as a cheap sanity check; the powered
   curve is run when the floor is being set or questioned.

### What I would do, and the argument against

I lean to **option 3**, because the project already set this precedent one round earlier and wrote
the reason into the job script: *"a control must not sit on the critical path of a run someone is
waiting for."* The floor is a settled parameter; re-measuring it at full power on every run is
re-deriving a constant. Option 2 is the pragmatic middle.

The argument against option 3 is real, though: a control that is opt-in is a control that stops being
run, and the floor has already been wrong twice (8 bp from the wrong instrument, 10 bp from the right
instrument at the wrong resolution). If the decoy curve only runs when someone remembers, the next
threshold change will be argued rather than measured. **That is the trade-off I cannot settle:
whether this control is a one-off calibration or a standing regression test.**

---

## Decision 8 — `min_assign = 6` sits exactly at the 1% tolerance

### The numbers

Cost per clip length, exact matching, chance computed against the **15 real** flanks (what a spurious
assignment must match) and validated against the decoys:

| clip length | molecules | chance vs real flanks | expected spurious | as a share of that band |
|---|---|---|---|---|
| 4 bp | 17 | 11.10% | 1.887 | 11.1% |
| 5 bp | 9 | 2.89% | 0.260 | 2.9% |
| **6 bp** | **35** | **0.730%** | **0.255** | **0.73%** |
| 7 bp | 29 | 0.183% | 0.053 | 0.18% |
| 8 bp | 8 | 0.046% | 0.004 | 0.05% |

Cumulative, by floor:

| floor | molecules admitted | expected spurious | share |
|---|---|---|---|
| 5 bp | 940 | 0.574 | 0.06% |
| **6 bp (current)** | **931** | **0.314** | **0.03%** |
| 7 bp | 896 | 0.059 | 0.01% |
| 8 bp | 867 | 0.006 | 0.00% |

### The honest reading

Two framings give different answers, and this is the crux.

**Across the whole accepted set**, 6 bp is comfortably safe: 0.31 expected spurious among 931
molecules is 0.03%, far inside a 1% tolerance.

**Within the band it admits**, 6 bp sits *at* the tolerance: the 6-bp molecules carry a 0.73% chance
each of matching a wrong real flank, and the 25 assignments that band contributes come with 0.26
expected spurious — about **1 spurious per 96 real**, i.e. ~1.0% of the assignments it adds. 7 bp is
0.18%, five-fold inside.

So the review's T1 observation reappears in a different form: the 6-bp band is the only one operating
at the limit rather than beneath it. Whether that matters depends on whether the tolerance applies to
the assignment set as a whole or to each increment of it — and that is a choice about what the
threshold means, not a measurement.

### Options

1. **Keep 6 bp.** +25 assignments over ~20 previously unlabelled cells, ~2% of the cell table, at
   0.26 expected spurious.
2. **Take 7 bp.** Costs those 25 assignments; expected spurious across the accepted set drops
   0.31 → 0.059, a five-fold margin improvement.

### Additional facts that bear on it

- The **mono-replicate control agrees with 6 bp**: the 96 molecules at 6–9 bp contribute **zero**
  cross-mouse errors, and the pooled cross-mouse rate *improved* from 0.62% to 0.50% when the floor
  moved to 6 with exact matching, because two of the five rejected 10–12 bp assignments were errors.
  That is an independent instrument and it does not object.
- The **wild-type tubes add nothing** at 6 bp: still one called cell carrying a junction, the Rx_1
  event, unchanged by the floor.
- `min_assign = 6` is **only defensible with `exact_below = 15`**. Under the `L %/% 10` budget a 6-bp
  clip is allowed one mismatch and 21 of those 35 molecules hit a decoy. The two constants are one
  decision; please do not treat the floor as adjustable on its own.

### What I would do

I lean to **keeping 6** on the strength of the mono-replicate control agreeing and the whole-set
figure being 0.03%. But I would not argue against 7: the 25 assignments are ~2% of the cell table and
buy nothing that any conclusion in the report rests on, and 7 bp is the only setting where every band
is comfortably inside the tolerance rather than one of them sitting on it.

**The question for you:** does the 1% tolerance apply to the assignment set as a whole, or to each
band admitted? On the first reading 6 is fine; on the second, 7 is the defensible floor.
