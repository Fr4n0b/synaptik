# Interpreting the output

A short guide to reading results without over-reading them. Aimed at someone
using this to decide what experiment to run next.

## The one-sentence version

The tool tells you **how much the field has looked, and how consistently it has
found the same thing** — not whether the thing is true.

## Start with the tier, not the ranking

The ranked list is a prioritisation aid. The tier is the actual finding.

| Tier | What to do with it |
|---|---|
| **strong** | Multiple groups, multiple methods, consistent direction. Reasonable to build on. Still worth reading two or three primary sources. |
| **moderate** | Real but thinner. Read the supporting papers before designing around it. |
| **weak** | Treat as a hypothesis. Often means directions conflict or all evidence uses one method. |
| **single-lab** | One group's finding, however well published. Not independently replicated. High-value to test yourself precisely *because* nobody else has. |
| **insufficient** | Evidence exists but weighs almost nothing — usually retraction or review-only assertions. |

## Numbers that should change your behaviour

**Independent labs = 1.** Everything else is noise until this is ≥2. Nine
papers from one group is one finding. This is the most common way to be misled
by a literature search, and the number most worth verifying by hand.

**Assay families = 1.** The result may be a property of the measurement rather
than the biology. If everything is imaging, an antibody or segmentation
artifact reproduces perfectly across labs. Look for an orthogonal readout
before committing.

**Leave-one-lab-out FAILS.** Delete the most prolific group and the association
evaporates. That is one group's narrative, not a field consensus.

**Direction consistency < 70%.** The literature disagrees about the *sign*.
This is genuinely useful information — it usually means the effect depends on
something nobody has pinned down (DIV, line background, culture conditions).
Worth reading the conflicting reports directly; the discrepancy often points at
the real variable.

**Null results = 0 across many reports.** Either an unusually robust effect or
a filtered literature. Assume the latter until shown otherwise.

## A worked example

```
pathway                        tier          enrich       q  labs  assays
neurotransmitter release       strong          16.9  1.6e-06     9       2
exotic kinase cascade          single-lab     225.2  3.5e-09     1       1
```

Row two has 13× the enrichment and a smaller q-value. It is much less
trustworthy: one group, one assay, no nulls across nine reports.

If you sorted by enrichment or by p-value — which is what most literature-mining
tools do — you would chase row two.

The correct read: row one is background you can build on; row two is an
interesting, untested claim. Both are useful, for opposite reasons.

## Reading experiment lookup

**EXACT** means same system, same maturity, same perturbation class, same
readout. If exact matches exist, read them before running anything.

**CLOSE** differs in one dimension — check the `differs:` line, since which
dimension matters enormously. A different DIV bin is usually disqualifying for
synaptic phenotypes; a different marker antibody usually is not.

**Consensus is disputed.** Do not average. Read the conflicting reports and
work out what differs between them — that difference is often the real finding.

**"All matches come from a single group."** The experiment has been done once,
by one lab. Your replication has genuine value.

## Reading splicing output

**Regulator inference** is a hypothesis about *why* the events occur. Cryptic
exon inclusion across multiple genes is the canonical TDP-43 nuclear-loss
signature and is inferred even when no paper names it — treat it as a strong
prior, not a measurement.

**`sign_confidence: uncertain`** means the reported ΔPSI conflicted with the
prose, or the detection tool has no standard sign convention. Verify against
the source figure before believing the direction.

**NMD risk: high** changes the therapeutic question. If the aberrant transcript
is degraded, your defect is *loss of protein*, and measuring transcript alone
will mislead you. Measure protein.

**ASO feasibility** is about the event class, not your specific event. It tells
you which strategy applies; it does not mean an ASO will work.

## Reading drug candidates

Only direction-matched compounds are shown by default. A compound is listed
because its mechanism **opposes** the observed change — an activator for a
down-regulated pathway.

**`bbb_plausible: false` does not disqualify a compound for dish work.** It
matters for in vivo follow-up only. This is why the two flags are separate.

**Efflux-substrate warnings override the CNS score.** Physicochemistry says the
compound should reach the brain; a transporter says otherwise. Loperamide is
the textbook case.

**"already tested" annotations are the highest-value line in the output.**
Somebody has put this compound on neurons and reported what happened. Read that
paper first.

## Questions this tool cannot answer

- Does pathway W *cause* phenotype P? (Co-occurrence cannot establish this.)
- What is the true effect size? (Most papers do not report variance.)
- Has anyone tried this and failed? (Unpublished nulls are unrecoverable.)
- Is this pathway relevant in *my* specific line and background? (No.)

If a result matters to your next six months of work, read the primary sources.
This tool tells you which ones, and how sceptical to be.
