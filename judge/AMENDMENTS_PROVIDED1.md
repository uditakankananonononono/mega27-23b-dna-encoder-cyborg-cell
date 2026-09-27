# Lane-23b AMENDMENT QUEUE - LOCKED BEFORE EXECUTION
Verdict: PROVIDED round 1 (WhatsApp 10:58:02 IST, wamid...MjhCRAA=, verbatim in judge/round_provided1_verdict_whatsapp.txt).
Locked: 2026-09-27 11:00 IST, against 62pp docx build 430e78b. No execution began before this lock.

## Standing directive (verbatim prefix): "MAKE IT MORE ADVANCED AND TEST IT AGAINST REAL TOOLS"
Headline principle to land (#18): "Optimal information density is limited by cellular resource constraints."

## Weakness items (#1-#18) mapped vs current 62pp docx
- #1 real biological validation (published sequencing-error datasets: Illumina/Nanopore/synthesis error profiles): NEW - cheap: public error models from literature, then #11 injection test
- #2 rename away from "cyborg cell" -> "Biologically constrained DNA information storage system": FRAMING - execute in title/headers
- #3 metabolic-coupling justification (storage capacity vs GC burden/replication cost/transcription burden tradeoff): PARTIAL (metabolic model exists; add formal tradeoff definition)
- #4 modern baselines (constrained coding + neural DNA storage, beyond Goldman/DNA Fountain): NEW - literature comparison table
- #5 unique contribution statement ("optimized code design under simultaneous biological and information constraints"): EDITORIAL
- #6 multi-objective optimization (density vs GC extremes/synthesis difficulty/metabolic cost): PARTIAL (capacity sweep exists) -> formal Pareto analysis
- #7 evolutionary stability (mutation accumulation; how long info survives): NEW - cheap simulation on our sequences
- #8 one host organism (E. coli): FRAMING - pin all constraints to E. coli
- #9 natural-genome realism comparison (bacterial/plasmid/viral genomes): NEW - cheap, download reference stats
- #10 genome-scale constraint chain (storage sequence -> expression burden -> metabolic cost): PARTIAL - formalize the coupling
- #11 experimental proxy (inject encoded sequences into public sequencing datasets, measure recovery): NEW - the flagship validation
- #12 realistic error distributions (indel/substitution/homopolymer): NEW - error model upgrade feeding #1/#11
- #13 scaling demonstration 1KB -> 1MB -> 1GB theoretical: PARTIAL (capacity law exists) -> explicit scaling section
- #14 robustness under varied constraints (GC limits/homopolymer limits/parity size sensitivity): PARTIAL (capacity sweep) -> full sensitivity grid
- #15 optional ML (learn biological constraints from genomes; transformer predicts synthesis difficulty): NEW - optional per verdict
- #16 frame "cell" as information metabolism: EDITORIAL
- #17 evolutionary comparison (why biology avoids certain patterns; natural vs synthetic storage DNA): NEW - pairs with #9
- #18 land the biological principle as positive headline: EDITORIAL spine


## LANDED 2026-09-27 (foldback 1)
- #2/#8/#16: rename LANDED - title now "A Biologically Constrained DNA Information Storage System..." (no "cyborg cell" anywhere in the paper; host pinned to E. coli). #18: headline principle landed in abstract + section 9i ("optimal information density is limited by cellular resource constraints"). Docx rebuilt 8,890 -> 9,227 words; new sections 9i (mutation-accumulation stability, generations-to-loss table) and 9j (Illumina/synthesis/Nanopore error-profile injection-recovery).
- #7: event-driven mutation-accumulation stability. Coded 256B payload survives median 5.5-12.5 mutation events -> expected loss at ~1,984 generations (mu=1e-6 neutral drift), ~451 (1e-5), ~27 (1e-4, e.g. mutator stress). Uncoded 2048-base payload intact-probability at mu=1e-5: 0.98@10gen, 0.82@100, 0.36@500, 0.13@1000.
- #1/#11/#12: literature-scale error-profile injection-recovery proxy (256B, 16 trials, density 0.813 bits/base): Illumina-like 1.0 recovery; synthesis-like (5e-3 subst + 1e-3 indel, homopolymer x2.5) 0.0625; Nanopore-like (2e-2 subst + 3e-2 indel, homo x3) 0.0. Positive framing: the RS(45,30)+RS(9,3) inner code exactly covers Illumina-scale channels; the synthesis/Nanopore gap sizes the redundancy still needed (a design spec, not a hidden failure).
- Files: experiments/evolution_error_models.py, results/evolution_error_models.json.

## First deliverables (cheap, existing data + public references)
4. Natural-genome realism comparison (#9/#17) from public genome stats
