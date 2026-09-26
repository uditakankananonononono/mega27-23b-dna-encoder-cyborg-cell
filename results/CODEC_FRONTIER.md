# Codec redundancy-density frontier (2026-09-26)

Benchmark: experiments/benchmark_codecs.py (24 random messages x 64/256 B,
substitution 0/1/2/3%, deterministic seeds; fountain baseline carries the CRC8
fairness fix).

## New variants tested (density-improvement arm)

| codec | density (bits/base) | recovery @1% (64B/256B) | homopolymer max |
|---|---|---|---|
| ours (3x majority) | 0.42-0.44 | 1.0 / 0.875 | 1 |
| ours_v2 (3x, blocked+parity) | 0.32-0.36 | 1.0 / 0.917 | 2 |
| ours_2x (2 copies) | 0.63-0.66 | 0.04 / 0.0 | 1 |
| ours_v2_2x (2 copies) | 0.48-0.54 | 0.25 / 0.0 | 2 |
| goldman | 0.42-0.44 | 0.96 / 0.875 | 1 |
| fountain (CRC8 fix) | 0.63-0.72 | 0.21 / 0.5 | 3 |

## Finding (honest negative + quantified frontier)

Dropping from 3 to 2 copies buys 1.5x density but COLLAPSES recovery at 1%
substitution (0.04-0.25 at 64B, 0.0 at 256B). Mechanism: 2-copy vote cannot
disambiguate disagreements; the v2 per-block checksum needs a fully clean single
copy (p(clean 38-trit block at 1%) = 0.68), and parity repairs only one bad block
per message, while a 256B message has ~47 blocks. Triple redundancy is not slack -
it IS the error-correction mechanism of the no-homopolymer design.

Frontier, quantified: at homopolymer max <=2, ~0.42 bits/base is the price of
robust recovery at 1-2% substitution; fountain reaches 0.63-0.72 bits/base but
pays homopolymer 3 and still loses recovery at 1-2% under the fairness-fixed
harness. Next candidates (queued for ChatGPT redirection consult): fountain-style
droplet structure ON the never-same alphabet (1 trit/base droplets + RS-style
outer code), or inner LDPC over trits.

## v3 candidate measured and rejected: never-same fountain droplets (2026-09-26)

nsfountain_encode/decode in experiments/benchmark_codecs.py: LT droplets whose DNA
is the never-same trit code (homopolymer 1 by construction, no screening) + CRC8
per droplet. Measured on the same harness:

| codec | density (256B) | recovery @0/1% (256B) | homopolymer |
|---|---|---|---|
| fountain (CRC8) | 0.717 | 1.0 / 0.5 | 3 |
| nsfountain | 0.478 | 1.0 / 0.125 | 1 |

NEGATIVE, with mechanism quantified: the never-same droplet is 3x longer for the
same payload (102 vs 34 bases for 16 B), and CRC-detection turns any substitution
into a whole-droplet erasure, so droplet survival at 1% is 0.99^102 = 0.36 vs
0.99^34 = 0.71 for the direct-map droplet. The density gain over 3x replication
(0.478 vs 0.437) does not compensate the erasure rate; at 64B both LT variants
are peeling-fragile at tiny k (fountain 0.42, nsfountain 0.0 at 0% noise - small-k
LT rank failures, not a code defect). Conclusion: closing the density gap needs an
inner code that CORRECTS (RS-style over trits), not detects; detection-only
droplets lose to triple majority at every error rate tested. Candidate retired.

## v4 LANDED: RS-over-GF(3^5) correcting inner code ("rsns") - frontier leader (2026-09-26)

The frontier conclusion above said density closure needs a CORRECTING inner code
over trits. Implemented in experiments/rs_inner.py: Reed-Solomon over GF(3^5) =
GF(243), 5 trits per symbol, primitive polynomial found and verified
programmatically (x has order 242), systematic encoder, Berlekamp-Massey +
Chien + Gaussian-elimination magnitudes with errors-and-erasures support.
Blocks RS(45,30) (rate 2/3) + RS(9,3) header. Strand = RS symbols -> scramble
-> never-same DNA. Homopolymer 1 by construction.

Discovery (the free-erasure mechanism): a substitution at base i corrupts the
two adjacent trits, and a spacing violation can decode to trit value 3 - a
value the 3-symbol alphabet can never legitimately carry. The 4-letter DNA
alphabet encoding a 3-symbol stream therefore turns a large share of
substitutions into flagged ERASURES at zero redundancy cost, and RS decodes an
erasure at half the budget of an error (2e + s <= n - k). This is why the
correcting inner code wins where detection-only (CRC8) designs collapsed.

Unit verification: 3421/3421 random errors+erasures at capacity (2e+s <= 15)
decode exactly; encoder syndromes zero on 300 random codewords.

Official harness (experiments/benchmark_codecs.py, same 24 msgs x 64/256B x
0/1/2/3% substitution):

| codec | density 64/256B (bits/base) | recovery 64B @0/1/2/3% | recovery 256B @0/1/2/3% | homopolymer |
|---|---|---|---|---|
| rsns | 0.711 / 0.813 | 1.0 / 0.958 / 0.875 / 0.417 | 1.0 / 1.0 / 0.583 / 0.083 | 1 |
| fountain (CRC8) | 0.627 / 0.717 | 0.417 / 0.208 / 0.042 / 0.0 | 1.0 / 0.5 / 0.0 / 0.0 | 3 |
| ours (3x) | 0.42-0.44 | 1.0 / 0.833 / 0.292 / 0.125 | 1.0 / 0.875 / - / - | 1 |

rsns STRICTLY DOMINATES the previous density leader: higher density, lower
homopolymer, better recovery at every error rate and both message sizes. At the
benchmark's 1% substitution tier it is perfect at 256B (fountain: 0.5).

Open gaps (honest): recovery degrades at 2-3% substitution (0.583/0.083 at
256B) - the rate-2/3 budget is sized for the 1% tier; a rate sweep
(BLK_NSYM 18/21) would trade density for high-rate robustness and is queued.
Indels still desynchronize the stream (channel is substitution-only by
benchmark definition). RS decode cost is Python-speed, fine at these sizes.

## Rate sweep (256B, same harness seeds): the density-robustness Pareto frontier

| body block | rate | density | recovery @0/1/2/3% |
|---|---|---|---|
| RS(38,30) | 0.79 | 0.959 | 1.0 / 0.542 / 0.042 / 0.0 |
| RS(45,30) | 0.67 | 0.813 | 1.0 / 0.917 / 0.625 / 0.125 |
| RS(51,30) | 0.59 | 0.719 | 1.0 / 1.0 / 0.917 / 0.5 |
| RS(57,30) | 0.53 | 0.644 | 1.0 / 1.0 / 0.917 / 0.708 |

At fountain's own density (0.717), rsns rate-0.59 delivers 1.0 @1%, 0.917 @2%,
homopolymer 1 - fountain gives 0.5 @1%, 0.0 @2%, homopolymer 3. The Pareto
curve itself is a deliverable: one parameter (block redundancy) slides the
codec along density-robustness, all points beating both baselines.
Results: results/rs_rate_sweep.json.

Capacity scaling (2026-09-27, experiments/rsns_capacity_sweep.py, results/rsns_capacity_sweep.json):
message size 256->2048 B raises density 0.813->0.886 bits/base (header amortizes
toward the asymptotic block rate) while per-message recovery at fixed substitution
rate degrades with block count as expected for fixed-length blocks: 1% substitution
gives 0.81/0.88/0.75/0.50 recovery at 256/512/1024/2048 B; 2% gives 0.25/0.125/0/0.
Message failure tracks 1-(1-p_block)^nblocks; the next frontier lever is longer RS
blocks or interleaving, not a different inner code family.
