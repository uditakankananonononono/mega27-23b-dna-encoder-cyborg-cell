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
