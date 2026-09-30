# Benchmark status: NOT validated against published codecs
September 30 audit withdraws catastrophic downstream SUBSTITUTION propagation for the implemented differential mapping; only adjacent trits change. Full-strand rejection is a detection effect.
Historical nominal substitution rates are attempts with expected realized rate 0.75*p. Strict-channel rerun is separate.
Fountain CRC8/uniform-low-degree reimplementation fails 8/32 zero-noise 128-byte messages. No published-codec superiority is established. Retain old data as local-harness history only.
See channel_semantics_audit.json and strict_substitution_benchmark.json, their executable scripts, and prominent paper correction.

External port audit (f97c1b8a81f5c5b819209d5b5e26c9c8f4439495): clean hashes pass for four 2KB payloads. Fixed-budget strict substitution tests are retained in fountain_port_audit.json. This does not restore superiority: RSNS uses a single 18,495-nt strand versus 180 Fountain oligos of 152 nt. A segmented RSNS control uses 64 oligos of 495 nt with ideal known order, no indexing charged, and no outer redundancy. At 1% it recovers 7/12 files versus long RSNS 9/12; at 2% 0/12 versus 1/12. Segmented density .517 versus external port .599 bits/encoded nt. Both comparisons exclude primers, read consensus and molecular-copy budgets and use only four independent payloads. Repeats are noise replicates. No matched-geometry, optimized-resource or biological verdict exists.

Near-matched clean geometry gate: 63x496-nt external oligos, 31,248 encoded nt (1.36% below local segmented budget). Four payloads x four parity settings x three homopolymer screens = 48 clean configurations. At max homopolymer 3, port peeling is 0/4 at every parity setting; exact GF(2) clean rank/solve recovers 3/4,3/4,2/4,4/4 for RS 2/8/16/24. At max homopolymer 4 peeling is 4/4,3/4,4/4,4/4; at max 5 it is 4/4 everywhere. Screen relaxation changes constraints. Exact Gaussian solve is only a mathematical clean diagnostic, not published decoder or noise performance. See fountain_clean_graph_audit.json and per-message files. Failed message-zero strict sweep retained but excluded from superiority reasoning. Frontier remains unresolved.

Clean-valid HP5 parity/noise sweep: near-matched 63x496-nt budget and shared four2KB payloads, three noise repeats. Whole-file recovery at0/.5/1/2/3% strict changed-base probability: RS2 12/0/0/0/0, RS8 12/12/2/0/0, RS16 12/12/12/1/0, RS24 12/12/12/12/0 of12. Local segmented HP1 control12/12/7/0/0. This refutes an unqualified superiority reading of the first external RS2 pilot. HP5 versus HP1, ideal local order, unoptimized budgets and excluded primers/copies prevent universal frontier or biological claims. Clean gates and hashes verified; see fountain_valid_noise_audit.json and per-message trial records/FASTA. Preserve all previous negatives.

## Source-level rate/checksum/claim correction (2026-09-30)
Exact v1 total length c*(6L+24); v2 c*38*(ceil((6L+24)/32)+1).
Three-copy limiting payload rates 0.444444 and 0.374269 bits/encoded nt,
excluding primers/index metadata. Byte packing is 6 trits per byte, not capacity.
Block suffix checksum keeps only low Fletcher byte s1 modulo 255. An explicit
compensating two-trit error passes it while changing the full Fletcher checksum.
No uniform 1/729 collision or universal exponential worst-case bound follows.
The old union-bound expression evaluates about 0.2036, not 1e-4, at its stated
parameters; it is not a justified decoder failure model. Independent-copy
strict-majority events do not equal four-letter plurality failure events.
43 reviewed paragraph corrections are logged verbatim, retaining numeric trials.
This is a targeted consistency pass, not proof that all manuscript claims have
been audited. Biological calibration and a fully matched codec frontier remain
open. Test suite: 33 passed; rendered changed pages inspected.
