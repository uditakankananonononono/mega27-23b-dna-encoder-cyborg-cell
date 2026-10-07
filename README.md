# mega27-23b-dna-encoder-cyborg-cell
- **DNA data encoder**: bytes -> base-3 trits -> SHA-256 keystream scrambler -> never-same rotating base code. Homopolymer runs > 1 impossible by construction; GC ~0.5; a Fletcher-based integrity suffix (the block variant keeps only the low Fletcher byte s1 mod 255, so an explicit compensating two-trit error passes it while changing the full Fletcher checksum; see results/BENCHMARK_STATUS.md) + 3-fold majority-vote error correction; failure detection is incomplete: strict-noise whole-file audits retain silent wrong recoveries. No never-silent-corruption guarantee.
- **Cyborg-cell metabolic logic loop**: chemostat ODE with replication-investment gate; analytic steady state verified against simulation; sensing-efficiency reallocation bound derived and verified monotone in sensor efficiency alpha.
Run: `pip install -e . && pytest`

Current manuscript and superseded draft status: see `paper/DRAFT_STATUS.md`. The historical `50p` filename is not an audited body-page completion claim.

### Mapping-free bounded-prefix rank ceiling
ConditionalK=67088chunks:4595receiver-predicatepassinguniqueseeds give coefficientrank<=4595,nullity>=62493;all4677RSexactuniqueseeds give rank<=4677,nullity>=62411. No seedmapping,actualrank,coveredcoordinatecount,peeling or payloaddecode is claimed. This countbound needsno Python3-for-Python2 PRNG substitution. See `results/dna_fountain_rank_bound_audit.json`.

### Pinned CPython2.7 sample compatibility
An actualcompiled2.7.18 oracle and explicitfloat-sample mirror agree on all4677unique-prefixseed vectors. ModernPython3sample differs on all4677. This is pinnednativebranch/sharedCDF compatibility,not originalunpinnedruntime,NumPybranch,originalauthororacle orpayloadidentity. Nopeeling/decodingperformed. See `results/python27_sample_parity_audit.json` and archivedvectors.

### Conditional graph coverage, not decoding
Pinned2.7.18/sharedCDF vectors give51819/52323 coveredcoordinates for4595screened/4677unscreeneduniqueseeds;bothhave9initialsingletons and10structurallypeeledcoordinates. Coverage/connectivitydoesnot implyrank or authenticatedchunkrecovery. Payloadbytesnotread. See `results/python27_graph_coverage_audit.json`.
