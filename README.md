# mega27-23b-dna-encoder-cyborg-cell
- **DNA data encoder**: bytes -> base-3 trits -> SHA-256 keystream scrambler -> never-same rotating base code. Homopolymer runs > 1 impossible by construction; GC ~0.5; a Fletcher-based integrity suffix (the block variant keeps only the low Fletcher byte s1 mod 255, so an explicit compensating two-trit error passes it while changing the full Fletcher checksum; see results/BENCHMARK_STATUS.md) + 3-fold majority-vote error correction; failure detection is incomplete: strict-noise whole-file audits retain silent wrong recoveries. No never-silent-corruption guarantee.
- **Cyborg-cell metabolic logic loop**: chemostat ODE with replication-investment gate; analytic steady state verified against simulation; sensing-efficiency reallocation bound derived and verified monotone in sensor efficiency alpha.
Run: `pip install -e . && pytest`

Current manuscript and superseded draft status: see `paper/DRAFT_STATUS.md`. The historical `50p` filename is not an audited body-page completion claim.
