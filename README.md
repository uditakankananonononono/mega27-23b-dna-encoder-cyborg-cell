# mega27-23b-dna-encoder-cyborg-cell
- **DNA data encoder**: bytes -> base-3 trits -> SHA-256 keystream scrambler -> never-same rotating base code. Homopolymer runs > 1 impossible by construction; GC ~0.5; Fletcher-16 integrity + 3-fold majority-vote error correction; loud failure (never silent corruption).
- **Cyborg-cell metabolic logic loop**: chemostat ODE with replication-investment gate; analytic steady state verified against simulation; sensing-efficiency reallocation bound derived and verified monotone in sensor efficiency alpha.
Run: `pip install -e . && pytest`
