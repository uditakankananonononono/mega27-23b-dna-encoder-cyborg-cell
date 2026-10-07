# DNA-Fountain original-read source qualification, 7 October 2026

Ingest-only. Plan SHA256 8cebcea88dfdc4c23fdf267cb9d0a2aa2bcf7c1aa76ca77d2a2bfec2f72d7e3a.

16 paired runs, 95,237,393,879 compressed FASTQ bytes. ENA read_count sum 479,469,649, not unique molecules or receiver successes. No reads downloaded or decoded. Per-run sizes, identifiers, checksums and exact observed ENA portal report URLs are in results/dna_fountain_source_audit.json; source TSV snapshots replay the tally.

Original source https://github.com/TeamErlich/dna-fountain pinned8ee2777aa5e9e101e5d756f7e2449f6672b66f1f. INSTALL https://github.com/TeamErlich/dna-fountain/blob/master/INSTALL specifies NumPy,Cython,reedsolo,tqdm,logging,SciPy without versions. README specifies PEAR stitched152nt filters and Python2 receiver. glass.pyx compares byte differences after RS re-encoding, not nt Hamming distance.

Original archive URL http://files.teamerlich.org/dna_fountain/dna-fountain-input-files.tar.gz fails DNS for HTTP and HTTPS here. Fork https://github.com/jdbrody/dna-fountain pinnedf97c1b8a81f5c5b819209d5b5e26c9c8f4439495 holds an input-file container, not verified encoded payload. Size2145290, MD54d21a7b72451302c0afc8c8f0c1eb1d6, SHA25614f42ef9a2c1d70d3a14199b18a97c9ded93866cc2912c5be35fb8ff592d7be5. Member list only inspected, no extraction or nested archive execution. Contains zipbomb-labelled test input. Zero-padding this container to2146816 bytes does not yield expectedMD58651e90d3a013178b816b63fdbb94b9b. Do not confuse the container with the encoded archive.

Recipe discrepancy: declared padding2116608 vs decoder67088*32=2146816, difference30208. Retain, do not silently repair. Conditional for bounded physical preprocessing, blocked for checksum-certified source reconstruction until exact payload/convention recovered. Neither outcome validates local bound-nonce codec.

Supporting context sources inspected: publication https://www.science.org/doi/10.1126/science.aaj2038 and index https://www.omicsdi.org/dataset/project/PRJEB19307. These are context; executable mechanics and tallies come from pinned code and primary ENA snapshots. No physical success-rate claims copied from publication.
