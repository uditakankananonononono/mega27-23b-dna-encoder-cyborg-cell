# Optional ML audit, September 30
The 130 natural inputs are accession windows from E. coli K-12 strain C3 WGS records, not MG1655 genes. FASTA headers and the fresh manifest preserve provenance.

The model discriminates natural accession windows from synthetic codec windows. It does not predict synthesis difficulty. All natural accessions share one strain/project; there is no independent-genome validation. Excluding the homopolymer feature does not exclude its signature contained in k-mers. No experimental synthesis labels are used. Biological-validation extension remains open.

Reproduce with `OPENBLAS_NUM_THREADS=1 python3 experiments/constraint_ml.py`. JSON records checksums, window starts, seeds, fold scores, ablations and conditional bootstrap intervals.
