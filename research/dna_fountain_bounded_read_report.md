# Original archive recovery and bounded physical-read preprocessing

Original encoded archive was not recovered. The mirrored container's README explicitly says the gift card was removed. Original git history/release page and Zenodo article deposit supply no encoded archive. Primary bioRxiv body specifies2,146,816bytes and67,088chunks, but README truncate value remains discrepant. Container truncation/padding candidates all fail the expected digest. No nested archive extracted or executed. Recovery is not certified.

A separately frozen fallback downloaded exactly8MiB from each ERR1816980 mate through verified HTTP206 byte ranges. Prefix SHA256s and full-file sizes are in the result ledger; these hashes are not ENA fullfile MD5s. Incomplete gzip streams yielded92,414 and86,355 complete FASTQ records, with terminal partial records excluded. First10,000 synchronized151nt read pairs were checked. No randomization or population-rate inference.

Strict exact-overlap stitching, minimum20nt, yielded6368pairs and rejected3632. This is not PEAR reproduction.5278stitched sequences were152nt;10 contained non-ACGT symbols,395 failed exactRS2 codeword check,4873 passed. Those4873 records contain4677unique seed IDs and4677unique sequence hashes. RS2 decode/re-encode equality uses reedsolo1.7.0, not original unpinned environment. No GC/homopolymer screening, soliton mapping, peeling or wholepayload decoding performed. RS-valid does not prove original source identity or full receiver acceptance.

Sources: https://www.biorxiv.org/content/10.1101/074237v4.full-text ; https://github.com/TeamErlich/dna-fountain/releases ; https://zenodo.org/records/889697 ; https://github.com/jdbrody/dna-fountain . Exact physical file URLs in result prefix ledger:
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_2.fastq.gz

No reconstruction or local bound-nonce physical validation claim. Page gate remains uncertified.
