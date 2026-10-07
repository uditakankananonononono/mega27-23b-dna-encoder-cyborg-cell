# Audit-local decoder input boundaries

The protected-nonce file decoder previously raised IndexError for empty or
undersized geometry, TypeError for nonpositive stripe counts and KeyError for
non-ACGT input. These are interface crashes, not accepted wrong scientific
outcomes. Checks now reject unknown detectors, invalid stripe types/counts,
empty or unequal strands, non-string strands, unaligned lengths, insufficient
frames for data plus parity, and non-ACGT alphabet with ValueError.

Twelve malformed-input regressions retain those boundaries. A clean encoded
empty byte payload is distinct from an empty strand and still decodes. The
production codec is untouched. All saved geometry/recovery negatives remain;
these checks do not add error correction, physical validation or guarantees
for arbitrary consensus/outer-parity corruption. Full suite: 68 passed.
