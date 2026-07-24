# Test fixtures

## `synthetic_hacyp1_destination.gb`

A synthetic circular destination plasmid used by `tests/test_vectors.py` to run the
`.gb` read/locate/splice/emit path on a real file rather than an in-memory string. It
is not a lab construct.

The backbone is random BsaI-clean filler. The public hAcyP1 CDS (`ACYP1_HUMAN`) sits on
the minus strand as a `CDS` feature, and a stretch of backbone is annotated `msGFP2` so
the tests can confirm backbone features carry over onto the emitted per-tile maps.

BioPython reads SnapGene `.dna` files but cannot write them, so this fixture is GenBank
only. To regenerate it, rebuild a circular record with the CDS on the reverse strand and
a BsaI-clean backbone, then write it with `Bio.SeqIO.write(..., "genbank")`.
