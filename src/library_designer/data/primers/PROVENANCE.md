# Bundled orthogonal primer sets

Primer sequences are written 5' to 3' as ordered. On load (`primers.load_primer_set`) any primer carrying the Golden Gate enzyme's recognition site is dropped, because a set validated for orthogonality is not necessarily free of your enzyme's site.

## `subramanian2018.csv` (default)

165 experimentally validated, mutually orthogonal 20-mer amplification primers. It is a flat pool: any two are orthogonal, so the tiler draws two per tile. Used combinatorially they specify 13,695 gene pairs.

- Source: Subramanian SK, Russ WP, Ranganathan R. "A set of experimentally validated, mutually orthogonal primers for combinatorially specifying genetic components." *Synthetic Biology* 3(1):ysx008 (2018). PMC7445780. doi:10.1093/synbio/ysx008
- Retrieved from the article's Supplementary Table S1 (`ysx008_supp_st_1.xlsx`), keeping the 165 rows flagged "Keep primer in orthogonal set? = Yes".
- License: CC BY-NC (facts and sequences; cite the paper).
- One of the 165 (`SUB151`, `TATAACAGGCTGCTGAGACC`) contains a BsaI site and is dropped automatically when `enzyme="BsaI"`, leaving 164 usable, or about 82 tiles.

Selectively amplifying sub-pools from a shared oligo pool with orthogonal flanking primers goes back to Kosuri et al., *Nat. Biotechnol.* 28:1295 (2010) (PMC3139991), whose primers were themselves seeded from Xu et al.'s 240,000 orthogonal 25-mers (*PNAS* 2009). The Subramanian 2018 set is the reusable, validated resource.

## `gck800.csv` (worked example)

7 synthetic forward/reverse primer pairs, included to show the paired-set format and to give the glucokinase tiled example a concrete primer source. They are illustrative sequences, not an experimentally validated orthogonal set, so screen your own primers before ordering. Capacity is 7 tiles.
