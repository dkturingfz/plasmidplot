# plasmidplot 0.1.1

## Changed behavior

* `pp_find_sites()` now numbers a restriction site by the last base before
  the cut on the top strand, the convention 'SnapGene' and NEB maps use.
  Positions are one lower than in 0.1.0, which reported the first base after
  the cut: `EcoRI` on pUC19 is now 396, not 397.

## Bug fixes

* `pp_find_sites()` finds overlapping sites. Two `BssHII` sites in
  `GCGCGCGC` were counted as one, so a double cutter could be reported as a
  unique one.
* `pp_find_sites()` searches both strands for a custom recognition sequence
  that is not a palindrome, such as `c(BsaI = "GGTCTC")`.
* `read_genbank()` and `read_embl()` read
  `join(complement(300..400),complement(100..200))` as the 100..400 feature
  it is, instead of one wrapping the origin, and no longer take the digits of
  an accession in a remote segment (`J00194.1:100..202`) for a position.
* The scale of a molecule a few dozen bases long stays on whole base pairs.

## Layout

* Crowded labels on a circular map are spread about the middle of their
  group and kept on the page. They used to be pushed one way only, off the
  bottom of the page or across the features.
* Crowded labels on a linear map are slanted and fanned out sideways. They
  used to stack into a tower with every leader running through the labels
  beneath it. Maps whose labels fit level, in up to three rows with no leader
  crossing a label, are drawn as before.
* A feature too short to see, such as one of a single base, is drawn at a
  small minimum width instead of vanishing.

## Checks

* `pp_marker()`, `pp_site()`, `pp_style()` and `plot()` reject a bad label,
  color, offset, width, background or shape number where it is given, rather
  than failing inside 'grid' at draw time.
* `pp_style()` reports an unnamed override instead of ignoring it.
