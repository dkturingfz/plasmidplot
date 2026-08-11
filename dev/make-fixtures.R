# Generates the test fixtures. Run from the package root; the sequence is
# built here rather than hand-typed so its length always matches the header
# and the enzyme sites are exactly the ones the tests expect.
setwd("i:/plasmidplot")
set.seed(42)

LEN <- 1200
implants <- list(
  EcoRI   = list(site = "GAATTC",   at = 100),
  BamHI   = list(site = "GGATCC",   at = 400),
  HindIII = list(site = "AAGCTT",   at = 760),
  NotI    = list(site = "GCGGCCGC", at = 1010)
)

repeat {
  seq <- paste(sample(c("A", "C", "G", "T"), LEN, replace = TRUE), collapse = "")
  for (nm in names(implants)) {
    s <- implants[[nm]]
    substr(seq, s$at, s$at + nchar(s$site) - 1L) <- s$site
  }
  # Accept only a draft where each implanted enzyme cuts exactly once, so the
  # fixture cannot drift into a second accidental site.
  counts <- vapply(implants, function(s) {
    hay <- paste0(seq, substr(seq, 1L, nchar(s$site) - 1L))
    length(gregexpr(s$site, hay)[[1]][gregexpr(s$site, hay)[[1]] > 0])
  }, integer(1))
  if (all(counts == 1L)) break
}

wrap_gb <- function(seq) {
  out <- character(0)
  for (i in seq(1, nchar(seq), by = 60)) {
    chunk <- substr(seq, i, min(i + 59, nchar(seq)))
    words <- substring(chunk, seq(1, nchar(chunk), 10),
                       pmin(seq(10, nchar(chunk) + 9, 10), nchar(chunk)))
    out <- c(out, sprintf("%9d %s", i, tolower(paste(words, collapse = " "))))
  }
  out
}

gb <- c(
  sprintf("LOCUS       pDemo                   %d bp ds-DNA     circular SYN 01-JAN-2026", LEN),
  "DEFINITION  Synthetic demo plasmid for plasmidplot tests.",
  "ACCESSION   .",
  "KEYWORDS    .",
  "FEATURES             Location/Qualifiers",
  sprintf("     source          1..%d", LEN),
  '                     /organism="synthetic construct"',
  "     CDS             120..380",
  '                     /gene="gfp"',
  "                     /label=GFP",
  '                     /ApEinfo_fwdcolor="#1baf7a"',
  "     CDS             complement(430..740)",
  '                     /gene="kan"',
  "                     /label=KanR",
  "     rep_origin      800..1000",
  "                     /label=ori",
  "                     /note=\"high-copy-number origin of",
  "                     replication\"",
  "     promoter        60..110",
  "                     /label=T7 promoter",
  sprintf("     misc_feature    join(%d..%d,1..40)", LEN - 60, LEN),
  "                     /label=origin-spanning region",
  "ORIGIN",
  wrap_gb(seq),
  "//"
)
writeLines(gb, "inst/extdata/pDemo.gb")

wrap_embl <- function(seq) {
  out <- character(0)
  for (i in seq(1, nchar(seq), by = 60)) {
    chunk <- substr(seq, i, min(i + 59, nchar(seq)))
    words <- substring(chunk, seq(1, nchar(chunk), 10),
                       pmin(seq(10, nchar(chunk) + 9, 10), nchar(chunk)))
    out <- c(out, sprintf("     %-66s%6d",
                          tolower(paste(words, collapse = " ")),
                          min(i + 59, nchar(seq))))
  }
  out
}

embl <- c(
  sprintf("ID   pDemo; SV 1; circular; DNA; STD; SYN; %d BP.", LEN),
  "XX",
  "DE   Synthetic demo plasmid for plasmidplot tests.",
  "XX",
  "FH   Key             Location/Qualifiers",
  "FH",
  sprintf("FT   source          1..%d", LEN),
  'FT                   /organism="synthetic construct"',
  "FT   CDS             120..380",
  "FT                   /label=GFP",
  "FT   CDS             complement(430..740)",
  "FT                   /label=KanR",
  "FT   rep_origin      800..1000",
  "FT                   /label=ori",
  "XX",
  sprintf("SQ   Sequence %d BP;", LEN),
  wrap_embl(seq),
  "//"
)
writeLines(embl, "inst/extdata/pDemo.embl")

writeLines(c(">pDemo synthetic demo plasmid",
             substring(seq, seq(1, nchar(seq), 70),
                       pmin(seq(70, nchar(seq) + 69, 70), nchar(seq)))),
           "inst/extdata/pDemo.fa")

cat("wrote pDemo.gb / .embl / .fa,", LEN, "bp\n")
cat("cut positions:",
    paste(sprintf("%s=%d", names(implants),
                  vapply(implants, function(s) s$at, numeric(1))),
          collapse = " "), "\n")
