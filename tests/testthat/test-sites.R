demo_gb <- function() system.file("extdata", "pDemo.gb", package = "plasmidplot")

test_that("pp_site validates and appends", {
  p <- plasmid("pX", 1000)
  p <- pp_site(p, 500, "EcoRI")
  expect_length(p$sites, 1L)
  expect_equal(p$sites[[1]]$position, 500)
  expect_error(pp_site(p, 0, "x"), "Position must lie")
  expect_error(pp_site(p, 1001, "x"), "Position must lie")
  expect_error(pp_site(p, c(1, 2)), "single number")
})

test_that("pp_features adds a whole data frame", {
  df <- data.frame(
    start = c(100, 400), end = c(300, 600),
    label = c("a", "b"), arrow = c("end", "none"),
    stringsAsFactors = FALSE
  )
  p <- pp_features(plasmid("pX", 1000), df)
  expect_length(p$markers, 2L)
  expect_equal(p$markers[[1]]$label, "a")
  expect_equal(p$markers[[1]]$arrow, "end")
  expect_equal(p$markers[[2]]$arrow, "none")
  # Unknown columns are ignored rather than fatal.
  df$notes <- c("x", "y")
  expect_no_error(pp_features(plasmid("pX", 1000), df))
  expect_error(pp_features(plasmid("pX", 1000), data.frame(start = 1)),
               "needs column")
})

test_that("iupac_regex expands ambiguity codes and rejects junk", {
  expect_equal(plasmidplot:::iupac_regex("GAATTC"), "GAATTC")
  expect_equal(plasmidplot:::iupac_regex("GGWCC"), "GG[AT]CC")
  expect_equal(plasmidplot:::iupac_regex("GGNCC"), "GG.CC")
  expect_error(plasmidplot:::iupac_regex("GGZCC"), "Z")
})

test_that("pp_enzymes recognition sequences and cut offsets are consistent", {
  tbl <- pp_enzymes()
  expect_true(all(nchar(tbl$site) >= 4L))
  # A cut offset must fall inside its own recognition sequence.
  expect_true(all(tbl$cut >= 0L & tbl$cut <= nchar(tbl$site)))
  expect_true(all(grepl("^[ACGTRYSWKMBDHVN]+$", tbl$site)))
  expect_false(anyDuplicated(tbl$enzyme) > 0)
})

test_that("pp_find_sites locates the implanted sites at the right positions", {
  p <- read_genbank(demo_gb())
  expect_equal(nchar(p$sequence), p$length)

  q <- pp_find_sites(p, c("EcoRI", "BamHI", "HindIII", "NotI"),
                     show_position = FALSE)
  labs <- vapply(q$sites, function(s) s$label, character(1))
  pos <- vapply(q$sites, function(s) s$position, numeric(1))
  expect_setequal(labs, c("EcoRI", "BamHI", "HindIII", "NotI"))
  # The fixture implants the recognition sequences at these starts. A site
  # is numbered by the last base before its top-strand cut.
  expect_equal(pos[labs == "EcoRI"], 100)          # G^AATTC
  expect_equal(pos[labs == "BamHI"], 400)          # G^GATCC
  expect_equal(pos[labs == "HindIII"], 760)        # A^AGCTT
  expect_equal(pos[labs == "NotI"], 1010 + 1)      # GC^GGCCGC
})

test_that("site numbering matches published pUC19 maps", {
  # The pUC19 polylinker from EcoRI to KpnI, placed so EcoRI starts at 396 as
  # it does in the real plasmid. SnapGene and NEB both print EcoRI (396),
  # SacI (406) and KpnI (412) for it.
  seq <- paste0(strrep("A", 395), "GAATTCGAGCTCGGTACC", strrep("A", 200))
  p <- plasmid("mcs", nchar(seq))
  p$sequence <- seq
  q <- pp_find_sites(p, c("EcoRI", "SacI", "KpnI"))
  expect_setequal(vapply(q$sites, function(s) s$label, character(1)),
                  c("EcoRI (396)", "SacI (406)", "KpnI (412)"))
  # A cut after the last base of the molecule wraps to it, not past it.
  seq <- paste0("C", strrep("A", 40), "GGTAC")
  p <- plasmid("wrapcut", nchar(seq))
  p$sequence <- seq
  q <- pp_find_sites(p, "KpnI", show_position = FALSE)
  expect_equal(q$sites[[1]]$position, nchar(seq))
})

test_that("pp_find_sites labels carry the position when asked", {
  p <- pp_find_sites(read_genbank(demo_gb()), "EcoRI")
  expect_equal(p$sites[[1]]$label, "EcoRI (100)")
})

test_that("pp_find_sites accepts custom recognition sequences", {
  p <- read_genbank(demo_gb())
  q <- pp_find_sites(p, c(MyEnz = "GAATTC"), show_position = FALSE)
  expect_equal(vapply(q$sites, function(s) s$label, character(1)), "MyEnz")
  # No cut offset is known for a custom site, so it reports the site start.
  expect_equal(q$sites[[1]]$position, 100)
})

test_that("pp_find_sites finds sites that straddle the origin", {
  # Put half of a GAATTC at the end and half at the start.
  seq <- paste0("ATTC", strrep("GGGG", 40), "GA")
  p <- plasmid("wrap", nchar(seq))
  p$sequence <- seq
  q <- pp_find_sites(p, "EcoRI", show_position = FALSE)
  expect_length(q$sites, 1L)
  # Recognition starts at the last two bases; the base before the G^AATTC
  # cut is the G itself.
  expect_equal(q$sites[[1]]$position, nchar(seq) - 1)

  # A linear molecule has no wrap, so the same sequence yields nothing.
  p$topology <- "linear"
  expect_length(pp_find_sites(p, "EcoRI")$sites, 0L)
})

test_that("unique_only keeps single cutters and max_sites bounds the rest", {
  seq <- paste0(strrep("A", 50), "GAATTC", strrep("A", 50), "GAATTC",
                strrep("A", 50))
  p <- plasmid("twice", nchar(seq))
  p$sequence <- seq
  expect_length(pp_find_sites(p, "EcoRI")$sites, 0L)
  expect_length(pp_find_sites(p, "EcoRI", unique_only = FALSE)$sites, 2L)
  expect_length(
    pp_find_sites(p, "EcoRI", unique_only = FALSE, max_sites = 1L)$sites, 0L)
})

test_that("overlapping recognition sites are all found", {
  # GCGCGCGC holds BssHII (GCGCGC) twice, two bases apart. Counting only the
  # first would report a double cutter as a unique one.
  seq <- paste0(strrep("A", 20), "GCGCGCGC", strrep("A", 20))
  p <- plasmid("overlap", nchar(seq))
  p$sequence <- seq
  expect_length(pp_find_sites(p, "BssHII")$sites, 0L)
  q <- pp_find_sites(p, "BssHII", unique_only = FALSE, show_position = FALSE)
  expect_equal(vapply(q$sites, function(s) s$position, numeric(1)),
               c(21, 23))                                    # G^CGCGC
})

test_that("custom non-palindromic sites are found on both strands", {
  # BsaI is GGTCTC; on the other strand the top strand reads GAGACC.
  fwd <- paste0(strrep("A", 10), "GGTCTC", strrep("A", 30))
  rev <- paste0(strrep("A", 30), "GAGACC", strrep("A", 10))
  for (case in list(list(seq = fwd, at = 11), list(seq = rev, at = 31))) {
    p <- plasmid("iis", nchar(case$seq))
    p$sequence <- case$seq
    q <- pp_find_sites(p, c(BsaI = "GGTCTC"), show_position = FALSE)
    expect_length(q$sites, 1L)
    expect_equal(q$sites[[1]]$position, case$at)
  }
  # One on each strand is two cuts, not a unique cutter.
  both <- paste0(fwd, rev)
  p <- plasmid("iis2", nchar(both))
  p$sequence <- both
  expect_length(pp_find_sites(p, c(BsaI = "GGTCTC"))$sites, 0L)
  expect_length(
    pp_find_sites(p, c(BsaI = "GGTCTC"), unique_only = FALSE)$sites, 2L)
  # A palindrome is the same site on either strand and is counted once.
  pal <- paste0(strrep("A", 10), "GAATTC", strrep("A", 30))
  p <- plasmid("pal", nchar(pal))
  p$sequence <- pal
  expect_length(pp_find_sites(p, c(MyEnz = "gaattc"))$sites, 1L)
})

test_that("revcomp_site complements ambiguity codes", {
  rc <- plasmidplot:::revcomp_site
  expect_equal(rc("GGTCTC"), "GAGACC")
  expect_equal(rc("GGWCC"), "GGWCC")
  expect_equal(rc("ACRYKMBDHVN"), "NBDHVKMRYGT")
  # The built-in table records a top-strand cut only, which is enough because
  # every entry is a palindrome. An asymmetric enzyme added here would need a
  # second cut offset before its reverse-strand sites could be placed.
  tbl <- pp_enzymes()
  expect_equal(vapply(tbl$site, rc, character(1), USE.NAMES = FALSE), tbl$site)
})

test_that("pp_find_sites reports missing sequences and unknown enzymes", {
  p <- read_genbank(demo_gb(), sequence = FALSE)
  expect_error(pp_find_sites(p), "no sequence")
  expect_error(pp_find_sites(read_genbank(demo_gb()), "NotAnEnzyme"),
               "Unknown enzyme")
})

test_that("a plasmid with no sequence prints and reports consistently", {
  # NA, NULL and "" all mean the same thing to a reader; none of them may
  # turn printing the object into an error.
  for (empty in list(NA_character_, NULL, "")) {
    p <- plasmid("x", 100)
    p$sequence <- empty
    expect_false(plasmidplot:::has_sequence(p))
    expect_output(print(p), "<plasmid> x")
    expect_error(pp_find_sites(p), "no sequence")
  }
  p <- plasmid("x", 6)
  p$sequence <- "GAATTC"
  expect_true(plasmidplot:::has_sequence(p))
  expect_output(print(p), "with sequence")
})

test_that("site search is case-insensitive", {
  p <- plasmid("lc", 60)
  p$sequence <- tolower(paste0(strrep("a", 27), "gaattc", strrep("a", 27)))
  expect_length(pp_find_sites(p, "EcoRI")$sites, 1L)
})

test_that("site labels share the layout with feature labels", {
  p <- read_genbank(demo_gb())
  p <- pp_find_sites(p, c("EcoRI", "BamHI"))
  s <- pp_style("angular")
  entries <- plasmidplot:::circ_label_entries(p, s)
  # Both kinds are present, and the site leaders start further out than the
  # feature leaders so their ticks are not crossed.
  expect_gt(length(entries), length(p$sites))
  anchors <- vapply(entries, function(e) e$anchor, numeric(1))
  texts <- vapply(entries, function(e) e$text, character(1))
  site_rows <- grepl("^(EcoRI|BamHI)", texts)
  expect_gt(min(anchors[site_rows]), max(anchors[!site_rows]))
})

test_that("plots with sites render in every style", {
  p <- pp_find_sites(read_genbank(demo_gb()))
  expect_gt(length(p$sites), 0L)
  tf <- tempfile(fileext = ".png")
  png(tf, width = 700, height = 700)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  for (nm in pp_style()) expect_no_error(plot(p, style = nm))
})
