gb_path <- function() {
  system.file("extdata", "pBR322.gb", package = "plasmidplot")
}

test_that("parse_location handles the GenBank location grammar", {
  pl <- plasmidplot:::parse_location
  expect_equal(pl("86..1276")[c("start", "end", "complement")],
               list(start = 86, end = 1276, complement = FALSE))
  expect_equal(pl("complement(3293..4153)")$complement, TRUE)
  expect_equal(pl("complement(3293..4153)")$start, 3293)
  # Partial boundaries are markers, not part of the number.
  expect_equal(pl("<1..>100")[c("start", "end")], list(start = 1, end = 100))
  # A join is taken first-segment start to last-segment end, so a feature
  # written in transcription order across the origin stays a wrap.
  expect_equal(pl("join(4300..4361,1..40)")[c("start", "end")],
               list(start = 4300, end = 40))
  expect_equal(pl("complement(join(4300..4361,1..40))")$complement, TRUE)
  expect_equal(pl("55")[c("start", "end")], list(start = 55, end = 55))
  expect_null(pl(""))
  expect_null(pl("unknown"))
})

test_that("read_genbank reads the LOCUS header", {
  p <- read_genbank(gb_path())
  expect_s3_class(p, "plasmid")
  expect_equal(p$name, "pBR322")
  expect_equal(p$length, 4361)
  expect_equal(p$topology, "circular")
})

test_that("read_genbank reads features, strands and wraps", {
  p <- read_genbank(gb_path())
  labs <- vapply(p$markers, function(m) m$label, character(1))
  # "source" spans the whole plasmid and is dropped by default.
  expect_false("Cloning vector pBR322" %in% labs)
  expect_true(all(c("TcR", "rop", "ori", "AmpR") %in% labs))

  tcr <- p$markers[[which(labs == "TcR")]]
  expect_equal(c(tcr$start, tcr$end), c(86, 1276))
  expect_equal(tcr$arrow, "end")
  expect_equal(tcr$type, "CDS")

  amp <- p$markers[[which(labs == "AmpR")]]
  expect_equal(amp$arrow, "start")

  # rep_origin is a site, not a transcribed unit: no arrowhead.
  expect_equal(p$markers[[which(labs == "ori")]]$arrow, "none")

  wrap <- p$markers[[which(labs == "origin-spanning region")]]
  expect_equal(c(wrap$start, wrap$end), c(4300, 40))
})

test_that("read_genbank honors label priority and qualifier wrapping", {
  p <- read_genbank(gb_path(), label_from = "gene")
  labs <- vapply(p$markers, function(m) m$label, character(1))
  # /gene wins where present; the feature type is the fallback.
  expect_true("tet" %in% labs)
  expect_true("bla" %in% labs)
  expect_true("rep_origin" %in% labs)

  # A /note wrapped over two lines is rejoined rather than truncated.
  feats <- plasmidplot:::parse_gb_features(readLines(gb_path(), warn = FALSE))
  ori <- Filter(function(f) identical(f$type, "rep_origin"), feats)[[1]]
  expect_match(plasmidplot:::unquote(ori$quals$note), "origin of replication$")
})

test_that("read_genbank filters by type", {
  p <- read_genbank(gb_path(), types = "CDS")
  expect_true(all(vapply(p$markers, function(m) m$type, character(1)) == "CDS"))
  expect_length(p$markers, 3L)

  p2 <- read_genbank(gb_path(), skip_types = c("source", "CDS"))
  expect_false("CDS" %in% vapply(p2$markers, function(m) m$type, character(1)))
})

test_that("color_by chooses between per-feature and per-type coloring", {
  s <- pp_style("angular")
  # Default: the three CDS features get three different colors, so adjacent
  # genes stay tellable apart.
  p <- read_genbank(gb_path())
  labs <- vapply(p$markers, function(m) m$label, character(1))
  cols <- plasmidplot:::marker_colors(p, s)
  cds <- which(labs %in% c("TcR", "rop", "AmpR"))
  expect_length(unique(cols[cds]), 3L)

  # color_by = "type": all three CDS share one color, ori differs.
  q <- read_genbank(gb_path(), color_by = "type")
  qlabs <- vapply(q$markers, function(m) m$label, character(1))
  qcols <- plasmidplot:::marker_colors(q, s)
  qcds <- which(qlabs %in% c("TcR", "rop", "AmpR"))
  expect_length(unique(qcols[qcds]), 1L)
  expect_false(qcols[qcds[1]] == qcols[which(qlabs == "ori")])
})

test_that("colors = 'file' picks up ApEinfo colors, 'style' ignores them", {
  p <- read_genbank(gb_path(), colors = "file")
  labs <- vapply(p$markers, function(m) m$label, character(1))
  expect_equal(p$markers[[which(labs == "TcR")]]$color, "#ff9ccd")
  # AmpR is on the complement strand, so its /ApEinfo_revcolor is used.
  expect_equal(p$markers[[which(labs == "AmpR")]]$color, "#ffef86")
  # A feature with no color qualifier still falls back to the palette.
  expect_true(is.na(p$markers[[which(labs == "ori")]]$color))

  q <- read_genbank(gb_path(), colors = "style")
  expect_true(all(is.na(vapply(q$markers, function(m) m$color, character(1)))))
})

test_that("read_genbank rejects a file with no LOCUS line", {
  tf <- tempfile(fileext = ".gb")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("not a genbank file", "at all"), tf)
  expect_error(read_genbank(tf), "No LOCUS line")
})

# --- SnapGene ---------------------------------------------------------------

# Build a minimal but structurally real .dna file: segments of one type byte,
# four big-endian length bytes, then the payload.
write_dna <- function(path, seq_len = 2000, circular = TRUE, xml = NULL) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  put <- function(type, payload) {
    writeBin(as.raw(type), con)
    writeBin(length(payload), con, size = 4L, endian = "big")
    if (length(payload)) writeBin(payload, con)
  }
  header <- c(charToRaw("SnapGene"), as.raw(c(0, 1, 0, 15, 0, 15)))
  put(9L, header)
  put(0L, c(as.raw(if (circular) 1L else 0L),
            charToRaw(paste(rep("a", seq_len), collapse = ""))))
  if (is.null(xml)) {
    xml <- paste0(
      '<Features nextValidID="3">',
      '<Feature recentID="0" name="AmpR" type="CDS" directionality="2">',
      '<Segment range="1200-1900" color="#ccffcc" type="standard"/>',
      '<Q name="gene"><V text="bla"/></Q>',
      '</Feature>',
      '<Feature recentID="1" name="ori" type="rep_origin">',
      '<Segment range="200-700" color="#ffff00" type="standard"/>',
      '</Feature>',
      '<Feature recentID="2" name="hub" type="promoter" directionality="3">',
      '<Segment range="800-900" color="#31b6c9" type="standard"/>',
      '</Feature>',
      '</Features>')
  }
  put(10L, charToRaw(xml))
  path
}

test_that("read_snapgene reads the stored .dna fixture", {
  skip_if_not_installed("xml2")
  # A file on disk, not one the test just built, so the binary reader is
  # exercised against bytes that were written once and left alone.
  dna <- system.file("extdata", "pDemo.dna", package = "plasmidplot")
  p <- read_snapgene(dna)
  expect_equal(p$length, 1200)
  expect_equal(p$topology, "circular")
  expect_equal(nchar(p$sequence), 1200)

  labs <- vapply(p$markers, function(m) m$label, character(1))
  expect_setequal(labs, c("GFP", "KanR", "ori", "T7 promoter"))
  expect_equal(p$markers[[which(labs == "GFP")]]$arrow, "end")
  expect_equal(p$markers[[which(labs == "KanR")]]$arrow, "start")
  expect_equal(p$markers[[which(labs == "T7 promoter")]]$arrow, "both")
  expect_equal(p$markers[[which(labs == "ori")]]$arrow, "none")

  # The same molecule as the GenBank fixture, so the two readers must agree.
  gb <- read_genbank(system.file("extdata", "pDemo.gb", package = "plasmidplot"))
  expect_equal(p$sequence, gb$sequence)
  for (lab in c("GFP", "KanR", "ori")) {
    dm <- p$markers[[which(labs == lab)]]
    glabs <- vapply(gb$markers, function(m) m$label, character(1))
    gm <- gb$markers[[which(glabs == lab)]]
    expect_equal(c(dm$start, dm$end), c(gm$start, gm$end), info = lab)
  }
})

test_that("read_snapgene parses segments, topology and features", {
  skip_if_not_installed("xml2")
  tf <- tempfile(fileext = ".dna")
  on.exit(unlink(tf), add = TRUE)
  write_dna(tf, seq_len = 2000, circular = TRUE)

  p <- read_snapgene(tf)
  expect_s3_class(p, "plasmid")
  expect_equal(p$length, 2000)
  expect_equal(p$topology, "circular")
  expect_equal(p$name, tools::file_path_sans_ext(basename(tf)))

  labs <- vapply(p$markers, function(m) m$label, character(1))
  expect_setequal(labs, c("AmpR", "ori", "hub"))

  amp <- p$markers[[which(labs == "AmpR")]]
  expect_equal(c(amp$start, amp$end), c(1200, 1900))
  expect_equal(amp$arrow, "start")   # directionality 2 = reverse strand
  expect_equal(amp$type, "CDS")

  # directionality 3 = bidirectional
  expect_equal(p$markers[[which(labs == "hub")]]$arrow, "both")
  # rep_origin is non-directional regardless of the file
  expect_equal(p$markers[[which(labs == "ori")]]$arrow, "none")
})

test_that("read_snapgene honors linear topology, colors and name override", {
  skip_if_not_installed("xml2")
  tf <- tempfile(fileext = ".dna")
  on.exit(unlink(tf), add = TRUE)
  write_dna(tf, seq_len = 2000, circular = FALSE)

  p <- read_snapgene(tf, name = "myPlasmid", colors = "file")
  expect_equal(p$topology, "linear")
  expect_equal(p$name, "myPlasmid")
  labs <- vapply(p$markers, function(m) m$label, character(1))
  expect_equal(p$markers[[which(labs == "ori")]]$color, "#ffff00")

  q <- read_snapgene(tf, colors = "style")
  expect_true(all(is.na(vapply(q$markers, function(m) m$color, character(1)))))
})

test_that("features outside the sequence are dropped, not clamped", {
  skip_if_not_installed("xml2")
  tf <- tempfile(fileext = ".dna")
  on.exit(unlink(tf), add = TRUE)
  # Features sit at 200-700, 800-900 and 1200-1900; only the first two fit.
  write_dna(tf, seq_len = 1000)
  p <- read_snapgene(tf)
  expect_setequal(vapply(p$markers, function(m) m$label, character(1)),
                  c("ori", "hub"))
})

test_that("read_snapgene rejects files without a SnapGene header", {
  skip_if_not_installed("xml2")
  tf <- tempfile(fileext = ".dna")
  on.exit(unlink(tf), add = TRUE)
  writeBin(charToRaw("this is not a SnapGene file at all"), tf)
  expect_error(read_snapgene(tf), "Not a SnapGene")
  expect_error(read_snapgene(tempfile()), "File not found")
})

test_that("read_plasmid dispatches on extension", {
  p <- read_plasmid(gb_path())
  expect_equal(p$name, "pBR322")
})

test_that("an imported plasmid plots in every style", {
  p <- read_genbank(gb_path())
  tf <- tempfile(fileext = ".png")
  png(tf, width = 600, height = 600)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  for (nm in pp_style()) expect_no_error(plot(p, style = nm))
})
