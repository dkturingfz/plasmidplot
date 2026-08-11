ext <- function(f) system.file("extdata", f, package = "plasmidplot")

test_that("detect_format reads content, not the extension", {
  d <- plasmidplot:::detect_format
  expect_equal(d(ext("pDemo.gb")), "genbank")
  expect_equal(d(ext("pDemo.embl")), "embl")
  expect_equal(d(ext("pDemo.fa")), "fasta")

  # A GenBank record saved under .dna is still GenBank.
  tf <- tempfile(fileext = ".dna")
  on.exit(unlink(tf), add = TRUE)
  file.copy(ext("pDemo.gb"), tf, overwrite = TRUE)
  expect_equal(d(tf), "genbank")

  # A bare sequence under .dna is plain sequence, not SnapGene.
  tf2 <- tempfile(fileext = ".dna")
  on.exit(unlink(tf2), add = TRUE)
  writeLines(c("ACGTACGTAC", "GTACGTACGT"), tf2)
  expect_equal(d(tf2), "raw")

  tf3 <- tempfile()
  on.exit(unlink(tf3), add = TRUE)
  writeLines("this is just prose, not a sequence", tf3)
  expect_error(d(tf3), "Unrecognized file format")
  expect_error(d(tempfile()), "File not found")
})

test_that("a binary file that is not SnapGene fails with a useful message", {
  # Reading arbitrary bytes as text throws an encoding error that says
  # nothing about the actual problem, so binary is rejected up front.
  tf <- tempfile(fileext = ".dna")
  on.exit(unlink(tf), add = TRUE)
  set.seed(1)
  writeBin(c(as.raw(0L), as.raw(sample(0:255, 400, TRUE))), tf)
  expect_error(plasmidplot:::detect_format(tf), "Unrecognized file format")
  expect_error(read_plasmid(tf), "Unrecognized file format")
})

test_that("a GenBank file with CRLF endings reads normally", {
  tf <- tempfile(fileext = ".gb")
  on.exit(unlink(tf), add = TRUE)
  con <- file(tf, "wb")
  writeLines(readLines(ext("pDemo.gb"), warn = FALSE), con, sep = "\r\n")
  close(con)
  p <- read_genbank(tf)
  expect_equal(p$length, 1200)
  expect_gt(length(p$markers), 0L)
})

test_that("detect_format spots a SnapGene file by its magic bytes", {
  tf <- tempfile(fileext = ".txt")
  on.exit(unlink(tf), add = TRUE)
  con <- file(tf, "wb")
  writeBin(as.raw(9L), con)
  writeBin(14L, con, size = 4L, endian = "big")
  writeBin(c(charToRaw("SnapGene"), as.raw(c(0, 1, 0, 15, 0, 15))), con)
  close(con)
  expect_equal(plasmidplot:::detect_format(tf), "snapgene")
})

test_that("read_fasta reads header, sequence and length", {
  p <- read_fasta(ext("pDemo.fa"))
  expect_equal(p$name, "pDemo")
  expect_equal(p$length, 1200)
  expect_equal(nchar(p$sequence), 1200)
  expect_equal(p$topology, "circular")
  expect_length(p$markers, 0L)
  expect_equal(read_fasta(ext("pDemo.fa"), circular = FALSE)$topology, "linear")
  expect_equal(read_fasta(ext("pDemo.fa"), name = "other")$name, "other")
})

test_that("read_fasta handles a headerless sequence and multiple records", {
  tf <- tempfile(fileext = ".seq")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("acgtacgtac", "gtacgtacgt"), tf)
  p <- read_fasta(tf)
  expect_equal(p$length, 20)
  expect_equal(p$sequence, "ACGTACGTACGTACGTACGT")
  expect_equal(p$name, tools::file_path_sans_ext(basename(tf)))

  tf2 <- tempfile(fileext = ".fa")
  on.exit(unlink(tf2), add = TRUE)
  writeLines(c(">first", "ACGTACGT", ">second", "TTTTTTTTTTTT"), tf2)
  expect_warning(p2 <- read_fasta(tf2), "2 records")
  expect_equal(p2$name, "first")
  expect_equal(p2$length, 8)
})

test_that("read_embl parses the header, features and sequence", {
  p <- read_embl(ext("pDemo.embl"))
  expect_equal(p$name, "pDemo")
  expect_equal(p$length, 1200)
  expect_equal(p$topology, "circular")
  expect_equal(nchar(p$sequence), 1200)

  labs <- vapply(p$markers, function(m) m$label, character(1))
  expect_setequal(labs, c("GFP", "KanR", "ori"))
  expect_equal(p$markers[[which(labs == "KanR")]]$arrow, "start")
  expect_equal(p$markers[[which(labs == "ori")]]$arrow, "none")
  expect_error(read_embl(ext("pDemo.fa")), "No ID line")
})

test_that("GenBank and EMBL of the same plasmid agree", {
  g <- read_genbank(ext("pDemo.gb"))
  e <- read_embl(ext("pDemo.embl"))
  expect_equal(g$length, e$length)
  expect_equal(g$sequence, e$sequence)
  for (lab in c("GFP", "KanR", "ori")) {
    gm <- g$markers[[which(vapply(g$markers, function(m) m$label,
                                  character(1)) == lab)]]
    em <- e$markers[[which(vapply(e$markers, function(m) m$label,
                                  character(1)) == lab)]]
    expect_equal(c(gm$start, gm$end), c(em$start, em$end))
    expect_equal(gm$arrow, em$arrow)
  }
})

test_that("read_plasmid dispatches every format", {
  expect_equal(read_plasmid(ext("pDemo.gb"))$length, 1200)
  expect_equal(read_plasmid(ext("pDemo.embl"))$length, 1200)
  expect_equal(read_plasmid(ext("pDemo.fa"))$length, 1200)
  expect_equal(read_plasmid(ext("pBR322.gb"))$name, "pBR322")
})

test_that("a sequence that contradicts the declared length is refused", {
  # Truncating the record would otherwise leave every computed position
  # silently wrong, so the sequence is dropped with a warning.
  lines <- readLines(ext("pDemo.gb"), warn = FALSE)
  origin <- grep("^ORIGIN", lines)
  tf <- tempfile(fileext = ".gb")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c(lines[seq_len(origin + 3L)], "//"), tf)

  expect_warning(p <- read_genbank(tf), "declares")
  expect_true(is.na(p$sequence))
  expect_equal(p$length, 1200)
})

test_that("a features-only record with no ORIGIN is fine", {
  p <- read_genbank(ext("pBR322.gb"))
  expect_true(is.na(p$sequence))
  expect_equal(p$length, 4361)
  expect_gt(length(p$markers), 0L)
})

test_that("sequence = FALSE skips the sequence", {
  expect_true(is.na(read_genbank(ext("pDemo.gb"), sequence = FALSE)$sequence))
  expect_true(is.na(read_embl(ext("pDemo.embl"), sequence = FALSE)$sequence))
})
