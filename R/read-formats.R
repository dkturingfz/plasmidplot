# Extensions lie: ".dna" is SnapGene's, but plain-sequence exports and even
# GenBank records get saved under it too. Sniff the content instead.
detect_format <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  con <- file(file, "rb")
  head <- readBin(con, "raw", n = 64L)
  close(con)
  if (length(head) >= 13L && head[1] == as.raw(9L) &&
      identical(head[6:13], charToRaw("SnapGene"))) {
    return("snapgene")
  }
  # A NUL byte means this is binary and not any of the text formats. Bail out
  # here: reading it as text throws an encoding error that tells the user
  # nothing about what actually went wrong.
  if (any(head == as.raw(0L))) stop(unrecognized_msg(file), call. = FALSE)

  lines <- tryCatch(readLines(file, warn = FALSE, n = 200L),
                    error = function(e) character(0))
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("Empty file: ", file, call. = FALSE)
  # useBytes throughout: an undeclared or mis-declared encoding must not turn
  # format detection into an error.
  if (any(grepl("^LOCUS", lines, useBytes = TRUE))) return("genbank")
  if (any(grepl("^ID   ", lines, useBytes = TRUE)) &&
      any(grepl("^(FT|SQ)", lines, useBytes = TRUE))) {
    return("embl")
  }
  if (grepl("^>", lines[1], useBytes = TRUE)) return("fasta")
  # No header of any kind: accept it if it reads as bare sequence.
  body <- paste(lines, collapse = "")
  letters_only <- gsub("[^A-Za-z]", "", body, useBytes = TRUE)
  if (nchar(letters_only) > 0 &&
      nchar(gsub("[ACGTUNRYSWKMBDHVacgtunryswkmbdhv]", "", letters_only,
                 useBytes = TRUE)) == 0) {
    return("raw")
  }
  stop(unrecognized_msg(file), call. = FALSE)
}

unrecognized_msg <- function(file) {
  paste0("Unrecognized file format: ", file,
         " (expected GenBank, EMBL, FASTA, plain sequence, or SnapGene)")
}

#' Read a plasmid from any supported file
#'
#' Detects the format from the file's content rather than its extension, so a
#' `.dna` file is read as SnapGene or as plain sequence depending on what it
#' actually contains, and a GenBank record saved as `.txt` still works.
#'
#' Recognized formats: SnapGene (`.dna` binary), GenBank, EMBL, FASTA, and a
#' bare sequence with no header at all.
#'
#' @param file Path to the file.
#' @param ... Passed to the underlying reader.
#' @return A `plasmid` object.
#' @seealso [read_genbank()], [read_snapgene()], [read_embl()], [read_fasta()]
#' @examples
#' gb <- system.file("extdata", "pBR322.gb", package = "plasmidplot")
#' plot(read_plasmid(gb))
#' @export
read_plasmid <- function(file, ...) {
  switch(detect_format(file),
    snapgene = read_snapgene(file, ...),
    genbank = read_genbank(file, ...),
    embl = read_embl(file, ...),
    fasta = ,
    raw = read_fasta(file, ...)
  )
}

#' Read a plasmid from a FASTA file or a bare sequence
#'
#' FASTA carries no feature table, so the result is a bare backbone: name,
#' length and sequence, with no markers. Add features with [pp_marker()], or
#' restriction sites with [pp_find_sites()], which needs exactly this.
#'
#' @param file Path to a FASTA file, or a file holding only sequence letters.
#' @param name Plasmid name. Defaults to the FASTA header (up to its first
#'   space), or the file name when there is no header.
#' @param circular Whether the sequence is circular (default `TRUE`). FASTA
#'   records no topology.
#' @param ... Ignored, so [read_plasmid()] can pass reader arguments through.
#' @return A `plasmid` object.
#' @examples
#' fa <- tempfile(fileext = ".fa")
#' writeLines(c(">pMini test plasmid", strrep("ATGC", 500)), fa)
#' p <- read_fasta(fa)
#' p
#' unlink(fa)
#' @export
read_fasta <- function(file, name = NULL, circular = TRUE, ...) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("Empty file: ", file, call. = FALSE)
  hdr <- grep("^>", lines)
  if (length(hdr) > 1L) {
    warning("FASTA holds ", length(hdr),
            " records; using the first. Split the file to plot the others.",
            call. = FALSE)
    lines <- lines[seq(hdr[1], hdr[2] - 1L)]
    hdr <- 1L
  }
  if (length(hdr)) {
    title <- trimws(sub("^>", "", lines[hdr[1]]))
    body <- lines[-seq_len(hdr[1])]
  } else {
    title <- ""
    body <- lines
  }
  seq <- toupper(gsub("[^A-Za-z]", "", paste(body, collapse = "")))
  if (!nchar(seq)) stop("No sequence found in ", file, call. = FALSE)

  nm <- name
  if (is.null(nm)) {
    nm <- sub("\\s.*$", "", title)
    if (!nzchar(nm)) nm <- tools::file_path_sans_ext(basename(file))
  }
  p <- plasmid(nm, nchar(seq))
  p$topology <- if (circular) "circular" else "linear"
  p$sequence <- seq
  p
}

#' Read a plasmid from an EMBL file
#'
#' EMBL's feature table uses the same location grammar as GenBank, so
#' `complement(...)`, `join(...)` and partial boundaries behave identically;
#' only the line prefixes differ.
#'
#' @inheritParams read_genbank
#' @return A `plasmid` object.
#' @seealso [read_genbank()], [read_plasmid()]
#' @examples
#' embl <- system.file("extdata", "pDemo.embl", package = "plasmidplot")
#' p <- read_embl(embl)
#' p
#' @export
read_embl <- function(file, name = NULL, types = NULL,
                      skip_types = "source",
                      label_from = c("label", "gene", "product",
                                     "standard_name", "note", "locus_tag"),
                      colors = c("style", "file"),
                      color_by = c("feature", "type"), arrows = TRUE,
                      sequence = TRUE) {
  colors <- match.arg(colors)
  color_by <- match.arg(color_by)
  lines <- readLines(file, warn = FALSE)
  id <- grep("^ID   ", lines)
  if (!length(id)) {
    stop("No ID line found; is this an EMBL file? ", file, call. = FALSE)
  }
  fields <- trimws(strsplit(sub("^ID   ", "", lines[id[1]]), ";")[[1]])
  nm <- sub("\\s.*$", "", fields[1])
  bp <- regmatches(lines[id[1]], regexpr("([0-9]+)\\s*BP", lines[id[1]],
                                         ignore.case = TRUE))
  len <- if (length(bp)) {
    as.numeric(sub("\\s*BP", "", bp, ignore.case = TRUE))
  } else NA_real_

  seq <- embl_sequence(lines)
  if (is.na(len) || len <= 0) len <- if (is.na(seq)) NA_real_ else nchar(seq)
  if (is.na(len) || len <= 0) {
    stop("Could not determine sequence length from ", file, call. = FALSE)
  }

  p <- plasmid(if (is.null(name)) nm else name, len)
  p$topology <- if (any(grepl("circular", fields, ignore.case = TRUE))) {
    "circular"
  } else "linear"
  if (sequence) p$sequence <- checked_sequence(seq, len)

  # Strip the "FT" prefix and the feature table becomes GenBank-shaped.
  ft <- grep("^FT", lines, value = TRUE)
  gb_like <- paste0("     ", substring(ft, 6))
  add_parsed_features(p, parse_gb_features(c("FEATURES", gb_like)),
                      types, skip_types, label_from, colors, arrows, color_by)
}

# EMBL sequence lives between the SQ line and the terminating "//", with a
# residue counter at the end of every line.
embl_sequence <- function(lines) {
  sq <- grep("^SQ ", lines)
  if (!length(sq)) return(NA_character_)
  body <- lines[(sq[1] + 1L):length(lines)]
  end <- grep("^//", body)
  if (length(end)) body <- body[seq_len(end[1] - 1L)]
  seq <- toupper(gsub("[^A-Za-z]", "", paste(body, collapse = "")))
  if (nzchar(seq)) seq else NA_character_
}
