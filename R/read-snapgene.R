# A SnapGene .dna file is a flat run of segments, each one byte of type, four
# big-endian bytes of length, then that many bytes of payload. The segments we
# need are 0 (the sequence, whose first byte holds the topology flags), 9 (the
# "SnapGene" header) and 10 (the feature table, which is XML).
snapgene_segments <- function(file) {
  con <- file(file, "rb")
  on.exit(close(con), add = TRUE)
  segs <- list()
  repeat {
    type <- readBin(con, "integer", n = 1L, size = 1L, signed = FALSE)
    if (!length(type)) break
    n <- readBin(con, "integer", n = 1L, size = 4L, endian = "big")
    if (!length(n) || is.na(n) || n < 0) break
    payload <- if (n > 0) readBin(con, "raw", n = n) else raw(0)
    if (length(payload) < n) break
    segs[[length(segs) + 1L]] <- list(type = type, data = payload)
  }
  segs
}

seg_of_type <- function(segs, type) {
  for (s in segs) if (s$type == type) return(s)
  NULL
}

# The header segment is "SnapGene" followed by binary version fields, so it
# must be matched as raw bytes: converting it to a string hits embedded NULs.
snapgene_header_ok <- function(seg) {
  !is.null(seg) && length(seg$data) >= 8L &&
    identical(seg$data[1:8], charToRaw("SnapGene"))
}

# XML payloads are sometimes NUL-padded to a block boundary; trim before
# converting, or rawToChar refuses the string.
seg_string <- function(seg) {
  if (is.null(seg)) return(NULL)
  bytes <- seg$data[seg$data != as.raw(0L)]
  if (!length(bytes)) return("")
  txt <- rawToChar(bytes)
  Encoding(txt) <- "UTF-8"
  txt
}

#' Read a plasmid from a SnapGene file
#'
#' Parses a SnapGene `.dna` file: sequence length and topology from the
#' sequence segment, and features (name, type, range, strand, color) from the
#' embedded feature table.
#'
#' SnapGene records carry no plasmid name, so the file's base name is used
#' unless `name` is given.
#'
#' @param file Path to a `.dna` file.
#' @param name Plasmid name. Defaults to the file name without its extension.
#' @param types Feature types to keep. `NULL` keeps everything not in
#'   `skip_types`.
#' @param skip_types Feature types to drop.
#' @param label_from Qualifier names to try, in order, when labeling a
#'   feature. `"name"` is the feature's own name attribute.
#' @param colors `"style"` takes colors from the plot style's palette;
#'   `"file"` keeps the colors stored in the SnapGene file.
#' @param color_by `"feature"` gives every feature its own palette color;
#'   `"type"` gives all features of one type the same color.
#' @param arrows Draw arrowheads for stranded features (default `TRUE`).
#'   SnapGene's bidirectional features become double-headed arrows.
#' @param sequence Keep the sequence on the result (default `TRUE`).
#'   [pp_find_sites()] needs it.
#' @return A `plasmid` object.
#' @note Requires the `xml2` package, which is only needed for this reader.
#' @seealso [read_genbank()], [read_plasmid()]
#' @examples
#' \dontrun{
#' p <- read_snapgene("pUC19.dna")
#' plot(p, style = "snapgene")
#' plot(read_snapgene("pUC19.dna", colors = "file"))
#' }
#' @export
read_snapgene <- function(file, name = NULL, types = NULL,
                          skip_types = character(),
                          label_from = c("label", "name", "gene", "product",
                                         "note"),
                          colors = c("style", "file"),
                          color_by = c("feature", "type"), arrows = TRUE,
                          sequence = TRUE) {
  colors <- match.arg(colors)
  color_by <- match.arg(color_by)
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("read_snapgene() needs the 'xml2' package. ",
         'Install it with install.packages("xml2").', call. = FALSE)
  }
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  segs <- snapgene_segments(file)
  if (!snapgene_header_ok(seg_of_type(segs, 9L))) {
    stop("Not a SnapGene .dna file (missing SnapGene header): ", file,
         call. = FALSE)
  }
  dna <- seg_of_type(segs, 0L)
  if (is.null(dna) || length(dna$data) < 2L) {
    stop("SnapGene file has no sequence segment: ", file, call. = FALSE)
  }
  len <- length(dna$data) - 1L
  circular <- bitwAnd(as.integer(dna$data[1]), 1L) == 1L

  nm <- if (!is.null(name)) name else tools::file_path_sans_ext(basename(file))
  p <- plasmid(nm, len)
  p$topology <- if (circular) "circular" else "linear"
  if (sequence) {
    p$sequence <- toupper(rawToChar(dna$data[-1]))
  }

  feats <- parse_snapgene_features(seg_string(seg_of_type(segs, 10L)))
  add_parsed_features(p, feats, types, skip_types, label_from, colors, arrows,
                      color_by)
}

parse_snapgene_features <- function(xml) {
  if (is.null(xml) || !nzchar(xml)) return(list())
  doc <- xml2::read_xml(xml)
  out <- list()
  for (node in xml2::xml_find_all(doc, "//Feature")) {
    segs <- xml2::xml_find_all(node, "./Segment")
    ranges <- xml2::xml_attr(segs, "range")
    ranges <- ranges[!is.na(ranges)]
    if (!length(ranges)) next
    bounds <- lapply(strsplit(ranges, "-", fixed = TRUE), as.numeric)
    first <- bounds[[1]]
    last <- bounds[[length(bounds)]]

    dir <- xml2::xml_attr(node, "directionality")
    arrow <- if (is.na(dir)) "none" else
      switch(dir, "1" = "end", "2" = "start", "3" = "both", "none")
    # Feed the shared marker builder a GenBank-shaped location so both
    # readers go through one code path.
    loc <- sprintf("%d..%d", as.integer(first[1]),
                   as.integer(last[length(last)]))
    if (identical(arrow, "start")) loc <- sprintf("complement(%s)", loc)

    quals <- list(name = xml2::xml_attr(node, "name"))
    for (q in xml2::xml_find_all(node, "./Q")) {
      key <- xml2::xml_attr(q, "name")
      val <- xml2::xml_attr(xml2::xml_find_first(q, "./V"), "text")
      if (!is.na(key) && !is.na(val)) quals[[key]] <- val
    }
    quals <- quals[!vapply(quals, function(v) is.null(v) || is.na(v),
                           logical(1))]

    type <- xml2::xml_attr(node, "type")
    out[[length(out) + 1L]] <- list(
      type = if (is.na(type)) "misc_feature" else type,
      location = loc,
      arrow = arrow,
      color = xml2::xml_attr(segs[[1]], "color"),
      quals = quals
    )
  }
  out
}

