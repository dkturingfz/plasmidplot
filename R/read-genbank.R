# Feature types drawn without an arrowhead even when the record gives them a
# strand: these are sites and regions rather than transcribed units.
pp_nondirectional <- c(
  "source", "rep_origin", "misc_feature", "protein_bind", "repeat_region",
  "misc_recomb", "oriT", "stem_loop", "variation", "misc_binding",
  "enhancer", "LTR", "centromere", "telomere", "STS", "gap"
)

#' Read a plasmid from a GenBank file
#'
#' Parses the `LOCUS` line for the name and length and the `FEATURES` block
#' for markers. Handles `complement(...)`, `join(...)`, and the `<`/`>`
#' partial-boundary markers. A `join()` whose first segment starts after its
#' last segment ends is kept as an origin-crossing feature.
#'
#' @param file Path to a `.gb`, `.gbk`, or `.genbank` file.
#' @param name Plasmid name. Defaults to the `LOCUS` identifier.
#' @param types Feature types to keep (e.g. `c("CDS", "rep_origin")`).
#'   `NULL` keeps everything not in `skip_types`.
#' @param skip_types Feature types to drop. `"source"` spans the whole
#'   plasmid and is dropped by default.
#' @param label_from Qualifier names to try, in order, when labeling a
#'   feature. Falls back to the feature type.
#' @param colors `"style"` takes colors from the plot style's palette.
#'   `"file"` honors `/ApEinfo_fwdcolor` qualifiers when the record has them
#'   (SnapGene and ApE write these), falling back to the style palette for
#'   features without one.
#' @param color_by `"feature"` gives every feature its own palette color, so
#'   neighbors stay distinguishable. `"type"` gives all features of one type
#'   the same color, which reads as a legend of feature classes instead.
#' @param arrows Draw arrowheads for stranded features (default `TRUE`).
#' @param sequence Keep the `ORIGIN` sequence on the result (default `TRUE`).
#'   [pp_find_sites()] needs it.
#' @return A `plasmid` object.
#' @seealso [read_snapgene()], [read_plasmid()], [pp_find_sites()]
#' @examples
#' gb <- system.file("extdata", "pBR322.gb", package = "plasmidplot")
#' p <- read_genbank(gb)
#' p
#' plot(p)
#' @export
read_genbank <- function(file, name = NULL, types = NULL,
                         skip_types = "source",
                         label_from = c("label", "gene", "product",
                                        "standard_name", "note", "locus_tag"),
                         colors = c("style", "file"),
                         color_by = c("feature", "type"), arrows = TRUE,
                         sequence = TRUE) {
  colors <- match.arg(colors)
  color_by <- match.arg(color_by)
  lines <- readLines(file, warn = FALSE)
  if (!length(lines)) stop("Empty file: ", file, call. = FALSE)

  hdr <- parse_locus(lines, file)
  feats <- parse_gb_features(lines)

  p <- plasmid(if (is.null(name)) hdr$name else name, hdr$length)
  p$topology <- hdr$topology
  if (sequence) p$sequence <- checked_sequence(origin_sequence(lines), hdr$length)
  add_parsed_features(p, feats, types, skip_types, label_from, colors, arrows,
                      color_by)
}

parse_locus <- function(lines, file) {
  i <- grep("^LOCUS", lines)
  if (!length(i)) stop("No LOCUS line found; is this a GenBank file? ", file,
                       call. = FALSE)
  ln <- lines[i[1]]
  rest <- trimws(sub("^LOCUS\\s*", "", ln))
  nm <- sub("\\s.*$", "", rest)
  bp <- regmatches(ln, regexpr("([0-9]+)\\s+bp", ln))
  len <- if (length(bp)) as.numeric(sub("\\s+bp", "", bp)) else NA_real_
  if (is.na(len)) len <- origin_length(lines)
  if (is.na(len) || len <= 0) {
    stop("Could not determine sequence length from ", file, call. = FALSE)
  }
  list(name = if (nzchar(nm)) nm else "plasmid", length = len,
       topology = if (grepl("circular", ln, ignore.case = TRUE)) "circular"
                  else "linear")
}

# The ORIGIN block interleaves residues with position counters; keep letters.
origin_sequence <- function(lines) {
  o <- grep("^ORIGIN", lines)
  if (!length(o) || o[1] == length(lines)) return(NA_character_)
  body <- lines[(o[1] + 1L):length(lines)]
  end <- grep("^//", body)
  if (length(end)) body <- body[seq_len(end[1] - 1L)]
  seq <- toupper(gsub("[^A-Za-z]", "", paste(body, collapse = "")))
  if (nzchar(seq)) seq else NA_character_
}

# Fallback when LOCUS carries no length.
origin_length <- function(lines) {
  seq <- origin_sequence(lines)
  if (is.na(seq)) NA_real_ else nchar(seq)
}

# A sequence shorter than the record claims (a truncated or excerpted file)
# would silently corrupt every position computed from it, so drop it rather
# than let pp_find_sites() report positions against the wrong length.
checked_sequence <- function(seq, len) {
  if (is.na(seq)) return(NA_character_)
  if (nchar(seq) != len) {
    warning("Sequence is ", nchar(seq), " bp but the record declares ",
            format(round(len)), " bp; ignoring the sequence. ",
            "Position-based functions will not work on this file.",
            call. = FALSE)
    return(NA_character_)
  }
  seq
}

# Collect one entry per feature: type, raw location string, qualifier list.
parse_gb_features <- function(lines) {
  start <- grep("^FEATURES", lines)
  if (!length(start)) return(list())
  body <- lines[(start[1] + 1L):length(lines)]
  stop_at <- grep("^(ORIGIN|CONTIG|//)", body)
  if (length(stop_at)) body <- body[seq_len(stop_at[1] - 1L)]

  feats <- list()
  cur <- NULL
  qual <- NULL
  flush_qual <- function() {
    if (!is.null(cur) && !is.null(qual)) {
      cur$quals[[qual$key]] <<- qual$value
    }
    qual <<- NULL
  }
  for (ln in body) {
    if (!nzchar(trimws(ln))) next
    indent <- nchar(ln) - nchar(sub("^ +", "", ln))
    if (indent > 0 && indent < 21) {
      # New feature: key at column 6, location from column 22.
      flush_qual()
      if (!is.null(cur)) feats[[length(feats) + 1L]] <- cur
      parts <- strsplit(trimws(ln), "\\s+")[[1]]
      cur <- list(type = parts[1],
                  location = paste(parts[-1], collapse = ""),
                  quals = list())
    } else if (!is.null(cur)) {
      txt <- trimws(ln)
      if (startsWith(txt, "/")) {
        flush_qual()
        kv <- sub("^/", "", txt)
        eq <- regexpr("=", kv, fixed = TRUE)
        if (eq > 0) {
          qual <- list(key = substr(kv, 1, eq - 1L),
                       value = substr(kv, eq + 1L, nchar(kv)))
        } else {
          qual <- list(key = kv, value = "")
        }
      } else if (!is.null(qual)) {
        # Wrapped qualifier value.
        qual$value <- paste0(qual$value, " ", txt)
      } else {
        # Wrapped location.
        cur$location <- paste0(cur$location, txt)
      }
    }
  }
  flush_qual()
  if (!is.null(cur)) feats[[length(feats) + 1L]] <- cur
  feats
}

# "complement(join(4000..4361,1..100))" -> start 4000, end 100, reverse strand.
parse_location <- function(loc) {
  if (!nzchar(loc)) return(NULL)
  complement <- grepl("complement", loc, fixed = TRUE)
  clean <- gsub("[<>]", "", loc)
  segs <- regmatches(clean, gregexpr("[0-9]+(\\.\\.[0-9]+)?", clean))[[1]]
  if (!length(segs)) return(NULL)
  nums_of <- function(s) as.numeric(strsplit(s, "..", fixed = TRUE)[[1]])
  first <- nums_of(segs[1])
  last <- nums_of(segs[length(segs)])
  list(start = first[1], end = last[length(last)], complement = complement)
}

unquote <- function(x) {
  x <- trimws(x)
  gsub('^"|"$', "", x)
}

# SnapGene stores a color on the feature itself; GenBank records written by
# SnapGene or ApE carry it as an /ApEinfo_*color qualifier instead.
file_color <- function(f, complement) {
  hex <- f$color
  if (is.null(hex)) {
    key <- if (complement) "ApEinfo_revcolor" else "ApEinfo_fwdcolor"
    hex <- f$quals[[key]]
    if (is.null(hex)) hex <- f$quals[["ApEinfo_fwdcolor"]]
  }
  if (is.null(hex) || is.na(hex)) return(NULL)
  hex <- unquote(hex)
  if (grepl("^#[0-9A-Fa-f]{6}$", hex)) hex else NULL
}

pick_label <- function(quals, label_from, type) {
  for (k in label_from) {
    if (!is.null(quals[[k]])) {
      v <- unquote(quals[[k]])
      if (nzchar(v)) return(v)
    }
  }
  type
}

# Shared by the GenBank and SnapGene readers: turn parsed records into markers.
add_parsed_features <- function(p, feats, types, skip_types, label_from,
                                colors, arrows, color_by = "feature") {
  for (f in feats) {
    if (!is.null(types) && !f$type %in% types) next
    if (f$type %in% skip_types) next
    lp <- parse_location(f$location)
    if (is.null(lp)) next
    if (lp$start < 1 || lp$end < 1 ||
        lp$start > p$length || lp$end > p$length) next

    col <- if (identical(colors, "file")) file_color(f, lp$complement) else NULL
    arrow <- "none"
    if (arrows && !f$type %in% pp_nondirectional) {
      arrow <- if (!is.null(f$arrow)) f$arrow
               else if (lp$complement) "start" else "end"
    }
    p <- pp_marker(p, lp$start, lp$end,
                   label = pick_label(f$quals, label_from, f$type),
                   color = col,
                   group = if (identical(color_by, "type")) f$type else NULL,
                   arrow = arrow)
    # Keep the type on the marker either way, so it stays available for
    # filtering and inspection even when it is not driving the color.
    p$markers[[length(p$markers)]]$type <- f$type
  }
  p
}
