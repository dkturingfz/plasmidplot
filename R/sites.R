#' Mark a single position on the plasmid
#'
#' Sites are drawn as a radial tick outside the feature ring with a label,
#' the way restriction sites are shown on a classic plasmid map. Unlike
#' [pp_marker()] they occupy one position rather than a span.
#'
#' @param p A `plasmid` object.
#' @param position Position in base pairs.
#' @param label Site label, e.g. an enzyme name.
#' @param color Tick color. Defaults to the style's `site_col`.
#' @return The updated `plasmid` object.
#' @seealso [pp_find_sites()] to add restriction sites from the sequence.
#' @examples
#' p <- plasmid("pUC19", 2686) |>
#'   pp_marker(146, 469, label = "lacZ", arrow = "end") |>
#'   pp_site(396, "EcoRI") |>
#'   pp_site(447, "BamHI")
#' plot(p)
#' @export
pp_site <- function(p, position, label = NULL, color = NULL) {
  stopifnot(inherits(p, "plasmid"))
  if (!is.numeric(position) || length(position) != 1L || is.na(position)) {
    stop("`position` must be a single number.", call. = FALSE)
  }
  if (position < 1 || position > p$length) {
    stop(sprintf("Position must lie in [1, %d].", round(p$length)),
         call. = FALSE)
  }
  p$sites[[length(p$sites) + 1L]] <- list(
    position = as.numeric(position),
    label = if (is.null(label)) NA_character_ else as.character(label),
    color = if (is.null(color)) NA_character_ else as.character(color)
  )
  p
}

#' Add many features at once from a data frame
#'
#' A vectorized [pp_marker()]: one row per feature. Only `start` and `end`
#' are required; any other column named after a [pp_marker()] argument is
#' used, and unknown columns are ignored.
#'
#' @param p A `plasmid` object.
#' @param df A data frame with columns `start` and `end`, and optionally
#'   `label`, `color`, `group`, `arrow`, `offset`, `width`.
#' @return The updated `plasmid` object.
#' @examples
#' feats <- data.frame(
#'   start = c(86, 2535, 3293),
#'   end   = c(1276, 3122, 4153),
#'   label = c("TcR", "ori", "AmpR"),
#'   arrow = c("end", "none", "start")
#' )
#' plot(pp_features(plasmid("pBR322", 4361), feats))
#' @export
pp_features <- function(p, df) {
  stopifnot(inherits(p, "plasmid"))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  missing_cols <- setdiff(c("start", "end"), names(df))
  if (length(missing_cols)) {
    stop("`df` needs column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  opt <- c("label", "color", "group", "arrow", "offset", "width")
  for (i in seq_len(nrow(df))) {
    args <- list(p = p, start = df$start[i], end = df$end[i])
    for (nm in intersect(opt, names(df))) {
      v <- df[[nm]][i]
      if (!is.na(v)) args[[nm]] <- if (is.factor(v)) as.character(v) else v
    }
    p <- do.call(pp_marker, args)
  }
  p
}

#' Common restriction enzymes
#'
#' The recognition sequences [pp_find_sites()] searches for.
#'
#' @return A data frame with columns `enzyme`, `site` (recognition sequence,
#'   IUPAC codes allowed) and `cut` (0-based cut offset on the top strand
#'   within the recognition sequence).
#' @examples
#' head(pp_enzymes())
#' subset(pp_enzymes(), enzyme %in% c("EcoRI", "BamHI"))
#' @export
pp_enzymes <- function() {
  data.frame(
    enzyme = c("AatII", "AflII", "AgeI", "ApaI", "AscI", "AvrII", "BamHI",
               "BglII", "BsiWI", "BspEI", "BsrGI", "BssHII", "ClaI", "DraI",
               "EagI", "EcoRI", "EcoRV", "FseI", "HindIII", "HpaI", "KpnI",
               "MfeI", "MluI", "NarI", "NcoI", "NdeI", "NheI", "NotI",
               "NsiI", "PacI", "PmeI", "PstI", "PvuI", "PvuII", "SacI",
               "SacII", "SalI", "SbfI", "ScaI", "SmaI", "SpeI", "SphI",
               "StuI", "XbaI", "XhoI", "XmaI"),
    site = c("GACGTC", "CTTAAG", "ACCGGT", "GGGCCC", "GGCGCGCC", "CCTAGG",
             "GGATCC", "AGATCT", "CGTACG", "TCCGGA", "TGTACA", "GCGCGC",
             "ATCGAT", "TTTAAA", "CGGCCG", "GAATTC", "GATATC", "GGCCGGCC",
             "AAGCTT", "GTTAAC", "GGTACC", "CAATTG", "ACGCGT", "GGCGCC",
             "CCATGG", "CATATG", "GCTAGC", "GCGGCCGC", "ATGCAT", "TTAATTAA",
             "GTTTAAAC", "CTGCAG", "CGATCG", "CAGCTG", "GAGCTC", "CCGCGG",
             "GTCGAC", "CCTGCAGG", "AGTACT", "CCCGGG", "ACTAGT", "GCATGC",
             "AGGCCT", "TCTAGA", "CTCGAG", "CCCGGG"),
    cut = c(5L, 1L, 1L, 5L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 3L, 1L, 1L,
            3L, 6L, 1L, 3L, 5L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 5L, 5L, 4L, 5L,
            4L, 3L, 5L, 4L, 1L, 6L, 3L, 3L, 1L, 5L, 3L, 1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}

# IUPAC ambiguity codes to a regex character class.
iupac_regex <- function(site) {
  map <- c(A = "A", C = "C", G = "G", T = "T", U = "T",
           R = "[AG]", Y = "[CT]", S = "[GC]", W = "[AT]", K = "[GT]",
           M = "[AC]", B = "[CGT]", D = "[AGT]", H = "[ACT]", V = "[ACG]",
           N = ".")
  chars <- strsplit(toupper(site), "")[[1]]
  unknown <- setdiff(chars, names(map))
  if (length(unknown)) {
    stop("Not a nucleotide sequence: '", site, "' contains ",
         paste(unique(unknown), collapse = ", "), call. = FALSE)
  }
  paste(map[chars], collapse = "")
}

#' Find restriction sites in the plasmid sequence
#'
#' Searches the sequence kept by the readers and adds a [pp_site()] for each
#' cut position. On a circular plasmid the search wraps the origin, so a site
#' straddling position 1 is not missed.
#'
#' @param p A `plasmid` object carrying a sequence (the readers keep one;
#'   [read_fasta()] gives you a bare backbone with just the sequence).
#' @param enzymes Enzyme names from [pp_enzymes()], or a named character
#'   vector of recognition sequences (IUPAC codes allowed), e.g.
#'   `c(MyEnz = "GGWCC")`. `NULL` searches every enzyme in [pp_enzymes()].
#' @param unique_only Keep only enzymes that cut exactly once (default
#'   `TRUE`). This is the classic plasmid-map view: an enzyme cutting a dozen
#'   times adds no information and buries the map in labels.
#' @param max_sites Drop enzymes cutting more than this many times. Ignored
#'   when `unique_only` is `TRUE`.
#' @param show_position Append the cut position to each label, e.g.
#'   `"EcoRI (396)"`.
#' @return The updated `plasmid` object.
#' @seealso [pp_site()], [pp_enzymes()]
#' @examples
#' p <- read_genbank(system.file("extdata", "pDemo.gb", package = "plasmidplot"))
#' p <- pp_find_sites(p, c("EcoRI", "BamHI", "HindIII", "NotI"))
#' p
#' plot(p)
#'
#' # Every unique cutter the built-in table knows about:
#' plot(pp_find_sites(read_genbank(
#'   system.file("extdata", "pDemo.gb", package = "plasmidplot"))))
#' @export
pp_find_sites <- function(p, enzymes = NULL, unique_only = TRUE,
                          max_sites = 3L, show_position = TRUE) {
  stopifnot(inherits(p, "plasmid"))
  if (!has_sequence(p)) {
    stop("This plasmid carries no sequence to search. Read it with ",
         "read_genbank(), read_snapgene(), read_embl() or read_fasta(), ",
         "and do not set sequence = FALSE.", call. = FALSE)
  }
  tbl <- pp_enzymes()
  if (is.null(enzymes)) {
    patterns <- stats::setNames(tbl$site, tbl$enzyme)
    cuts <- stats::setNames(tbl$cut, tbl$enzyme)
  } else if (!is.null(names(enzymes))) {
    patterns <- enzymes
    cuts <- stats::setNames(rep(0L, length(enzymes)), names(enzymes))
  } else {
    unknown <- setdiff(enzymes, tbl$enzyme)
    if (length(unknown)) {
      stop("Unknown enzyme(s): ", paste(unknown, collapse = ", "),
           ". See pp_enzymes() for the built-in list, or pass a named ",
           "vector of recognition sequences.", call. = FALSE)
    }
    i <- match(enzymes, tbl$enzyme)
    patterns <- stats::setNames(tbl$site[i], tbl$enzyme[i])
    cuts <- stats::setNames(tbl$cut[i], tbl$enzyme[i])
  }

  seq <- toupper(p$sequence)
  circular <- !identical(p$topology, "linear")
  for (nm in names(patterns)) {
    rx <- iupac_regex(patterns[[nm]])
    w <- nchar(patterns[[nm]])
    # Extend by w-1 bases so a site spanning the origin is still matched.
    hay <- if (circular && nchar(seq) > w) {
      paste0(seq, substr(seq, 1L, w - 1L))
    } else seq
    hits <- gregexpr(rx, hay)[[1]]
    if (hits[1] == -1L) next
    pos <- ((hits - 1L + cuts[[nm]]) %% nchar(seq)) + 1L
    pos <- sort(unique(pos))
    if (unique_only && length(pos) != 1L) next
    if (!unique_only && length(pos) > max_sites) next
    for (bp in pos) {
      lab <- if (show_position) sprintf("%s (%d)", nm, as.integer(bp)) else nm
      p <- pp_site(p, bp, lab)
    }
  }
  p
}
