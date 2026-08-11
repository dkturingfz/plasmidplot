#' Create a plasmid object
#'
#' A `plasmid` object holds the plasmid name, its length in base pairs, and
#' the list of features (markers) added with [pp_marker()]. Plot it with
#' [plot.plasmid()].
#'
#' @param name Plasmid name, shown in the center of the map.
#' @param length Plasmid length in base pairs (a single positive number).
#' @return An object of class `plasmid`.
#' @examples
#' p <- plasmid("pBR322", 4361)
#' p <- pp_marker(p, 86, 1276, label = "TcR")
#' plot(p)
#' @export
plasmid <- function(name, length) {
  stopifnot(is.character(name), length(name) == 1L)
  if (!is.numeric(length) || length(length) != 1L || is.na(length) || length <= 0) {
    stop("`length` must be a single positive number of base pairs.", call. = FALSE)
  }
  structure(
    list(name = name, length = as.numeric(length), markers = list(),
         sites = list(), topology = "circular", sequence = NA_character_),
    class = "plasmid"
  )
}

#' Add a feature marker to a plasmid
#'
#' Markers are drawn as colored arcs on the plasmid backbone. A marker whose
#' `end` is smaller than its `start` is taken to wrap across the origin
#' (position 1). Calls are pipe-friendly: each returns the updated object.
#'
#' @param p A `plasmid` object created by [plasmid()].
#' @param start,end Feature boundaries in base pairs (1-based, inclusive).
#' @param label Optional text label, drawn outside the map with a leader line.
#' @param color Arc color. Defaults to the next color of the style palette.
#' @param group Optional grouping key. Markers sharing a group take the same
#'   palette color, so e.g. every CDS can be drawn in one color without
#'   naming the color here. Ignored when `color` is given.
#' @param arrow Arrowhead placement: `"none"`, `"end"` (points clockwise, for
#'   a forward-strand feature), `"start"` (counter-clockwise, reverse strand),
#'   or `"both"` (bidirectional).
#' @param offset Radial offset from the backbone in npc units. Positive moves
#'   the marker outward; useful to separate overlapping features.
#' @param width Arc thickness in npc units. Defaults to the style's
#'   `arc`.
#' @return The updated `plasmid` object.
#' @examples
#' p <- plasmid("pUC19", 2686) |>
#'   pp_marker(146, 469, label = "lacZ-alpha", arrow = "end") |>
#'   pp_marker(1629, 2486, label = "AmpR", arrow = "start")
#' plot(p)
#' @export
pp_marker <- function(p, start, end, label = NULL, color = NULL, group = NULL,
                      arrow = c("none", "end", "start", "both"),
                      offset = 0, width = NULL) {
  stopifnot(inherits(p, "plasmid"))
  arrow <- match.arg(arrow)
  for (v in list(start = start, end = end)) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop("`start` and `end` must be single numbers.", call. = FALSE)
    }
  }
  if (start < 1 || start > p$length || end < 1 || end > p$length) {
    stop(sprintf("Positions must lie in [1, %d].", round(p$length)), call. = FALSE)
  }
  p$markers[[length(p$markers) + 1L]] <- list(
    start = start, end = end,
    label = if (is.null(label)) NA_character_ else as.character(label),
    color = if (is.null(color)) NA_character_ else as.character(color),
    group = if (is.null(group)) NA_character_ else as.character(group),
    arrow = arrow, offset = offset, width = width
  )
  p
}

# A plasmid may carry no sequence at all, either because the file had none or
# because the reader was told to drop it.
has_sequence <- function(x) {
  s <- x$sequence
  length(s) == 1L && !is.na(s) && nzchar(s)
}

#' @export
print.plasmid <- function(x, ...) {
  cat(sprintf("<plasmid> %s (%s bp, %s)%s\n",
              x$name, format(round(x$length), big.mark = ","),
              if (is.null(x$topology)) "circular" else x$topology,
              if (has_sequence(x)) " with sequence" else ""))
  cat(sprintf("  %d marker(s), %d site(s)\n",
              length(x$markers), length(x$sites)))
  for (m in x$markers) {
    cat(sprintf("  %7d..%-7d %s%s\n", round(m$start), round(m$end),
                if (is.na(m$label)) "(unlabeled)" else m$label,
                switch(m$arrow, none = "", end = " ->", start = " <-",
                       both = " <->")))
  }
  for (st in x$sites) {
    cat(sprintf("  %7d          %s\n", round(st$position),
                if (is.na(st$label)) "(unlabeled)" else st$label))
  }
  invisible(x)
}
