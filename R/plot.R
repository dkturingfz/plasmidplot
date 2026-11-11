#' Plot a plasmid map
#'
#' Renders the plasmid as a circular or linear map: a backbone with a
#' base-pair scale, colored feature arcs (with optional arrowheads), callout
#' labels with leader lines, restriction sites, and the plasmid name and size.
#'
#' The layout follows the plasmid's `topology` unless the style says
#' otherwise; see the `layout` shape parameter in [pp_style()].
#'
#' Labels are placed so that they do not overlap. On a circular map they are
#' spread up and down each side of the ring. On a linear map they sit level
#' in up to three rows while that stays tidy, and are otherwise slanted and
#' fanned out sideways, which is what keeps a dozen restriction sites a few
#' bases apart readable. A feature too short to see at the map's scale, such
#' as a single-base variant, is drawn at a small minimum width.
#'
#' @param x A `plasmid` object.
#' @param style A preset name or a [pp_style()] object. Call `pp_style()` for
#'   the available presets: `"angular"`, `"classic"`, `"dark"`, `"snapgene"`,
#'   `"soft"`, `"neon"`, `"minimal"`, `"blueprint"`.
#' @param bg Background to paint behind the map. `NULL` (the default) paints
#'   nothing, so the device or surrounding viewport shows through -- styles
#'   never impose a background of their own. Pass a color, or
#'   `pp_canvas(style)` for the ground the style was designed against, which
#'   the dark presets need to be legible.
#' @param newpage Start a new grid page first (default `TRUE`). Set to
#'   `FALSE` to draw into an existing viewport, e.g. for panel layouts.
#' @param ... Ignored.
#' @return Invisibly, `x`.
#' @examples
#' p <- plasmid("pBR322", 4361) |>
#'   pp_marker(86, 1276, label = "TcR", arrow = "end") |>
#'   pp_marker(2535, 3122, label = "ori") |>
#'   pp_marker(3293, 4153, label = "AmpR", arrow = "start")
#' plot(p, style = "angular")
#' plot(p, style = "dark", bg = pp_canvas("dark"))
#' plot(p, style = pp_style("angular", layout = "linear"))
#' @export
plot.plasmid <- function(x, style = "angular", bg = NULL, newpage = TRUE, ...) {
  s <- resolve_style(style)
  if (!is.null(bg) && (length(bg) != 1L || (!is.na(bg) && !is_color(bg)))) {
    stop("`bg` must be a single color, or NULL to paint nothing.",
         call. = FALSE)
  }
  if (newpage) grid::grid.newpage()
  ground <- if (is.null(bg)) s$bg else bg
  if (!is.null(ground) && !is.na(ground)) {
    grid::grid.rect(gp = grid::gpar(fill = ground, col = NA))
  }
  cols <- marker_colors(x, s)
  if (identical(resolve_layout(s, x), "linear")) {
    draw_linear(x, s, cols)
  } else {
    # A circle needs a square viewport; a linear map wants the full width.
    grid::pushViewport(grid::viewport(
      width = grid::unit(1, "snpc"), height = grid::unit(1, "snpc")
    ))
    on.exit(grid::popViewport(), add = TRUE)
    draw_circular(x, s, cols)
  }
  invisible(x)
}

# Palette colors are handed out in order of first appearance. Markers sharing
# a `group` share the color assigned to that group's first member.
marker_colors <- function(x, s) {
  n <- length(x$markers)
  cols <- character(n)
  slot <- integer(n)          # 0 marks a marker that named its own color
  group_slot <- integer(0)
  used <- 0L
  # First pass assigns slot numbers, so the palette is asked for exactly the
  # number of colors the map needs before any color is looked up.
  for (i in seq_len(n)) {
    m <- x$markers[[i]]
    if (!is.na(m$color)) next
    g <- if (is.null(m$group)) NA_character_ else m$group
    if (!is.na(g) && g %in% names(group_slot)) {
      slot[i] <- group_slot[[g]]
      next
    }
    used <- used + 1L
    slot[i] <- used
    if (!is.na(g)) group_slot[[g]] <- used
  }
  pal <- palette_colors(s$palette, used)
  for (i in seq_len(n)) {
    cols[i] <- if (slot[i] == 0L) {
      x$markers[[i]]$color
    } else {
      pal[(slot[i] - 1L) %% length(pal) + 1L]
    }
  }
  cols
}

marker_gp <- function(s, fill) {
  grid::gpar(fill = fill, col = resolve_border(s$arc_border, fill),
             lwd = s$arc_lwd, linejoin = "mitre")
}

site_color <- function(s, st) {
  if (!is.na(st$color)) return(st$color)
  if (is.na(s$site_col)) s$label_col else s$site_col
}

marker_thickness <- function(s, m) {
  if (is.null(m$width) || is.na(m$width)) s$arc else m$width
}

# ---------------------------------------------------------------------------
# Circular layout
# ---------------------------------------------------------------------------

draw_circular <- function(x, s, cols) {
  circ_backbone(x, s)
  circ_ticks(x, s)
  circ_markers(x, s, cols)
  circ_sites(x, s)
  draw_labels(circ_label_entries(x, s), s, len = x$length)
  circ_title(x, s)
}

circ_backbone <- function(x, s) {
  if (identical(s$backbone, "none")) return(invisible())
  if (identical(s$backbone, "ring")) {
    poly <- sector_polygon(0, x$length, x$length,
                           s$radius - s$track / 2, s$radius + s$track / 2)
    grid::grid.polygon(poly$x, poly$y, default.units = "npc",
                       gp = grid::gpar(fill = s$track_fill, col = s$track_col))
  } else {
    grid::grid.circle(0.5, 0.5, r = s$radius, default.units = "npc",
                      gp = grid::gpar(col = s$backbone_col,
                                      lwd = s$backbone_lwd, fill = NA))
  }
}

circ_ticks <- function(x, s) {
  interval <- if (is.na(s$tick_every)) nice_interval(x$length) else s$tick_every
  r_in <- band_extent(s)[["inner"]]
  for (bp in tick_positions(x$length, interval, "circular")) {
    th <- pp_theta(bp, x$length)
    grid::grid.lines(
      pp_x(th, c(r_in, r_in - s$tick_len)), pp_y(th, c(r_in, r_in - s$tick_len)),
      default.units = "npc", gp = grid::gpar(col = s$tick_col, lwd = 0.8))
    if (s$tick_labels && bp > 0) {
      lab_r <- r_in - s$tick_len - 0.010
      grid::grid.text(
        format(bp, big.mark = ",", scientific = FALSE),
        pp_x(th, lab_r), pp_y(th, lab_r), default.units = "npc",
        rot = tick_rot(th),
        gp = grid::gpar(col = s$tick_col, fontsize = s$tick_fontsize))
    }
  }
}

# Rotate tick labels along the circle, flipped on the lower half so text
# never reads upside down.
tick_rot <- function(theta) {
  deg <- (theta * 180 / pi - 90) %% 360
  if (deg > 90 && deg < 270) deg <- deg + 180
  deg %% 360
}

circ_markers <- function(x, s, cols) {
  for (i in seq_along(x$markers)) {
    m <- x$markers[[i]]
    w <- marker_thickness(s, m)
    R <- marker_radius(s, m$offset, w)
    sp <- drawn_span(m, x$length, R)
    hbp <- head_length(sp[2] - sp[1], x$length, m$arrow, s$arrow_deg)
    poly <- marker_polygon(sp[1], sp[2], x$length, R - w / 2, R + w / 2,
                           m$arrow, hbp, s$arrow_flare)
    grid::grid.polygon(poly$x, poly$y, default.units = "npc",
                       gp = marker_gp(s, cols[i]))
  }
}

circ_sites <- function(x, s) {
  if (!length(x$sites)) return(invisible())
  r0 <- band_extent(s)[["outer"]]
  for (st in x$sites) {
    th <- pp_theta(st$position, x$length)
    grid::grid.lines(
      pp_x(th, c(r0, r0 + s$site_len)), pp_y(th, c(r0, r0 + s$site_len)),
      default.units = "npc",
      gp = grid::gpar(col = site_color(s, st), lwd = s$site_lwd))
  }
}

circ_label_entries <- function(x, s) {
  out <- list()
  for (m in x$markers) {
    if (is.na(m$label)) next
    sp <- span_bp(m, x$length)
    w <- marker_thickness(s, m)
    out[[length(out) + 1L]] <- list(
      bp = mean(sp) %% x$length, text = m$label,
      anchor = marker_radius(s, m$offset, w) + w / 2 + 0.006,
      col = s$label_col, fontsize = s$label_fontsize,
      fontface = s$label_fontface)
  }
  for (st in x$sites) {
    if (is.na(st$label)) next
    out[[length(out) + 1L]] <- list(
      bp = st$position, text = st$label,
      anchor = band_extent(s)[["outer"]] + s$site_len + 0.004,
      col = site_color(s, st), fontsize = s$site_fontsize, fontface = "plain")
  }
  out
}

circ_title <- function(x, s) {
  if (!s$show_title) return(invisible())
  grid::grid.text(x$name, 0.5, 0.52, default.units = "npc",
                  gp = grid::gpar(col = s$title_col, fontsize = s$title_fontsize,
                                  fontface = "bold"))
  grid::grid.text(paste(format(round(x$length), big.mark = ","), "bp"),
                  0.5, 0.475, default.units = "npc",
                  gp = grid::gpar(col = s$subtitle_col,
                                  fontsize = s$subtitle_fontsize))
}

# ---------------------------------------------------------------------------
# Linear layout
# ---------------------------------------------------------------------------

draw_linear <- function(x, s, cols) {
  lin_backbone(x, s)
  lin_ticks(x, s)
  lin_markers(x, s, cols)
  lin_sites(x, s)
  lin_labels(x, s)
  lin_title(x, s)
}

lin_backbone <- function(x, s) {
  if (identical(s$backbone, "none")) return(invisible())
  y <- lin_baseline()
  if (identical(s$backbone, "ring")) {
    grid::grid.rect(x = (lin_x0(s) + lin_x1(s)) / 2, y = y,
                    width = lin_x1(s) - lin_x0(s), height = s$track,
                    default.units = "npc",
                    gp = grid::gpar(fill = s$track_fill, col = s$track_col))
  } else {
    grid::grid.lines(c(lin_x0(s), lin_x1(s)), c(y, y), default.units = "npc",
                     gp = grid::gpar(col = s$backbone_col,
                                     lwd = s$backbone_lwd))
  }
}

lin_ticks <- function(x, s) {
  interval <- if (is.na(s$tick_every)) nice_interval(x$length) else s$tick_every
  y0 <- lin_baseline() + band_span(s)[["near"]]
  for (bp in tick_positions(x$length, interval, "linear")) {
    xx <- lin_x(bp, x$length, s)
    grid::grid.lines(c(xx, xx), c(y0, y0 - s$tick_len), default.units = "npc",
                     gp = grid::gpar(col = s$tick_col, lwd = 0.8))
    if (s$tick_labels) {
      grid::grid.text(format(round(bp), big.mark = ",", scientific = FALSE),
                      xx, y0 - s$tick_len - 0.012, default.units = "npc",
                      gp = grid::gpar(col = s$tick_col,
                                      fontsize = s$tick_fontsize))
    }
  }
}

lin_markers <- function(x, s, cols) {
  for (i in seq_along(x$markers)) {
    m <- x$markers[[i]]
    w <- marker_thickness(s, m)
    y <- lin_baseline() + anchor_offset(s, w) + m$offset
    pieces <- linear_pieces(m, x$length)
    for (k in seq_along(pieces)) {
      bp <- pieces[[k]]
      # An origin-crossing feature keeps its head only on the piece that
      # actually carries the end it points at.
      arrow <- m$arrow
      if (length(pieces) > 1L) {
        arrow <- if (k == 1L) {
          if (m$arrow %in% c("start", "both")) "start" else "none"
        } else {
          if (m$arrow %in% c("end", "both")) "end" else "none"
        }
      }
      poly <- lin_marker_polygon(bp[1], bp[2], x$length, y, w, arrow, s)
      grid::grid.polygon(poly$x, poly$y, default.units = "npc",
                         gp = marker_gp(s, cols[i]))
    }
  }
}

# The linear twin of marker_polygon(): one closed outline so the body and the
# arrowhead never show a seam between them.
lin_marker_polygon <- function(bp1, bp2, len, y, w, arrow, s) {
  xx <- at_least(lin_x(bp1, len, s), lin_x(bp2, len, s), pp_min_extent)
  x1 <- xx[1]
  x2 <- xx[2]
  n_heads <- switch(arrow, none = 0L, both = 2L, 1L)
  head <- if (n_heads == 0L) 0 else {
    min((x2 - x1) / (n_heads * 2),
        (lin_x1(s) - lin_x0(s)) * s$arrow_deg / 360)
  }
  ext <- s$arrow_flare * w
  head_end <- arrow %in% c("end", "both")
  head_start <- arrow %in% c("start", "both")
  a <- if (head_start) x1 + head else x1
  b <- if (head_end) x2 - head else x2

  # Top edge forward, end head, bottom edge back, start head. Every corner of
  # the body has to be named or the top edge slopes into the head and the
  # feature renders as a wedge.
  xs <- c(a, b)
  ys <- c(y + w / 2, y + w / 2)
  if (head_end) {
    xs <- c(xs, b, x2, b)
    ys <- c(ys, y + w / 2 + ext, y, y - w / 2 - ext)
  }
  xs <- c(xs, b, a)
  ys <- c(ys, y - w / 2, y - w / 2)
  if (head_start) {
    xs <- c(xs, a, x1, a)
    ys <- c(ys, y - w / 2 - ext, y, y + w / 2 + ext)
  }
  list(x = xs, y = ys)
}

lin_sites <- function(x, s) {
  if (!length(x$sites)) return(invisible())
  y0 <- lin_baseline() + band_span(s)[["far"]]
  for (st in x$sites) {
    xx <- lin_x(st$position, x$length, s)
    grid::grid.lines(c(xx, xx), c(y0, y0 + s$site_len), default.units = "npc",
                     gp = grid::gpar(col = site_color(s, st), lwd = s$site_lwd))
  }
}

lin_label_entries <- function(x, s) {
  out <- list()
  for (m in x$markers) {
    if (is.na(m$label)) next
    w <- marker_thickness(s, m)
    pieces <- linear_pieces(m, x$length)
    # Label the longer piece, so a wrapped feature is named once.
    widest <- pieces[[which.max(vapply(pieces, function(b) b[2] - b[1],
                                       numeric(1)))]]
    out[[length(out) + 1L]] <- list(
      x = lin_x(mean(widest), x$length, s), text = m$label,
      anchor = lin_baseline() + anchor_offset(s, w) + m$offset + w / 2 + 0.004,
      col = s$label_col, fontsize = s$label_fontsize,
      fontface = s$label_fontface)
  }
  for (st in x$sites) {
    if (is.na(st$label)) next
    out[[length(out) + 1L]] <- list(
      x = lin_x(st$position, x$length, s), text = st$label,
      anchor = lin_baseline() + band_span(s)[["far"]] + s$site_len + 0.004,
      col = site_color(s, st), fontsize = s$site_fontsize, fontface = "plain")
  }
  out
}

# Width as drawn: a bold label is wider than the same text set plain, and a
# label measured plain would be let run off the edge of the page.
text_width_npc <- function(txt, fontsize, fontface = "plain") {
  grid::convertWidth(
    grid::grobWidth(grid::textGrob(txt, gp = grid::gpar(fontsize = fontsize,
                                                        fontface = fontface))),
    "npc", valueOnly = TRUE)
}

# Level text in stacked rows reads best, so it is kept for as long as it
# stays tidy: no more rows than this, clear of the title, and no leader
# running through a label on its way to the row above.
lin_max_rows <- 3L

# Slanted labels start on one line and rise at this angle.
lin_slant_deg <- 45

# Rows for level labels: each goes in the lowest row with room for it. On a
# linear map labels compete for horizontal room, so a single row would
# overlap as soon as two features sit close together.
lin_stack_rows <- function(lx, widths, pad = 0.012) {
  row_end <- numeric(0)      # right edge of the last label placed in each row
  row_of <- integer(length(lx))
  for (i in order(lx)) {
    r <- which(row_end < lx[i] - widths[i] / 2 - pad)
    r <- if (length(r)) r[1] else length(row_end) + 1L
    row_of[i] <- r
    row_end[r] <- lx[i] + widths[i] / 2
  }
  row_of
}

# How high labels may reach: clear of the title when there is one.
lin_label_ceiling <- function(s) {
  if (!s$show_title) return(0.99)
  0.94 - grid::convertHeight(grid::unit(s$title_fontsize * 0.9, "points"),
                             "npc", valueOnly = TRUE)
}

# Where every label goes, and where its leader ends. Returns one row per
# entry: text position, rotation and justification, and the leader's far end.
lin_label_layout <- function(entries, s) {
  xs <- vapply(entries, function(e) e$x, numeric(1))
  sizes <- vapply(entries, function(e) e$fontsize, numeric(1))
  widths <- vapply(entries, function(e) {
    text_width_npc(e$text, e$fontsize, e$fontface)
  }, numeric(1))
  row_h <- grid::convertHeight(grid::unit(max(sizes) * 1.7, "points"), "npc",
                               valueOnly = TRUE)
  anchors <- vapply(entries, function(e) e$anchor, numeric(1))
  y_first <- max(anchors) + 0.035
  ceiling <- lin_label_ceiling(s)
  margin <- 0.006

  # A label on the first or last feature would otherwise hang off the page,
  # so pull it inside and let its leader slant instead.
  lx <- pmin(pmax(xs, widths / 2 + margin), 1 - widths / 2 - margin)
  row_of <- lin_stack_rows(lx, widths)
  rows <- max(row_of)
  y <- y_first + (row_of - 1L) * row_h
  leader_y <- y - row_h * 0.32
  if (rows <= lin_max_rows && y_first + (rows - 0.5) * row_h <= ceiling &&
      !lin_leader_hits_label(xs, anchors, lx, leader_y, y, widths, row_of)) {
    return(data.frame(x = lx, y = y, rot = 0, hjust = 0.5,
                      leader_x = lx, leader_y = leader_y))
  }

  # Too crowded for rows. Slanted labels need only a line height of
  # horizontal room each, so a dozen sites a few bases apart fan out side by
  # side. All of them start on one line that every leader stops short of, so
  # no leader can cross a label however tight the cluster.
  #
  # The slant eases off if the longest label would otherwise reach the
  # title, but not so far that the labels lie back down on one another.
  inch <- function(w) grid::convertWidth(grid::unit(w, "npc"), "inches",
                                         valueOnly = TRUE)
  room <- grid::convertHeight(grid::unit(max(ceiling - y_first, 0), "npc"),
                              "inches", valueOnly = TRUE)
  deg <- min(lin_slant_deg,
             max(20, asin(min(1, room / max(inch(max(widths)), 1e-9))) *
                   180 / pi))
  rad <- deg * pi / 180
  # Baselines one line apart, measured across the slant.
  step <- grid::convertWidth(grid::unit(max(sizes) * 1.25 / sin(rad), "points"),
                             "npc", valueOnly = TRUE)
  n <- length(xs)
  lo <- margin + step / 2
  hi <- 1 - margin - widths * cos(rad)      # each label must end on the page
  if (n > 1L) step <- min(step, max(min(hi) - lo, 0) / (n - 1L))
  lx <- spread_positions(xs, step, lo, max(hi))
  # A long label near the right-hand edge has less room than its neighbors:
  # move it in, and whatever it then crowds along with it.
  ord <- order(lx)
  for (k in rev(seq_along(ord))) {
    i <- ord[k]
    cap <- if (k < n) lx[ord[k + 1L]] - step else Inf
    lx[i] <- max(lo, min(lx[i], hi[i], cap))
  }
  # Labels moved far sideways would be reached by leaders lying almost flat
  # along the map. Raise the whole line by what its longest reach needs to
  # keep a visible slope, as far as the room under the title allows.
  rise <- inch(max(widths)) * sin(rad)
  lift <- min(inch(max(abs(lx - xs))) * tan(25 * pi / 180),
              max(room - rise, 0))
  y0 <- y_first + grid::convertHeight(grid::unit(lift, "inches"), "npc",
                                      valueOnly = TRUE)
  data.frame(x = lx, y = y0, rot = deg, hjust = 0,
             leader_x = lx, leader_y = y0 - row_h * 0.25)
}

# Does any leader, on its way up to a label in a higher row, run through a
# label in a row beneath it? The stack is only worth keeping when none does.
lin_leader_hits_label <- function(xs, anchors, lx, leader_y, y, widths,
                                  row_of, pad = 0.002) {
  for (i in which(row_of > 1L)) {
    for (j in which(row_of < row_of[i])) {
      at <- xs[i] + (y[j] - anchors[i]) / (leader_y[i] - anchors[i]) *
        (lx[i] - xs[i])
      if (abs(at - lx[j]) < widths[j] / 2 + pad) return(TRUE)
    }
  }
  FALSE
}

lin_labels <- function(x, s) {
  entries <- lin_label_entries(x, s)
  if (!length(entries)) return(invisible())
  lay <- lin_label_layout(entries, s)
  for (i in seq_along(entries)) {
    e <- entries[[i]]
    grid::grid.lines(c(e$x, lay$leader_x[i]), c(e$anchor, lay$leader_y[i]),
                     default.units = "npc",
                     gp = grid::gpar(col = s$leader_col, lwd = s$leader_lwd))
    grid::grid.text(e$text, lay$x[i], lay$y[i], default.units = "npc",
                    rot = lay$rot[i], hjust = lay$hjust[i], vjust = 0.5,
                    gp = grid::gpar(col = e$col, fontsize = e$fontsize,
                                    fontface = e$fontface))
  }
}

lin_title <- function(x, s) {
  if (!s$show_title) return(invisible())
  grid::grid.text(x$name, lin_x0(s), 0.94, default.units = "npc", hjust = 0,
                  gp = grid::gpar(col = s$title_col, fontsize = s$title_fontsize,
                                  fontface = "bold"))
  grid::grid.text(paste(format(round(x$length), big.mark = ","), "bp",
                        if (identical(x$topology, "linear")) "linear" else
                          "circular, shown linear"),
                  lin_x1(s), 0.94, default.units = "npc", hjust = 1,
                  gp = grid::gpar(col = s$subtitle_col,
                                  fontsize = s$subtitle_fontsize))
}

# ---------------------------------------------------------------------------
# Shared circular label placement
# ---------------------------------------------------------------------------

draw_labels <- function(entries, s, len) {
  if (!length(entries)) return(invisible())
  bp <- vapply(entries, function(e) e$bp, numeric(1))
  sizes <- vapply(entries, function(e) e$fontsize, numeric(1))
  anchors <- vapply(entries, function(e) e$anchor, numeric(1))

  min_gap <- grid::convertHeight(
    grid::unit(max(sizes) * 1.6, "points"), "npc", valueOnly = TRUE)
  # The label ring has to clear every leader's start, or a style that pushes
  # features outward would run its labels through its own arcs.
  r_label <- max(s$label_r, max(anchors) + 0.02)
  lay <- layout_labels(bp, len, r_label, min_gap)

  for (k in seq_along(entries)) {
    e <- entries[[k]]
    th <- pp_theta(e$bp, len)
    grid::grid.lines(
      c(pp_x(th, e$anchor), lay$x[k] - lay$side[k] * 0.008),
      c(pp_y(th, e$anchor), lay$y[k]), default.units = "npc",
      gp = grid::gpar(col = s$leader_col, lwd = s$leader_lwd))
    grid::grid.text(e$text, lay$x[k], lay$y[k], default.units = "npc",
                    hjust = lay$hjust[k], vjust = 0.5,
                    gp = grid::gpar(col = e$col, fontsize = e$fontsize,
                                    fontface = e$fontface))
  }
}
