# Internal geometry helpers. All coordinates live in a square "snpc" viewport
# with center (0.5, 0.5); radii are fractions of the viewport (max usable 0.5).

# Angle (radians) for a base-pair position: origin at 12 o'clock, clockwise.
pp_theta <- function(pos, len) pi / 2 - 2 * pi * pos / len

pp_x <- function(theta, r) 0.5 + r * cos(theta)
pp_y <- function(theta, r) 0.5 + r * sin(theta)

# Base-pair sample points along an arc, dense enough to look smooth.
arc_seq <- function(bp1, bp2, len, step_deg = 1) {
  n <- max(2L, ceiling(abs(bp2 - bp1) / len * 360 / step_deg))
  seq(bp1, bp2, length.out = n + 1L)
}

# Annular sector polygon between bp1..bp2 spanning radii r0..r1.
sector_polygon <- function(bp1, bp2, len, r0, r1) {
  theta <- pp_theta(arc_seq(bp1, bp2, len), len)
  list(
    x = c(pp_x(theta, r1), rev(pp_x(theta, r0))),
    y = c(pp_y(theta, r1), rev(pp_y(theta, r0)))
  )
}

# One closed outline for a whole feature, arrowheads included. Body and head
# must be a single polygon: drawing them as two overlapping shapes leaves an
# antialiasing seam along the join where the background bleeds through.
#
# The path runs: outer arc forward, out to the head flare, to the tip, back to
# the inner flare, inner arc back, then the start-side head if there is one.
marker_polygon <- function(bp1, bp2, len, r0, r1, arrow = "none",
                           head_bp = 0, flare = 0.35) {
  ext <- flare * (r1 - r0)
  rm <- (r0 + r1) / 2
  head_end <- arrow %in% c("end", "both")
  head_start <- arrow %in% c("start", "both")
  b1 <- if (head_start) bp1 + head_bp else bp1
  b2 <- if (head_end) bp2 - head_bp else bp2
  th <- pp_theta(arc_seq(b1, b2, len), len)

  x <- pp_x(th, r1)
  y <- pp_y(th, r1)
  if (head_end) {
    tb <- pp_theta(b2, len)
    tt <- pp_theta(bp2, len)
    x <- c(x, pp_x(tb, r1 + ext), pp_x(tt, rm), pp_x(tb, r0 - ext))
    y <- c(y, pp_y(tb, r1 + ext), pp_y(tt, rm), pp_y(tb, r0 - ext))
  }
  x <- c(x, rev(pp_x(th, r0)))
  y <- c(y, rev(pp_y(th, r0)))
  if (head_start) {
    tb <- pp_theta(b1, len)
    tt <- pp_theta(bp1, len)
    x <- c(x, pp_x(tb, r0 - ext), pp_x(tt, rm), pp_x(tb, r1 + ext))
    y <- c(y, pp_y(tb, r0 - ext), pp_y(tt, rm), pp_y(tb, r1 + ext))
  }
  list(x = x, y = y)
}

# Base-pair length of one arrowhead, capped so the head never eats the body.
head_length <- function(span, len, arrow, arrow_deg) {
  n <- switch(arrow, none = 0L, both = 2L, 1L)
  if (n == 0L) return(0)
  min(span / (n * 2), len * arrow_deg / 360)
}

# Half-thickness of the backbone itself: a ring has width, a line does not.
backbone_half <- function(s) {
  if (identical(s$backbone, "ring")) s$track / 2 else 0
}

# Signed displacement of a feature band from the backbone. `anchor` decides
# whether features sit on the backbone, ride clear of it on the outward side,
# or run on the inward side -- the same rule in both layouts, which is why a
# style keeps its character when switched from circular to linear.
anchor_offset <- function(s, w = NULL) {
  if (is.null(w) || is.na(w)) w <- s$arc
  switch(s$anchor,
    center = 0,
    outside = backbone_half(s) + s$gap + w / 2,
    inside = -(backbone_half(s) + s$gap + w / 2)
  )
}

# Radial center of a feature band on a circular map.
marker_radius <- function(s, offset = 0, w = NULL) {
  s$radius + anchor_offset(s, w) + offset
}

# Nearest and furthest extent anything on the backbone occupies, measured from
# the backbone's own center line. Ticks stay clear on the near side; sites and
# labels stay clear on the far side. Works for both layouts.
band_span <- function(s) {
  off <- anchor_offset(s)
  c(near = min(-backbone_half(s), off - s$arc / 2),
    far = max(backbone_half(s), off + s$arc / 2))
}

# Same thing expressed as circular radii.
band_extent <- function(s) {
  sp <- band_span(s)
  c(inner = s$radius + sp[["near"]], outer = s$radius + sp[["far"]])
}

# --- linear layout ----------------------------------------------------------
# A linear map runs left to right across the viewport, with labels above it
# and the scale below. `radius` means the same thing in both layouts -- how
# much room the map takes -- so it sets the half-width here, keeping all seven
# shape parameters live whichever way the map is drawn.
lin_half <- function(s) min(s$radius * 1.40, 0.44)
lin_x0 <- function(s) 0.5 - lin_half(s)
lin_x1 <- function(s) 0.5 + lin_half(s)
lin_baseline <- function() 0.42

lin_x <- function(bp, len, s) {
  lin_x0(s) + (bp - 1) / max(len - 1, 1) * (lin_x1(s) - lin_x0(s))
}

# A feature that wraps the origin cannot exist on a linear molecule, so it is
# drawn as the two pieces it would occupy, rather than silently spanning the
# whole map.
linear_pieces <- function(m, len) {
  if (m$end >= m$start) return(list(c(m$start, m$end)))
  list(c(m$start, len), c(1, m$end))
}

# A tick interval giving roughly `target` major ticks around the circle.
nice_interval <- function(len, target = 8) {
  raw <- len / target
  base <- 10^floor(log10(raw))
  cand <- c(1, 2, 2.5, 5, 10) * base
  cand[which.min(abs(len / cand - target))]
}

# Where the scale is marked. Shared by both layouts so they cannot disagree
# about the numbers, and guarded: an interval wider than the molecule leaves
# no major ticks at all rather than asking seq() for a backwards sequence.
tick_positions <- function(len, interval, layout = "circular") {
  majors <- if (is.finite(interval) && interval > 0 &&
                len - interval / 2 >= interval) {
    seq(interval, len - interval / 2, by = interval)
  } else {
    numeric(0)
  }
  # A circle has no ends, so it marks the origin; a linear molecule has two
  # real endpoints and marks both.
  if (identical(layout, "linear")) {
    unique(c(1, majors, len))
  } else {
    unique(c(0, majors))
  }
}

# Normalize a marker span so end >= start (features may wrap the origin).
span_bp <- function(m, len) {
  end <- if (m$end < m$start) m$end + len else m$end
  c(m$start, end)
}

# Greedy vertical spread of labels on each side of the map so they don't
# overlap. Returns anchor x/y and hjust per label.
layout_labels <- function(mid_bp, len, label_r, min_gap) {
  theta <- pp_theta(mid_bp, len)
  x <- pp_x(theta, label_r)
  y <- pp_y(theta, label_r)
  side <- ifelse(cos(theta) >= 0, 1, -1)
  for (sd in c(1, -1)) {
    idx <- which(side == sd)
    if (length(idx) < 2L) next
    ord <- idx[order(y[idx], decreasing = TRUE)]
    for (k in 2:length(ord)) {
      prev <- ord[k - 1L]
      cur <- ord[k]
      if (y[prev] - y[cur] < min_gap) y[cur] <- y[prev] - min_gap
    }
    # Re-anchor pushed labels on the label circle where possible.
    for (i in idx) {
      dy <- y[i] - 0.5
      x[i] <- 0.5 + sd * sqrt(max(label_r^2 - dy^2, 0.03^2))
    }
  }
  list(x = x, y = y, hjust = ifelse(side == 1, 0, 1), side = side)
}
