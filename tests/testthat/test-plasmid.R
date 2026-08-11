test_that("plasmid() validates input", {
  expect_s3_class(plasmid("pX", 5000), "plasmid")
  expect_error(plasmid("pX", -1), "positive")
  expect_error(plasmid("pX", c(1, 2)))
})

test_that("pp_marker() appends markers and validates positions", {
  p <- plasmid("pX", 5000)
  p <- pp_marker(p, 100, 900, label = "a")
  p <- pp_marker(p, 4800, 200, label = "wraps")
  expect_length(p$markers, 2L)
  expect_error(pp_marker(p, 0, 100), "Positions")
  expect_error(pp_marker(p, 100, 6000), "Positions")
})

test_that("span_bp unwraps origin-crossing features", {
  expect_equal(plasmidplot:::span_bp(list(start = 4800, end = 200), 5000),
               c(4800, 5200))
  expect_equal(plasmidplot:::span_bp(list(start = 100, end = 900), 5000),
               c(100, 900))
})

test_that("nice_interval gives a sensible number of ticks", {
  for (len in c(2686, 4361, 9200, 48502)) {
    k <- len / plasmidplot:::nice_interval(len)
    expect_gte(k, 4)
    expect_lte(k, 14)
  }
})

test_that("marker_polygon is one closed ring, so body and head never seam", {
  # The whole point of the single-polygon build: an arrowed feature must be
  # ONE path. Two shapes would leave an antialiasing gap along the join.
  for (arrow in c("none", "end", "start", "both")) {
    poly <- plasmidplot:::marker_polygon(100, 900, 5000, 0.27, 0.33,
                                         arrow = arrow, head_bp = 80)
    expect_length(poly$x, length(poly$y))
    expect_true(all(is.finite(poly$x)), info = arrow)
    expect_true(all(is.finite(poly$y)), info = arrow)
    # Every vertex stays inside the viewport.
    expect_true(all(poly$x > 0 & poly$x < 1), info = arrow)
    expect_true(all(poly$y > 0 & poly$y < 1), info = arrow)
  }
})

test_that("arrowheads flare past the body and tip at the feature boundary", {
  r0 <- 0.27
  r1 <- 0.33
  rmid <- (r0 + r1) / 2
  poly_of <- function(arrow, head) {
    plasmidplot:::marker_polygon(100, 900, 5000, r0, r1, arrow, head)
  }
  radius <- function(p) sqrt((p$x - 0.5)^2 + (p$y - 0.5)^2)

  plain <- poly_of("none", 0)
  expect_equal(max(radius(plain)), r1, tolerance = 1e-9)
  expect_equal(min(radius(plain)), r0, tolerance = 1e-9)

  # A head is wider than the body it grows out of, on both edges.
  for (arrow in c("end", "start", "both")) {
    p <- poly_of(arrow, 80)
    expect_gt(max(radius(p)), r1)
    expect_lt(min(radius(p)), r0)
  }

  # The tip lands on the mid-radius at exactly the feature boundary, so an
  # arrowed feature still ends where the user said it ends.
  tip <- function(bp) {
    th <- plasmidplot:::pp_theta(bp, 5000)
    c(plasmidplot:::pp_x(th, rmid), plasmidplot:::pp_y(th, rmid))
  }
  near_tip <- function(p, bp) {
    t <- tip(bp)
    min(sqrt((p$x - t[1])^2 + (p$y - t[2])^2))
  }
  expect_lt(near_tip(poly_of("end", 80), 900), 1e-9)
  expect_lt(near_tip(poly_of("start", 80), 100), 1e-9)
  expect_lt(near_tip(poly_of("both", 80), 900), 1e-9)
  expect_lt(near_tip(poly_of("both", 80), 100), 1e-9)
})

test_that("head_length never lets arrowheads exceed the feature", {
  expect_equal(plasmidplot:::head_length(1000, 5000, "none", 6), 0)
  # A tiny feature gets a head no longer than half its span.
  expect_lte(plasmidplot:::head_length(10, 5000, "end", 6), 5)
  # Two heads on a tiny feature each get at most a quarter.
  expect_lte(plasmidplot:::head_length(10, 5000, "both", 6), 2.5)
  # A long feature is capped by arrow_deg, not by its own length.
  expect_equal(plasmidplot:::head_length(4000, 3600, "end", 10), 100)
})

test_that("tick_positions survives an interval wider than the molecule", {
  # seq(from, to, by) throws when `to` falls behind `from`, so an interval
  # bigger than the plasmid used to abort the whole plot.
  for (layout in c("circular", "linear")) {
    expect_no_error(plasmidplot:::tick_positions(1000, 5000, layout))
    expect_no_error(plasmidplot:::tick_positions(1000, 700, layout))
    expect_no_error(plasmidplot:::tick_positions(1, 1, layout))
  }
  # A linear molecule still marks both of its real ends.
  expect_equal(plasmidplot:::tick_positions(1000, 5000, "linear"), c(1, 1000))
  # A circle has no ends, so only the origin survives.
  expect_equal(plasmidplot:::tick_positions(1000, 5000, "circular"), 0)
})

test_that("both layouts mark the scale at the same numbers", {
  # The two layouts drawing 200/400/600 versus 201/401/601 for one plasmid
  # reads as an off-by-one bug to anyone comparing them.
  circ <- plasmidplot:::tick_positions(1200, 200, "circular")
  lin <- plasmidplot:::tick_positions(1200, 200, "linear")
  expect_equal(circ, c(0, 200, 400, 600, 800, 1000))
  expect_equal(lin, c(1, 200, 400, 600, 800, 1000, 1200))
  # Every interior mark agrees; the layouts differ only at the ends.
  expect_true(all(setdiff(circ, 0) %in% lin))
})

test_that("a tick interval wider than the plasmid still plots", {
  p <- pp_marker(plasmid("p", 1000), 100, 400, label = "x")
  tf <- tempfile(fileext = ".png")
  png(tf, width = 500, height = 500)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  for (layout in c("circular", "linear")) {
    expect_no_error(
      plot(p, style = pp_style("angular", layout = layout, tick_every = 5000)))
  }
})

test_that("tick_every rejects values that cannot be an interval", {
  expect_error(pp_style("angular", tick_every = 0), "positive number")
  expect_error(pp_style("angular", tick_every = -100), "positive number")
  expect_error(pp_style("angular", tick_every = c(100, 200)), "positive number")
  expect_no_error(pp_style("angular", tick_every = NA))
  expect_no_error(pp_style("angular", tick_every = 250))
})

test_that("pp_style() presets resolve and reject unknown fields", {
  for (nm in pp_style()) expect_s3_class(pp_style(nm), "pp_style")
  expect_error(pp_style("nope"), "Unknown preset")
  expect_equal(pp_style("angular", label_fontsize = 12)$label_fontsize, 12)
  expect_error(pp_style("angular", nope = 1), "Unknown style field")
})

test_that("no preset paints a background, but each names its canvas", {
  # "Styles must not change the background" is the contract: pp_style() only
  # recommends a ground, plot() paints one only when asked.
  for (nm in pp_style()) {
    expect_true(is.na(pp_style(nm)$bg), info = nm)
    expect_true(plasmidplot:::is_color(pp_canvas(nm)), info = nm)
  }
  expect_equal(pp_canvas("dark"), "#1a1a19")
  expect_equal(pp_canvas(pp_style("neon")), "#0b0f19")
})

test_that("presets are structurally distinct, not recolors", {
  # A style is a construction, not a color scheme. If two presets share every
  # geometric field they are the same style twice and one of them should go.
  geometry_of <- function(nm) {
    s <- pp_style(nm)
    paste(s$backbone, s$radius, s$track,
          s$anchor, s$gap, s$arc, sep = "/")
  }
  nms <- pp_style()
  g <- stats::setNames(vapply(nms, geometry_of, character(1)), nms)
  # `dark` is deliberately angular's construction restepped for a dark
  # ground; every other preset must stand on its own geometry.
  expect_equal(unname(g["dark"]), unname(g["angular"]))
  rest <- g[setdiff(nms, "dark")]
  expect_equal(length(unique(rest)), length(rest))
  # And the three placements are all actually used by some preset.
  anchors <- vapply(nms, function(nm) pp_style(nm)$anchor, character(1))
  expect_setequal(unique(anchors), c("center", "outside", "inside"))
})

test_that("anchor moves features clear of the backbone", {
  base <- pp_style("angular")
  centered <- plasmidplot:::marker_radius(base)
  outside <- plasmidplot:::marker_radius(
    pp_style("angular", anchor = "outside"))
  inside <- plasmidplot:::marker_radius(
    pp_style("angular", anchor = "inside"))
  expect_equal(centered, base$radius)
  expect_gt(outside, centered)
  expect_lt(inside, centered)
  # Outside features clear the ring: their inner edge sits past its outer one.
  expect_gte(outside - base$arc / 2,
             base$radius + base$track / 2)
  expect_error(plasmidplot:::marker_radius(
    pp_style("angular", anchor = "sideways")), "should be one of")
})

test_that("band_extent covers backbone and features whatever the anchor", {
  # The bounds are reached exactly whenever one part is the outermost, so
  # compare with a tolerance rather than strictly.
  eps <- 1e-9
  for (a in c("center", "outside", "inside")) {
    s <- pp_style("angular", anchor = a)
    ext <- plasmidplot:::band_extent(s)
    mid <- plasmidplot:::marker_radius(s)
    expect_lte(ext[["inner"]], mid - s$arc / 2 + eps)
    expect_gte(ext[["outer"]], mid + s$arc / 2 - eps)
    expect_lte(ext[["inner"]], s$radius - s$track / 2 + eps)
    expect_gte(ext[["outer"]], s$radius + s$track / 2 - eps)
  }
})

test_that("style names and palette names never collide", {
  # A style is a look; a palette is colors. Sharing a word between the two
  # namespaces makes pp_style("x") and pp_style(palette = "x") read as the
  # same request when they are not, so keep them disjoint.
  expect_length(intersect(pp_style(), pp_palette()), 0L)
})

test_that("styles for dark grounds use a palette named for it", {
  # The _dark suffix is the signal that a palette needs a dark canvas, so a
  # dark preset must carry one and a light preset must not.
  dark_styles <- c("dark", "neon", "blueprint")
  palette_name_of <- function(colors) {
    for (nm in pp_palette()) {
      if (identical(colors, pp_palette(nm))) return(nm)
    }
    NA_character_
  }
  for (nm in pp_style()) {
    pal <- palette_name_of(pp_style(nm)$palette)
    expect_false(is.na(pal), info = nm)
    expect_equal(grepl("_dark$", pal), nm %in% dark_styles, info = nm)
  }
})

test_that("palettes resolve by name and every preset carries real colors", {
  expect_true(all(c("default", "default_dark", "muted", "vivid", "jewel",
                    "candy", "neon_dark") %in% pp_palette()))
  expect_equal(pp_style("angular", palette = "vivid")$palette,
               pp_palette("vivid"))
  expect_error(pp_palette("nope"), "Unknown palette")
  # A literal color vector is fine; a typo'd palette name is caught here
  # rather than failing later inside grid.
  expect_equal(pp_style("angular", palette = c("red", "#00ff00"))$palette,
               c("red", "#00ff00"))
  expect_error(pp_style("angular", palette = "ocean"), "known palette name")
  for (nm in pp_style()) {
    pal <- pp_style(nm)$palette
    expect_gte(length(pal), 8L)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)), info = nm)
  }
})

test_that("grouped markers share a palette color", {
  p <- plasmid("pX", 5000) |>
    pp_marker(100, 500, group = "CDS") |>
    pp_marker(600, 900, group = "ori") |>
    pp_marker(1000, 1400, group = "CDS") |>
    pp_marker(1500, 1900, color = "#123456")
  cols <- plasmidplot:::marker_colors(p, pp_style("angular"))
  expect_equal(cols[1], cols[3])
  expect_false(cols[1] == cols[2])
  expect_equal(cols[4], "#123456")
})

test_that("pp_shade lightens and darkens", {
  expect_equal(pp_shade("#808080", -1), "#000000")
  expect_equal(pp_shade("#808080", 1), "#FFFFFF")
  expect_equal(plasmidplot:::resolve_border(NA, "#2a78d6"), NA)
  expect_equal(plasmidplot:::resolve_border("auto", "#ffffff"), "#999999")
})

test_that("plotting runs without error in every style", {
  p <- plasmid("pX", 5000) |>
    pp_marker(100, 1500, label = "gene A", arrow = "end") |>
    pp_marker(2000, 2100, label = "tiny") |>
    pp_marker(2500, 3000, label = "two-way", arrow = "both") |>
    pp_marker(4800, 400, label = "wraps", arrow = "start")
  tf <- tempfile(fileext = ".png")
  png(tf, width = 600, height = 600)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  for (nm in pp_style()) expect_no_error(plot(p, style = nm))
  expect_no_error(plot(p, style = pp_style("minimal", backbone = "none")))
})
