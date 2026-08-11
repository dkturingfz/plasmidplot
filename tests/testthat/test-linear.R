test_that("layout follows the plasmid's topology unless overridden", {
  circ <- plasmid("c", 1000)
  lin <- plasmid("l", 1000)
  lin$topology <- "linear"
  auto <- pp_style("angular")
  expect_equal(plasmidplot:::resolve_layout(auto, circ), "circular")
  expect_equal(plasmidplot:::resolve_layout(auto, lin), "linear")
  # An explicit layout wins in both directions.
  expect_equal(
    plasmidplot:::resolve_layout(pp_style("angular", layout = "linear"), circ),
    "linear")
  expect_equal(
    plasmidplot:::resolve_layout(pp_style("angular", layout = "circular"), lin),
    "circular")
  expect_error(pp_style("angular", layout = "spiral"), "should be one of")
})

test_that("lin_x maps the sequence across the drawing area", {
  s <- pp_style("angular")
  expect_equal(plasmidplot:::lin_x(1, 1000, s), plasmidplot:::lin_x0(s))
  expect_equal(plasmidplot:::lin_x(1000, 1000, s), plasmidplot:::lin_x1(s))
  mid <- plasmidplot:::lin_x(500.5, 1000, s)
  expect_equal(mid, (plasmidplot:::lin_x0(s) + plasmidplot:::lin_x1(s)) / 2)
  # Monotone, so features never cross each other.
  xs <- plasmidplot:::lin_x(c(1, 250, 500, 750, 1000), 1000, s)
  expect_false(is.unsorted(xs))
  # Centered whatever the extent.
  wide <- pp_style("angular", radius = 0.33)
  expect_equal(plasmidplot:::lin_x0(wide) + plasmidplot:::lin_x1(wide), 1)
})

test_that("radius sets how much room the map takes in BOTH layouts", {
  # radius would otherwise be inert on a linear map, leaving one of the seven
  # shape parameters with no effect and the styles looking alike.
  small <- pp_style("angular", radius = 0.24)
  big <- pp_style("angular", radius = 0.33)
  span <- function(s) plasmidplot:::lin_x1(s) - plasmidplot:::lin_x0(s)
  expect_gt(span(big), span(small))
  expect_gt(plasmidplot:::marker_radius(big), plasmidplot:::marker_radius(small))
  # It stays inside the viewport however large radius gets.
  huge <- pp_style("angular", radius = 2)
  expect_gte(plasmidplot:::lin_x0(huge), 0)
  expect_lte(plasmidplot:::lin_x1(huge), 1)
})

test_that("a wrapping feature becomes two pieces on a linear map", {
  one <- plasmidplot:::linear_pieces(list(start = 100, end = 900), 1000)
  expect_length(one, 1L)
  expect_equal(one[[1]], c(100, 900))

  two <- plasmidplot:::linear_pieces(list(start = 900, end = 100), 1000)
  expect_length(two, 2L)
  expect_equal(two[[1]], c(900, 1000))
  expect_equal(two[[2]], c(1, 100))
})

test_that("anchor_offset places features the same way in both layouts", {
  # The anchor rule is shared, which is what keeps a style recognizable when
  # it is switched from circular to linear.
  for (a in c("center", "outside", "inside")) {
    s <- pp_style("angular", anchor = a)
    off <- plasmidplot:::anchor_offset(s)
    expect_equal(plasmidplot:::marker_radius(s), s$radius + off)
  }
  expect_equal(plasmidplot:::anchor_offset(pp_style("angular",
                                                    anchor = "center")), 0)
  expect_gt(plasmidplot:::anchor_offset(pp_style("angular",
                                                 anchor = "outside")), 0)
  expect_lt(plasmidplot:::anchor_offset(pp_style("angular",
                                                 anchor = "inside")), 0)
})

test_that("lin_marker_polygon is one closed outline with a flared head", {
  s <- pp_style("angular")
  y <- 0.4
  w <- 0.06
  plain <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, "none", s)
  expect_equal(max(plain$y), y + w / 2, tolerance = 1e-9)
  expect_equal(min(plain$y), y - w / 2, tolerance = 1e-9)

  for (arrow in c("end", "start", "both")) {
    poly <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, arrow, s)
    expect_length(poly$x, length(poly$y))
    expect_true(all(is.finite(poly$x)), info = arrow)
    # The head flares past the body on both edges.
    expect_gt(max(poly$y), y + w / 2)
    expect_lt(min(poly$y), y - w / 2)
  }

  # The tip reaches the feature boundary and stops there.
  e <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, "end", s)
  expect_equal(max(e$x), plasmidplot:::lin_x(900, 1000, s), tolerance = 1e-9)
  st <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, "start", s)
  expect_equal(min(st$x), plasmidplot:::lin_x(100, 1000, s), tolerance = 1e-9)
})

test_that("a linear feature body is a rectangle, not a wedge", {
  # Naming only some of the body corners lets the top edge slope into the
  # head, which turns every feature into a triangle.
  s <- pp_style("angular")
  y <- 0.4
  w <- 0.06
  flat <- function(poly, level) sum(abs(poly$y - level) < 1e-9)

  plain <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, "none", s)
  expect_length(plain$x, 4L)
  expect_equal(flat(plain, y + w / 2), 2L)
  expect_equal(flat(plain, y - w / 2), 2L)

  # With a head, the body still contributes two vertices to each long edge:
  # its far corner and the base of the head.
  for (arrow in c("end", "start", "both")) {
    poly <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, arrow, s)
    expect_equal(flat(poly, y + w / 2), 2L, info = arrow)
    expect_equal(flat(poly, y - w / 2), 2L, info = arrow)
  }

  # The head takes only a slice: the body keeps most of the feature's width.
  e <- plasmidplot:::lin_marker_polygon(100, 900, 1000, y, w, "end", s)
  body <- max(e$x[abs(e$y - (y + w / 2)) < 1e-9]) -
    min(e$x[abs(e$y - (y + w / 2)) < 1e-9])
  full <- plasmidplot:::lin_x(900, 1000, s) - plasmidplot:::lin_x(100, 1000, s)
  expect_gt(body, full * 0.75)
})

test_that("linear labels are pulled inside the viewport at the ends", {
  # A long label on the very first or last feature would hang off the page.
  long <- strrep("wide label ", 4)
  p <- plasmid("edges", 1000) |>
    pp_marker(1, 60, label = long) |>
    pp_marker(940, 1000, label = long)
  tf <- tempfile(fileext = ".png")
  png(tf, width = 500, height = 400)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  s <- pp_style("angular", layout = "linear")
  entries <- plasmidplot:::lin_label_entries(p, s)
  widths <- vapply(entries, function(e) {
    plasmidplot:::text_width_npc(e$text, e$fontsize)
  }, numeric(1))
  xs <- vapply(entries, function(e) e$x, numeric(1))
  clamped <- pmin(pmax(xs, widths / 2 + 0.006), 1 - widths / 2 - 0.006)
  expect_true(all(clamped - widths / 2 >= 0))
  expect_true(all(clamped + widths / 2 <= 1))
  # Labels this wide at both ends must actually have been moved, so the
  # clamp is doing work rather than the assertion passing vacuously.
  expect_true(any(abs(clamped - xs) > 1e-9))
  expect_no_error(plot(p, style = s))
})

test_that("linear labels stack into rows instead of overlapping", {
  # Ten features crowded into a short stretch cannot share one row.
  p <- plasmid("crowded", 2000)
  for (i in 1:10) {
    p <- pp_marker(p, i * 40, i * 40 + 25, label = paste0("feature", i))
  }
  tf <- tempfile(fileext = ".png")
  png(tf, width = 900, height = 500)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  expect_no_error(plot(p, style = pp_style("angular", layout = "linear")))
})

test_that("linear maps render in every style, with features and sites", {
  p <- pp_find_sites(
    read_genbank(system.file("extdata", "pDemo.gb", package = "plasmidplot")))
  tf <- tempfile(fileext = ".png")
  png(tf, width = 1000, height = 500)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  for (nm in pp_style()) {
    expect_no_error(plot(p, style = pp_style(nm, layout = "linear")))
  }
})

test_that("a linear plasmid read from file plots linear by default", {
  fa <- tempfile(fileext = ".fa")
  on.exit(unlink(fa), add = TRUE)
  writeLines(c(">pLin", strrep("ACGT", 300)), fa)
  p <- read_fasta(fa, circular = FALSE)
  expect_equal(p$topology, "linear")
  expect_equal(plasmidplot:::resolve_layout(pp_style("angular"), p), "linear")

  tf <- tempfile(fileext = ".png")
  png(tf, width = 900, height = 500)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  expect_no_error(plot(pp_marker(p, 100, 500, label = "gene", arrow = "end")))
})
