test_that("a palette function is called with the number of colors needed", {
  asked <- integer(0)
  fake <- function(n) {
    asked <<- c(asked, n)
    grDevices::hcl.colors(n, "Dark 3")
  }
  p <- plasmid("pX", 5000) |>
    pp_marker(100, 900) |>
    pp_marker(1000, 1900) |>
    pp_marker(2000, 2900)
  cols <- plasmidplot:::marker_colors(p, pp_style("angular", palette = fake))
  expect_equal(asked, 3L)
  expect_length(unique(cols), 3L)
})

test_that("a palette function is asked only for slots that need a color", {
  asked <- integer(0)
  fake <- function(n) {
    asked <<- c(asked, n)
    grDevices::hcl.colors(max(n, 2L), "Dark 3")
  }
  # Two markers share a group and one names its own color, so only two
  # palette slots are actually needed.
  p <- plasmid("pX", 5000) |>
    pp_marker(100, 900, group = "CDS") |>
    pp_marker(1000, 1900, group = "CDS") |>
    pp_marker(2000, 2900, group = "ori") |>
    pp_marker(3000, 3900, color = "#123456")
  cols <- plasmidplot:::marker_colors(p, pp_style("angular", palette = fake))
  expect_equal(asked, 2L)
  expect_equal(cols[1], cols[2])
  expect_equal(cols[4], "#123456")
})

test_that("palette functions that run short warn and cycle", {
  short <- function(n) c("#2a78d6", "#eb6834")
  p <- plasmid("pX", 5000)
  for (i in 1:4) p <- pp_marker(p, i * 1000 - 900, i * 1000 - 600)
  expect_warning(
    cols <- plasmidplot:::marker_colors(p, pp_style("angular", palette = short)),
    "colors will repeat"
  )
  expect_equal(cols[1], cols[3])
  expect_equal(cols[2], cols[4])
})

test_that("palette functions returning NA padding are tolerated", {
  # ggsci-style: caps out and pads with NA rather than erroring.
  padded <- function(n) c("#2a78d6", "#eb6834", "#1baf7a", rep(NA_character_,
                                                               max(0, n - 3L)))
  expect_equal(plasmidplot:::palette_colors(padded, 3L),
               c("#2a78d6", "#eb6834", "#1baf7a"))
  expect_warning(plasmidplot:::palette_colors(padded, 6L), "colors will repeat")
})

test_that("bad palette functions are reported clearly", {
  expect_error(plasmidplot:::palette_colors(function(n) 1:n, 3L),
               "character vector of colors")
  expect_error(plasmidplot:::palette_colors(function(n) rep("nope", n), 3L),
               "not colors")
  expect_error(plasmidplot:::palette_colors(function(n) character(0), 3L),
               "no usable colors")
})

test_that("a plasmid plots with an external-style palette function", {
  p <- pp_demo_plasmid()
  tf <- tempfile(fileext = ".png")
  png(tf, width = 600, height = 600)
  on.exit({
    dev.off()
    unlink(tf)
  }, add = TRUE)
  expect_no_error(plot(p, style = pp_style("angular",
                                           palette = grDevices::rainbow)))
  expect_no_error(plot(p, style = pp_style("angular",
                                           palette = function(n) {
                                             grDevices::hcl.colors(n, "Set2")
                                           })))
})

# --- pp_check_palette -------------------------------------------------------

test_that("pp_check_palette agrees with the reference numbers", {
  # These are the values the built-in palettes were validated against; they
  # pin the color math (Machado simulation, OKLab dE, WCAG contrast).
  chk <- pp_check_palette("default")
  expect_equal(chk$status[chk$check == "Lightness band"], "PASS")
  expect_equal(chk$status[chk$check == "Chroma floor"], "PASS")
  expect_match(chk$detail[chk$check == "CVD separation"], "dE 9\\.1 \\(protan\\)")
  expect_match(chk$detail[chk$check == "Normal-vision floor"], "dE 19\\.6")
  expect_equal(chk$status[chk$check == "Contrast vs surface"], "WARN")

  dark <- pp_check_palette("default_dark", mode = "dark")
  expect_true(all(dark$status[dark$check != "Contrast vs surface"] == "PASS"))
  expect_match(dark$detail[dark$check == "CVD separation"], "dE 8\\.4")
})

test_that("pp_check_palette flags a palette that is genuinely unsafe", {
  browns <- c("#8c5a3c", "#c9a227", "#6b8e5a", "#b5651d",
              "#7d6b5d", "#a3623c", "#556b4f", "#8a7250")
  chk <- pp_check_palette(browns)
  expect_false(attr(chk, "ok"))
  expect_equal(chk$status[chk$check == "CVD separation"], "FAIL")
  expect_equal(chk$status[chk$check == "Normal-vision floor"], "FAIL")
  expect_match(chk$detail[chk$check == "CVD separation"], "dE 3\\.3")
  expect_output(print(chk), "NOT SAFE")
})

test_that("pp_check_palette separates legibility from house style", {
  # neon_dark is readable but deliberately brighter than the lightness band.
  chk <- pp_check_palette("neon_dark", mode = "dark", surface = "#0b0f19")
  expect_equal(chk$status[chk$check == "Lightness band"], "FAIL")
  expect_equal(chk$status[chk$check == "Normal-vision floor"], "PASS")
  expect_output(print(chk), "off-house-style")
  expect_output(print(pp_check_palette("jewel")), "good")
})

test_that("pp_check_palette accepts names, vectors and functions", {
  expect_s3_class(pp_check_palette("vivid"), "pp_palette_check")
  expect_s3_class(pp_check_palette(c("#E64B35", "#4DBBD5", "#00A087")),
                  "pp_palette_check")
  chk <- pp_check_palette(function(n) grDevices::hcl.colors(n, "Dark 3"), n = 5L)
  expect_equal(length(attr(chk, "palette")), 5L)
  # 8-digit hex with alpha, as ggsci emits, is accepted.
  expect_s3_class(pp_check_palette(c("#E64B35FF", "#4DBBD5FF", "#3C5488FF")),
                  "pp_palette_check")
  expect_error(pp_check_palette(c("#2a78d6", "notacolor")), "Not colors")
  expect_error(pp_check_palette("#2a78d6"), "at least two")
})

test_that("pp_check_palette pair modes and surfaces are honored", {
  pal <- pp_palette("default")
  adj <- pp_check_palette(pal, pairs = "adjacent")
  all <- pp_check_palette(pal, pairs = "all")
  expect_equal(attr(adj, "pairs"), "adjacent")
  # Checking every pair can only find a worse worst-case than neighbors alone.
  worst <- function(chk) {
    as.numeric(sub(".*dE ([0-9.]+).*", "\\1",
                   chk$detail[chk$check == "Normal-vision floor"]))
  }
  expect_lte(worst(all), worst(adj))
  expect_equal(attr(pp_check_palette(pal, surface = "#ffffff"), "surface"),
               "#ffffff")
})
