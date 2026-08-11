# Thresholds and transforms below mirror the reference validator the built-in
# palettes were checked with, so a palette you bring is held to the same bar.
# ΔE is Euclidean distance in OKLab, x100.
pp_thresholds <- list(
  band = list(light = c(0.43, 0.77), dark = c(0.48, 0.67)),  # OKLCH L
  chroma_floor = 0.10,
  cvd_target = 8.0,     # min(protan, deutan) on the active pair list
  cvd_floor = 6.0,      # below this is a hard failure
  normal_floor = 15.0,  # unsimulated vision, worst pair
  contrast_min = 3.0,   # WCAG vs the surface
  surface = list(light = "#fcfcfb", dark = "#1a1a19")
)

# Machado, Oliveira & Fernandes (2009) transforms at severity 1.0, in linear
# RGB. The thresholds above are calibrated to this model specifically.
pp_machado <- list(
  protan = matrix(c( 0.152286,  1.052583, -0.204868,
                     0.114503,  0.786281,  0.099216,
                    -0.003882, -0.048116,  1.051998),
                  nrow = 3, byrow = TRUE),
  deutan = matrix(c( 0.367322,  0.860646, -0.227968,
                     0.280085,  0.672501,  0.047413,
                    -0.011820,  0.042940,  0.968881),
                  nrow = 3, byrow = TRUE),
  tritan = matrix(c( 1.255528, -0.076749, -0.178779,
                    -0.078411,  0.930809,  0.147602,
                     0.004733,  0.691367,  0.303900),
                  nrow = 3, byrow = TRUE)
)

# 3 x n matrix of linear-light RGB.
lin_rgb <- function(colors) {
  srgb <- grDevices::col2rgb(colors) / 255
  ifelse(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055)^2.4)
}

oklab_from_lin <- function(rgb) {
  l <- (0.4122214708 * rgb[1, ] + 0.5363325363 * rgb[2, ] + 0.0514459929 * rgb[3, ])^(1 / 3)
  m <- (0.2119034982 * rgb[1, ] + 0.6806995451 * rgb[2, ] + 0.1073969566 * rgb[3, ])^(1 / 3)
  s <- (0.0883024619 * rgb[1, ] + 0.2817188376 * rgb[2, ] + 0.6299787005 * rgb[3, ])^(1 / 3)
  rbind(
    L = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
    a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
    b = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s
  )
}

simulate_cvd <- function(rgb, kind) {
  pmin(pmax(pp_machado[[kind]] %*% rgb, 0), 1)
}

relative_luminance <- function(colors) {
  rgb <- lin_rgb(colors)
  0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
}

wcag_contrast <- function(a, b) {
  la <- relative_luminance(a)
  lb <- relative_luminance(b)
  (pmax(la, lb) + 0.05) / (pmin(la, lb) + 0.05)
}

# Worst-case pair distances under one vision model.
pair_delta_e <- function(colors, idx, kind = NULL) {
  rgb <- lin_rgb(colors)
  if (!is.null(kind)) rgb <- simulate_cvd(rgb, kind)
  lab <- oklab_from_lin(rgb)
  apply(idx, 1, function(p) {
    100 * sqrt(sum((lab[, p[1]] - lab[, p[2]])^2))
  })
}

#' Check a palette for colorblind safety and contrast
#'
#' Runs the measurable checks a categorical palette has to pass, using the
#' same transforms and thresholds the built-in palettes were validated with.
#' Use it on any palette you bring from another package — `ggsci`,
#' `RColorBrewer`, `viridis`, `ggpubr` — before relying on it.
#'
#' @details
#' Five checks are computed. Two of them — **CVD separation** and the
#' **normal-vision floor** — decide whether a reader can tell two features
#' apart at all; failing either means the palette is unsafe as it stands. The
#' other three govern how the palette looks and how it sits on the surface,
#' and failing them means the colors are legible but sit outside the band the
#' built-in palettes hold to. The printed summary says which kind failed,
#' because most palettes from other packages miss the cosmetic bands while
#' remaining perfectly readable.
#'
#' * **Lightness band** — OKLCH L inside the band for the mode.
#' * **Chroma floor** — OKLCH C at or above 0.10, below which a hue reads gray.
#' * **CVD separation** — OKLab ΔE between slots under simulated protanopia
#'   and deuteranopia (tritanopia reported alongside). At or above 8 passes;
#'   6 to 8 is a warning band that is sound only with a second channel — on a
#'   plasmid map every feature carries a text label, which supplies it. Below
#'   6 fails.
#' * **Normal-vision floor** — worst pair at or above ΔE 15 under ordinary
#'   vision. This one is a hard gate: below it, full-color readers cannot tell
#'   the two apart either, and a text label does not excuse it.
#' * **Contrast vs surface** — WCAG ratio of at least 3:1. A warning here
#'   obliges visible labels, which the maps always draw.
#'
#' @param colors A vector of colors, a [pp_palette()] name, or a palette
#'   function, in which case `n` colors are drawn from it.
#' @param mode `"light"` or `"dark"`; selects the lightness band and the
#'   default surface.
#' @param surface The background the marks sit on. Defaults to the mode's
#'   standard surface. Pass [pp_canvas()] of the style you will plot with.
#' @param pairs `"adjacent"` checks neighboring slots, which is what matters
#'   when features sit side by side. `"all"` checks every pair, which is the
#'   stricter bar and rarely passes past three or four colors.
#' @param n Number of colors to draw when `colors` is a function.
#' @return An object of class `pp_palette_check`: a data frame of one row per
#'   check with columns `check`, `status` (`"PASS"`, `"WARN"` or `"FAIL"`) and
#'   `detail`, carrying an `ok` attribute that is `FALSE` if anything failed.
#' @examples
#' pp_check_palette("default")
#' pp_check_palette("default_dark", mode = "dark")
#'
#' # A palette brought from elsewhere:
#' pp_check_palette(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))
#'
#' \dontrun{
#' pp_check_palette(ggsci::pal_npg("nrc"), n = 8)
#' pp_check_palette(RColorBrewer::brewer.pal(8, "Dark2"))
#' }
#' @export
pp_check_palette <- function(colors, mode = c("light", "dark"), surface = NULL,
                             pairs = c("adjacent", "all"), n = 8L) {
  mode <- match.arg(mode)
  pairs <- match.arg(pairs)
  if (is.character(colors) && length(colors) == 1L &&
      colors %in% names(pp_palettes)) {
    colors <- pp_palettes[[colors]]
  }
  if (is.function(colors)) colors <- palette_colors(colors, n)
  bad <- colors[!vapply(colors, is_color, logical(1))]
  if (length(bad)) {
    stop("Not colors: ", paste(unique(bad), collapse = ", "), call. = FALSE)
  }
  k <- length(colors)
  if (k < 2L) stop("Need at least two colors to check.", call. = FALSE)
  if (is.null(surface)) surface <- pp_thresholds$surface[[mode]]

  idx <- if (identical(pairs, "adjacent")) {
    cbind(seq_len(k - 1L), seq(2L, k))
  } else {
    t(utils::combn(k, 2L))
  }

  lab <- oklab_from_lin(lin_rgb(colors))
  L <- lab["L", ]
  C <- sqrt(lab["a", ]^2 + lab["b", ]^2)
  band <- pp_thresholds$band[[mode]]

  rows <- list()
  add <- function(check, status, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, status = status, detail = detail,
      stringsAsFactors = FALSE
    )
  }
  off <- which(L < band[1] | L > band[2])
  add("Lightness band",
      if (length(off)) "FAIL" else "PASS",
      if (length(off)) {
        paste0("outside L ", band[1], "-", band[2], ": ",
               paste(sprintf("%s (%.3f)", colors[off], L[off]), collapse = ", "))
      } else sprintf("all %d inside L %.2f-%.2f", k, band[1], band[2]))

  gray <- which(C < pp_thresholds$chroma_floor)
  add("Chroma floor",
      if (length(gray)) "FAIL" else "PASS",
      if (length(gray)) {
        paste0("reads gray: ",
               paste(sprintf("%s (%.3f)", colors[gray], C[gray]), collapse = ", "))
      } else sprintf("all %d at or above %.2f", k, pp_thresholds$chroma_floor))

  d_pro <- pair_delta_e(colors, idx, "protan")
  d_deu <- pair_delta_e(colors, idx, "deutan")
  d_tri <- pair_delta_e(colors, idx, "tritan")
  d_cvd <- pmin(d_pro, d_deu)
  w <- which.min(d_cvd)
  worst_kind <- if (d_pro[w] <= d_deu[w]) "protan" else "deutan"
  add("CVD separation",
      if (d_cvd[w] >= pp_thresholds$cvd_target) "PASS"
      else if (d_cvd[w] >= pp_thresholds$cvd_floor) "WARN" else "FAIL",
      sprintf("worst %s vs %s dE %.1f (%s), tritan %.1f",
              colors[idx[w, 1]], colors[idx[w, 2]], d_cvd[w], worst_kind,
              d_tri[w]))

  d_norm <- pair_delta_e(colors, idx)
  wn <- which.min(d_norm)
  add("Normal-vision floor",
      if (d_norm[wn] >= pp_thresholds$normal_floor) "PASS" else "FAIL",
      sprintf("worst %s vs %s dE %.1f%s", colors[idx[wn, 1]],
              colors[idx[wn, 2]], d_norm[wn],
              if (d_norm[wn] >= pp_thresholds$normal_floor) "" else
                " - below 15, hard to tell apart even with full color vision"))

  cr <- wcag_contrast(colors, surface)
  low <- which(cr < pp_thresholds$contrast_min)
  add("Contrast vs surface",
      if (length(low)) "WARN" else "PASS",
      if (length(low)) {
        paste0("below 3:1, needs visible labels: ",
               paste(sprintf("%s (%.2f)", colors[low], cr[low]), collapse = ", "))
      } else sprintf("all %d at or above 3:1 on %s", k, surface))

  out <- do.call(rbind, rows)
  structure(out, class = c("pp_palette_check", "data.frame"),
            ok = !any(out$status == "FAIL"), palette = colors,
            mode = mode, surface = surface, pairs = pairs)
}

# The two separation checks decide whether a reader can tell two features
# apart; the other three shape how a palette looks and sits on the surface.
# Reporting them apart keeps a cosmetic miss from reading like a blocker.
pp_legibility_checks <- c("CVD separation", "Normal-vision floor")

#' @export
print.pp_palette_check <- function(x, ...) {
  cat(sprintf("Palette check: %d colors, %s mode, surface %s, %s pairs\n",
              length(attr(x, "palette")), attr(x, "mode"),
              attr(x, "surface"), attr(x, "pairs")))
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  [%-4s] %-21s %s\n", x$status[i], x$check[i], x$detail[i]))
  }
  legible <- x[x$check %in% pp_legibility_checks, ]
  cosmetic <- x[!x$check %in% pp_legibility_checks, ]
  bad_legible <- legible$check[legible$status == "FAIL"]
  bad_cosmetic <- cosmetic$check[cosmetic$status == "FAIL"]

  if (length(bad_legible)) {
    cat("  -> NOT SAFE: ", paste(bad_legible, collapse = ", "),
        " failed. Some features will be indistinguishable; ",
        "cut colors, reorder slots, or facet.\n", sep = "")
  } else if (length(bad_cosmetic)) {
    cat("  -> readable, but off-house-style: ",
        paste(bad_cosmetic, collapse = ", "),
        " failed. Features stay tellable apart; the palette just sits ",
        "outside the lightness/chroma band the built-ins hold to.\n", sep = "")
  } else {
    cat("  -> good (WARN items rely on the labels the maps already draw)\n")
  }
  invisible(x)
}
