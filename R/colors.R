# Lighten (amount > 0) or darken (amount < 0) a color toward white or black.
pp_shade <- function(col, amount) {
  m <- grDevices::col2rgb(col, alpha = FALSE) / 255
  m <- if (amount < 0) m * (1 + amount) else m + (1 - m) * amount
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

#' Built-in categorical palettes
#'
#' Named color sets that can be passed to [pp_style()] as `palette = "<name>"`.
#'
#' A palette is only the colors. The geometry, chrome and typography of a map
#' are a *style* — see [pp_style()] — and the two are chosen independently:
#' any palette can be paired with any style. Names ending in `_dark` are
#' stepped for a dark ground and will be washed out on a light one.
#'
#' @details
#' Every palette here was checked pair-by-pair under protanopia, deuteranopia
#' and tritanopia simulation against the surface it is meant for, and against
#' a normal-vision separation floor:
#'
#' * `"default"` (light), `"default_dark"` (dark) and `"jewel"` (light) clear
#'   every gate, contrast included.
#' * `"vivid"` and `"candy"` (light) clear every separation gate.
#' * `"muted"` (light) clears the separation floors, with its closest
#'   colorblind pair in the band that is sound only because each feature also
#'   carries a text label — identity here never rests on hue alone.
#' * `"neon_dark"` clears contrast and both separation floors on a dark
#'   ground. It sits deliberately brighter than the others; that is the look.
#'
#' Slot order carries the safety, so take the colors in the order given.
#' Palettes of near-neighbor hues (all-brown, all-blue sets) were dropped
#' during development because adjacent features became indistinguishable
#' under red-green colorblindness.
#'
#' @param name Palette name. If missing, returns the names of all palettes.
#' @return A character vector of hex colors, or of palette names.
#' @examples
#' pp_palette()
#' pp_palette("muted")
#' plot(pp_demo_plasmid(), style = pp_style("angular", palette = "vivid"))
#' @export
pp_palette <- function(name) {
  if (missing(name)) return(names(pp_palettes))
  if (!is.character(name) || length(name) != 1L || !name %in% names(pp_palettes)) {
    stop("Unknown palette. Available: ",
         paste(names(pp_palettes), collapse = ", "), call. = FALSE)
  }
  pp_palettes[[name]]
}

# Slot order is the colorblind-safety mechanism, not decoration: each of these
# was checked pair-by-pair in this order against protan/deutan/tritan
# simulation. Reordering a palette invalidates that; re-check if you do.
pp_palettes <- list(
  default = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
              "#e87ba4", "#008300", "#4a3aa7", "#e34948"),
  # The same eight hues re-stepped for a dark ground, not a separate palette.
  default_dark = c("#3987e5", "#d95926", "#199e70", "#c98500",
                   "#d55181", "#008300", "#9085e9", "#e66767"),
  muted = c("#5b93d1", "#e08461", "#3fae8b", "#d9a327",
            "#dd7fa8", "#7aa653", "#8f83c9", "#dd6f6f"),
  vivid = c("#0b6fd8", "#ff5a1f", "#00a878", "#e09000",
            "#e5399b", "#00934a", "#5a2ec4", "#e01b24"),
  jewel = c("#1f6fb2", "#c1442f", "#118a6e", "#b8860b",
            "#8e3b8e", "#2f7d32", "#4b3f9e", "#b02a4a"),
  candy = c("#3d8bfd", "#ff7043", "#26c6a2", "#d9aa05",
            "#ec5f9e", "#66bb44", "#7e6cf0", "#ef5350"),
  neon_dark = c("#4cc9f0", "#f72585", "#4ade80", "#fbbf24",
                "#c77dff", "#22d3ee", "#a3e635", "#fb7185")
)

# Accept a palette name, a literal vector of colors, or a palette *function*
# taking n and returning n colors -- which is the shape ggsci, RColorBrewer,
# scales and viridis all expose, so those work with no adapter. Colors are
# checked here rather than at draw time, so a typo'd palette name is reported
# against the pp_style() call that caused it instead of failing inside grid.
resolve_palette <- function(p) {
  if (is.function(p)) return(p)
  if (is.character(p) && length(p) == 1L && p %in% names(pp_palettes)) {
    return(pp_palettes[[p]])
  }
  if (!is.character(p) || !length(p)) {
    stop("`palette` must be a palette name, a vector of colors, or a ",
         "function of n returning n colors. Available palettes: ",
         paste(names(pp_palettes), collapse = ", "), call. = FALSE)
  }
  bad <- p[!vapply(p, is_color, logical(1))]
  if (length(bad)) {
    stop("`palette` is neither a known palette name nor a vector of colors ",
         "(invalid: ", paste(unique(bad), collapse = ", "), "). ",
         "Available palettes: ", paste(names(pp_palettes), collapse = ", "),
         call. = FALSE)
  }
  p
}

# Draw n colors from whatever the style's palette turned out to be. Palette
# functions are called with the count actually needed, so a map with twelve
# features gets twelve distinct colors instead of a cycled eight.
palette_colors <- function(p, n) {
  n <- max(as.integer(n), 1L)
  if (is.function(p)) {
    out <- withCallingHandlers(
      p(n),
      warning = function(w) invokeRestart("muffleWarning")
    )
    if (!is.character(out)) {
      stop("The `palette` function must return a character vector of colors.",
           call. = FALSE)
    }
    out <- out[!is.na(out)]
    bad <- out[!vapply(out, is_color, logical(1))]
    if (length(bad)) {
      stop("The `palette` function returned values that are not colors: ",
           paste(unique(bad), collapse = ", "), call. = FALSE)
    }
    if (!length(out)) {
      stop("The `palette` function returned no usable colors.", call. = FALSE)
    }
    # Many palette functions cap out (ggsci's top out around 10) and pad with
    # NA. Say so once rather than silently repeating colors.
    if (length(out) < n) {
      warning("The palette function supplied ", length(out), " colors for ",
              n, " features; colors will repeat. Pass a longer palette, or ",
              "group features with pp_marker(group=).", call. = FALSE)
    }
    return(out)
  }
  p
}

is_color <- function(x) {
  !is.na(x) && !inherits(try(grDevices::col2rgb(x), silent = TRUE), "try-error")
}

# `NA` = no border, "auto" = a darker step of the fill itself.
resolve_border <- function(spec, fill) {
  if (is.null(spec)) return(NA)
  if (length(spec) == 1L && is.na(spec)) return(NA)
  if (identical(spec, "auto")) return(pp_shade(fill, -0.4))
  spec
}
