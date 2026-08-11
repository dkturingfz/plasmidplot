# A style is six shape parameters plus secondary ink and type settings. The
# shape parameters are the design; the presets below are nothing more than
# named combinations of them, so anything a preset can express you can write
# out by hand.
pp_shape_fields <- c("layout", "anchor", "backbone", "radius", "track", "arc",
                     "gap")

pp_style_base <- function() {
  list(
    # --- shape -------------------------------------------------------------
    layout = "auto",      # "auto" follows the plasmid, else circular/linear
    anchor = "center",    # features sit on the backbone, "outside", "inside"
    backbone = "ring",    # "ring" (a band), "line" (a stroke), "none"
    radius = 0.30,        # circular only: backbone radius, npc
    track = 0.045,        # backbone thickness
    arc = 0.058,          # feature thickness
    gap = 0.012,          # clearance when anchor is outside or inside

    # --- ink ---------------------------------------------------------------
    # Styles never paint a background: whatever the device or the surrounding
    # viewport already has shows through. `canvas` only *recommends* a ground
    # for the dark presets; pass it to plot(bg=) to actually paint it.
    bg = NA,
    canvas = "#fcfcfb",
    track_fill = "#e1e0d9",
    track_col = NA,
    backbone_col = "#0b0b0b",
    backbone_lwd = 1.2,
    palette = "default",
    arc_border = NA,
    arc_lwd = 0.8,
    arrow_deg = 6,
    arrow_flare = 0.35,

    # --- labels, sites, scale, title ---------------------------------------
    label_r = 0.365,
    label_fontsize = 9,
    label_fontface = "plain",
    label_col = "#0b0b0b",
    leader_col = "#898781",
    leader_lwd = 0.8,
    site_col = NA,
    site_lwd = 0.9,
    site_len = 0.022,
    site_fontsize = 8,
    tick_every = NA,
    tick_len = 0.012,
    tick_col = "#898781",
    tick_fontsize = 7,
    tick_labels = TRUE,
    show_title = TRUE,
    title_col = "#0b0b0b",
    title_fontsize = 14,
    subtitle_col = "#52514e",
    subtitle_fontsize = 9
  )
}

# Each preset is a combination of the shape parameters (plus the ink a dark
# ground needs). Nothing here is reachable only through a preset.
pp_presets <- list(
  angular = list(),
  classic = list(
    backbone = "line", radius = 0.26, anchor = "outside",
    gap = 0.008, arc = 0.050,
    canvas = "#ffffff",
    track_fill = NA,
    arc_border = "#0b0b0b",
    leader_col = "#52514e",
    tick_col = "#52514e"
  ),
  dark = list(
    canvas = "#1a1a19",
    track_fill = "#383835",
    backbone_col = "#ffffff",
    palette = "default_dark",
    label_col = "#ffffff",
    title_col = "#ffffff",
    subtitle_col = "#c3c2b7"
  ),
  snapgene = list(
    track = 0.012, arc = 0.042,
    canvas = "#ffffff",
    track_fill = "#c8c8c8",
    arc_border = "auto",
    arc_lwd = 0.7,
    arrow_deg = 5,
    arrow_flare = 0.3,
    leader_col = "#b0b0b0",
    tick_col = "#8c8c8c",
    tick_fontsize = 6.5
  ),
  soft = list(
    radius = 0.28, track = 0.09, arc = 0.09,
    canvas = "#fdfcfa",
    track_fill = "#ece9e2",
    palette = "muted",
    arrow_deg = 8,
    arrow_flare = 0.25,
    label_fontsize = 10,
    label_col = "#3b3a37",
    leader_col = "#b8b5ad",
    tick_labels = FALSE,
    tick_len = 0.02,
    tick_col = "#c4c1b8",
    title_col = "#2b2a28",
    subtitle_col = "#6d6b66"
  ),
  neon = list(
    radius = 0.25, track = 0.03, anchor = "outside", gap = 0.018, arc = 0.050,
    canvas = "#0b0f19",
    track_fill = "#1c2438",
    backbone_col = "#e6edf7",
    palette = "neon_dark",
    label_col = "#e6edf7",
    label_fontsize = 9.5,
    leader_col = "#5b6b86",
    tick_col = "#6e7f9c",
    title_col = "#ffffff",
    subtitle_col = "#9fb0cc"
  ),
  minimal = list(
    backbone = "line", arc = 0.022,
    canvas = "#ffffff",
    backbone_col = "#c3c2b7",
    backbone_lwd = 0.8,
    track_fill = NA,
    arrow_deg = 4,
    label_r = 0.345,
    label_fontsize = 8.5,
    label_col = "#3b3a37",
    leader_col = "#c8c6bf",
    leader_lwd = 0.6,
    tick_len = 0.008,
    tick_fontsize = 6.5,
    tick_col = "#a8a59d",
    title_fontsize = 12,
    subtitle_fontsize = 8
  ),
  blueprint = list(
    backbone = "line", radius = 0.33, anchor = "inside",
    gap = 0.014, arc = 0.050,
    canvas = "#0d2137",
    backbone_col = "#7fb2d9",
    backbone_lwd = 1,
    track_fill = NA,
    palette = "default_dark",
    arc_border = "#cfe3f2",
    arc_lwd = 0.9,
    label_col = "#e8f2fa",
    label_fontface = "bold",
    label_r = 0.385,
    leader_col = "#5b8bb0",
    tick_col = "#7fb2d9",
    title_col = "#ffffff",
    subtitle_col = "#a8cbe4"
  )
)

#' Create or customize a plasmid map style
#'
#' A style is built from a handful of **shape** parameters plus secondary ink
#' and type settings. Set the shape parameters directly; the presets are only
#' named combinations of them, and everything a preset does can be written out
#' by hand.
#'
#' @section Shape:
#' These decide the construction of the map. They are the arguments worth
#' reaching for first:
#'
#' \describe{
#'   \item{`layout`}{`"auto"` draws a circular map for a circular plasmid and
#'     a linear one for a linear plasmid, following the object's `topology`.
#'     Force it with `"circular"` or `"linear"`.}
#'   \item{`anchor`}{Where features sit relative to the backbone: `"center"`
#'     on it, `"outside"` clear of it (outward on a circle, above the line on
#'     a linear map), or `"inside"`.}
#'   \item{`backbone`}{`"ring"` draws the backbone as a filled band, `"line"`
#'     as a single stroke, `"none"` omits it.}
#'   \item{`radius`}{Circular only: the backbone's radius, in npc units of a
#'     square viewport, so the usable maximum is 0.5.}
#'   \item{`track`}{Backbone thickness, used when `backbone = "ring"`.}
#'   \item{`arc`}{Feature thickness.}
#'   \item{`gap`}{Clearance between backbone and features when `anchor` is
#'     `"outside"` or `"inside"`.}
#' }
#'
#' @section Presets:
#' Every preset is exactly the shape below, plus colors. `dark` is `angular`
#' restepped for a dark ground, which is why their shapes match.
#'
#' | preset | layout | anchor | backbone | radius | track | arc |
#' | --- | --- | --- | --- | --- | --- | --- |
#' | `angular` | auto | center | ring | 0.30 | 0.045 | 0.058 |
#' | `classic` | auto | outside | line | 0.26 | 0.045 | 0.050 |
#' | `dark` | auto | center | ring | 0.30 | 0.045 | 0.058 |
#' | `snapgene` | auto | center | ring | 0.30 | 0.012 | 0.042 |
#' | `soft` | auto | center | ring | 0.28 | 0.090 | 0.090 |
#' | `neon` | auto | outside | ring | 0.25 | 0.030 | 0.050 |
#' | `minimal` | auto | center | line | 0.30 | 0.045 | 0.022 |
#' | `blueprint` | auto | inside | line | 0.33 | 0.045 | 0.050 |
#'
#' Their default palettes are `default`, except `soft` (`muted`), `dark` and
#' `blueprint` (`default_dark`) and `neon` (`neon_dark`). A palette is an
#' independent choice — see [pp_palette()].
#'
#' @param preset Base preset to start from; see the table above. Call
#'   `pp_style()` with no arguments for the names.
#' @param ... Named overrides. The shape parameters are `layout`, `anchor`,
#'   `backbone`, `radius`, `track`, `arc` and `gap`. The rest are:
#'   `bg` (`NA` in every preset — styles do not paint a background; see
#'   [pp_canvas()]), `canvas`, `track_fill`, `track_col`, `backbone_col`,
#'   `backbone_lwd`, `palette` (a [pp_palette()] name, a vector of colors, or
#'   a function of `n` returning `n` colors — the shape `ggsci`,
#'   `RColorBrewer`, `scales` and `viridis` all expose, so those drop straight
#'   in; check any of them with [pp_check_palette()]), `arc_border` (`NA`,
#'   `"auto"` for a darker step of the fill, or a color), `arc_lwd`,
#'   `arrow_deg`, `arrow_flare`, `label_r`, `label_fontsize`,
#'   `label_fontface`, `label_col`, `leader_col`, `leader_lwd`, `site_col`
#'   (`NA` follows `label_col`), `site_lwd`, `site_len`, `site_fontsize`,
#'   `tick_every` (`NA` = automatic), `tick_len`, `tick_col`, `tick_fontsize`,
#'   `tick_labels`, `show_title`, `title_col`, `title_fontsize`,
#'   `subtitle_col`, `subtitle_fontsize`.
#' @return An object of class `pp_style`, or a character vector of preset
#'   names when called with no arguments.
#' @examples
#' pp_style()
#' p <- pp_demo_plasmid()
#'
#' # Compose a style from the shape parameters, no preset involved.
#' plot(p, style = pp_style(anchor = "outside", backbone = "line",
#'                          radius = 0.26, arc = 0.05))
#'
#' # Or start from a preset and change one thing.
#' plot(p, style = pp_style("angular", arc = 0.09))
#' plot(p, style = pp_style("minimal", layout = "linear"))
#'
#' # Any palette vector or palette function works, e.g.
#' # pp_style("angular", palette = ggsci::pal_npg("nrc"))
#' plot(p, style = pp_style("angular", palette = grDevices::hcl.colors(8, "Set2")))
#' @export
pp_style <- function(preset, ...) {
  if (missing(preset) && !...length()) return(names(pp_presets))
  if (missing(preset)) preset <- "angular"
  if (!is.character(preset) || length(preset) != 1L ||
      !preset %in% names(pp_presets)) {
    stop("Unknown preset. Available: ",
         paste(names(pp_presets), collapse = ", "), call. = FALSE)
  }
  s <- utils::modifyList(pp_style_base(), pp_presets[[preset]])
  dots <- list(...)
  if (length(dots)) {
    unknown <- setdiff(names(dots), names(s))
    if (length(unknown)) {
      stop("Unknown style field(s): ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    s[names(dots)] <- dots
  }
  s$palette <- resolve_palette(s$palette)
  if (length(s$tick_every) != 1L ||
      (!is.na(s$tick_every) && (!is.numeric(s$tick_every) || s$tick_every <= 0))) {
    stop("`tick_every` must be a single positive number, or NA for automatic.",
         call. = FALSE)
  }
  s$layout <- match.arg(s$layout, c("auto", "circular", "linear"))
  s$anchor <- match.arg(s$anchor, c("center", "outside", "inside"))
  s$backbone <- match.arg(s$backbone, c("ring", "line", "none"))
  structure(s, class = "pp_style", preset = preset)
}

#' The background a style is designed for
#'
#' Styles do not paint a background; the device or surrounding viewport shows
#' through. This returns the ground a style was designed against, so a dark
#' preset can be given one explicitly.
#'
#' @param style A preset name or [pp_style()] object.
#' @return A color string.
#' @examples
#' pp_canvas("dark")
#' p <- pp_demo_plasmid()
#' plot(p, style = "neon", bg = pp_canvas("neon"))
#' @export
pp_canvas <- function(style) resolve_style(style)$canvas

#' @export
print.pp_style <- function(x, ...) {
  cat(sprintf("<pp_style> preset '%s'\n", attr(x, "preset")))
  cat("  shape: ",
      paste(sprintf("%s=%s", pp_shape_fields,
                    vapply(x[pp_shape_fields], format, character(1))),
            collapse = "  "), "\n", sep = "")
  cat(sprintf("  palette: %d colors\n",
              if (is.function(x$palette)) NA_integer_ else length(x$palette)))
  invisible(x)
}

# Accept either a preset name or a ready pp_style object.
resolve_style <- function(style) {
  if (inherits(style, "pp_style")) return(style)
  if (is.character(style) && length(style) == 1L) return(pp_style(style))
  stop("`style` must be a preset name or a pp_style() object.", call. = FALSE)
}

# "auto" defers to the plasmid; anything else is the user overriding it.
resolve_layout <- function(s, p) {
  if (!identical(s$layout, "auto")) return(s$layout)
  if (identical(p$topology, "linear")) "linear" else "circular"
}
