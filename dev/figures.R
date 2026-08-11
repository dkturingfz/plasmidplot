# Generates everything in man/figures. One script so the sheets stay
# consistent with each other: each varies exactly ONE dimension and says so.
setwd("i:/plasmidplot")
pkgload::load_all(quiet = TRUE)
unlink(list.files("man/figures", full.names = TRUE))
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

DARK_STYLES <- c("dark", "neon", "blueprint")
DARK_PALETTES <- c("default_dark", "neon_dark")

# What actually differs between the presets, stated outright on each panel.
STYLE_NOTES <- c(
  angular = "features on a broad ring",
  classic = "features outside a thin circle",
  dark = "angular, dark ground",
  snapgene = "hairline track, edged features",
  soft = "fat features, unnumbered scale",
  neon = "slim features outside the track",
  minimal = "hairline everything",
  blueprint = "features inside the circle"
)

demo <- pp_demo_plasmid()
mapped <- pp_find_sites(
  read_genbank(system.file("extdata", "pDemo.gb", package = "plasmidplot")))

sheet <- function(items, path, draw, ncol = 4, w = 3200, h = 1680) {
  nrow <- ceiling(length(items) / ncol)
  ragg::agg_png(path, width = w, height = h, res = 150)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "#ffffff", col = NA))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
  for (i in seq_along(items)) {
    grid::pushViewport(grid::viewport(
      layout.pos.row = (i - 1L) %/% ncol + 1L,
      layout.pos.col = (i - 1L) %% ncol + 1L))
    draw(items[i])
    grid::popViewport()
  }
  grid::popViewport()
  dev.off()
  cat("wrote", path, "\n")
}

caption <- function(text, on_dark, note = NULL, y = 0.055) {
  grid::grid.text(text, y = y, gp = grid::gpar(
    fontsize = 15, fontface = "bold",
    col = if (on_dark) "#ffffff" else "#0b0b0b"))
  if (!is.null(note)) {
    grid::grid.text(note, y = y - 0.048, gp = grid::gpar(
      fontsize = 10, col = if (on_dark) "#9aa4b2" else "#6d6b66"))
  }
}

# --- 01 styles, circular: colors held fixed, only the look varies -----------
sheet(pp_style(), "man/figures/styles-circular.png", function(st) {
  on_dark <- st %in% DARK_STYLES
  plot(demo, bg = pp_canvas(st), newpage = FALSE,
       style = pp_style(st, show_title = FALSE,
                        palette = if (on_dark) "default_dark" else "default"))
  caption(st, on_dark, STYLE_NOTES[[st]])
})

# --- 02 styles, linear: same styles, same fixed colors, other layout --------
sheet(pp_style(), "man/figures/styles-linear.png", function(st) {
  on_dark <- st %in% DARK_STYLES
  plot(mapped, bg = pp_canvas(st), newpage = FALSE,
       style = pp_style(st, layout = "linear", show_title = FALSE,
                        palette = if (on_dark) "default_dark" else "default"))
  caption(st, on_dark, STYLE_NOTES[[st]], y = 0.08)
}, ncol = 2, w = 3000, h = 2000)

# --- 03 palettes: style held fixed, only the colors vary --------------------
sheet(pp_palette(), "man/figures/palettes.png", function(pal) {
  on_dark <- pal %in% DARK_PALETTES
  st <- if (on_dark) {
    pp_style("angular", palette = pal, show_title = FALSE,
             track_fill = "#383835", label_col = "#ffffff",
             tick_col = "#898781")
  } else {
    pp_style("angular", palette = pal, show_title = FALSE)
  }
  plot(demo, style = st, newpage = FALSE,
       bg = if (on_dark) "#1a1a19" else "#fcfcfb")
  caption(pal, on_dark)
})

# --- 04 palettes from other packages ----------------------------------------
externals <- list(
  "ggsci npg" = ggsci::pal_npg("nrc")(8),
  "ggsci aaas" = ggsci::pal_aaas()(8),
  "ggsci lancet" = ggsci::pal_lancet()(8),
  "ggsci jco" = ggsci::pal_jco()(8),
  "ggpubr npg" = ggpubr::get_palette("npg", 8),
  "brewer Dark2" = RColorBrewer::brewer.pal(8, "Dark2"),
  "brewer Set2" = RColorBrewer::brewer.pal(8, "Set2"),
  "viridis" = viridisLite::viridis(8)
)
sheet(names(externals), "man/figures/palettes-external.png", function(nm) {
  plot(demo, newpage = FALSE, bg = "#fcfcfb",
       style = pp_style("angular", palette = externals[[nm]],
                        show_title = FALSE))
  caption(nm, FALSE)
})

# --- 05 the two layouts, same plasmid, same style ---------------------------
ragg::agg_png("man/figures/layouts.png", width = 2400, height = 1500, res = 160)
grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "#ffffff", col = NA))
grid::pushViewport(grid::viewport(layout = grid::grid.layout(
  2, 1, heights = grid::unit(c(1.35, 1), "null"))))
grid::pushViewport(grid::viewport(layout.pos.row = 1))
plot(mapped, style = pp_style("angular", layout = "circular"), newpage = FALSE)
grid::popViewport()
grid::pushViewport(grid::viewport(layout.pos.row = 2))
plot(mapped, style = pp_style("angular", layout = "linear"), newpage = FALSE)
grid::popViewport()
grid::popViewport()
dev.off()
cat("wrote man/figures/layouts.png\n")

# --- 06 restriction sites ----------------------------------------------------
ragg::agg_png("man/figures/sites.png", width = 1500, height = 1500, res = 220)
plot(mapped, bg = "#ffffff")
dev.off()
cat("wrote man/figures/sites.png\n")

# --- 07 imported straight from GenBank --------------------------------------
ragg::agg_png("man/figures/import.png", width = 1400, height = 1400, res = 220)
plot(read_genbank(system.file("extdata", "pBR322.gb", package = "plasmidplot")),
     style = "snapgene", bg = "#ffffff")
dev.off()
cat("wrote man/figures/import.png\n")

cat("\n", length(list.files("man/figures")), "figures\n")
