library(plasmidplot)
cat(sprintf("%-10s %-8s %-8s %-8s %-6s %-6s %-6s\n",
            "style", "layout", "anchor", "backbone", "radius", "track", "arc"))
for (s in pp_style()) {
  st <- pp_style(s)
  cat(sprintf("%-10s %-8s %-8s %-8s %-6.2f %-6.3f %-6.3f\n",
              s, st$layout, st$anchor, st$backbone, st$radius, st$track, st$arc))
}
cat("\npalettes:", paste(pp_palette(), collapse = " "), "\n")
cat("style/palette name overlap:",
    length(intersect(pp_style(), pp_palette())), "\n")
