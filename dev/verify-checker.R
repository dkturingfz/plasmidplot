setwd("i:/plasmidplot")
pkgload::load_all(quiet = TRUE)
for (nm in pp_palette()) {
  mode <- if (grepl("_dark$", nm)) "dark" else "light"
  surf <- if (nm == "neon_dark") "#0b0f19" else NULL
  cat("=====", nm, "(", mode, ")\n")
  print(pp_check_palette(nm, mode = mode, surface = surf))
}
cat("===== a deliberately bad one (all browns)\n")
print(pp_check_palette(c("#8c5a3c", "#c9a227", "#6b8e5a", "#b5651d",
                         "#7d6b5d", "#a3623c", "#556b4f", "#8a7250")))
