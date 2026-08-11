#' Example plasmid used by the demos
#'
#' Builds a small pBR322 map (positions approximate) without plotting it.
#'
#' @return A `plasmid` object.
#' @examples
#' pp_demo_plasmid()
#' @export
pp_demo_plasmid <- function() {
  p <- plasmid("pBR322", 4361)
  p <- pp_marker(p, 86, 1276, label = "TcR", arrow = "end")
  p <- pp_marker(p, 1915, 2106, label = "rop", arrow = "start")
  p <- pp_marker(p, 2535, 3122, label = "ori")
  p <- pp_marker(p, 3293, 4153, label = "AmpR", arrow = "start")
  p
}

#' Draw a demo plasmid map
#'
#' Plots [pp_demo_plasmid()] in the given style. Handy for previewing styles;
#' see [pp_style()] for the list of presets.
#'
#' The demo paints the style's own [pp_canvas()] so the dark presets are
#' legible on a white device; ordinary plotting leaves the background alone.
#'
#' @param style A preset name or [pp_style()] object; see [plot.plasmid()].
#' @return Invisibly, the demo `plasmid` object.
#' @examples
#' pp_demo("angular")
#' pp_demo("neon")
#' @export
pp_demo <- function(style = "angular") {
  plot(pp_demo_plasmid(), style = style, bg = pp_canvas(style))
}
