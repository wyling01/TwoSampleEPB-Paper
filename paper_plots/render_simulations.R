#!/usr/bin/env Rscript
# Render the main-text simulation figures from replicate-level result CSVs.
# Base R/Cairo only. All artwork is vector, at the manuscript's 5.7-inch width.

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script) != 1L) stop("Run this file with Rscript.")
script <- gsub("~+~", " ", script, fixed = TRUE)
script_dir <- dirname(normalizePath(script))
result_root <- file.path(dirname(script_dir), "Simulation", "Simulation_result")
out <- file.path(script_dir, "figures")
if (!dir.exists(result_root)) stop("Simulation result folder not found: ", result_root)
is_macos <- identical(unname(Sys.info()["sysname"]), "Darwin")
if (!is_macos && !capabilities("cairo")) stop("An R installation with Cairo support is required.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

methods <- c("EV_test", "Welch", "B_F", "VREPB", "DVEPB")
source_names <- c(EV_test = "ev", Welch = "welch", B_F = "b_f", VREPB = "vrepb", DVEPB = "dvepb")
labels <- c("EV", "Welch", "BF", "VREPB", "DVEPB")
colors <- c(EV_test = "#BD5B64", Welch = "#A77B21", B_F = "#777777",
            VREPB = "#008875", DVEPB = "#315B8A")
shapes <- c(EV_test = 3, Welch = 17, B_F = 5, VREPB = 0, DVEPB = 16)
sizes <- c(EV_test = 0.72, Welch = 0.67, B_F = 0.65, VREPB = 0.82, DVEPB = 0.60)
line_types <- c(EV_test = "44", Welch = "13", B_F = "4313", VREPB = "solid", DVEPB = "solid")
ink <- "#252525"
axis_ink <- "#888888"

scenario_files <- function(setting, ka) {
  kb <- ka:(ka + 6L)
  file.path(result_root, tolower(setting), paste0("KA_", ka, "_KB_", kb, "_result.csv"))
}

read_scenario <- function(setting, ka) {
  kb_values <- ka:(ka + 6L)
  paths <- scenario_files(setting, ka)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing result file(s):\n", paste(missing, collapse = "\n"))
  result <- lapply(c(FDR = "fdp", Power = "power"), function(suffix) {
    rows <- lapply(seq_along(paths), function(i) {
      dat <- read.csv(paths[i], check.names = FALSE)
      required <- paste0(source_names, "_", suffix)
      if (!all(required %in% names(dat))) stop("Unexpected columns in: ", paths[i])
      values <- vapply(required, function(column) mean(dat[[column]]), numeric(1))
      if (any(!is.finite(values)) || any(values < 0) || any(values > 1))
        stop("Invalid simulation results in: ", paths[i])
      row <- data.frame(n1 = ka, n2 = kb_values[i])
      row[methods] <- as.list(unname(values))
      row
    })
    do.call(rbind, rows)
  })
  result
}

intro_available <- all(file.exists(scenario_files("Unequal", 3L)))
intro <- if (intro_available) read_scenario("Unequal", 3L) else NULL
overview <- list(read_scenario("Unequal", 5L), read_scenario("Equal", 3L),
                 read_scenario("Equal", 5L), read_scenario("Diffuse", 5L))

# Place plotting regions in inches so labels retain their physical size.
region <- function(x, y, width, height) {
  device <- par("din")
  par(fig = c(x, x + width, y, y + height) / rep(device, each = 2),
      mar = c(0, 0, 0, 0), new = TRUE)
}

panel <- function(dat, metric, x, y, width, height, y_axis = TRUE, x_axis = TRUE,
                  upper = if (metric == "FDR") 0.4 else 0.7, label_every = 2L) {
  region(x, y, width, height)
  ticks <- if (metric == "FDR") seq(0, upper, 0.1) else seq(0, 0.6, 0.2)
  stopifnot(all(is.finite(as.matrix(dat[methods]))),
            all(as.matrix(dat[methods]) >= 0), all(as.matrix(dat[methods]) <= 1))
  plot(NA_real_, NA_real_, xlim = range(dat$n2) + c(-0.24, 0.24),
       ylim = c(-0.035 * upper, 1.035 * upper), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  abline(h = ticks, col = "#ECECEC", lwd = 0.5)
  if (metric == "FDR") abline(h = 0.1, col = "#555555", lwd = 0.85, lty = "66")
  for (method in methods)
    lines(dat$n2, dat[[method]], col = colors[method], lwd = 1.35, lty = line_types[method])
  for (method in c("B_F", "Welch", "DVEPB", "VREPB", "EV_test"))
    points(dat$n2, dat[[method]], col = colors[method], pch = shapes[method],
           cex = sizes[method], lwd = 0.95)
  if (y_axis) {
    axis(2, at = ticks, labels = sprintf("%.1f", ticks), las = 1,
         lwd = 0, lwd.ticks = 0.55, col.ticks = axis_ink, tcl = -0.18,
         mgp = c(0, 0.34, 0), cex.axis = 0.94, col.axis = ink)
    mtext(metric, side = 2, line = 2.25, cex = 1.05, col = ink)
  }
  tick_labels <- ifelse((seq_along(dat$n2) - 1L) %% label_every == 0L, dat$n2, "")
  axis(1, at = dat$n2, labels = if (x_axis) tick_labels else FALSE,
       lwd = 0.55, col = axis_ink, col.ticks = axis_ink, tcl = -0.18,
       mgp = c(0, 0.34, 0), cex.axis = 0.94, col.axis = ink, gap.axis = -1)
}

canvas <- function() {
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, par("din")[1]), ylim = c(0, par("din")[2]), xaxs = "i", yaxs = "i")
}

shared_legend <- function(y) {
  canvas()
  segment <- 0.24
  gap <- 0.19
  widths <- segment + 0.07 + strwidth(labels)
  x <- (par("din")[1] - sum(widths) - gap * (length(methods) - 1L)) / 2
  for (j in seq_along(methods)) {
    segments(x, y, x + segment, y, col = colors[j], lwd = 1.35, lty = line_types[j])
    points(x + segment / 2, y, col = colors[j], pch = shapes[j], cex = sizes[j], lwd = 0.95)
    text(x + segment + 0.07, y, labels[j], adj = c(0, 0.5), col = ink)
    x <- x + widths[j] + gap
  }
}

draw_intro <- function() {
  panel(intro$FDR, "FDR", 0.50, 0.62, 2.20, 1.38, upper = 0.3, label_every = 1L)
  panel(intro$Power, "Power", 3.39, 0.62, 2.20, 1.38, label_every = 1L)
  canvas()
  text(c(1.60, 4.49), 0.32, expression(italic(K)[B]), cex = 1.05, col = ink)
  shared_legend(0.11)
}

draw_overview <- function() {
  x <- 0.50 + (0:3) * 1.31
  for (j in seq_along(overview)) {
    panel(overview[[j]]$FDR, "FDR", x[j], 2.12, 1.17, 1.22, y_axis = j == 1L, x_axis = FALSE)
    panel(overview[[j]]$Power, "Power", x[j], 0.62, 1.17, 1.22, y_axis = j == 1L)
  }
  canvas()
  text(x + 1.17 / 2, 3.81,
       c("(a) Unequal\nvariances", "(b) Equal\nvariances", "(c) Equal\nvariances", "(d) Diffuse\nvariance ratios"),
       cex = 1.02, col = ink)
  text(x + 1.17 / 2, 3.49,
       expression(italic(K)[A] == 5, italic(K)[A] == 3, italic(K)[A] == 5, italic(K)[A] == 5),
       cex = 1, col = ink)
  text(mean(x + 1.17 / 2), 0.32, expression(italic(K)[B]), cex = 1.05, col = ink)
  shared_legend(0.11)
}

render <- function(name, height, draw) {
  path <- file.path(out, paste0(name, ".pdf"))
  if (is_macos) {
    quartz(file = path, type = "pdf", width = 5.7, height = height,
           family = "Helvetica", pointsize = 9, bg = "white")
  } else {
    cairo_pdf(path, width = 5.7, height = height, family = "Helvetica",
              pointsize = 9, bg = "white", onefile = FALSE)
  }
  on.exit(dev.off())
  par(mar = rep(0, 4), family = "Helvetica", ps = 9, cex = 1,
      lend = "round", ljoin = "round", fg = ink, col = ink)
  plot.new()
  draw()
}

outputs <- character()
if (intro_available) {
  render("unequal_ka3", 2.12, draw_intro)
  outputs <- c(outputs, "unequal_ka3.pdf")
} else {
  message("Skipping unequal_ka3.pdf: unequal KA=3 result files are not present.")
}
render("overview", 4.02, draw_overview)
outputs <- c(outputs, "overview.pdf")
message("Wrote ", paste(outputs, collapse = " and "), " in ", normalizePath(out))
