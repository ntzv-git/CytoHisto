#!/usr/bin/env Rscript
# CytoHisto (R command-line version): publication-ready flow cytometry histograms.
# Reads one or more .fcs files, draws one curve per file (overlaid) and saves the figure.
#
# Usage:
#   Rscript cytohisto.R [shared options] -f file1.fcs [file options] -f file2.fcs [...]
#
# Run "Rscript cytohisto.R -h" for the list of options.

VERSION <- "1.3.1"

suppressPackageStartupMessages(library(ggplot2))

help_text <- paste0("CytoHisto ", VERSION, " (R command-line version)\n", '
Usage:
  Rscript cytohisto.R [shared options] -f A.fcs [options for A] -f B.fcs [options for B] ...

Tip: drag a file from your file manager into the terminal to paste its path.

File options (placed AFTER their -f; placed BEFORE the first -f they become the
default for every file):
  -f, --file PATH            .fcs file (repeat for as many files as needed)
  -n, --name TEXT            curve name (legend or label; default: file name)
  -c, --color COLOUR         HTML/CSS colour: "#D55E00", "#d50", "rgb(213,94,0)", crimson, teal...
                             (quote it: otherwise the shell reads # as a comment)
  -w, --linewidth MM         line width (default 0.5)
  -t, --linetype TYPE        solid, dashed, dotted, dotdash, longdash (default solid)
  -s, --smooth N             number of classes merged (default 1); the histogram always has
                             1024 classes over the full scale, so neither smoothing nor the
                             output channels change the Y axis
  -g, --gate MIN,MAX         highlighted region in the curve colour + %, mean, CV in the console
  -p, --parameter NAME       parameter to plot (default: DNA if present, else the first non-time one)
  -i, --in-channels N        input range of this file (default: read from the file, $PnR);
                             files with different ranges are each converted with their own
      --label-pos X,Y        position of the curve label (with -L; default: above the highest
                             visible peak)

Shared options:
  -o, --output FILE          figure.pdf, .png, .tiff or .svg; no extension -> pdf + png + tiff
                             (default histogram)
  -x, --xlab TEXT            X axis title (default "Fluorescence intensity (channels)")
  -y, --ylab TEXT            Y axis title (default "Count")
  -C, --channels N           output channels = X axis resolution, same for all files (default 1024)
  -X, --xlim MIN,MAX         X zoom (e.g. 50,250), limited to 0 - output channels; the Y axis
                             fits the visible peaks
  -Y, --ylim MIN,MAX         Y zoom (e.g. 0,300); default: 0 to just above the highest peak
  -L, --labels               write each curve name above its peak (instead of a legend)
  -P, --pixels WxH           image size in pixels (default 2000x1000); text, lines and ticks
                             scale with the width
  -W, --width MM             PDF/SVG width in mm (default 170); the height follows the pixel ratio
      --font-size PT         all text sizes at once (default 9; values 7)
      --font-labels PT       curve names (labels or legend)     (default 9)
      --font-xtitle PT       X axis title                       (default 9)
      --font-ytitle PT       Y axis title                       (default 9)
      --font-xvalues PT      X axis values                      (default 7)
      --font-yvalues PT      Y axis values                      (default 7)
                             Sizes are for the reference layout and scale with the image width.
      --x-steps L[,G[,T]]    X steps: labels, grid + long ticks, short ticks; missing values are
                             derived (G = L/2, T = G/5). Default: auto from the visible range
      --y-steps L[,T]        Y steps: labels + grid, short ticks (T = L/5 if omitted). Default: auto
  -v, --version              show the version
  -h, --help                 show this help

Example:
  Rscript cytohisto.R -s 2 -L -X 50,250 -o comparison.png -P 3000 \\
      -f control.fcs -n "Diploid" -c grey40 -t dashed \\
      -f 313.fcs     -n "313" -c "#D55E00" -w 0.8 -g 85,105
')

# ---- HTML/CSS colours ----------------------------------------------------------
# R already knows #RRGGBB and most CSS names; add the missing ones, #RGB and rgb(r, g, b).
css_missing <- c(aqua = "#00FFFF", fuchsia = "#FF00FF", lime = "#00FF00", olive = "#808000",
                 silver = "#C0C0C0", teal = "#008080", crimson = "#DC143C", indigo = "#4B0082",
                 rebeccapurple = "#663399", darkcyan = "#008B8B", darkmagenta = "#8B008B")
css_color <- function(x) {
  y <- tolower(gsub(" ", "", x))
  if (grepl("^#[0-9a-f]{3}$", y)) y <- paste0("#", paste(rep(strsplit(substring(y, 2), "")[[1]], each = 2), collapse = ""))
  if (grepl("^rgb\\(", y)) {
    v <- as.numeric(strsplit(gsub("^rgb\\(|\\)$", "", y), ",")[[1]])
    y <- rgb(v[1], v[2], v[3], maxColorValue = 255)
  }
  if (y %in% names(css_missing)) y <- css_missing[[y]]
  ok <- tryCatch({ col2rgb(y); TRUE }, error = function(e) FALSE)
  if (!ok) stop("Unknown colour: ", x)
  y
}

# ---- Command-line arguments ------------------------------------------------------
parse_args <- function(args) {
  g <- list(output = "histogram", xlab = "Fluorescence intensity (channels)", ylab = "Count",
            channels = 1024, x_steps = NULL, y_steps = NULL,
            xlim = NULL, labels = FALSE, width = 170, height = 85, dpi = 600,
            fonts = c(labels = 9, xtitle = 9, ytitle = 9, xvalues = 7, yvalues = 7),
            pixels = c(2000, 1000))
  defaults <- list(name = NULL, color = NULL, linewidth = 0.5, linetype = "solid", smooth = 1,
                   gate = NULL, parameter = NULL, label_pos = NULL, in_range = NULL)
  files <- list()
  nums <- function(x, n, opt) {
    v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
    if (!length(v) %in% n || anyNA(v))
      stop(opt, " expects ", paste(n, collapse = " or "), " comma-separated number(s): ", x)
    v
  }
  i <- 1
  while (i <= length(args)) {
    opt <- args[i]
    if (opt %in% c("-h", "--help")) { cat(help_text); quit(status = 0) }
    if (opt %in% c("-v", "--version")) { cat("CytoHisto", VERSION, "\n"); quit(status = 0) }
    if (opt %in% c("-L", "--labels")) { g$labels <- TRUE; i <- i + 1; next }
    if (i == length(args)) stop("Option ", opt, " expects a value.")
    val <- args[i + 1]
    # file option: applies to the last -f, or to every file if no -f was given yet
    set <- function(key, v) {
      if (length(files)) files[[length(files)]][key] <<- list(v) else defaults[key] <<- list(v)
    }
    switch(opt,
      "-f" = , "--file"        = { files[[length(files) + 1]] <- c(list(path = val), defaults) },
      "-n" = , "--name"        = set("name", val),
      "-c" = , "--color"       = set("color", css_color(val)),
      "-w" = , "--linewidth"   = set("linewidth", nums(val, 1, opt)),
      "-t" = , "--linetype"    = set("linetype", val),
      "-s" = , "--smooth"      = set("smooth", nums(val, 1, opt)),
      "-g" = , "--gate"        = set("gate", sort(nums(val, 2, opt))),
      "-p" = , "--parameter"   = set("parameter", val),
      "--label-pos"            = set("label_pos", nums(val, 2, opt)),
      "-o" = , "--output"      = { g$output <- val },
      "-x" = , "--xlab"        = { g$xlab <- val },
      "-y" = , "--ylab"        = { g$ylab <- val },
      "-i" = , "--in-channels" = set("in_range", nums(val, 1, opt)),
      "-C" = , "--channels"    = { g$channels <- nums(val, 1, opt) },
      "-X" = , "--xlim"        = { g$xlim <- sort(nums(val, 2, opt)) },
      "-Y" = , "--ylim"        = { g$ylim <- sort(nums(val, 2, opt)) },
      "-W" = , "--width"       = { g$width <- nums(val, 1, opt) },
      "-P" = , "--pixels"      = { g$pixels <- suppressWarnings(as.numeric(strsplit(tolower(val), "x")[[1]]))
                                   if (length(g$pixels) != 2 || anyNA(g$pixels) || any(g$pixels < 100))
                                     stop("-P expects width x height in pixels, e.g. 2000x1000: ", val) },
      "--font-size"            = { v <- nums(val, 1, opt)
                                   g$fonts[] <- c(v, v, v, v * 7 / 9, v * 7 / 9) },
      "--font-labels"          = { g$fonts[["labels"]] <- nums(val, 1, opt) },
      "--font-xtitle"          = { g$fonts[["xtitle"]] <- nums(val, 1, opt) },
      "--font-ytitle"          = { g$fonts[["ytitle"]] <- nums(val, 1, opt) },
      "--font-xvalues"         = { g$fonts[["xvalues"]] <- nums(val, 1, opt) },
      "--font-yvalues"         = { g$fonts[["yvalues"]] <- nums(val, 1, opt) },
      "--x-steps"              = { g$x_steps <- nums(val, 1:3, opt) },
      "--y-steps"              = { g$y_steps <- nums(val, 1:2, opt) },
      stop("Unknown option: ", opt, " (see -h)"))
    i <- i + 2
  }
  if (!length(files)) { cat(help_text); stop("No .fcs file given (option -f).") }
  if (g$channels <= 1 || any(vapply(files, function(f) isTRUE(f$in_range <= 1), TRUE)))
    stop("Channel numbers must be positive.")
  list(g = g, files = files)
}

# ---- Minimal FCS reader ------------------------------------------------------------
# Handles parameters stored with different bit widths ($PnB), which flowCore cannot read
# ("object 'dat' not found").
read_fcs <- function(f) {
  con <- file(f, "rb"); on.exit(close(con))
  header <- rawToChar(readBin(con, "raw", 58))
  pos <- as.numeric(substring(header, c(11, 19, 27, 35), c(18, 26, 34, 42)))

  seek(con, pos[1])
  text <- readBin(con, "raw", pos[2] - pos[1] + 1)
  text <- rawToChar(text[text != as.raw(0)])
  delim <- substr(text, 1, 1)
  text <- gsub(paste0(delim, delim), "\001", text, fixed = TRUE, useBytes = TRUE)  # escaped delimiter
  fields <- gsub("\001", delim, strsplit(substring(text, 2), delim, fixed = TRUE, useBytes = TRUE)[[1]],
                 fixed = TRUE, useBytes = TRUE)
  fields <- fields[seq_len(length(fields) %/% 2 * 2)]
  kw <- setNames(trimws(fields[c(FALSE, TRUE)]), toupper(trimws(fields[c(TRUE, FALSE)])))

  if (is.na(pos[3]) || pos[3] == 0) pos[3:4] <- as.numeric(kw[c("$BEGINDATA", "$ENDDATA")])
  npar <- as.integer(kw[["$PAR"]]); ntot <- as.integer(kw[["$TOT"]])
  p      <- seq_len(npar)
  bits   <- as.integer(kw[paste0("$P", p, "B")])
  ranges <- as.numeric(kw[paste0("$P", p, "R")])
  names  <- kw[paste0("$P", p, "N")]
  endian <- if (kw[["$BYTEORD"]] %in% c("1,2,3,4", "1,2")) "little" else "big"
  type   <- kw[["$DATATYPE"]]

  seek(con, pos[3])
  if (type %in% c("F", "D")) {
    v <- readBin(con, "double", n = npar * ntot, size = if (type == "F") 4 else 8, endian = endian)
    m <- matrix(v, ncol = npar, byrow = TRUE)
  } else if (type == "I") {
    nbytes <- bits %/% 8
    raw <- readBin(con, "raw", n = sum(nbytes) * ntot)
    raw <- matrix(as.integer(raw), ncol = sum(nbytes), byrow = TRUE)
    end <- cumsum(nbytes); start <- end - nbytes + 1
    m <- sapply(p, function(i) {
      cols <- start[i]:end[i]
      if (endian == "big") cols <- rev(cols)           # least significant byte first
      v <- as.vector(raw[, cols, drop = FALSE] %*% (256^(seq_along(cols) - 1)))
      used <- if (ranges[i] > 1) ceiling(log2(ranges[i])) else bits[i]   # bits actually used ($PnR)
      if (used < bits[i]) v <- v %% 2^used
      v
    })
    m <- matrix(m, ncol = npar)
  } else stop("$DATATYPE '", type, "' is not supported")

  list(data = setNames(as.data.frame(m), names), ranges = setNames(ranges, names))
}

default_parameter <- function(names) {
  if (any(toupper(names) == "DNA")) return(names[toupper(names) == "DNA"][1])
  rest <- names[!grepl("TIME", toupper(names))]
  if (length(rest)) rest[1] else names[1]
}

# Range used to look for a peak: ignore 5 % at each end of the full scale
# (debris / threshold at the low end, saturation channel at the high end)
inside_edges <- function(x) x > 0.05 * out & x < 0.95 * out

# Round step (1, 2, 2.5 or 5 x 10^n) giving at most n_max intervals over a span
nice_step <- function(span, n_max) {
  k <- 10^floor(log10(span / n_max))
  cand <- k * c(1, 2, 2.5, 5, 10)
  cand[which(span / cand <= n_max)[1]]
}

# ---- Read files and count events ---------------------------------------------------
a <- parse_args(commandArgs(trailingOnly = TRUE))
g <- a$g
n <- length(a$files)
palette <- c("black", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
out <- g$channels
xlim <- if (!is.null(g$xlim)) g$xlim else c(0, out)
if (xlim[1] < 0 || xlim[2] > out) {                 # never beyond the channel range
  message("X zoom limited to 0-", out)
  xlim <- c(max(xlim[1], 0), min(xlim[2], out))
}
if (xlim[2] <= xlim[1]) stop("Empty X zoom after limiting it to 0-", out)
BASE_BINS <- 1024   # histogram classes over the full scale, whatever the output channels
g$height <- g$width * g$pixels[2] / g$pixels[1]   # proportions given by the pixel size
k <- g$width / 170   # scale factor: everything is designed for 170 mm, then proportional

curves <- list(); stats <- list()
for (i in seq_len(n)) {
  f <- a$files[[i]]
  if (!file.exists(f$path)) stop("File not found: ", f$path)
  fcs <- read_fcs(f$path)
  d <- fcs$data
  if (is.null(f$parameter)) f$parameter <- default_parameter(names(d))
  if (!f$parameter %in% names(d))
    stop("Parameter '", f$parameter, "' not found in ", basename(f$path), ". Available: ",
         paste(names(d), collapse = ", "))
  # each file is converted with its own input range (set with -i, otherwise read from the file)
  in_range <- if (!is.null(f$in_range)) f$in_range else fcs$ranges[[f$parameter]]
  v <- d[[f$parameter]] * out / in_range                     # input range -> output channels
  if (is.null(f$name)) f$name <- sub("\\.fcs$", "", basename(f$path), ignore.case = TRUE)
  if (is.null(f$color)) f$color <- palette[(i - 1) %% length(palette) + 1]

  # histogram with a fixed resolution (BASE_BINS classes over the full scale): the output
  # channels only set the X unit, so the Y axis and the curve shape do not depend on them
  base <- out / BASE_BINS                               # class width in output channels
  w <- base * f$smooth                                  # smoothing merges classes
  nb <- ceiling(BASE_BINS / f$smooth)
  cnt <- tabulate(pmin(pmax(floor(v / w) + 1, 1), nb), nb)
  x <- (seq_len(nb) - 1) * w + (w - base) / 2           # at 1024 channels: x = channel number
  y <- cnt / f$smooth                                   # events per class
  dc <- data.frame(x = x, y = y)
  curves[[i]] <- dc[dc$x >= xlim[1] & dc$x <= xlim[2], ]   # only the visible range is drawn

  ok <- inside_edges(x)
  s <- data.frame(file = basename(f$path), label = f$name, events = length(v),
                  peak = if (any(ok)) round(x[ok][which.max(y[ok])], 1) else NA,
                  gate = NA, pct_in_gate = NA, mean = NA, cv_pct = NA)
  if (!is.null(f$gate)) {
    inside <- v[v >= f$gate[1] & v <= f$gate[2]]
    s$gate <- paste(f$gate, collapse = "-")
    s$pct_in_gate <- round(100 * length(inside) / length(v), 1)
    s$mean <- round(mean(inside), 1)
    s$cv_pct <- round(100 * sd(inside) / mean(inside), 2)
  }
  stats[[i]] <- s
  a$files[[i]] <- f
}
labels <- make.unique(vapply(a$files, `[[`, "", "name"))
for (i in seq_len(n)) curves[[i]]$name <- rep(labels[i], nrow(curves[[i]]))

# ---- Axes ----------------------------------------------------------------------------
ymax <- max(vapply(curves, function(d) max(c(d$y, 1)), 0))   # highest visible peak, all files
y_steps <- if (is.null(g$y_steps)) {
  p <- nice_step(ymax, 8)
  c(major = p, minor = p / if (p / 10^floor(log10(p)) == 2) 4 else 5)
} else c(major = g$y_steps[1], minor = if (length(g$y_steps) > 1) g$y_steps[2] else g$y_steps[1] / 5)
if (ymax / y_steps[["minor"]] > 400 || ymax / y_steps[["major"]] > 60) {   # would draw thousands of ticks
  message("Y steps too small for a peak of ", round(ymax), ": automatic steps used")
  p <- nice_step(ymax, 8)
  y_steps <- c(major = p, minor = p / if (p / 10^floor(log10(p)) == 2) 4 else 5)
}
if (!is.null(g$ylim)) {                              # Y zoom: fixed range, steps from its span
  bottom <- g$ylim[1]; top <- g$ylim[2]
  if (is.null(g$y_steps) || diff(g$ylim) / y_steps[["minor"]] > 400) {
    p <- nice_step(diff(g$ylim), 8)
    y_steps <- c(major = p, minor = p / if (p / 10^floor(log10(p)) == 2) 4 else 5)
  }
} else {
  bottom <- 0
  top <- ceiling(ymax * (if (g$labels) 1.12 else 1.05) / y_steps[["minor"]]) * y_steps[["minor"]]
}
y_major <- seq(ceiling(bottom / y_steps[["major"]] - 1e-9) * y_steps[["major"]], top, by = y_steps[["major"]])
y_minor <- setdiff(round(seq(ceiling(bottom / y_steps[["minor"]] - 1e-9) * y_steps[["minor"]], top,
                             by = y_steps[["minor"]]), 9), round(y_major, 9))
for (i in seq_len(n)) curves[[i]]$y <- pmin(pmax(curves[[i]]$y, bottom), top)   # stay in the frame
message(sprintf("Y axis: highest peak %g -> steps %g / %g", ymax, y_steps[["major"]], y_steps[["minor"]]))

# labels every L (about 6 over the visible range), grid + long ticks every L/2,
# short ticks every L/10. E.g. 0-1024 -> 200 / 100 / 20 ; zoom 50-250 -> 50 / 25 / 5
x_steps <- if (is.null(g$x_steps)) { p <- nice_step(diff(xlim), 6); c(p, p / 2, p / 10) } else g$x_steps
if (length(x_steps) == 1) x_steps <- c(x_steps, x_steps / 2)
if (length(x_steps) == 2) x_steps <- c(x_steps, x_steps[2] / 5)
x_steps <- setNames(x_steps, c("labels", "grid", "ticks"))
if (diff(xlim) / x_steps[["ticks"]] > 400 || diff(xlim) / x_steps[["labels"]] > 60) {
  message("X steps too small for ", diff(xlim), " channels: automatic steps used")
  p <- nice_step(diff(xlim), 6)
  x_steps <- c(labels = p, grid = p / 2, ticks = p / 10)
}
on_axis <- function(step) seq(ceiling(xlim[1] / step - 1e-9) * step, xlim[2], by = step)
x_grid <- on_axis(x_steps[["grid"]])
x_ticks <- setdiff(round(on_axis(x_steps[["ticks"]]), 9), round(x_grid, 9))
message(sprintf("X axis: labels %g, grid %g, ticks %g", x_steps[["labels"]], x_steps[["grid"]],
                x_steps[["ticks"]]))

# ---- Figure --------------------------------------------------------------------------
pg <- ggplot()     # grid first, so that it stays under the curves
grid_y <- y_major[y_major > bottom]
grid_x <- x_grid[x_grid > xlim[1] & x_grid <= xlim[2]]
if (length(grid_y)) pg <- pg + geom_hline(yintercept = grid_y, colour = "grey85", linewidth = 0.25 * k)
if (length(grid_x)) pg <- pg + geom_vline(xintercept = grid_x, colour = "grey85", linewidth = 0.25 * k)

# gates: region highlighted in the curve colour, under the curves
for (f in a$files) {
  if (is.null(f$gate) || f$gate[2] <= xlim[1] || f$gate[1] >= xlim[2]) next
  pg <- pg + annotate("rect", xmin = max(f$gate[1], xlim[1]), xmax = min(f$gate[2], xlim[2]),
                      ymin = bottom, ymax = top, fill = f$color, alpha = 0.15)
}

for (i in seq_len(n))
  pg <- pg + geom_line(data = curves[[i]], aes(x = x, y = y, colour = name),
                       linewidth = a$files[[i]]$linewidth * k, linetype = a$files[[i]]$linetype)

pg <- pg +
  scale_colour_manual(values = setNames(vapply(a$files, `[[`, "", "color"), labels),
                      breaks = labels, name = NULL) +
  guides(colour = if (n > 1 && !g$labels) guide_legend(override.aes = list(
    linewidth = vapply(a$files, `[[`, 0, "linewidth") * k,
    linetype = vapply(a$files, `[[`, "", "linetype"))) else "none")

# labels: each curve name above its highest visible peak (away from the edges)
if (g$labels) {
  for (i in seq_len(n)) {
    f <- a$files[[i]]; dc <- curves[[i]]
    if (!is.null(f$label_pos)) { px <- f$label_pos[1]; py <- f$label_pos[2] } else {
      zone <- inside_edges(dc$x)
      if (!any(zone)) next
      j <- which(zone)[which.max(dc$y[zone])]
      px <- dc$x[j]; py <- min(dc$y[j] + 0.02 * (top - bottom), top - 0.08 * (top - bottom))
    }
    pg <- pg + annotate("text", x = px, y = py, label = labels[i], colour = f$color,
                        vjust = 0, size = g$fonts[["labels"]] * k / .pt)
  }
}

# short Y ticks, drawn by hand so that it works with any ggplot2 version
tick_y <- grid::segmentsGrob(x0 = unit(0, "npc"), x1 = unit(-1.2 * k, "mm"),
                             y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
                             gp = grid::gpar(col = "grey20", lwd = 9 * k / 22 * .pt))
for (y in y_minor) pg <- pg + annotation_custom(tick_y, xmin = -Inf, xmax = -Inf, ymin = y, ymax = y)

pg <- pg + scale_x_continuous(breaks = x_grid, labels = function(x)
  ifelse(abs(x / x_steps[["labels"]] - round(x / x_steps[["labels"]])) < 1e-9, x, ""))
tick_x <- grid::segmentsGrob(x0 = unit(0.5, "npc"), x1 = unit(0.5, "npc"),
                             y0 = unit(0, "npc"), y1 = unit(-0.8 * k, "mm"),
                             gp = grid::gpar(col = "grey20", lwd = 9 * k / 22 * .pt))
for (xm in x_ticks) pg <- pg + annotation_custom(tick_x, xmin = xm, xmax = xm, ymin = -Inf, ymax = -Inf)

pg <- pg +
  scale_y_continuous(breaks = y_major) +
  coord_cartesian(xlim = xlim, ylim = c(bottom, top), expand = FALSE, clip = "off") +
  labs(x = g$xlab, y = g$ylab) +
  theme_classic(base_size = 9 * k) +   # line widths follow the image width, not the fonts
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = g$fonts[["xvalues"]] * k),
        axis.text.y = element_text(size = g$fonts[["yvalues"]] * k),
        axis.title.x = element_text(size = g$fonts[["xtitle"]] * k),
        axis.title.y = element_text(size = g$fonts[["ytitle"]] * k),
        legend.text = element_text(size = g$fonts[["labels"]] * k),
        axis.ticks.length.x = unit(1.5 * k, "mm"),
        axis.ticks.length.y = unit(2.5 * k, "mm"),
        legend.justification = c(1, 1),
        legend.key.width = unit(8 * k, "mm"),
        legend.background = element_blank(), legend.key = element_blank(),
        plot.margin = margin(6 * k, 10 * k, 4 * k, 4 * k))
pg <- pg + if (packageVersion("ggplot2") < "3.5.0") theme(legend.position = c(0.99, 0.99)) else
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99))

# ---- Save ----------------------------------------------------------------------------
dpi <- g$pixels[1] / (g$width / 25.4)
ext <- tolower(tools::file_ext(g$output))
targets <- if (ext %in% c("pdf", "png", "tiff", "tif", "svg")) g$output else
  paste0(g$output, c(".pdf", ".png", ".tiff"))
for (target in targets) {
  e <- tolower(tools::file_ext(target))
  if (e == "png") ggsave(target, pg, width = g$width, height = g$height, units = "mm", dpi = dpi)
  else if (e %in% c("tiff", "tif"))
    ggsave(target, pg, width = g$width, height = g$height, units = "mm", dpi = dpi, compression = "lzw")
  else ggsave(target, pg, width = g$width, height = g$height, units = "mm")
}

print(do.call(rbind, stats), row.names = FALSE)
message("Figure saved: ", paste(targets, collapse = ", "),
        if (any(grepl("\\.(png|tiff?)$", targets)))
          sprintf(" (%d x %d px)", round(g$width / 25.4 * dpi), round(g$height / 25.4 * dpi)) else "")
