#' Open a publication-quality graphics device
#'
#' Enforces the conventions expected by Nature-family journals: sans-serif
#' font family, high resolution (default 600 dpi, exceeding Nature's
#' minimum of 300 dpi for combination art), RGB color space, and physical
#' sizing in inches/mm rather than pixels so figures scale correctly at
#' print size. All `plot_*()` functions in this file call this internally;
#' call it directly only if you are assembling a custom multi-panel figure
#' with base graphics (e.g. `par(mfrow = c(1,2))`) around calls to the
#' individual `plot_*()` "panel" functions.
#'
#' @param path Output file path. Extension determines format if `format`
#'   is not supplied (`.png`, `.tiff`/`.tif`, `.pdf`).
#' @param width_mm,height_mm Physical figure size in millimeters. Defaults
#'   (89 x 89 mm) match Nature's single-column width; use `183` mm for
#'   `width_mm` for a double-column (full-page-width) figure.
#' @param dpi Resolution in dots per inch for raster formats. Default 600.
#' @param format One of `"png"`, `"tiff"`, `"pdf"`. Default inferred from
#'   `path`'s extension. TIFF (LZW-compressed) is offered because it
#'   remains the most widely-accepted raster submission format across
#'   Nature-family journals; PDF is offered for vector figures (e.g. the
#'   mediation path diagram) that should not be rasterized at all.
#' @param font_family Font family passed to the graphics device. Default
#'   `"sans"` (maps to Helvetica/Arial on most systems), matching Nature's
#'   typography guidelines.
#' @return Invisibly, the resolved output `path`. The device is left open;
#'   call [grDevices::dev.off()] after drawing (or use [with_nature_device()]
#'   to handle this automatically).
#' @export
open_nature_device <- function(path, width_mm = 89, height_mm = 89, dpi = 600,
                                format = NULL, font_family = "sans") {
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(path))
    format <- switch(ext, png = "png", tif = "tiff", tiff = "tiff", pdf = "pdf",
                      stop("Could not infer format from path extension '", ext,
                           "'; supply `format` explicitly.", call. = FALSE))
  }
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  switch(format,
    png = grDevices::png(path, width = width_in, height = height_in, units = "in",
                          res = dpi, family = font_family, type = "cairo",
                          bg = "white"),
    tiff = grDevices::tiff(path, width = width_in, height = height_in, units = "in",
                            res = dpi, family = font_family, compression = "lzw",
                            bg = "white"),
    pdf = grDevices::pdf(path, width = width_in, height = height_in,
                          family = font_family, useDingbats = FALSE)
  )
  invisible(path)
}

#' Run plotting code inside a publication-quality device, closing it safely
#'
#' @param path,width_mm,height_mm,dpi,format,font_family As in
#'   [open_nature_device()].
#' @param expr An expression (typically a `{ ... }` block) that draws the
#'   figure using base graphics.
#' @return Invisibly, the output `path`.
#' @export
with_nature_device <- function(path, expr, width_mm = 89, height_mm = 89,
                                dpi = 600, format = NULL, font_family = "sans") {
  open_nature_device(path, width_mm, height_mm, dpi, format, font_family)
  on.exit(grDevices::dev.off())
  force(expr)
  invisible(path)
}

#' @keywords internal
.nature_par <- function(mar = c(4, 4.5, 2, 1.5)) {
  graphics::par(mar = mar, mgp = c(2.4, 0.7, 0), tcl = -0.3, las = 1,
                cex.axis = 0.85, cex.lab = 0.95, font.lab = 1, lend = 1,
                family = "sans", xpd = FALSE)
}

#' Nature-style color palette (colorblind-safe, categorical)
#' @export
nature_palette <- c(
  blue   = "#0072B2", vermillion = "#D55E00", teal = "#009E73",
  amber  = "#E69F00", purple = "#CC79A7", gray = "#666666"
)

#' Figure: nested-model discrimination (AUC by model, with 95% CI)
#'
#' Renders the nested-model AUC comparison (M1-M5 or whichever subset was
#' fit, see [compare_nested_models()]) as a Nature-style panel: filled
#' points with 95% CI whiskers, model labels naming the *added predictor*
#' at each step (not "M1"/"M2") so the figure is interpretable without the
#' Methods text, and horizontal reference gridlines only (no chart-junk
#' vertical gridlines).
#'
#' @param nested_comparison Output of [compare_nested_models()], or its
#'   `$comparison` data frame directly.
#' @param auc_ci Optional data frame with columns `model`, `ci_low`,
#'   `ci_high` (e.g. from bootstrap or DeLong CIs on each model's AUC via
#'   `pROC::ci.auc()`); if omitted, points are drawn without whiskers.
#' @param path Output file path (png/tiff/pdf).
#' @param ... Passed to [with_nature_device()] (`width_mm`, `height_mm`,
#'   `dpi`, `format`).
#' @return Invisibly, the output path.
#' @export
plot_nested_auc <- function(nested_comparison, auc_ci = NULL, path, ...) {
  cmp <- if (inherits(nested_comparison, "traceAsthma_nested")) nested_comparison$comparison else nested_comparison
  label_map <- c(
    M1 = "Asthma PRS", M2 = "Multi-trait PRS", M3 = "+ eQTM score",
    M4 = "+ TRACE score", M5 = "+ clinical covariates"
  )
  labels <- ifelse(cmp$model %in% names(label_map), label_map[cmp$model], cmp$predictors)

  with_nature_device(path, {
    .nature_par(mar = c(6, 4.5, 1.5, 1.5))
    n <- nrow(cmp)
    graphics::plot(seq_len(n), cmp$auc, type = "n", xaxt = "n",
                   xlab = "", ylab = "Area under the ROC curve",
                   xlim = c(0.5, n + 0.5),
                   ylim = c(max(0, min(cmp$auc, na.rm = TRUE) - 0.08),
                            min(1, max(cmp$auc, na.rm = TRUE) + 0.08)),
                   bty = "l")
    graphics::grid(nx = NA, ny = NULL, col = "grey88", lty = 1)
    if (!is.null(auc_ci)) {
      ord <- match(cmp$model, auc_ci$model)
      graphics::segments(seq_len(n), auc_ci$ci_low[ord], seq_len(n), auc_ci$ci_high[ord],
                          col = nature_palette["gray"], lwd = 1.3)
    }
    graphics::points(seq_len(n), cmp$auc, pch = 21, bg = nature_palette["blue"],
                      col = "black", cex = 1.4, lwd = 1)
    graphics::axis(1, at = seq_len(n), labels = FALSE)
    graphics::text(seq_len(n), graphics::par("usr")[3] - 0.02 * diff(graphics::par("usr")[3:4]),
                   labels = labels, srt = 35, adj = c(1, 1), xpd = TRUE, cex = 0.8)
  }, ...)
}

#' Figure: calibration plot (observed vs. predicted asthma probability)
#'
#' @param validation Output of [validate_model()] (uses `$calibration$plot_data`).
#' @param path Output file path.
#' @param model_label Character, used in the legend to name the model
#'   being calibrated (e.g. `"MPRS + eQTM + TRACE"`), not a generic "Model".
#' @param ... Passed to [with_nature_device()].
#' @export
plot_calibration <- function(validation, path, model_label = "Integrated model", ...) {
  pd <- validation$calibration$plot_data
  with_nature_device(path, {
    .nature_par()
    graphics::plot(pd$mean_predicted, pd$mean_observed, type = "n",
                   xlim = c(0, 1), ylim = c(0, 1), bty = "l",
                   xlab = "Mean predicted probability of asthma",
                   ylab = "Observed proportion with asthma")
    graphics::abline(0, 1, col = "grey60", lty = 2, lwd = 1.2)
    graphics::points(pd$mean_predicted, pd$mean_observed, pch = 21,
                      bg = nature_palette["vermillion"], col = "black", cex = 1.3)
    graphics::lines(pd$mean_predicted, pd$mean_observed, col = nature_palette["vermillion"], lwd = 1.2)
    graphics::legend("topleft", legend = c(model_label, "Perfect calibration"),
                      col = c(nature_palette["vermillion"], "grey60"),
                      pch = c(21, NA), pt.bg = c(nature_palette["vermillion"], NA),
                      lty = c(1, 2), bty = "n", cex = 0.8, seg.len = 1.6)
  }, ...)
}

#' Figure: ROC curve(s) for one or more nested models overlaid
#'
#' @param roc_list Named list of `pROC::roc` objects (or of `list(observed=,
#'   predicted=)` pairs), one per model to overlay; names are used
#'   directly as legend labels, so name them descriptively (e.g.
#'   `list("Asthma PRS only" = roc1, "+ eQTM + TRACE" = roc2)`), not
#'   generically ("Model 1").
#' @param path Output file path.
#' @param ... Passed to [with_nature_device()].
#' @export
plot_roc_curves <- function(roc_list, path, ...) {
  requirePkg("pROC", "ROC curve plotting")
  cols <- unname(nature_palette)[seq_along(roc_list)]

  with_nature_device(path, {
    .nature_par()
    graphics::plot(c(0, 1), c(0, 1), type = "n", bty = "l",
                   xlab = "1 - Specificity", ylab = "Sensitivity")
    graphics::abline(0, 1, col = "grey75", lty = 3)
    aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
    for (i in seq_along(roc_list)) {
      r <- roc_list[[i]]
      graphics::lines(1 - r$specificities[order(r$sensitivities)],
                       r$sensitivities[order(r$sensitivities)],
                       col = cols[i], lwd = 1.8)
    }
    legend_labels <- sprintf("%s (AUC = %.2f)", names(roc_list), aucs)
    graphics::legend("bottomright", legend = legend_labels, col = cols, lwd = 1.8,
                      bty = "n", cex = 0.8, seg.len = 1.6)
  }, ...)
}

#' Figure: TRACE transcription-factor weight profile
#'
#' Horizontal bar chart of the elastic-net weights contributing to a
#' fitted TRACE score (see [summarize_trace_weights()]), colored by
#' direction of association, with TF gene symbols as direct axis labels
#' (never "TF1", "TF2").
#'
#' @param trace_fit Output of [compute_trace_score()].
#' @param path Output file path.
#' @param top_n Number of top TFs (by |weight|) to display. Default all.
#' @param ... Passed to [with_nature_device()].
#' @export
plot_trace_weights <- function(trace_fit, path, top_n = NULL, ...) {
  df <- summarize_trace_weights(trace_fit, top_n = top_n)
  df <- df[order(df$weight), ]  # ascending for horizontal bar plot top-to-bottom

  with_nature_device(path, {
    .nature_par(mar = c(4, 7, 1.5, 1.5))
    cols <- ifelse(df$weight > 0, nature_palette["vermillion"], nature_palette["blue"])
    bp <- graphics::barplot(df$weight, horiz = TRUE, col = cols, border = NA,
                             names.arg = df$tf, las = 1, cex.names = 0.85,
                             xlab = "Elastic-net weight (contribution to TRACE score)")
    graphics::abline(v = 0, col = "black", lwd = 0.8)
    graphics::legend("bottomright", legend = c("Activating (raises risk)", "Repressive (lowers risk)"),
                      fill = c(nature_palette["vermillion"], nature_palette["blue"]),
                      border = NA, bty = "n", cex = 0.75)
  }, ...)
}

#' Figure: mediation pathway diagram with path coefficients
#'
#' Draws the MPRS -> eQTM -> TRACE -> Asthma structural path diagram with
#' each arrow directly labeled by its estimated path coefficient and
#' significance, plus the total indirect (mediated) effect summarized in
#' a caption line -- avoiding a generic unlabeled box-and-arrow figure.
#' Recommended as a vector PDF (`format = "pdf"`) since it contains no
#' raster elements.
#'
#' @param sem_result Output of [fit_mediation_sem()].
#' @param path Output file path.
#' @param ... Passed to [with_nature_device()].
#' @export
plot_mediation_path <- function(sem_result, path, ...) {
  pe <- lavaan::parameterEstimates(sem_result$fit)
  get_est <- function(lbl) pe[pe$label == lbl, c("est", "pvalue")]
  a1 <- get_est("a1"); b1 <- get_est("b1"); c3 <- get_est("c3"); c1 <- get_est("c1")
  sig <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""

  nodes <- data.frame(
    label = c("MPRS", "eQTM", "TRACE", "Asthma"),
    x = c(0, 1, 2, 3), y = c(0, 0, 0, 0)
  )

  with_nature_device(path, {
    .nature_par(mar = c(3, 1, 1, 1))
    graphics::plot(NA, xlim = c(-0.3, 3.3), ylim = c(-0.6, 0.9), axes = FALSE,
                   xlab = "", ylab = "")
    for (i in seq_len(nrow(nodes) - 1)) {
      graphics::arrows(nodes$x[i] + 0.28, nodes$y[i], nodes$x[i + 1] - 0.28, nodes$y[i + 1],
                        length = 0.08, lwd = 1.6, col = "black")
    }
    graphics::text(0.5, 0.18, sprintf("a = %.2f%s", a1$est, sig(a1$pvalue)), cex = 0.8)
    graphics::text(1.5, 0.18, sprintf("b = %.2f%s", b1$est, sig(b1$pvalue)), cex = 0.8)
    graphics::text(2.5, 0.18, sprintf("c3 = %.2f%s", c3$est, sig(c3$pvalue)), cex = 0.8)
    graphics::arrows(0, -0.35, 3, -0.35, length = 0.08, lwd = 1, lty = 2, col = "grey40")
    graphics::text(1.5, -0.48, sprintf("Direct effect c1 = %.2f%s", c1$est, sig(c1$pvalue)),
                   cex = 0.75, col = "grey40")
    for (i in seq_len(nrow(nodes))) {
      graphics::rect(nodes$x[i] - 0.28, -0.15, nodes$x[i] + 0.28, 0.15, col = "white", border = "black")
      graphics::text(nodes$x[i], 0, nodes$label[i], cex = 0.9, font = 2)
    }
    graphics::mtext(sprintf("%.1f%% of the total MPRS effect on asthma is mediated through eQTM \u2192 TRACE",
                             100 * sem_result$proportion_mediated),
                     side = 1, line = 1.5, cex = 0.75)
  }, ...)
}
