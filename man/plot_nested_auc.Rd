\name{plot_nested_auc}
\alias{plot_nested_auc}
\alias{plot_calibration}
\alias{plot_roc_curves}
\alias{plot_trace_weights}
\alias{plot_mediation_path}
\alias{open_nature_device}
\alias{with_nature_device}
\alias{nature_palette}
\title{Publication-quality (Nature-style) figure export}
\usage{
open_nature_device(path, width_mm = 89, height_mm = 89, dpi = 600,
  format = NULL, font_family = "sans")

with_nature_device(path, expr, width_mm = 89, height_mm = 89, dpi = 600,
  format = NULL, font_family = "sans")

plot_nested_auc(nested_comparison, auc_ci = NULL, path, ...)

plot_calibration(validation, path, model_label = "Integrated model", ...)

plot_roc_curves(roc_list, path, ...)

plot_trace_weights(trace_fit, path, top_n = NULL, ...)

plot_mediation_path(sem_result, path, ...)
}
\arguments{
\item{path}{Output file path; extension (\code{.png}, \code{.tiff}/\code{.tif},
\code{.pdf}) determines format unless \code{format} is supplied.}
\item{width_mm, height_mm}{Physical figure size in millimeters. Defaults
(89 x 89 mm) match a Nature single-column figure; use \code{width_mm = 183}
for a double-column figure.}
\item{dpi}{Resolution in dots per inch. Default 600 (exceeds Nature's
300 dpi minimum for combination art).}
\item{format}{One of \code{"png"}, \code{"tiff"}, \code{"pdf"}.}
\item{font_family}{Graphics device font family. Default \code{"sans"}
(Helvetica/Arial), per Nature typography guidelines.}
\item{expr}{An expression (typically a \code{\{ ... \}} block) that draws
a figure using base graphics.}
\item{nested_comparison}{Output of \code{compare_nested_models()}.}
\item{auc_ci}{Optional data frame with columns \code{model}, \code{ci_low},
\code{ci_high} for AUC confidence-interval whiskers.}
\item{validation}{Output of \code{validate_model()}.}
\item{model_label}{Character; names the model in the legend (used
directly, not a generic "Model" label).}
\item{roc_list}{Named list of \code{pROC::roc} objects; names are used
directly as legend labels -- name them descriptively (e.g.
\code{"Asthma PRS only"}), not generically (\code{"Model 1"}).}
\item{trace_fit}{Output of \code{compute_trace_score()}.}
\item{top_n}{Number of top transcription factors (by \code{|weight|}) to display.}
\item{sem_result}{Output of \code{fit_mediation_sem()}.}
\item{...}{Passed to \code{with_nature_device()} (\code{width_mm},
\code{height_mm}, \code{dpi}, \code{format}).}
}
\value{
Invisibly, the output file path.
}
\description{
A family of functions that render every major TRACE-Asthma analysis
output (nested-model discrimination, calibration, ROC curves, TRACE
transcription-factor weight profile, mediation path diagram) as
publication-ready figures: 600 dpi by default, sans-serif typography,
colorblind-safe categorical palette (\code{nature_palette}), and direct,
non-generic axis/legend labeling (gene symbols and named predictors, not
"TF1"/"Model A") throughout. \code{open_nature_device()} and
\code{with_nature_device()} are the underlying device helpers, exported
for building custom multi-panel figures around individual \code{plot_*()}
"panel" functions with base graphics (\code{par(mfrow = ...)}).
}
\examples{
\dontrun{
plot_trace_weights(trace_fit, path = "figure3_trace_weights.tiff",
                    width_mm = 89, height_mm = 89, dpi = 600)
}
}
