
#' Plot MCMC diagnostics.
#' 
#' Plot MCMC diagnostics (traceplot and autocorrelation) of the log-posterior probability
#' for a \code{ouija_fit} object.
#' 
#' Further assessment of convergence can be done using \code{rstan} functions on the
#' underlying STAN object (accessed through \code{oui$fit}).
#' 
#' @param oui A \code{ouija_fit} object
#' @param arrange How to arrange the plots. If "vertical", traceplot and autocorrelation are 
#' arranged in one column, while if "horizontal" traceplot and autocorrelation are arranged
#' in one row.
#' @export
#' @importFrom cowplot plot_grid
#' 
#' @return A \code{ggplot2} object
#' @examples 
#' data(oui)
#' plot_diagnostics(oui)
plot_diagnostics <- function(oui, arrange = c("vertical", "horizontal")) {
  stopifnot(is(oui, "ouija_fit"))

  if (oui$inference_type == "vb") {
    warning("Diagnostic plots only make sense for HMC inference")
  }

  arrange <- match.arg(arrange)
  nrow <- switch(arrange,
                 vertical = 2,
                 horizontal = 1)
  plt <- cowplot::plot_grid(stan_trace(oui$fit, "lp__"),
                            stan_ac(oui$fit, "lp__"),
                            nrow = nrow)
  return(plt)
}


#' @name consistency
#' @param cmo An optional ordered consistency matrix
#' @param interpolate Passed to \code{geom_raster}
#' @export
#' @examples
#' data(oui)
#' plot_consistency(oui)
plot_consistency <- function(oui, cmo = NULL, interpolate = FALSE) {
  x <- value <- y <- NULL
  if (is.null(cmo)) cmo <- consistency_matrix_ordered(oui)
  diag(cmo) <- NA
  cmo_df <- as_data_frame(cmo)
  names(cmo_df) <- seq_len(oui$N)
  cmo_df$x <- as.numeric(names(cmo_df))
  cmo_df_tidy <- gather(cmo_df, y, value, -x)
  cmo_df_tidy$y <- as.numeric(cmo_df_tidy$y)
  ggplot(cmo_df_tidy, aes(x = x, y = y, fill = value)) +
    geom_raster(interpolate = interpolate) +
    viridis::scale_fill_viridis(name = "P(i>j)") +
    xlab(sprintf("Pseudotime order \u2192")) +
    ylab(sprintf("Pseudotime order \u2192")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.background = element_blank())
}

#' @name plotexprs
#' @export
#' @import viridis
#' @importFrom rstan extract
#' @importFrom coda mcmc
plot_switch_times <- function(oui) {
  Gene <- t0_mean <- lower <- upper <- NULL
  vpal <- viridis_pal()(8)
  k_trace <- extract(oui$fit, "k")$k
  kmean <- colMeans(k_trace)
  t0 <- extract(oui$fit, "t0")$t0
  t0_means <- colMeans(t0)
  t0_interval <- coda::HPDinterval(mcmc(t0))
  t0_df <- data_frame(t0_mean = t0_means,
                      lower = t0_interval[, 1],
                      upper = t0_interval[, 2],
                      kmean = kmean)
  t0_df$Gene <- colnames(oui$Y[, oui$response_type == "switch"])
  t0_df$Gene <- factor(t0_df$Gene, t0_df$Gene[order(t0_means)])
  ggplot(t0_df, aes(x = Gene, y = t0_mean, fill = kmean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey60",
                  width = 0.5, alpha = 0.5) +
    coord_flip() +
    geom_point(color = "grey50", shape = 21, size = 3) +
    ylab("Switch point") +
    scale_fill_gradient2(name = "Regulation", low = vpal[1], high = vpal[5]) +
    scale_color_gradient2(name = "Regulation", low = vpal[1], high = vpal[5]) +
    theme(legend.position = "top")
}

#' @name plotexprs
#' @export
plot_peak_times <- function(oui) {
  if (sum(oui$response_type == "transient") == 0) {
    stop("Fit must contain transient genes to plot peak times")
  }
  Gene <- p_mean <- lower <- upper <- NULL

  p_trace <- extract(oui$fit, "p")$p
  p_means <- colMeans(p_trace)
  p_interval <- coda::HPDinterval(mcmc(p_trace))
  p_df <- data_frame(p_mean = p_means,
                     lower = p_interval[, 1],
                     upper = p_interval[, 2])

  p_df$Gene <- colnames(oui$Y[, oui$response_type == "transient"])
  p_df$Gene <- factor(p_df$Gene, p_df$Gene[order(p_means)])
  ggplot(p_df, aes(x = Gene, y = p_mean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey60",
                  width = 0.5, alpha = 0.5) +
    coord_flip() +
    geom_point(color = "grey50", shape = 21, size = 3, fill = "grey10") +
    ylab("Peak time") +
    theme(legend.position = "top") +
    ylim(c(0, 1))
}

#' Plot fitted expression, switch times, and peak times
#' 
#' Plot the expression along pseudotime with the maximum a posteriori (MAP)
#' mean function as fitted by Ouija, or the MAP estimates of the
#' switch and peak times along with 95% HPD credible interval error bars.
#' 
#' @param oui A \code{ouija_fit} to plot
#' @param ncol The number of columns for the expression plot
#' @param nrow The number of rows for the expression plot
#' 
#' @name plotexprs
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr gather
#' @export
#' 
#' @examples 
#' data(oui)
#' plot_expression(oui)
#' plot_switch_times(oui)
#' \dontrun{
#' plot_peak_times(oui)
#' }
#' 
#' @return 
#' A \code{ggplot2} object plotting the expression, switch times, or peak times
plot_expression <- function(oui, ncol = 2, nrow = NULL) {
  gene <- ouija_pseudotime <- NULL

  expr_df <- as_data_frame(oui$Y) %>%
    mutate(ouija_pseudotime = map_pseudotime(oui)) %>%
    gather(gene, expression, -ouija_pseudotime)


  mu_df <- as_data_frame(predicted_expression(oui)) %>%
    mutate(ouija_pseudotime = map_pseudotime(oui)) %>%
    gather(gene, expression, -ouija_pseudotime) %>%
    arrange(ouija_pseudotime)


  ggplot(expr_df, aes(x = ouija_pseudotime)) +
    geom_point(aes(y = expression), alpha = 0.65, color = "grey30") +
    facet_wrap(~ gene, ncol = ncol, nrow = nrow,
               strip.position = "right") +
    geom_line(data = mu_df, aes(y = expression), size = 1.2, alpha = 0.7,
              color = "red") +
    scale_color_brewer(palette = "Set2", name = "Cell type") +
    theme(legend.position = "top",
          strip.background = element_rect(fill = "white")) +
    ylab("Normalised log expression") +
    xlab("Ouija pseudotime") +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5))
}
