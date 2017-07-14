




#' Plot a \code{ouija_fit}
#' 
#' Plot a \code{ouija_fit} object. Returns a plot of either
#' \itemize{
#' \item{diagnostic} Trace and autocorrelation plots of the log-posterior
#' probability. Underlying call is to \code{\link{plot_ouija_fit_diagnostics}}
#' \item{behaviour}  Gene expression as a function of the MAP pseudotime with a red
#' line denoting the mean sigmoid trend
#' \item{heatmap} A heatmap of gene expression as a function of pseudotime
#' across different pseudotime samples.
#' \item{pp} Density plots comparing prior to posterior distributions for the
#' activation parameters.
#' \item{dropout} The relationship between latent expression
#' value and dropout probability. 
#' }
#' 
#' @param x An object of class \code{ouija_fit}
#' @param what One of
#' \itemize{
#' \item \code{heatmap} This produces a heatmap of gene expression as a function of pseudotime
#' across different pseudotime samples. Underlying call is to \code{\link{plot_ouija_fit_heatmap}}.
#' \item \code{behaviour} This plots gene expression as a function of the MAP pseudotime with a red
#' line denoting the mean sigmoid trend. Underlying call is to
#' \code{\link{plot_ouija_fit_behaviour}}
#' \item \code{diagnostic} This returns trace and autocorrelation plots of the log-posterior
#' probability. Underlying call is to \code{\link{plot_ouija_fit_diagnostics}}
#' \item \code{pp} This returns density plots of posterior distributions for either activation
#' strength parameters \code{k} or activation time parameters \code{t0}.
#' \item \code{dropout} Returns a plot showing the relationship between latent expression
#' value and dropout probability. 
#' Underlying call is to \code{\link{plot_ouija_fit_dropout_probability}}
#' }
#' @param ... Additional arguments passed to the corresponding functions
#' 
#' @return A \code{ggplot2} plot.
#' @method plot ouija_fit
#' @export
#' 
#' @examples 
#' data(oui)
#' plot(oui)
plot.ouija_fit <- function(x, what = c("behaviour", "behavior", "diagnostic", 
                                       "heatmap", "pp", "dropout"), ...) {
  what <- match.arg(what)
  plt <- switch(what,
                heatmap = plot_ouija_fit_heatmap(x, ...),
                behaviour = plot_ouija_fit_behaviour(x, ...),
                behavior = plot_ouija_fit_behaviour(x, ...),
                diagnostic = plot_ouija_fit_diagnostics(x, ...),
                pp = plot_ouija_fit_pp(x, ...),
                dropout = plot_ouija_fit_dropout_probability(x, ...))
  return(plt)
}

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
#' plot_ouija_fit_diagnostics(oui)
plot_ouija_fit_diagnostics <- function(oui, arrange = c("vertical", "horizontal")) {
  stopifnot(is(oui, "ouija_fit"))
  
  if(oui$inference_type == "vb") {
    warning("Diagnostic plots only make sense for HMC inference")
  }
  
  arrange <- match.arg(arrange)
  nrow <- switch(arrange,
                 vertical = 2,
                 horizontal = 1)
  plt <- cowplot::plot_grid(stan_trace(oui$fit, "lp__"), stan_ac(oui$fit, "lp__"), nrow = nrow)
  return(plt)
}

#' Plot heatmaps of gene expression changes
#' 
#' Produces a heatmap of gene expression as a function of pseudotime
#' across different pseudotime samples.
#' 
#' @param oui An object of class \code{ouija_fit}
#' @param samples Number of posterior pseudotime samples to use (number of rows of heatmap)
#' @param genes A vector that subsets the gene expression matrix. Defaults to the first \code{g}
#' genes, where \code{g} is either 4 or the number of genes in the model if less than 4.
#' @param output If \code{grid} then \code{cowplot::plot_grid} is called and a grid plot
#' of all genes is returned. If \code{plotlist} then a list of \code{ggplot2} objects is returned 
#' for the user to customise further
#' @param show_legend Logical. If \code{TRUE} then the legend (ie gene expression magnitude) is
#' displayed for each heatmap.
#' @param ... Additional arguments passed to \code{cowplot::plot_grid}
#' 
#' @importFrom rstan extract
#' @importFrom cowplot plot_grid
#' @importFrom viridis scale_fill_viridis
#' @importFrom reshape2 melt
#' @importFrom methods is
#' 
#' @export
#' 
#' @return A \code{ggplot2} object.
#' @examples 
#' data(oui)
#' plot_ouija_fit_heatmap(oui)
plot_ouija_fit_heatmap <- function(oui, samples = 50, genes = seq_len(min(oui$G, 6)),
                                   output = c("grid", "plotlist"), 
                                   show_legend = FALSE, ...) {
  stopifnot(is(oui, "ouija_fit"))
  output <- match.arg(output)
  
  ttrace <- extract(oui$fit, "t")$t
  to_sample <- sample(seq_len(min(nrow(ttrace), samples)))
  ttrace <- ttrace[to_sample, ]
  cell_orders <- apply(ttrace, 1, order)
  # apply over genes
  plts <- lapply(genes, function(g) {
    yg <- oui$Y[,g]
    Xg <- apply(cell_orders, 2, function(or) yg[or])
    Xg <- data.frame(Xg)
    names(Xg) <- 1:ncol(Xg)
    Xg$x <- 1:nrow(Xg)
    Xm <- melt(Xg, id.vars = "x", variable.name = "y", value.name = "Expression")
    
    x <- y <- Expression <- NULL
    plt <- ggplot(Xm, aes(x = x, y = y, fill = Expression)) + geom_tile() +
      xlab("Cell") + ylab(expression("Pseudotime\n sample")) + 
      scale_fill_viridis() + 
      theme_bw() + 
      theme(axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank())
    if(!show_legend) plt <- plt + theme(legend.position = "none")
    return( plt )
  })
  
  if(output == "grid") {
    return(plot_grid(plotlist = plts, scale = 0.9, ...))
  } else {
    return( plts )
  }
}


#' Plot gene expression as a function of MAP pseudotime
#' 
#' Plot gene expression as a function of the MAP pseudotime with a red
#' line denoting a LOESS fit (showing the overall trend). Genes are plotted with
#' one per grid square (using a call to \code{facet_wrap(~ gene)}).
#' 
#' @param oui An object of class \code{ouija_fit}
#' @param genes A vector that subsets the gene expression matrix. Defaults to the first \code{g}
#' genes, where \code{g} is either 4 or the number of genes in the model if less than 4.
#' @param expression_units The label for the y-axis of the plot. Defaults to
#' "log2(TPM + 1)".
#' 
#' @importFrom reshape2 melt
#' @importFrom rstan extract
#' @importFrom MCMCglmm posterior.mode
#' @importFrom coda mcmc
#' @importFrom dplyr inner_join
#' @import ggplot2
#' @importFrom methods is
#' 
#' @export
#' 
#' @return An object of class \code{ggplot2}
#' @examples 
#' data(oui)
#' plot_ouija_fit_behaviour(oui)
plot_ouija_fit_behaviour <- function(oui, genes = seq_len(min(oui$G, 6)),
                                     expression_units = "log2(TPM+1)") {
  stopifnot(is(oui, "ouija_fit"))
  tmap <- map_pseudotime(oui)
  Y <- oui$Y
  
  ## want to plot sigmoid function so need MAP estimates
  sig_map <- predicted_expression(oui)
  
  
  ## Create data frame for gene expression values
  Y <- Y[,genes]
  dy <- data.frame(Y, pseudotime = tmap)
  dm <- melt(dy, id.vars = "pseudotime", 
             variable.name = "gene", 
             value.name = "expression")
  
  S <- sig_map[,genes]
  ds <- data.frame(S, pseudotime = tmap)
  dm2 <- melt(ds, id.vars = "pseudotime",
              variable.name = "gene",
              value.name = "predicted_expression")
  
  dm_joined <- inner_join(dm, dm2, by = c("pseudotime", "gene"))
  pseudotime <- predicted_expression <- NULL
  
  plt <- ggplot(dm_joined, aes(x = pseudotime, y = expression, colour = "Measured")) + 
    geom_point() +
    facet_wrap(~ gene, scales = "free_y") +
    xlab("MAP pseudotime") + ylab(paste("Expression", expression_units)) + 
    cowplot::theme_cowplot() +
    scale_x_continuous(breaks = c(0, 0.5, 1))
  plt <- plt + 
    geom_line(aes(x = pseudotime, y = predicted_expression, color = 'Predicted'), 
              size = 2, alpha = 0.7) +
    scale_colour_manual(values = c("Predicted" = "red", "Measured" = "black"), name = element_blank()) +
    theme(legend.position = "bottom") 
  return( plt )
}

#' Plot heatmaps showing comparisons of measured data and imputed
#' 
#' @param oui An object of class \code{ouija_fit}.
#' @param return_plotlist If TRUE then the list of \code{ggplot}s is returned
#' instead of being plotted with \code{cowplot::plot_grid}
#' 
#' 
#' @export
#' @import ggplot2
#' @importFrom methods is
#' 
#' @return Either a list of plots of class \code{ggplot} or a single 
#' \code{ggplot} showing them
#' 
#' @examples 
#' data(oui)
#' plot_ouija_fit_comparison(oui)
plot_ouija_fit_comparison <- function(oui, return_plotlist = FALSE) {
  stopifnot(is(oui, "ouija_fit"))
  X <- oui$Y
  Z <- rexprs(oui)
  tmap <- map_pseudotime(oui)
  
  make_tile_plot <- function(X, tmap) {
    Xp <- apply(X, 2, function(x) (x - min(x)) / (max(x) - min(x)))
    
    dfx <- data.frame(Xp, pseudotime = rank(tmap)) %>%
      melt(id.vars = "pseudotime", variable.name = "gene", value.name = "expression")
    
    pseudotime <- gene <- NULL
    
    ggplot(dfx, aes(x = pseudotime, y = gene, fill = expression)) + 
      geom_tile() + scale_fill_viridis(name = expression("Relative\nexpression")) + 
      theme_bw() + 
      theme(axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.title.y = element_blank()) +
      xlab("Pseudotime order")
  }
  plts <- lapply(list(X, Z), make_tile_plot, tmap)
  if(return_plotlist) return(plts)
  cowplot::plot_grid(plotlist = plts, nrow = 2, 
                     labels = c("Measured", "Reconstructed"))
}


#' Plot dropout probability
#' 
#' Plot the probability of a dropout as a function of latent expression value. In ouija
#' this is implemented via logit regression, so the probability of dropout is related
#' to latent expression \eqn{\mu_{ig}} via
#' \deqn{\frac{1}{1 + \exp(-(\beta_0 + \beta_1 \mu_{ig}))}}
#' The red curve shows the MAP estimate of the relationship, while the grey lines show 
#' posterior samples of the relationship.
#' 
#' @param oui An object of class \code{ouija_fit}
#' @param posterior_samples Number of posterior samples to add to the plot. If 0, only
#' the MAP estimate is plotted.
#' 
#' @importFrom rstan extract
#' @importFrom MCMCglmm posterior.mode
#' @importFrom coda mcmc
#' @importFrom methods is
#' @export
#' 
#' @return An object of type \code{ggplot}
#' 
#' @examples 
#' data(oui)
#' plot_ouija_fit_dropout_probability(oui)
plot_ouija_fit_dropout_probability <- function(oui, posterior_samples = 40) {
  stopifnot(is(oui, "ouija_fit"))
  x_range <- range(as.vector(oui$Y))
  dsig <- function(x, beta0, beta1) 1 / (1 + exp(-(beta0 + beta1 * x)))
  ext <- extract(oui$fit, "beta")
  beta_map <- posterior.mode(mcmc(ext$beta))
  
  Mean_expression <- NULL
  plt <- ggplot(data.frame(Mean_expression = x_range), aes(x = Mean_expression)) 
  
  
  if(posterior_samples > 0) {
    total_samples <- dim(ext$beta)[1]
    to_sample <- sample(total_samples, posterior_samples)
    for(i in seq_len(posterior_samples)) {
      beta <- ext$beta[to_sample[i], ]
      plt <- plt + stat_function(fun = dsig, 
                                 args = list(beta0 = beta[1], beta1 = beta[2]),
                                 alpha = 0.3)
    }
  }
  plt <- plt + stat_function(fun = dsig, 
                             args = list(beta0 = beta_map[1], beta1 = beta_map[2]),
                             colour = "red") +
    ylab("Dropout probability") + xlab("Latent expression")
  return( plt )
}

#' Prior-posterior density plots for activation parameters
#' 
#' @param oui An object of class \code{ouija_fit}
#' @param genes A numeric vector indicating indices of genes to plot
#' @param param Either plot the activation strength \code{k} parameter or
#' activation time \code{t0} parameter
#' 
#' @return An object of class \code{ggplot2}
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats rnorm
#' @importFrom methods is
#' 
#' @examples 
#' data(oui)
#' plot_ouija_fit_pp(oui)
plot_ouija_fit_pp <- function(oui, genes = seq_len(ncol(oui$Y)), param = c("k", "t0")) {
  stopifnot(is(oui, "ouija_fit"))
  ngenes <- ncol(oui$Y)
  if(any(genes > ngenes)) stop("Genes outside possible range")
  
  param <- match.arg(param)
  
  label <- switch(param,
                  k = "activation strength",
                  t0 = "activation time")
  label <- paste("Prior-posterior for", label)
  
  mu <- sd <- NULL
  if(param == "k") {
    mu <- oui$strengths[genes]
    sd <- oui$strength_sd[genes]
  } else {
    mu <- oui$times[genes]
    sd <- oui$time_sd[genes]
  }
  dfs <- lapply(genes, function(g) {
    posterior <- extract(oui$fit, param)[[param]][,g]
    prior <- rnorm(length(posterior), mu[g], sd[g])
    
    dm <- data.frame(prior = prior, posterior = posterior) %>%
      melt(variable.name = "Prob", value.name = "Sample") %>%
      dplyr::mutate(Gene = colnames(oui$Y)[g])
    return(dm)
  })
  dfs <- do.call(rbind, dfs)
  
  Sample <- Prob <- NULL
  
  ggplot(dfs, aes(x = Sample, fill = Prob)) + geom_density(colour = "Black", alpha = 0.8) +
    cowplot::theme_cowplot() +
    ylab("Density") + xlab("Parameter value") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ Gene, scales = "free") +
    ggtitle(label)
}

#' Plot confusion matrix
#' @export
plot_confusion <- function(oui, cmo = NULL, interpolate = FALSE) {
  if(is.null(cmo)) cmo <- confusion_matrix_ordered(oui)
  diag(cmo) <- NA
  cmo_df <- as_data_frame(cmo)
  names(cmo_df) <- seq_len(oui$N)
  cmo_df$x <- as.numeric(names(cmo_df))
  cmo_df_tidy <- gather(cmo_df, y, value, -x)
  cmo_df_tidy$y <- as.numeric(cmo_df_tidy$y)
  ggplot(cmo_df_tidy, aes(x = x, y = y, fill = value)) +
    geom_raster(interpolate = interpolate) + 
    viridis::scale_fill_viridis(name = "P(i>j)") +
    xlab(sprintf("Pseudotime order \u2192")) + ylab(sprintf("Pseudotime order \u2192")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank())
}

#' Plot switch times
#' @import viridis
#' @importFrom rstan extract
#' @importFrom coda mcmc
plot_switch_times <- function(oui) {
  vpal <- viridis_pal()(8)
  k_trace <- extract(oui$fit, "k")$k
  kmean <- colMeans(k_trace)
  t0 <- extract(oui$fit, "t0")$t0
  t0_means <- colMeans(t0)
  t0_interval <- HPDinterval(mcmc(t0))
  t0_df <- data_frame(t0_mean = t0_means, lower = t0_interval[,1], upper = t0_interval[,2],
                      kmean = kmean)
  t0_df$Gene <- colnames(oui$Y[, oui$response_type == "switch"])
  t0_df$Gene <- factor(t0_df$Gene, t0_df$Gene[order(t0_means)])
  ggplot(t0_df, aes(x = Gene, y = t0_mean, fill = kmean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey60", 
                  width = 0.5, alpha = 0.5) +
    coord_flip() +
    geom_point(color = 'grey50', shape = 21, size = 3) +
    ylab("Switch point") +
    scale_fill_gradient2(name = "Regulation", low = vpal[1], high = vpal[5]) +
    scale_color_gradient2(name = "Regulation", low = vpal[1], high = vpal[5]) +
    theme(legend.position = "top")
}

#' Plot transient times
#' @export
plot_transient_times <- function(oui) {
  p_trace <- extract(oui$fit, "p")$p
  p_means <- colMeans(p_trace)
  p_interval <- HPDinterval(mcmc(p_trace))
  p_df <- data_frame(p_mean = p_means, lower = p_interval[,1], upper = p_interval[,2])
  
  p_df$Gene <- colnames(oui$Y[, oui$response_type == "transient"])
  p_df$Gene <- factor(p_df$Gene, p_df$Gene[order(p_means)])
  ggplot(p_df, aes(x = Gene, y = p_mean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey60", 
                  width = 0.5, alpha = 0.5) +
    coord_flip() +
    geom_point(color = 'grey50', shape = 21, size = 3, fill = 'grey10') +
    ylab("Peak time") +
    theme(legend.position = "top") +
    ylim(c(0,1))
}

#' Plot fitted expression
#' @import dplyr
#' @importFrom tidyr gather
#' @export
plot_expression <- function(oui) {
  
  expr_df <- as_data_frame(oui$Y) %>% 
    mutate(ouija_pseudotime = map_pseudotime(oui)) %>% 
    gather(gene, expression, -ouija_pseudotime)
  
  
  mu_df <- as_data_frame(predicted_expression(oui)) %>% 
    mutate(ouija_pseudotime = map_pseudotime(oui)) %>% 
    gather(gene, expression, -ouija_pseudotime) %>% 
    arrange(ouija_pseudotime)
  
  
  ggplot(expr_df, aes(x = ouija_pseudotime)) +
    geom_point(aes(y = expression), alpha = 0.65, color = 'grey30') + 
    facet_wrap(~ gene, ncol = 2, 
               strip.position = 'right') +
    geom_line(data = mu_df, aes(y = expression), size = 1.2, alpha = 0.7,
              color = 'red') +
    scale_color_brewer(palette = "Set2", name = "Cell type") +
    theme(legend.position = "top",
          strip.background = element_rect(fill="white")) +
    ylab("Normalised log expression") +
    xlab("Ouija pseudotime") +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5)) 
  
}


plot_regulations <- function(reg_df) {
   ggplot(reg_df, aes(x = param_diffs)) +
     geom_histogram() + facet_wrap(~ label) +
     labs(x = "Difference in switch times", 
          subtitle = "Posterior difference in regulation time")
}