# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(stats)))

# Generación de los datos ------------------------------------------------------

n           <- 100
lambda_true <- 3
y_max       <- 10

y <- rpois(
     n      = n,
     lambda = lambda_true
)

values <- sort(unique(y))

counts <- tabulate(
     bin   = match(y, values),
     nbins = length(values)
)

cat("Media muestral:", round(mean(y), 3), "\n")
cat("Varianza muestral:", round(var(y), 3), "\n")

print(table(y))

# Log-verosimilitud marginal ---------------------------------------------------

log_marginal_likelihood <- function(
          alpha,
          lambda,
          values,
          counts,
          n
) {
     log_numerator <- 0
     
     for (j in seq_along(values)) {
          log_q <- dpois(
               x      = values[j],
               lambda = lambda,
               log    = TRUE
          )
          
          if (!is.finite(log_q)) {
               return(-Inf)
          }
          
          log_numerator <- log_numerator + log(alpha) + log_q
          
          if (counts[j] > 1) {
               alpha_q <- alpha * exp(log_q)
               
               log_numerator <- log_numerator +
                    sum(log(alpha_q + seq_len(counts[j] - 1)))
          }
     }
     
     log_denominator <- sum(log(alpha + 0:(n - 1)))
     
     log_numerator - log_denominator
}

# Log-posterior en la escala transformada --------------------------------------

log_posterior_log_scale <- function(
          log_alpha,
          log_lambda,
          values,
          counts,
          n,
          a_alpha,
          b_alpha,
          a_lambda,
          b_lambda
) {
     alpha  <- exp(log_alpha)
     lambda <- exp(log_lambda)
     
     out <- log_marginal_likelihood(
          alpha  = alpha,
          lambda = lambda,
          values = values,
          counts = counts,
          n      = n
     ) +
          dgamma(
               x     = alpha,
               shape = a_alpha,
               rate  = b_alpha,
               log   = TRUE
          ) +
          dgamma(
               x     = lambda,
               shape = a_lambda,
               rate  = b_lambda,
               log   = TRUE
          ) +
          log_alpha +
          log_lambda
     
     if (!is.finite(out)) {
          return(-Inf)
     }
     
     out
}

# Algoritmo MCMC ---------------------------------------------------------------

mcmc <- function(
          y,
          a_alpha,
          b_alpha,
          a_lambda,
          b_lambda,
          n_burn,
          n_save,
          n_thin,
          adapt_interval = 50,
          target_accept = 0.44
) {
     n <- length(y)
     
     values <- sort(unique(y))
     
     counts <- tabulate(
          bin   = match(y, values),
          nbins = length(values)
     )
     
     B <- n_burn + n_save * n_thin
     
     log_alpha  <- log(a_alpha / b_alpha)
     log_lambda <- log(max(mean(y), 0.1))
     
     log_sd_alpha  <- log(0.30)
     log_sd_lambda <- log(0.20)
     
     current_log_post <- log_posterior_log_scale(
          log_alpha  = log_alpha,
          log_lambda = log_lambda,
          values     = values,
          counts     = counts,
          n          = n,
          a_alpha    = a_alpha,
          b_alpha    = b_alpha,
          a_lambda   = a_lambda,
          b_lambda   = b_lambda
     )
     
     alpha_chain  <- numeric(n_save)
     lambda_chain <- numeric(n_save)
     
     batch_accept_alpha  <- 0
     batch_accept_lambda <- 0
     
     post_accept_alpha  <- 0
     post_accept_lambda <- 0
     
     save_index <- 0
     
     for (b in seq_len(B)) {
          # Actualización de alpha
          proposal_log_alpha <- log_alpha + rnorm(
               n    = 1,
               mean = 0,
               sd   = exp(log_sd_alpha)
          )
          
          proposal_log_post <- log_posterior_log_scale(
               log_alpha  = proposal_log_alpha,
               log_lambda = log_lambda,
               values     = values,
               counts     = counts,
               n          = n,
               a_alpha    = a_alpha,
               b_alpha    = b_alpha,
               a_lambda   = a_lambda,
               b_lambda   = b_lambda
          )
          
          if (log(runif(1)) < proposal_log_post - current_log_post) {
               log_alpha        <- proposal_log_alpha
               current_log_post <- proposal_log_post
               
               if (b <= n_burn) {
                    batch_accept_alpha <- batch_accept_alpha + 1
               } else {
                    post_accept_alpha <- post_accept_alpha + 1
               }
          }
          
          # Actualización de lambda
          proposal_log_lambda <- log_lambda + rnorm(
               n    = 1,
               mean = 0,
               sd   = exp(log_sd_lambda)
          )
          
          proposal_log_post <- log_posterior_log_scale(
               log_alpha  = log_alpha,
               log_lambda = proposal_log_lambda,
               values     = values,
               counts     = counts,
               n          = n,
               a_alpha    = a_alpha,
               b_alpha    = b_alpha,
               a_lambda   = a_lambda,
               b_lambda   = b_lambda
          )
          
          if (log(runif(1)) < proposal_log_post - current_log_post) {
               log_lambda       <- proposal_log_lambda
               current_log_post <- proposal_log_post
               
               if (b <= n_burn) {
                    batch_accept_lambda <- batch_accept_lambda + 1
               } else {
                    post_accept_lambda <- post_accept_lambda + 1
               }
          }
          
          # Adaptación durante el calentamiento
          if (b <= n_burn && b %% adapt_interval == 0) {
               batch_index     <- b / adapt_interval
               adaptation_rate <- 1 / sqrt(batch_index)
               
               acceptance_alpha  <- batch_accept_alpha / adapt_interval
               acceptance_lambda <- batch_accept_lambda / adapt_interval
               
               log_sd_alpha <- log_sd_alpha +
                    adaptation_rate * (acceptance_alpha - target_accept)
               
               log_sd_lambda <- log_sd_lambda +
                    adaptation_rate * (acceptance_lambda - target_accept)
               
               log_sd_alpha <- min(
                    max(log_sd_alpha, log(0.01)),
                    log(2)
               )
               
               log_sd_lambda <- min(
                    max(log_sd_lambda, log(0.01)),
                    log(2)
               )
               
               batch_accept_alpha  <- 0
               batch_accept_lambda <- 0
          }
          
          # Almacenamiento de las extracciones posteriores
          if (b > n_burn && (b - n_burn) %% n_thin == 0) {
               save_index <- save_index + 1
               
               alpha_chain[save_index]  <- exp(log_alpha)
               lambda_chain[save_index] <- exp(log_lambda)
          }
     }
     
     n_post_iterations <- n_save * n_thin
     
     list(
          alpha_chain      = alpha_chain,
          lambda_chain     = lambda_chain,
          acceptance_alpha = post_accept_alpha / n_post_iterations,
          acceptance_lambda = post_accept_lambda / n_post_iterations,
          proposal_sd_alpha = exp(log_sd_alpha),
          proposal_sd_lambda = exp(log_sd_lambda)
     )
}

# Simulación posterior ---------------------------------------------------------

a_alpha  <- 1
b_alpha  <- 1
a_lambda <- 1
b_lambda <- 1

n_burn         <- 10000
n_save         <- 10000
n_thin         <- 10
adapt_interval <- 50
target_accept  <- 0.44

set.seed(123)

samples <- mcmc(
     y              = y,
     a_alpha        = a_alpha,
     b_alpha        = b_alpha,
     a_lambda       = a_lambda,
     b_lambda       = b_lambda,
     n_burn         = n_burn,
     n_save         = n_save,
     n_thin         = n_thin,
     adapt_interval = adapt_interval,
     target_accept  = target_accept
)

cat(
     "Tasa de aceptación de alpha:",
     round(samples$acceptance_alpha, 3),
     "\n"
)

cat(
     "Tasa de aceptación de lambda:",
     round(samples$acceptance_lambda, 3),
     "\n"
)

# Resúmenes posteriores --------------------------------------------------------

summarize_chain <- function(
          x,
          parameter,
          acceptance
) {
     q <- quantile(
          x     = x,
          probs = c(0.025, 0.5, 0.975),
          names = FALSE
     )
     
     data.frame(
          parameter  = parameter,
          mean       = mean(x),
          sd         = sd(x),
          q025       = q[1],
          median     = q[2],
          q975       = q[3],
          acceptance = acceptance
     )
}

posterior_summary <- rbind(
     summarize_chain(
          x          = samples$alpha_chain,
          parameter  = "alpha",
          acceptance = samples$acceptance_alpha
     ),
     summarize_chain(
          x          = samples$lambda_chain,
          parameter  = "lambda",
          acceptance = samples$acceptance_lambda
     )
)

rownames(posterior_summary) <- c("alpha", "lambda")

print(
     round(
          posterior_summary[, -1],
          digits = 3
     )
)

# Diagnósticos de las cadenas --------------------------------------------------

par(
     mfrow = c(2, 2),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = seq_along(samples$alpha_chain),
     y    = samples$alpha_chain,
     type = "p",
     pch  = ".",
     col  = "royalblue3",
     xlab = "Iteración guardada",
     ylab = expression(alpha),
     main = ""
)

acf(
     x       = samples$alpha_chain,
     lag.max = 50,
     main    = "",
     xlab    = "Rezago",
     ylab    = "Autocorrelación",
     col     = "royalblue3"
)

plot(
     x    = seq_along(samples$lambda_chain),
     y    = samples$lambda_chain,
     type = "p",
     pch  = ".",
     col  = "firebrick3",
     xlab = "Iteración guardada",
     ylab = expression(lambda),
     main = ""
)

acf(
     x       = samples$lambda_chain,
     lag.max = 50,
     main    = "",
     xlab    = "Rezago",
     ylab    = "Autocorrelación",
     col     = "firebrick3"
)

dev.off()

# Simulación de las masas posteriores de G -------------------------------------

draw_dirichlet <- function(shape) {
     z <- rgamma(
          n     = length(shape),
          shape = shape,
          rate  = 1
     )
     
     z / sum(z)
}

build_mass_draws <- function(
          y,
          alpha_chain,
          lambda_chain,
          y_max
) {
     y_grid <- 0:y_max
     n      <- length(y)
     B      <- length(alpha_chain)
     
     observed_counts <- vapply(
          X         = y_grid,
          FUN       = function(k) sum(y == k),
          FUN.VALUE = numeric(1)
     )
     
     tail_count <- sum(y > y_max)
     
     G_mass_chain <- matrix(
          data = NA_real_,
          nrow = B,
          ncol = length(y_grid) + 1
     )
     
     predictive_mass_chain <- matrix(
          data = NA_real_,
          nrow = B,
          ncol = length(y_grid) + 1
     )
     
     for (b in seq_len(B)) {
          base_mass <- dpois(
               x      = y_grid,
               lambda = lambda_chain[b]
          )
          
          base_tail <- ppois(
               q          = y_max,
               lambda     = lambda_chain[b],
               lower.tail = FALSE
          )
          
          shape <- c(
               alpha_chain[b] * base_mass + observed_counts,
               alpha_chain[b] * base_tail + tail_count
          )
          
          predictive_mass_chain[b, ] <- shape / (alpha_chain[b] + n)
          
          G_mass_chain[b, ] <- draw_dirichlet(shape = shape)
     }
     
     list(
          y_grid               = y_grid,
          observed_counts      = observed_counts,
          tail_count           = tail_count,
          predictive_mass_chain = predictive_mass_chain,
          G_mass_chain         = G_mass_chain
     )
}

mass_draws <- build_mass_draws(
     y            = y,
     alpha_chain  = samples$alpha_chain,
     lambda_chain = samples$lambda_chain,
     y_max        = y_max
)

# Distribución predictiva posterior --------------------------------------------

y_grid <- mass_draws$y_grid

true_mass <- c(
     dpois(
          x      = y_grid,
          lambda = lambda_true
     ),
     ppois(
          q          = y_max,
          lambda     = lambda_true,
          lower.tail = FALSE
     )
)

empirical_mass <- c(
     mass_draws$observed_counts,
     mass_draws$tail_count
) / n

predictive_mean <- colMeans(
     mass_draws$predictive_mass_chain
)

credible_limits <- apply(
     X      = mass_draws$G_mass_chain,
     MARGIN = 2,
     FUN    = quantile,
     probs  = c(0.025, 0.975)
)

predictive_summary <- data.frame(
     category = c(
          as.character(y_grid),
          paste0("$>", y_max, "$")
     ),
     true            = true_mass,
     empirical       = empirical_mass,
     predictive_mean = predictive_mean,
     lower           = credible_limits[1, ],
     upper           = credible_limits[2, ]
)

print(
     predictive_summary,
     digits    = 4,
     row.names = FALSE
)

# Figura de la distribución predictiva -----------------------------------------

plot_index <- seq_along(y_grid)

true_plot       <- predictive_summary$true[plot_index]
empirical_plot  <- predictive_summary$empirical[plot_index]
predictive_plot <- predictive_summary$predictive_mean[plot_index]
lower_plot      <- predictive_summary$lower[plot_index]
upper_plot      <- predictive_summary$upper[plot_index]

pdf(
     file      = "mdp_poisson_predictive_pmf.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

y_upper <- 1.10 * max(
     true_plot,
     empirical_plot,
     upper_plot
)

plot(
     x    = y_grid,
     y    = true_plot,
     type = "b",
     pch  = 16,
     lwd  = 2,
     col  = "black",
     ylim = c(0, y_upper),
     xlab = "Conteo",
     ylab = "Probabilidad",
     main = ""
)

arrows(
     x0     = y_grid,
     y0     = lower_plot,
     x1     = y_grid,
     y1     = upper_plot,
     angle  = 90,
     code   = 3,
     length = 0.0,
     lwd    = 1.5,
     col    = "firebrick3"
)

lines(
     x    = y_grid,
     y    = predictive_plot,
     type = "b",
     pch  = 16,
     lwd  = 2,
     col  = "firebrick3"
)

legend(
     "topright",
     legend = c(
          "Poisson(3)",
          "Media e IC 95%"
     ),
     col = c("black", "firebrick3"),
     pch = c(16, 16),
     lty = c(1, 1),
     lwd = c(2, 2),
     bty = "n",
     cex = 0.85
)

dev.off()

# Fin --------------------------------------------------------------------------