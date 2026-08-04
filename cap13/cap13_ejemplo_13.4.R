# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Especificación ---------------------------------------------------------------

n <- 100
B <- 100000

scenario <- c(
     "Concentrada alrededor de 1",
     "Dispersa alrededor de 1",
     "Concentrada alrededor de 4"
)

a_alpha <- c(4, 1, 4)
b_alpha <- c(4, 1, 1)

prior_legend <- c(
     "G(4, 4)",
     "G(1, 1)",
     "G(4, 1)"
)

# Simulación de alpha ----------------------------------------------------------

set.seed(123)

alpha_prior <- lapply(
     X   = seq_along(scenario),
     FUN = function(j) {
          rgamma(
               n     = B,
               shape = a_alpha[j],
               rate  = b_alpha[j]
          )
     }
)

# Simulación de K_n ------------------------------------------------------------

set.seed(123)

K_prior <- lapply(
     X   = alpha_prior,
     FUN = function(alpha_j) {
          K_j <- rep(
               x     = 1L,
               times = B
          )
          
          for (i in 2:n) {
               K_j <- K_j + rbinom(
                    n    = B,
                    size = 1,
                    prob = alpha_j / (alpha_j + i - 1)
               )
          }
          
          K_j
     }
)

# Resúmenes --------------------------------------------------------------------

mean_K <- vapply(
     X         = K_prior,
     FUN       = mean,
     FUN.VALUE = numeric(1)
)

sd_K <- vapply(
     X         = K_prior,
     FUN       = sd,
     FUN.VALUE = numeric(1)
)

quantile_probs <- c(
     0.025,
     0.500,
     0.975
)

quantile_K <- t(
     vapply(
          X         = K_prior,
          FUN       = quantile,
          FUN.VALUE = numeric(length(quantile_probs)),
          probs     = quantile_probs,
          names     = FALSE
     )
)

tab_prior <- data.frame(
     Especificacion = scenario,
     a_alpha        = a_alpha,
     b_alpha        = b_alpha,
     Media_alpha    = a_alpha / b_alpha,
     DE_alpha       = sqrt(a_alpha) / b_alpha,
     Media_K        = mean_K,
     DE_K           = sd_K,
     Q_0025         = quantile_K[, 1],
     Q_0500         = quantile_K[, 2],
     Q_0975         = quantile_K[, 3]
)

print(tab_prior)

# Distribuciones empíricas de K_n ----------------------------------------------

max_K <- max(unlist(K_prior))

K_grid <- seq(
     from = 1,
     to   = max_K
)

prob_K <- sapply(
     X   = K_prior,
     FUN = function(K_j) {
          tabulate(
               bin   = K_j,
               nbins = max_K
          ) / length(K_j)
     }
)

# Figura -----------------------------------------------------------------------

scenario_col <- c(
     "#0072B2",
     "#009E73",
     "#D55E00"
)

pdf(
     file      = "dpm_prior_number_components.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

matplot(
     x    = K_grid,
     y    = prob_K,
     type = "o",
     lty  = 1,
     lwd  = 2,
     pch  = c(16, 17, 15),
     col  = scenario_col,
     xlab = expression(K[100]),
     ylab = "Probabilidad",
     ylim = c(0, 1.05 * max(prob_K))
)

legend(
     x      = "topright",
     legend = prior_legend,
     lty    = 1,
     lwd    = 2,
     pch    = c(16, 17, 15),
     col    = scenario_col,
     bty    = "n"
)

dev.off()

# Fin --------------------------------------------------------------------------