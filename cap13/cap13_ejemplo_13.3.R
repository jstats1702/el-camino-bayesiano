# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Especificación del modelo ----------------------------------------------------

m0       <- 0
sigma2   <- 1
s02      <- c(0.25, 4)
scenario <- c("Concentrada", "Dispersa")
B        <- 100000

# Resúmenes analíticos ---------------------------------------------------------

sd_theta      <- sqrt(s02)
sd_predictive <- sqrt(sigma2 + s02)

q_theta_lower <- qnorm(
     p    = 0.025,
     mean = m0,
     sd   = sd_theta
)

q_theta_upper <- qnorm(
     p    = 0.975,
     mean = m0,
     sd   = sd_theta
)

q_predictive_lower <- qnorm(
     p    = 0.025,
     mean = m0,
     sd   = sd_predictive
)

q_predictive_upper <- qnorm(
     p    = 0.975,
     mean = m0,
     sd   = sd_predictive
)

prob_tail <- pnorm(
     q    = -3,
     mean = m0,
     sd   = sd_predictive
) + 1 - pnorm(
     q    = 3,
     mean = m0,
     sd   = sd_predictive
)

# Simulación de Monte Carlo ----------------------------------------------------

set.seed(123)

theta_prior <- lapply(
     X   = s02,
     FUN = function(s02_j) {
          rnorm(
               n    = B,
               mean = m0,
               sd   = sqrt(s02_j)
          )
     }
)

y_prior <- lapply(
     X   = theta_prior,
     FUN = function(theta_j) {
          rnorm(
               n    = B,
               mean = theta_j,
               sd   = sqrt(sigma2)
          )
     }
)

# Tabla de resultados ----------------------------------------------------------

tab_prior <- data.frame(
     Especificacion             = scenario,
     s02                        = s02,
     SD_theta                   = sd_theta,
     Limite_theta_inferior      = q_theta_lower,
     Limite_theta_superior      = q_theta_upper,
     SD_predictiva              = sd_predictive,
     Limite_predictivo_inferior = q_predictive_lower,
     Limite_predictivo_superior = q_predictive_upper,
     Prob_abs_Y_mayor_3         = prob_tail
)

print(tab_prior)

# Densidades -------------------------------------------------------------------

# Mallas de evaluación
theta_grid <- seq(
     from       = -6,
     to         = 6,
     length.out = 1000
)

y_grid <- seq(
     from       = -8,
     to         = 8,
     length.out = 1000
)

# Densidades analíticas
density_base <- sapply(
     X   = s02,
     FUN = function(s02_j) {
          dnorm(
               x    = theta_grid,
               mean = m0,
               sd   = sqrt(s02_j)
          )
     }
)

density_predictive <- sapply(
     X   = s02,
     FUN = function(s02_j) {
          dnorm(
               x    = y_grid,
               mean = m0,
               sd   = sqrt(sigma2 + s02_j)
          )
     }
)

# Densidades aproximadas mediante simulación
density_predictive_mc <- lapply(
     X   = y_prior,
     FUN = function(y_j) {
          density(
               x    = y_j,
               from = min(y_grid),
               to   = max(y_grid),
               n    = length(y_grid)
          )
     }
)

# Configuración gráfica
scenario_col <- c("#0072B2", "#D55E00")

max_density_mc <- max(
     vapply(
          X         = density_predictive_mc,
          FUN       = function(density_j) max(density_j$y),
          FUN.VALUE = numeric(1)
     )
)

ylim_density <- c(
     0,
     max(
          density_base,
          density_predictive,
          max_density_mc
     )
)

# Distribuciones base ----------------------------------------------------------

pdf(
     file      = "dpm_base_distribution.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

matplot(
     x    = theta_grid,
     y    = density_base,
     type = "l",
     lty  = c(1, 2),
     lwd  = 3,
     col  = scenario_col,
     ylim = ylim_density,
     xlab = expression(theta),
     ylab = "Densidad",
     main = ""
)

legend(
     x      = "topright",
     legend = scenario,
     lty    = c(1, 2),
     lwd    = 3,
     col    = scenario_col,
     bty    = "n"
)

dev.off()

# Distribuciones predictivas previas -------------------------------------------

pdf(
     file      = "dpm_prior_predictive_density.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

matplot(
     x    = y_grid,
     y    = density_predictive,
     type = "l",
     lty  = 1,
     lwd  = 3,
     col  = scenario_col,
     ylim = ylim_density,
     xlab = expression(y),
     ylab = "Densidad",
     main = ""
)

lines(
     x   = density_predictive_mc[[1]]$x,
     y   = density_predictive_mc[[1]]$y,
     lty = 3,
     lwd = 2,
     col = scenario_col[1]
)

lines(
     x   = density_predictive_mc[[2]]$x,
     y   = density_predictive_mc[[2]]$y,
     lty = 3,
     lwd = 2,
     col = scenario_col[2]
)

legend(
     x = "topright",
     legend = c(
          "Concentrada: analítica",
          "Concentrada: simulada",
          "Dispersa: analítica",
          "Dispersa: simulada"
     ),
     lty = c(1, 3, 1, 3),
     lwd = c(3, 2, 3, 2),
     col = c(
          scenario_col[1],
          scenario_col[1],
          scenario_col[2],
          scenario_col[2]
     ),
     bty = "n",
     cex = 0.82
)

dev.off()

# Fin --------------------------------------------------------------------------