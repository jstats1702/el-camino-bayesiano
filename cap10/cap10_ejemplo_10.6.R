# Configuración ----------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

# Tamaño de muestra y valor verdadero del parámetro
n      <- 20
theta0 <- 2

# Datos simulados
y <- rexp(n, rate = theta0)

# Hiperparámetros de la distribución previa Gamma(a0, b0)
a0 <- 1
b0 <- 1

# Parámetros de la distribución posterior exacta
a_post <- a0 + n
b_post <- b0 + sum(y)

# Resumen de la posterior exacta
media_exacta <- a_post / b_post
sd_exacta    <- sqrt(a_post / b_post^2)

IC_exacto <- qgamma(
     p     = c(0.025, 0.975),
     shape = a_post,
     rate  = b_post
)

# Algoritmo de Metropolis--Hastings --------------------------------------------

# Logaritmo del núcleo posterior
log_posterior <- function(theta, y, a0, b0) {
     if (theta <= 0) {
          return(-Inf)
     }
     
     n <- length(y)
     
     (a0 + n - 1) * log(theta) - (b0 + sum(y)) * theta
}

# Logaritmo de la densidad de propuesta lognormal
log_propuesta <- function(theta_nuevo, theta_actual, delta2) {
     dlnorm(
          x       = theta_nuevo,
          meanlog = log(theta_actual),
          sdlog   = sqrt(delta2),
          log     = TRUE
     )
}

mcmc <- function(
          y,
          a0,
          b0,
          B,
          burn_in,
          delta2,
          theta_ini,
          seed = 123
) {
     set.seed(seed)
     
     # Almacenamiento
     theta_mcmc <- numeric(B)
     
     # Valor inicial
     theta <- theta_ini
     
     # Contador de aceptaciones
     acr <- 0
     
     for (b in seq_len(B)) {
          # Propuesta lognormal
          theta_p <- rlnorm(
               n       = 1,
               meanlog = log(theta),
               sdlog   = sqrt(delta2)
          )
          
          # Razón de aceptación en escala logarítmica
          log_post_star <- log_posterior(theta_p, y, a0, b0)
          log_post_curr <- log_posterior(theta, y, a0, b0)
          
          log_prop_curr <- log_propuesta(theta, theta_p, delta2)
          log_prop_star <- log_propuesta(theta_p, theta, delta2)
          
          log_r <- log_post_star - log_post_curr +
               log_prop_curr - log_prop_star
          
          # Aceptar o rechazar
          if (log(runif(1)) <= min(0, log_r)) {
               theta <- theta_p
               acr   <- acr + 1
          }
          
          # Almacenar
          theta_mcmc[b] <- theta
     }
     
     # Muestras posteriores después del período de calentamiento
     theta_post <- theta_mcmc[(burn_in + 1):B]
     
     # Tasa empírica de aceptación
     tasa_aceptacion <- acr / B
     
     list(
          theta_mcmc      = theta_mcmc,
          theta_post      = theta_post,
          tasa_aceptacion = tasa_aceptacion
     )
}

# Ajuste del algoritmo ---------------------------------------------------------

B         <- 11000
burn_in   <- 1000
delta2    <- 0.75
theta_ini <- 5

muestras <- mcmc(
     y         = y,
     a0        = a0,
     b0        = b0,
     B         = B,
     burn_in   = burn_in,
     delta2    = delta2,
     theta_ini = theta_ini,
     seed      = 123
)

# Extraer resultados
theta_mcmc <- muestras$theta_mcmc
theta_post <- muestras$theta_post
acr        <- muestras$tasa_aceptacion

round(acr, 3)

# Cadena -----------------------------------------------------------------------

# Traza de la cadena
pdf(
     file      = "mh_exponencial_juguete_traceplot.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = seq_along(theta_mcmc),
     y    = theta_mcmc,
     type = "p",
     pch  = 16,
     cex  = 0.25,
     col  = adjustcolor("black", 0.5),
     xlab = "Iteración posterior",
     ylab = expression(theta)
)

box()

dev.off()

# Histograma posterior ---------------------------------------------------------

pdf(
     file      = "mh_exponencial_juguete_histograma.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     x      = theta_post,
     freq   = FALSE,
     col    = "gray80",
     border = "white",
     main   = "",
     xlab   = expression(theta),
     ylab   = "Densidad"
)

# Posterior analítica
theta_grid <- seq(
     from       = min(theta_post),
     to         = max(theta_post),
     length.out = 1000
)

lines(
     x   = theta_grid,
     y   = dgamma(theta_grid, shape = a_post, rate = b_post),
     col = 2,
     lty = 1,
     lwd = 2
)

box()

dev.off()

# Sensibilidad a la escala de la propuesta -------------------------------------

delta2_grid <- seq(
     from = 0.01,
     to   = 5,
     by   = 0.01
)

resumen_delta2 <- data.frame(
     delta2     = delta2_grid,
     aceptacion = NA,
     media      = NA,
     sd         = NA,
     acf1       = NA
)

for (j in seq_along(delta2_grid)) {
     muestras_j <- mcmc(
          y         = y,
          a0        = a0,
          b0        = b0,
          B         = B,
          burn_in   = burn_in,
          delta2    = delta2_grid[j],
          theta_ini = theta_ini,
          seed      = 123
     )
     
     theta_post_j <- muestras_j$theta_post
     
     resumen_delta2$aceptacion[j] <- muestras_j$tasa_aceptacion
     resumen_delta2$media[j]      <- mean(theta_post_j)
     resumen_delta2$sd[j]         <- sd(theta_post_j)
     resumen_delta2$acf1[j]       <- acf(
          theta_post_j,
          plot = FALSE
     )$acf[2]
}

# Gráfico de sensibilidad: tasa de aceptación ----------------------------------

pdf(
     file      = "mh_exponencial_sensibilidad_delta2_aceptacion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = resumen_delta2$delta2,
     y    = resumen_delta2$aceptacion,
     type = "p",
     pch  = 16,
     cex  = 0.3,
     col  = 1,
     ylim = c(0, 1),
     xlab = expression(delta^2),
     ylab = "Tasa de aceptación",
     main = ""
)

box()

dev.off()

# Gráfico de sensibilidad: autocorrelación de orden 1 --------------------------

pdf(
     file      = "mh_exponencial_sensibilidad_delta2_acf1.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = resumen_delta2$delta2,
     y    = resumen_delta2$acf1,
     type = "p",
     pch  = 16,
     cex  = 0.3,
     col  = 1,
     ylim = c(0, 1),
     xlab = expression(delta^2),
     ylab = "Autocorrelación de orden 1",
     main = ""
)

box()

dev.off()

# Fin --------------------------------------------------------------------------