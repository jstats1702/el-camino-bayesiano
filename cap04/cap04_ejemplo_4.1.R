# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Parámetros de la distribución Gamma
a <- 3
b <- 2

# Tamaños de muestra Monte Carlo
B_values <- c(10, 30, 1000)

# Simulación -------------------------------------------------------------------

set.seed(123)

theta_mc <- vector("list", length(B_values))

for (j in seq_along(B_values)) {
     theta_mc[[j]] <- rgamma(B_values[j], shape = a, rate = b)
}

# Gráficos ---------------------------------------------------------------------

pdf(file = "gamma_monte_carlo.pdf", width = 7.5, height = 5, pointsize = 15)

par(
     mfrow = c(2, 3),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Histogramas con densidad teórica
for (j in seq_along(B_values)) {
     hist(
          x      = theta_mc[[j]],
          prob   = TRUE,
          xlim   = c(0, 6),
          ylim   = c(0, 0.8),
          xlab   = expression(theta),
          ylab   = "Densidad",
          main   = paste0("B = ", B_values[j]),
          col    = "gray90",
          border = "gray90"
     )
     
     curve(
          dgamma(x, shape = a, rate = b),
          col = "blue",
          lwd = 2,
          add = TRUE,
          n   = 1000
     )
}

# Funciones de distribución acumulada
for (j in seq_along(B_values)) {
     x <- sort(theta_mc[[j]])
     B <- length(x)
     
     plot(
          x    = x,
          y    = seq_len(B) / B,
          type = "s",
          col  = "gray60",
          lwd  = 2,
          xlim = c(0, 6),
          ylim = c(0, 1),
          xlab = expression(theta),
          ylab = "Distribución acumulada",
          main = paste0("B = ", B_values[j])
     )
     
     curve(
          pgamma(x, shape = a, rate = b),
          col = "blue",
          lwd = 2,
          add = TRUE,
          n   = 1000
     )
}

dev.off()

# Fin --------------------------------------------------------------------------