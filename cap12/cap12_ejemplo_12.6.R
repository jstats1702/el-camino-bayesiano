# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(stats)))

alpha   <- 5
n_small <- 10
n_large <- 100

x_grid <- seq(
     from       = -5,
     to         = 5,
     length.out = 2000
)

# Generación de las muestras ---------------------------------------------------

set.seed(123)

y_large <- rcauchy(
     n        = n_large,
     location = 0,
     scale    = 1
)

# La muestra pequeña corresponde a las primeras diez observaciones
y_small <- y_large[seq_len(n_small)]

# Funciones de distribución ----------------------------------------------------

# Función de distribución verdadera
F_true <- function(x) {
     pcauchy(
          q        = x,
          location = 0,
          scale    = 1
     )
}

# Distribución base previa
G0 <- function(x) {
     pnorm(
          q    = x,
          mean = 0,
          sd   = 1
     )
}

# Funciones de distribución empíricas
G_emp_small <- ecdf(y_small)
G_emp_large <- ecdf(y_large)

# Medias posteriores
G_post_small <- function(x) {
     (alpha * G0(x) + n_small * G_emp_small(x)) / (alpha + n_small)
}

G_post_large <- function(x) {
     (alpha * G0(x) + n_large * G_emp_large(x)) / (alpha + n_large)
}

# Pesos de la distribución previa y de la distribución empírica ----------------

weights <- data.frame(
     n                = c(n_small, n_large),
     weight_G0        = alpha / (alpha + c(n_small, n_large)),
     weight_empirical = c(n_small, n_large) /
          (alpha + c(n_small, n_large))
)

round(weights, 3)

# Función para construir cada gráfico ------------------------------------------

plot_cdf_comparison <- function(
          G_emp,
          G_post,
          show_legend = FALSE
) {
     plot(
          x    = x_grid,
          y    = F_true(x_grid),
          type = "l",
          lty  = 1,
          lwd  = 2.5,
          col  = "black",
          ylim = c(0, 1),
          xlab = expression(x),
          ylab = "Función de distribución",
          main = ""
     )
     
     lines(
          x    = x_grid,
          y    = G0(x_grid),
          type = "l",
          lty  = 1,
          lwd  = 2,
          col  = "gray60"
     )
     
     lines(
          x    = x_grid,
          y    = G_emp(x_grid),
          type = "s",
          lty  = 1,
          lwd  = 1.75,
          col  = "royalblue3"
     )
     
     lines(
          x    = x_grid,
          y    = G_post(x_grid),
          type = "l",
          lty  = 1,
          lwd  = 2.5,
          col  = "firebrick3"
     )
     
     if (show_legend) {
          legend(
               "topleft",
               legend = c(
                    "Cauchy(0, 1)",
                    "G0 = N(0, 1)",
                    "Distribución empírica",
                    "Media posterior"
               ),
               col = c(
                    "black",
                    "gray60",
                    "royalblue3",
                    "firebrick3"
               ),
               lty = 1,
               lwd = c(2.5, 2, 1.75, 2.5),
               bty = "n",
               cex = 0.85
          )
     }
}

# Figura para n = 10 -----------------------------------------------------------

pdf(
     file      = "dp_cauchy_cdf_posterior_n10.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot_cdf_comparison(
     G_emp       = G_emp_small,
     G_post      = G_post_small,
     show_legend = TRUE
)

dev.off()

# Figura para n = 100 ----------------------------------------------------------

pdf(
     file      = "dp_cauchy_cdf_posterior_n100.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot_cdf_comparison(
     G_emp       = G_emp_large,
     G_post      = G_post_large,
     show_legend = FALSE
)

dev.off()

# Fin --------------------------------------------------------------------------