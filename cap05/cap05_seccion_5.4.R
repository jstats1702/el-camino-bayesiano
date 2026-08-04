# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Cargar datos
dat <- read.csv(
     file             = file.path(dir_work, "Examen_Saber_11_20251.txt"),
     sep              = ";",
     stringsAsFactors = FALSE
)

dim(dat)

# Revisión inicial por departamento
table(dat$estu_depto_reside, useNA = "ifany")

round(100 * table(dat$estu_depto_reside, useNA = "ifany") / nrow(dat), 1)

# Limpiar variable de departamento
dat$estu_depto_reside <- trimws(dat$estu_depto_reside)

# Remover registros sin departamento y registros del exterior
dat <- dat[
     !is.na(dat$estu_depto_reside) &
          dat$estu_depto_reside != "" &
          !dat$estu_depto_reside %in% c("EXTRANJERO", "EXTRANGERO"),
     ,
     drop = FALSE
]

# Verificar departamentos restantes
table(dat$estu_depto_reside)
unique(dat$estu_depto_reside)

# Puntajes de Matemáticas para Bogotá
y_bogota <- dat[
     dat$estu_depto_reside == "BOGOTÁ",
     "punt_matematicas"
]

# Remover valores faltantes
y_bogota <- as.numeric(y_bogota)
y_bogota <- y_bogota[is.finite(y_bogota)]

# Tamaño de muestra y media
length(y_bogota)
round(mean(y_bogota), 3)

# Muestreo proporcional por departamento --------------------------------------

prop_muestra  <- 0.05
departamentos <- sort(unique(dat$estu_depto_reside))

set.seed(123)

idx_muestra <- integer(0)

for (depto in departamentos) {
     idx_depto     <- which(dat$estu_depto_reside == depto)
     n_depto_total <- length(idx_depto)
     
     if (n_depto_total < 100) {
          n_depto <- min(10, n_depto_total)
     } else {
          n_depto <- max(1, floor(n_depto_total * prop_muestra))
     }
     
     idx_muestra <- c(
          idx_muestra,
          sample(idx_depto, size = n_depto, replace = FALSE)
     )
}

# Base final
dat_muestra <- dat[idx_muestra, , drop = FALSE]

dim(dat_muestra)
table(dat_muestra$estu_depto_reside)

# Datos de Bogotá --------------------------------------------------------------

dat_bogota <- dat_muestra[
     dat_muestra$estu_depto_reside == "BOGOTÁ",
     ,
     drop = FALSE
]

# Puntaje de matemáticas
y <- dat_bogota$punt_matematicas
y <- as.numeric(y)
y <- y[!is.na(y)]

# Resumen de los datos
summary(y)

# Tamaño de muestra
(n <- length(y))

# Estadísticos suficientes
(yb <- mean(y))
(s2 <- var(y))

round(yb, 3)
round(s2, 3)
round(sqrt(s2), 3)

# Visualización ----------------------------------------------------------------

pdf(
     file      = "matematicas_bogota_histograma.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Densidades para ajustar los límites del gráfico
den_kernel <- density(y)
x_grid     <- seq(min(y), max(y), length.out = 10000)
den_normal <- dnorm(x_grid, mean = yb, sd = sqrt(s2))

ylim <- range(c(den_kernel$y, den_normal), finite = TRUE)

# Histograma en escala de densidad
hist(
     y,
     probability = TRUE,
     col         = "gray85",
     border      = "white",
     breaks      = "FD",
     main        = "",
     xlab        = "Puntaje de Matemáticas",
     ylab        = "Densidad",
     ylim        = ylim
)

# Estimación kernel de la densidad
lines(den_kernel, lwd = 2)

# Curva Normal con media y desviación estándar muestrales
lines(x_grid, den_normal, lwd = 2, lty = 2)

# Leyenda
legend(
     "topright",
     legend = c("Kernel", "Normal"),
     lwd    = c(2, 2),
     lty    = c(1, 2),
     bty    = "n"
)

dev.off()

# Hiperparámetros --------------------------------------------------------------

mu0 <- 50
k0  <- 1
nu0 <- 3
s20 <- ((nu0 - 2) / nu0) * 10^2

# Distribución posterior -------------------------------------------------------

kn  <- k0 + n
nun <- nu0 + n
mun <- (k0 / kn) * mu0 + (n / kn) * yb
s2n <- (nu0 * s20 + (n - 1) * s2 + (k0 * n / kn) * (yb - mu0)^2) / nun

round(kn, 4)
round(nun, 4)
round(mun, 4)
round(s2n, 4)

# Función de densidad de una Gamma-Inversa ------------------------------------

dinvgamma0 <- function(x, a, b, log = FALSE) {
     log_density         <- a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x
     log_density[x <= 0] <- -Inf
     
     if (log) {
          return(log_density)
     }
     
     exp(log_density)
}

# Función de densidad de una t -------------------------------------------------

dt0 <- function(x, nu, mu, sigma2, log = FALSE) {
     if (nu <= 0) stop("nu debe ser positivo.")
     if (sigma2 <= 0) stop("sigma2 debe ser positivo.")
     
     # Constante de normalización
     log_const <- lgamma((nu + 1) / 2) -
          lgamma(nu / 2) -
          0.5 * log(nu * pi * sigma2)
     
     # Término dependiente de x
     log_kernel <- -((nu + 1) / 2) *
          log1p((x - mu)^2 / (nu * sigma2))
     
     # Densidad en escala logarítmica
     log_density <- log_const + log_kernel
     
     if (log) {
          return(log_density)
     }
     
     exp(log_density)
}

# Muestras de la distribución posterior ---------------------------------------

# Número de muestras
B <- 10000

set.seed(123)

# Muestras de sigma^2
sigma2_mc <- 1 / rgamma(
     B,
     shape = nun / 2,
     rate  = nun * s2n / 2
)

# Muestras de theta
theta_mc <- rnorm(
     B,
     mean = mun,
     sd   = sqrt(sigma2_mc / kn)
)

# Visualización de muestras posteriores ---------------------------------------

pdf(
     file      = "matematicas_bogota_posterior_muestras.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Límites de los ejes
xlim_theta  <- range(theta_mc)
ylim_sigma2 <- range(sigma2_mc)

# Distribución conjunta
plot(
     theta_mc,
     sigma2_mc,
     pch  = 16,
     cex  = 0.3,
     col  = adjustcolor("black", alpha.f = 0.3),
     xlim = xlim_theta,
     ylim = ylim_sigma2,
     xlab = expression(theta),
     ylab = expression(sigma^2),
     main = ""
)

dev.off()

# Evaluación de la densidad posterior conjunta --------------------------------

g <- 25

theta_grid <- seq(
     from       = xlim_theta[1],
     to         = xlim_theta[2],
     length.out = g
)

sigma2_grid <- seq(
     from       = ylim_sigma2[1],
     to         = ylim_sigma2[2],
     length.out = g
)

# Evaluar y almacenar la distribución posterior en escala log
lp <- matrix(NA_real_, nrow = g, ncol = g)

for (i in seq_len(g)) {
     for (j in seq_len(g)) {
          lp[i, j] <- dnorm(
               theta_grid[i],
               mean = mun,
               sd   = sqrt(sigma2_grid[j] / kn),
               log  = TRUE
          ) +
               dinvgamma0(
                    sigma2_grid[j],
                    a   = nun / 2,
                    b   = nun * s2n / 2,
                    log = TRUE
               )
     }
}

# Visualización de la distribución posterior ----------------------------------

pdf(
     file      = "matematicas_bogota_posterior_densidad.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1, 1),
     mgp   = c(1.75, 0.75, 0)
)

# Gráfico 3D
persp(
     x        = theta_grid,
     y        = sigma2_grid,
     z        = exp(lp),
     theta    = 30,
     phi      = 30,
     expand   = 1,
     cex.axis = 0.5,
     cex.lab  = 0.5,
     xlab     = "Media",
     ylab     = "Varianza",
     zlab     = "Densidad",
     col      = "gray95",
     ticktype = "detailed",
     main     = ""
)

dev.off()

# Distribución marginal de theta -----------------------------------------------

pdf(
     file      = "matematicas_bogota_posterior_theta.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

x_theta <- seq(
     from       = min(theta_mc),
     to         = max(theta_mc),
     length.out = 10000
)

y_theta <- dt0(
     x      = x_theta,
     nu     = nun,
     mu     = mun,
     sigma2 = s2n / kn,
     log    = FALSE
)

# Distribución marginal de theta
hist(
     theta_mc,
     freq   = FALSE,
     col    = "gray85",
     border = "white",
     breaks = 25,
     xlab   = expression(theta),
     ylab   = "Densidad",
     main   = "",
     ylim   = c(0, max(y_theta, na.rm = TRUE))
)

lines(x_theta, y_theta, col = 4, lwd = 2)

legend(
     "topright",
     legend = c("MC", "Exacta"),
     bty    = "n",
     lwd    = c(NA, 2),
     pch    = c(15, NA),
     col    = c("gray85", 4),
     pt.cex = 2
)

dev.off()

# Distribución marginal de sigma^2 --------------------------------------------

pdf(
     file      = "matematicas_bogota_posterior_sigma2.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

x_sigma2 <- seq(
     from       = min(sigma2_mc),
     to         = max(sigma2_mc),
     length.out = 10000
)

y_sigma2 <- dinvgamma0(
     x_sigma2,
     a = nun / 2,
     b = nun * s2n / 2
)

# Distribución marginal de sigma^2
hist(
     sigma2_mc,
     freq   = FALSE,
     col    = "gray85",
     border = "white",
     breaks = 25,
     xlab   = expression(sigma^2),
     ylab   = "Densidad",
     main   = "",
     ylim   = c(0, max(y_sigma2, na.rm = TRUE))
)

lines(x_sigma2, y_sigma2, col = 4, lwd = 2)

legend(
     "topright",
     legend = c("MC", "Exacta"),
     bty    = "n",
     lwd    = c(NA, 2),
     pch    = c(15, NA),
     col    = c("gray85", 4),
     pt.cex = 2
)

dev.off()

# Distribución predictiva posterior -------------------------------------------

set.seed(1234)

y_new <- rnorm(
     n    = B,
     mean = theta_mc,
     sd   = sqrt(sigma2_mc)
)

# Función auxiliar para resumir simulaciones ----------------------------------

resumir_mc <- function(x) {
     c(
          estimacion = mean(x),
          CV         = sd(x) / mean(x),
          lim_inf_95 = quantile(x, 0.025),
          lim_sup_95 = quantile(x, 0.975)
     )
}

# Tabla de inferencia ----------------------------------------------------------

tab <- rbind(
     "Media"  = resumir_mc(theta_mc),
     "DE"     = resumir_mc(sqrt(sigma2_mc)),
     "y pred" = resumir_mc(y_new)
)

# Ajustar nombres y redondear
colnames(tab) <- c("Estimación", "CV", "L. Inf. 95%", "L. Sup. 95%")
tab <- round(tab, 3)
tab

# Funciones auxiliares para calcular la asimetría y la curtosis ---------------

asimetria <- function(x) {
     x <- x[is.finite(x)]
     mean((x - mean(x))^3) / sd(x)^3
}

curtosis <- function(x) {
     x <- x[is.finite(x)]
     mean((x - mean(x))^4) / sd(x)^4
}

# Estadísticos observados ------------------------------------------------------

t_obs <- c(
     media     = mean(y),
     de        = sd(y),
     mediana   = median(y),
     riq       = IQR(y),
     asimetria = asimetria(y),
     curtosis  = curtosis(y)
)

round(t_obs, 3)

# Inicializar matriz para almacenar estadísticas de prueba ---------------------

t_mc <- matrix(
     NA_real_,
     nrow = B,
     ncol = length(t_obs)
)

colnames(t_mc) <- names(t_obs)

# Distribución predictiva posterior -------------------------------------------

set.seed(1234)

for (i in seq_len(B)) {
     # Datos simulados
     y_rep <- rnorm(
          n    = n,
          mean = theta_mc[i],
          sd   = sqrt(sigma2_mc[i])
     )
     
     # Estadísticos de prueba
     t_mc[i, ] <- c(
          media     = mean(y_rep),
          de        = sd(y_rep),
          mediana   = median(y_rep),
          riq       = IQR(y_rep),
          asimetria = asimetria(y_rep),
          curtosis  = curtosis(y_rep)
     )
}

# ppp -------------------------------------------------------------------------

ppp <- colMeans(
     sweep(t_mc, 2, t_obs, FUN = "<")
)

round(ppp, 4)

# Visualización de estadísticas de prueba -------------------------------------

for (j in seq_along(t_obs)) {
     pdf(
          file      = paste0("matematicas_bogota_chequeo_", names(t_obs)[j], ".pdf"),
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     # Histograma predictivo posterior del estadístico de prueba
     hist(
          t_mc[, j],
          freq   = FALSE,
          col    = "gold2",
          border = "gold2",
          xlab   = "t",
          ylab   = expression(p(t ~ "|" ~ bold(y))),
          main   = ""
     )
     
     # Estadístico observado
     abline(v = t_obs[j], col = 1, lwd = 2, lty = 1)
     
     dev.off()
}

# Análisis de sensibilidad -----------------------------------------------------

# Número de muestras
B <- 10000

# Hiperparámetros
nu0_vec <- 3:6
mu0     <- 50
k0      <- 1

# Funciones auxiliares ---------------------------------------------------------

momentos_previos_sigma2 <- function(nu0, s20) {
     media <- nu0 * s20 / (nu0 - 2)
     
     varianza <- if (nu0 > 4) {
          2 * nu0^2 * s20^2 / ((nu0 - 2)^2 * (nu0 - 4))
     } else {
          Inf
     }
     
     c(
          media_prior_sigma2 = media,
          var_prior_sigma2   = varianza
     )
}

# Simulación para cada valor de nu0 -------------------------------------------

tab_list <- list()

set.seed(123)

for (nu0 in nu0_vec) {
     # Hiperparámetro que fija E(sigma^2) = 10^2
     s20 <- ((nu0 - 2) / nu0) * 10^2
     
     momentos_prior <- momentos_previos_sigma2(nu0, s20)
     
     # Hiperparámetros posteriores
     kn  <- k0 + n
     nun <- nu0 + n
     mun <- (k0 / kn) * mu0 + (n / kn) * yb
     s2n <- (nu0 * s20 + (n - 1) * s2 + (k0 * n / kn) * (yb - mu0)^2) / nun
     
     # Muestras posteriores
     sigma2_mc <- 1 / rgamma(
          B,
          shape = nun / 2,
          rate  = nun * s2n / 2
     )
     
     theta_mc <- rnorm(
          B,
          mean = mun,
          sd   = sqrt(sigma2_mc / kn)
     )
     
     sigma_mc <- sqrt(sigma2_mc)
     
     # Resúmenes posteriores
     tab_nu0 <- rbind(
          theta = resumir_mc(theta_mc),
          sigma = resumir_mc(sigma_mc)
     )
     
     tab_nu0 <- data.frame(
          nu0                = nu0,
          cantidad           = rownames(tab_nu0),
          media_prior_sigma2 = momentos_prior["media_prior_sigma2"],
          var_prior_sigma2   = momentos_prior["var_prior_sigma2"],
          tab_nu0,
          row.names          = NULL
     )
     
     tab_list[[as.character(nu0)]] <- tab_nu0
}

# Tabla final ------------------------------------------------------------------

tab_sensibilidad <- do.call(rbind, tab_list)

colnames(tab_sensibilidad) <- c(
     "nu0",
     "cantidad",
     "media_prior_sigma2",
     "var_prior_sigma2",
     "Estimación",
     "CV",
     "L. Inf. 95%",
     "L. Sup. 95%"
)

# Aproximar resultados a 3 decimales
cols_redondear <- c(
     "Estimación",
     "CV",
     "L. Inf. 95%",
     "L. Sup. 95%"
)

tab_sensibilidad[, cols_redondear] <- round(
     tab_sensibilidad[, cols_redondear],
     3
)

tab_sensibilidad

# Fin --------------------------------------------------------------------------