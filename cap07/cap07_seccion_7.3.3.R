# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Paquetes
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(sf)))
suppressMessages(suppressWarnings(library(stringi)))

# Cargar datos
dat <- read.csv(
     file             = file.path("Examen_Saber_11_20251.txt"),
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
dat <- dat[idx_muestra, , drop = FALSE]

dim(dat)
table(dat$estu_depto_reside)

# Estadísticos suficientes -----------------------------------------------------

# m : número de grupos (departamentos)
deptos <- sort(unique(dat$estu_cod_reside_depto))
m      <- length(deptos)

# n : número de individuos (estudiantes)
n <- nrow(dat)

# Tratamiento de datos
# y  : puntaje de los estudiantes (c)
# Y  : puntaje de los estudiantes por departamento (list)
# g  : identificador secuencial de los departamentos (c)
# nj : número de estudiantes por departamento (c)
# yb : promedios por departamento (c)
# s2 : varianzas por departamento (c)

y <- dat$punt_matematicas
Y <- vector(mode = "list", length = m)
g <- rep(NA_integer_, n)

for (j in seq_len(m)) {
     idx    <- dat$estu_cod_reside_depto == deptos[j]
     g[idx] <- j
     Y[[j]] <- y[idx]
}

# Tabla
estadisticos <- dat %>%
     group_by(estu_cod_reside_depto) %>%
     summarise(
          codigo  = first(estu_cod_reside_depto),
          nombre  = first(estu_depto_reside),
          nj      = dplyr::n(),
          yb      = mean(punt_matematicas, na.rm = TRUE),
          med     = median(punt_matematicas, na.rm = TRUE),
          s2      = var(punt_matematicas, na.rm = TRUE),
          s       = sd(punt_matematicas, na.rm = TRUE),
          min     = min(punt_matematicas, na.rm = TRUE),
          max     = max(punt_matematicas, na.rm = TRUE),
          .groups = "drop"
     ) %>%
     arrange(codigo) %>%
     select(codigo, nombre, nj, yb, med, s2, s, min, max)

as.data.frame(estadisticos[, c(2, 3, 4, 5, 7, 8, 9)])

# Vectores resumen
nj <- estadisticos$nj
yb <- estadisticos$yb
s2 <- estadisticos$s2

# Promedio global
round(mean(y), 3)

# Distribución previa modelo homocedástico -------------------------------------

# Hiperparámetros
mu0 <- 50
g20 <- 10^2

eta0 <- 1
t20  <- 10^2

nu0 <- 1
s20 <- 10^2

# Muestreador de Gibbs modelo homocedástico ------------------------------------

mcmc <- function(B, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20) {
     # Frecuencia para imprimir progreso
     ncat <- max(1, floor(B / 10))
     
     # Número total de observaciones y grupos
     n <- sum(nj)
     m <- length(nj)
     
     # Almacenamiento de la cadena
     THETA <- matrix(NA_real_, nrow = B, ncol = m + 4)
     
     # Valores iniciales
     theta <- yb
     sig2  <- mean(s2, na.rm = TRUE)
     mu    <- mean(theta, na.rm = TRUE)
     tau2  <- var(theta, na.rm = TRUE)
     
     # Cadena de MCMC
     for (b in seq_len(B)) {
          # Actualizar theta_j
          v_theta <- 1 / (1 / tau2 + nj / sig2)
          m_theta <- v_theta * (mu / tau2 + nj * yb / sig2)
          theta   <- rnorm(n = m, mean = m_theta, sd = sqrt(v_theta))
          
          # Actualizar sigma^2
          ss_sig2 <- sum((nj - 1) * s2 + nj * (yb - theta)^2)
          a_sig2  <- 0.5 * (nu0 + n)
          b_sig2  <- 0.5 * (nu0 * s20 + ss_sig2)
          sig2    <- 1 / rgamma(n = 1, shape = a_sig2, rate = b_sig2)
          
          # Actualizar mu
          v_mu <- 1 / (1 / g20 + m / tau2)
          m_mu <- v_mu * (mu0 / g20 + sum(theta) / tau2)
          mu   <- rnorm(n = 1, mean = m_mu, sd = sqrt(v_mu))
          
          # Actualizar tau^2
          ss_tau2 <- sum((theta - mu)^2)
          a_tau2  <- 0.5 * (eta0 + m)
          b_tau2  <- 0.5 * (eta0 * t20 + ss_tau2)
          tau2    <- 1 / rgamma(n = 1, shape = a_tau2, rate = b_tau2)
          
          # Log-verosimilitud usando estadísticos suficientes por grupo
          ll <- sum(
               -0.5 * nj * log(2 * pi * sig2) -
                    0.5 * ((nj - 1) * s2 + nj * (yb - theta)^2) / sig2
          )
          
          # Almacenar iteración
          THETA[b, ] <- c(theta, sig2, mu, tau2, ll)
          
          # Imprimir progreso
          if (b %% ncat == 0) {
               cat(100 * round(b / B, 1), "% completado ... \n", sep = "")
          }
     }
     
     # Salida
     colnames(THETA) <- c(paste0("theta", seq_len(m)), "sig2", "mu", "tau2", "ll")
     THETA <- as.data.frame(THETA)
     
     return(list(THETA = THETA))
}

# Ajuste del modelo homocedástico ----------------------------------------------

set.seed(123)
chain_1 <- mcmc(B = 10000, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20)

# Distribución previa ----------------------------------------------------------

# Hiperparámetros
mu0 <- 50
g20 <- 10^2

eta0 <- 1
t20  <- 10^2

lam0 <- 1

al0 <- 1
be0 <- 1 / 10^2

nus0 <- 1:50  # Rango para p(nu | rest)

# Muestreador de Gibbs modelo heterocedástico ----------------------------------

mcmc <- function(B, y, nj, yb, s2, mu0, g20, eta0, t20, lam0, al0, be0, nus0) {
     # Frecuencia para imprimir progreso
     ncat <- max(1, floor(B / 10))
     
     # Número total de observaciones y grupos
     n <- sum(nj)
     m <- length(nj)
     
     # Corregir varianzas no definidas o nulas
     s2 <- ifelse(is.na(s2), 0, s2)
     s2 <- pmax(s2, .Machine$double.eps)
     
     # Almacenamiento de la cadena
     THETA <- matrix(NA_real_, nrow = B, ncol = 2 * m + 5)
     
     # Valores iniciales
     theta <- yb
     sig2  <- s2
     mu    <- mean(theta, na.rm = TRUE)
     tau2  <- max(var(theta, na.rm = TRUE), .Machine$double.eps)
     nu    <- min(nus0)
     ups2  <- 100
     
     # Cadena MCMC
     for (b in seq_len(B)) {
          # Actualizar theta_j
          v_theta <- 1 / (1 / tau2 + nj / sig2)
          m_theta <- v_theta * (mu / tau2 + nj * yb / sig2)
          theta   <- rnorm(n = m, mean = m_theta, sd = sqrt(v_theta))
          
          # Actualizar sigma_j^2
          ss_sig2 <- (nj - 1) * s2 + nj * (yb - theta)^2
          a_sig2  <- 0.5 * (nu + nj)
          b_sig2  <- 0.5 * (nu * ups2 + ss_sig2)
          sig2    <- 1 / rgamma(n = m, shape = a_sig2, rate = b_sig2)
          
          # Actualizar mu
          v_mu <- 1 / (1 / g20 + m / tau2)
          m_mu <- v_mu * (mu0 / g20 + sum(theta) / tau2)
          mu   <- rnorm(n = 1, mean = m_mu, sd = sqrt(v_mu))
          
          # Actualizar tau^2
          ss_tau2 <- sum((theta - mu)^2)
          a_tau2  <- 0.5 * (eta0 + m)
          b_tau2  <- 0.5 * (eta0 * t20 + ss_tau2)
          tau2    <- 1 / rgamma(n = 1, shape = a_tau2, rate = b_tau2)
          
          # Actualizar nu
          lpnu <- 0.5 * m * nus0 * log(0.5 * nus0 * ups2) -
               m * lgamma(0.5 * nus0) -
               0.5 * nus0 * sum(log(sig2)) -
               nus0 * (lam0 + 0.5 * ups2 * sum(1 / sig2))
          
          pnu <- exp(lpnu - max(lpnu))
          pnu <- pnu / sum(pnu)
          
          nu <- sample(
               x    = nus0,
               size = 1,
               prob = pnu
          )
          
          # Actualizar ups2
          a_ups2 <- al0 + 0.5 * m * nu
          b_ups2 <- be0 + 0.5 * nu * sum(1 / sig2)
          ups2   <- rgamma(n = 1, shape = a_ups2, rate = b_ups2)
          
          # Log-verosimilitud usando estadísticos suficientes por grupo
          ll <- sum(
               -0.5 * nj * log(2 * pi * sig2) -
                    0.5 * ((nj - 1) * s2 + nj * (yb - theta)^2) / sig2
          )
          
          # Almacenar iteración
          THETA[b, ] <- c(theta, sig2, mu, tau2, nu, ups2, ll)
          
          # Imprimir progreso
          if (b %% ncat == 0) {
               cat(100 * round(b / B, 1), "% completado ... \n", sep = "")
          }
     }
     
     # Salida
     colnames(THETA) <- c(
          paste0("theta", seq_len(m)),
          paste0("sig2", seq_len(m)),
          "mu", "tau2", "nu", "ups2", "ll"
     )
     
     THETA <- as.data.frame(THETA)
     
     return(list(THETA = THETA))
}

# Ajuste del modelo heterocedástico --------------------------------------------

set.seed(123)
chain_2 <- mcmc(B = 10000, y, nj, yb, s2, mu0, g20, eta0, t20, lam0, al0, be0, nus0)

# Cadenas de la log-verosimilitud ----------------------------------------------

# Comparación de la log-verosimilitud de los modelos
round(mean(chain_1$THETA$ll), 3)
round(mean(chain_2$THETA$ll), 3)

# DIC modelo homocedástico -----------------------------------------------------

theta_hat  <- colMeans(chain_1$THETA[, seq_len(m), drop = FALSE])
sigma2_hat <- mean(chain_1$THETA$sig2)

lpy_m1 <- sum(
     dnorm(
          x    = y,
          mean = theta_hat[g],
          sd   = sqrt(sigma2_hat),
          log  = TRUE
     )
)

pDIC_m1 <- 2 * (lpy_m1 - mean(chain_1$THETA$ll))
dic_m1  <- -2 * lpy_m1 + 2 * pDIC_m1

tab_DIC_homo <- data.frame(
     Modelo   = "Homocedástico",
     lpy      = lpy_m1,
     pDIC     = pDIC_m1,
     DIC      = dic_m1,
     row.names = NULL
)

tab_DIC_homo

# DIC modelo heterocedástico ---------------------------------------------------

theta_hat  <- colMeans(chain_2$THETA[, seq_len(m), drop = FALSE])
sigma2_hat <- colMeans(chain_2$THETA[, m + seq_len(m), drop = FALSE])

lpy_m2 <- sum(
     dnorm(
          x    = y,
          mean = theta_hat[g],
          sd   = sqrt(sigma2_hat[g]),
          log  = TRUE
     )
)

pDIC_m2 <- 2 * (lpy_m2 - mean(chain_2$THETA$ll))
dic_m2  <- -2 * lpy_m2 + 2 * pDIC_m2

tab_DIC_hetero <- data.frame(
     Modelo   = "Heterocedástico",
     lpy      = lpy_m2,
     pDIC     = pDIC_m2,
     DIC      = dic_m2,
     row.names = NULL
)

tab_DIC_hetero

# WAIC modelo homocedástico ----------------------------------------------------

lppd_i_m1  <- numeric(n)
pWAIC_i_m1 <- numeric(n)

for (i in seq_len(n)) {
     log_lik_i <- dnorm(
          x    = y[i],
          mean = chain_1$THETA[, g[i]],
          sd   = sqrt(chain_1$THETA$sig2),
          log  = TRUE
     )
     
     # Cálculo estable de log(mean(exp(log_lik_i)))
     max_log_i <- max(log_lik_i)
     
     lppd_i_m1[i] <- max_log_i + log(mean(exp(log_lik_i - max_log_i)))
     pWAIC_i_m1[i] <- 2 * (lppd_i_m1[i] - mean(log_lik_i))
}

lppd_m1  <- sum(lppd_i_m1)
pWAIC_m1 <- sum(pWAIC_i_m1)
waic_m1  <- -2 * lppd_m1 + 2 * pWAIC_m1

tab_WAIC_homo <- data.frame(
     Modelo   = "Homocedástico",
     lppd     = lppd_m1,
     pWAIC    = pWAIC_m1,
     WAIC     = waic_m1,
     row.names = NULL
)

tab_WAIC_homo

# WAIC modelo heterocedástico --------------------------------------------------

lppd_i_m2  <- numeric(n)
pWAIC_i_m2 <- numeric(n)

for (i in seq_len(n)) {
     log_lik_i <- dnorm(
          x    = y[i],
          mean = chain_2$THETA[, g[i]],
          sd   = sqrt(chain_2$THETA[, m + g[i]]),
          log  = TRUE
     )
     
     # Cálculo estable de log(mean(exp(log_lik_i)))
     max_log_i <- max(log_lik_i)
     
     lppd_i_m2[i] <- max_log_i + log(mean(exp(log_lik_i - max_log_i)))
     pWAIC_i_m2[i] <- 2 * (lppd_i_m2[i] - mean(log_lik_i))
}

lppd_m2  <- sum(lppd_i_m2)
pWAIC_m2 <- sum(pWAIC_i_m2)
waic_m2  <- -2 * lppd_m2 + 2 * pWAIC_m2

tab_WAIC_hetero <- data.frame(
     Modelo   = "Heterocedástico",
     lppd     = lppd_m2,
     pWAIC    = pWAIC_m2,
     WAIC     = waic_m2,
     row.names = NULL
)

tab_WAIC_hetero

# Consolidación de los criterios de información -------------------------------

tab_modelos <- data.frame(
     Modelo   = c("Homocedástico", "Heterocedástico"),
     lp       = c(lpy_m1, lpy_m2),
     pDIC     = c(pDIC_m1, pDIC_m2),
     DIC      = c(dic_m1, dic_m2),
     lppd     = c(lppd_m1, lppd_m2),
     pWAIC    = c(pWAIC_m1, pWAIC_m2),
     WAIC     = c(waic_m1, waic_m2),
     row.names = NULL
)

tab_modelos

# Fin --------------------------------------------------------------------------