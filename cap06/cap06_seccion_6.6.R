# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
dir_work <- "~/Dropbox/UN/bayes_book"
setwd(dir_work)

# Librerías
suppressMessages(suppressWarnings(library(mvtnorm)))

# Datos ------------------------------------------------------------------------

# Puntajes de comprensión lectora antes y después de la instrucción
pre_test <- c(
     59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49,
     72, 34, 70, 34, 50, 41, 52, 60, 34, 28, 35
)

post_test <- c(
     77, 39, 46, 26, 38, 43, 68, 86, 77, 60, 50,
     59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33
)

# Matriz de datos
Y <- cbind(
     pre_test  = pre_test,
     post_test = post_test
)

# Dimensiones
n <- nrow(Y)
p <- ncol(Y)

# Inspección descriptiva
summary(Y)

# Estadísticos suficientes -----------------------------------------------------

# Media muestral
yb <- colMeans(Y)

# Matriz de covarianza muestral
SS <- cov(Y)

# Resultados
round(yb, 1)
round(SS, 1)

# Gráfico de dispersión --------------------------------------------------------

pdf(
     file      = "comprension_lectura_dispersion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = Y[, "pre_test"],
     y    = Y[, "post_test"],
     pch  = 16,
     col  = adjustcolor(4, alpha.f = 0.7),
     cex  = 1.1,
     xlab = "Antes",
     ylab = "Después",
     main = ""
)

# Líneas de referencia en las medias muestrales
abline(
     v   = mean(Y[, "pre_test"]),
     col = "gray70",
     lwd = 2,
     lty = 2
)

abline(
     h   = mean(Y[, "post_test"]),
     col = "gray70",
     lwd = 2,
     lty = 2
)

# Punto correspondiente al vector de medias muestral
points(
     x   = mean(Y[, "pre_test"]),
     y   = mean(Y[, "post_test"]),
     pch = 16,
     col = "gray70",
     cex = 1.3
)

dev.off()

# Correlación muestral
cor_y <- cor(Y[, "pre_test"], Y[, "post_test"])

round(cor_y, 3)

# Hiperparámetros --------------------------------------------------------------

mu0 <- c(50, 50)

L0 <- matrix(
     data  = c(
          278, 139,
          139, 278
     ),
     nrow  = 2,
     ncol  = 2,
     byrow = TRUE
)

nu0 <- 4

S0 <- matrix(
     data  = c(
          278, 139,
          139, 278
     ),
     nrow  = 2,
     ncol  = 2,
     byrow = TRUE
)

# Función auxiliar -------------------------------------------------------------

# Generar una matriz de WI(nu, S^{-1})
riwishart <- function(nu, S) {
     solve(stats::rWishart(n = 1, df = nu, Sigma = solve(S))[, , 1])
}

# Muestreador de Gibbs ---------------------------------------------------------

mcmc <- function(Y, mu0, L0, nu0, S0, B, seed = 123, verbose = TRUE) {
     # Cantidades de los datos
     n  <- nrow(Y)
     p  <- ncol(Y)
     yb <- colMeans(Y)
     SS <- cov(Y)
     
     # Cantidades previas fijas
     iL0   <- solve(L0)
     L0mu0 <- iL0 %*% mu0
     nun   <- nu0 + n
     SSn   <- S0 + (n - 1) * SS
     
     # Almacenamiento
     THETA <- matrix(data = NA_real_, nrow = B, ncol = p)
     colnames(THETA) <- paste0("theta_", 1:p)
     
     idx_sigma <- expand.grid(row = seq_len(p), col = seq_len(p))
     
     SIGMA <- matrix(data = NA_real_, nrow = B, ncol = p * p)
     colnames(SIGMA) <- paste0("sigma2_", idx_sigma$row, idx_sigma$col)
     
     YS <- matrix(data = NA_real_, nrow = B, ncol = p)
     colnames(YS) <- paste0("y_pred_", colnames(Y))
     
     LL <- numeric(B)
     
     # Inicialización desde la distribución previa
     Sigma <- riwishart(nu = nu0, S = S0)
     
     # Cantidad para mostrar progreso
     paso_progreso <- max(1, floor(B / 10))
     
     # Cadena
     set.seed(seed)
     
     for (b in seq_len(B)) {
          # Actualización de theta
          iSigma <- solve(Sigma)
          Ln     <- solve(iL0 + n * iSigma)
          mun    <- Ln %*% (L0mu0 + n * iSigma %*% yb)
          
          theta <- as.numeric(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = as.numeric(mun),
                    sigma = Ln
               )
          )
          
          # Actualización de Sigma
          Sn    <- SSn + n * tcrossprod(yb - theta)
          Sigma <- riwishart(nu = nun, S = Sn)
          
          # Distribución predictiva posterior
          YS[b, ] <- as.numeric(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = theta,
                    sigma = Sigma
               )
          )
          
          # Log-verosimilitud
          LL[b] <- sum(
               mvtnorm::dmvnorm(
                    x     = Y,
                    mean  = theta,
                    sigma = Sigma,
                    log   = TRUE
               )
          )
          
          # Almacenamiento
          THETA[b, ] <- theta
          SIGMA[b, ] <- as.vector(Sigma)
          
          if (verbose && b %% paso_progreso == 0) {
               cat(sprintf("%.0f%% completado...\n", 100 * b / B))
          }
     }
     
     list(
          THETA = as.data.frame(THETA),
          SIGMA = as.data.frame(SIGMA),
          LL    = LL,
          YS    = as.data.frame(YS)
     )
}

# Ajuste del modelo ------------------------------------------------------------

fit <- mcmc(
     Y       = Y,
     mu0     = mu0,
     L0      = L0,
     nu0     = nu0,
     S0      = S0,
     B       = 10000,
     seed    = 123,
     verbose = TRUE
)

# Extraer resultados -----------------------------------------------------------

THETA <- fit$THETA
SIGMA <- fit$SIGMA
LL    <- fit$LL
YS    <- fit$YS

# Función para graficar cadenas ------------------------------------------------

plot_cadena <- function(x, ylab, col, cex, file) {
     pdf(
          file      = file,
          width     = 6,
          height    = 4,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3.1, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = seq_along(x),
          y    = x,
          type = "p",
          pch  = 16,
          col  = adjustcolor(col, alpha.f = 0.1),
          cex  = cex,
          xlab = "Iteración",
          ylab = ylab,
          main = ""
     )
     
     dev.off()
}

# Métricas de precisión Monte Carlo -------------------------------------------

metricas_mcmc <- function(x) {
     # Convertir la cadena al formato mcmc
     x_mcmc <- coda::mcmc(x)
     
     # Calcular métricas
     neff <- round(as.numeric(coda::effectiveSize(x_mcmc)), 1)
     eemc <- round(stats::sd(x) / sqrt(neff), 4)
     cvmc <- round(eemc / abs(mean(x)), 4)
     
     c(
          neff = neff,
          eemc = eemc,
          cvmc = cvmc
     )
}

# Diagnóstico R-hat con múltiples cadenas -------------------------------------

calcular_rhat <- function(Y, mu0, L0, nu0, S0, B, M = 3, seed = 123) {
     # Ejecutar múltiples cadenas
     ajustes <- vector(mode = "list", length = M)
     
     for (m in seq_len(M)) {
          ajustes[[m]] <- mcmc(
               Y       = Y,
               mu0     = mu0,
               L0      = L0,
               nu0     = nu0,
               S0      = S0,
               B       = B,
               seed    = seed + m,
               verbose = FALSE
          )
     }
     
     # Convertir THETA a mcmc.list
     cadenas_theta <- coda::mcmc.list(
          lapply(
               X   = ajustes,
               FUN = function(x) coda::mcmc(as.matrix(x$THETA))
          )
     )
     
     # Convertir SIGMA a mcmc.list
     cadenas_sigma <- coda::mcmc.list(
          lapply(
               X   = ajustes,
               FUN = function(x) coda::mcmc(as.matrix(x$SIGMA))
          )
     )
     
     # Convertir LL a mcmc.list
     cadenas_ll <- coda::mcmc.list(
          lapply(
               X   = ajustes,
               FUN = function(x) {
                    coda::mcmc(matrix(x$LL, ncol = 1, dimnames = list(NULL, "LL")))
               }
          )
     )
     
     # Calcular R-hat
     rhat_theta <- coda::gelman.diag(
          x          = cadenas_theta,
          autoburnin = FALSE
     )$psrf
     
     rhat_sigma <- coda::gelman.diag(
          x          = cadenas_sigma,
          autoburnin = FALSE
     )$psrf
     
     rhat_ll <- coda::gelman.diag(
          x          = cadenas_ll,
          autoburnin = FALSE
     )$psrf
     
     # Organizar resultados
     tab_theta <- data.frame(
          parametro   = rownames(rhat_theta),
          rhat        = rhat_theta[, "Point est."],
          rhat_ic_sup = rhat_theta[, "Upper C.I."],
          tipo        = "theta",
          row.names   = NULL
     )
     
     tab_sigma <- data.frame(
          parametro   = rownames(rhat_sigma),
          rhat        = rhat_sigma[, "Point est."],
          rhat_ic_sup = rhat_sigma[, "Upper C.I."],
          tipo        = "sigma",
          row.names   = NULL
     )
     
     tab_ll <- data.frame(
          parametro   = rownames(rhat_ll),
          rhat        = rhat_ll[, "Point est."],
          rhat_ic_sup = rhat_ll[, "Upper C.I."],
          tipo        = "logverosimilitud",
          row.names   = NULL
     )
     
     tab_rhat <- rbind(tab_theta, tab_sigma, tab_ll)
     
     rownames(tab_rhat) <- tab_rhat$parametro
     tab_rhat$parametro <- NULL
     
     list(
          rhat    = tab_rhat,
          ajustes = ajustes
     )
}

# Resumen posterior ------------------------------------------------------------

resumen_posterior <- function(x, nivel_confianza = 0.95) {
     # Probabilidades para el intervalo de credibilidad
     alpha <- 1 - nivel_confianza
     probs <- c(alpha / 2, 1 - alpha / 2)
     
     # Resumen posterior
     media <- mean(x)
     cv    <- sd(x) / abs(media)
     ic    <- quantile(x, probs = probs, names = FALSE)
     
     round(
          c(
               media = media,
               cv    = cv,
               li    = ic[1],
               ls    = ic[2]
          ),
          digits = 3
     )
}

# Convergencia -----------------------------------------------------------------

# Cadenas para la log-verosimilitud
plot_cadena(
     x    = LL,
     ylab = "Log-verosimilitud",
     col  = 1,
     cex  = 0.5,
     file = "comprension_lectura_cadena_logverosimilitud.pdf"
)

# Cadenas para cada componente de theta
ylab_theta <- list(
     expression(theta[1]),
     expression(theta[2])
)

for (j in seq_along(THETA)) {
     plot_cadena(
          x    = THETA[[j]],
          ylab = ylab_theta[[j]],
          col  = 1,
          cex  = 0.5,
          file = paste0("comprension_lectura_cadena_", colnames(THETA)[j], ".pdf")
     )
}

# Cadenas para cada componente de Sigma
ylab_sigma <- list(
     expression(sigma[1]^2),
     expression(sigma[21]),
     expression(sigma[12]),
     expression(sigma[2]^2)
)

for (j in seq_along(SIGMA)) {
     plot_cadena(
          x    = SIGMA[[j]],
          ylab = ylab_sigma[[j]],
          col  = 1,
          cex  = 0.5,
          file = paste0("comprension_lectura_cadena_", colnames(SIGMA)[j], ".pdf")
     )
}

# Métricas de convergencia para log-verosimilitud
metricas_mcmc(LL)

# Métricas de convergencia para THETA
tab_theta <- t(
     sapply(
          X   = THETA,
          FUN = metricas_mcmc
     )
)

round(tab_theta, 5)

# Métricas de convergencia para SIGMA
tab_sigma <- t(
     sapply(
          X   = SIGMA,
          FUN = metricas_mcmc
     )
)

round(tab_sigma, 5)

# R-hat
diagnostico_rhat <- calcular_rhat(
     Y    = Y,
     mu0  = mu0,
     L0   = L0,
     nu0  = nu0,
     S0   = S0,
     B    = 10000,
     M    = 3,
     seed = 123
)

round(diagnostico_rhat$rhat[, c("rhat", "rhat_ic_sup")], 4)

# Inferencia -------------------------------------------------------------------

# Diferencia posterior entre medias
delta_theta <- THETA$theta_2 - THETA$theta_1

# Diferencia predictiva posterior
delta_pred <- YS$y_pred_post_test - YS$y_pred_pre_test

# Correlación posterior inducida por Sigma
rho_sigma <- SIGMA$sigma2_12 / sqrt(SIGMA$sigma2_11 * SIGMA$sigma2_22)

# Resúmenes posteriores
tab_resumen <- rbind(
     delta_theta = resumen_posterior(
          x               = delta_theta,
          nivel_confianza = 0.95
     ),
     delta_pred = resumen_posterior(
          x               = delta_pred,
          nivel_confianza = 0.95
     ),
     rho_sigma = resumen_posterior(
          x               = rho_sigma,
          nivel_confianza = 0.95
     )
)

tab_resumen

# Probabilidades posteriores
prob_delta_theta_pos <- mean(delta_theta > 0)
prob_delta_pred_pos  <- mean(delta_pred > 0)

round(prob_delta_theta_pos, 4)
round(prob_delta_pred_pos, 4)

# Gráfico posterior de theta_2 vs theta_1 --------------------------------------

lim_theta <- range(THETA[, 1], THETA[, 2])

pdf(
     file      = "comprension_lectura_posterior_delta_theta.pdf",
     width     = 5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = THETA[, 1],
     y    = THETA[, 2],
     pch  = 16,
     col  = adjustcolor(4, alpha.f = 0.1),
     cex  = 0.5,
     xlim = lim_theta,
     ylim = lim_theta,
     xlab = expression(theta[1]),
     ylab = expression(theta[2]),
     main = "",
     asp  = 1
)

abline(
     v   = mean(THETA[, 1]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     h   = mean(THETA[, 2]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     a   = 0,
     b   = 1,
     col = "gray30",
     lwd = 2
)

points(
     x   = mean(THETA[, 1]),
     y   = mean(THETA[, 2]),
     pch = 3,
     col = 2,
     lwd = 2,
     cex = 1.3
)

dev.off()

# Gráfico predictivo posterior de y_2^* vs y_1^* -------------------------------

lim_pred <- range(YS[, 1], YS[, 2])

pdf(
     file      = "comprension_lectura_posterior_delta_pred.pdf",
     width     = 5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = YS[, 1],
     y    = YS[, 2],
     pch  = 16,
     col  = adjustcolor(4, alpha.f = 0.1),
     cex  = 0.5,
     xlim = lim_pred,
     ylim = lim_pred,
     xlab = expression(tilde(y)[1]),
     ylab = expression(tilde(y)[2]),
     main = "",
     asp  = 1
)

abline(
     v   = mean(YS[, 1]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     h   = mean(YS[, 2]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     a   = 0,
     b   = 1,
     col = "gray30",
     lwd = 2
)

points(
     x   = mean(YS[, 1]),
     y   = mean(YS[, 2]),
     pch = 3,
     col = 2,
     lwd = 2,
     cex = 1.3
)

dev.off()

# Distribución posterior de la correlación -------------------------------------

RHO <- SIGMA[, 3] / sqrt(SIGMA[, 1] * SIGMA[, 4])

pdf(
     file      = "comprension_lectura_posterior_correlacion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     x      = RHO,
     freq   = FALSE,
     breaks = 30,
     col    = "gray90",
     border = "white",
     xlim   = c(0, 1),
     xlab   = expression(rho),
     ylab   = "Densidad",
     main   = ""
)

abline(
     v   = mean(RHO),
     col = 2,
     lwd = 2,
     lty = 2
)

abline(
     v   = quantile(RHO, probs = c(0.025, 0.975)),
     col = 4,
     lwd = 2,
     lty = 2
)

legend(
     x      = "topleft",
     legend = c("Media", "IC 95%"),
     col    = c(2, 4),
     fill   = c(2, 4),
     border = c(2, 4),
     bty    = "n"
)

dev.off()

# Fin --------------------------------------------------------------------------