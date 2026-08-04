# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

# Parámetros de la simulación
n      <- 100
H      <- 2
theta0 <- c(0, 5)
sig20  <- c(0.5, 1)
omega0 <- c(2 / 3, 1 / 3)

# Simulación de datos
set.seed(123)

xi <- sample(
     x       = 1:H,
     size    = n,
     replace = TRUE,
     prob    = omega0
)

y <- rnorm(
     n    = n,
     mean = theta0[xi],
     sd   = sqrt(sig20[xi])
)

# Resumen de los datos simulados
table(xi)

# Valores poblacionales verdaderos --------------------------------------------

xi_true    <- xi
theta_true <- theta0
sig2_true  <- sig20
omega_true <- omega0

# Función de densidad verdadera ------------------------------------------------

f_true <- function(x) {
     sum(
          omega_true * dnorm(
               x    = x,
               mean = theta_true,
               sd   = sqrt(sig2_true)
          )
     )
}

# Análisis exploratorio --------------------------------------------------------

# Resumen
round(mean(y), 3)

round(median(y), 3)

round(sd(y), 3)

# Histograma -------------------------------------------------------------------

pdf(
     file      = "simulacion_modelo_mezcla_finita_histograma.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

x0 <- seq(
     from       = min(y) - 2,
     to         = max(y) + 2,
     length.out = 1000
)

y0 <- NULL

for (i in 1:length(x0)) {
     y0[i] <- f_true(x0[i])
}

hist(
     x        = y,
     breaks   = 20,
     freq     = FALSE,
     col      = "gray90",
     border   = "white",
     main     = "",
     xlab     = "y",
     ylab     = "Densidad",
     cex.axis = 0.85,
     cex.lab  = 0.95
)

lines(
     x   = x0,
     y   = y0,
     lwd = 3,
     col = 2
)

box()

dev.off()

# Distribuciones condicionales completas ---------------------------------------

sample_theta <- function(nh, ybh, H, mu0, gam02, theta, sig2) {
     for (h in 1:H) {
          if (nh[h] == 0) {
               theta[h] <- rnorm(n = 1, mean = mu0, sd = sqrt(gam02))
          } else {
               v2 <- 1 / (1 / gam02 + nh[h] / sig2[h])
               m  <- v2 * (mu0 / gam02 + nh[h] * ybh[h] / sig2[h])
               
               theta[h] <- rnorm(n = 1, mean = m, sd = sqrt(v2))
          }
     }
     
     return(theta)
}

sample_sig2 <- function(nh, ybh, ssh, H, nu0, sig02, theta, sig2) {
     for (h in 1:H) {
          if (nh[h] == 0) {
               sig2[h] <- 1 / rgamma(
                    n     = 1,
                    shape = 0.5 * nu0,
                    rate  = 0.5 * nu0 * sig02
               )
          } else {
               a <- 0.5 * (nu0 + nh[h])
               b <- 0.5 * (nu0 * sig02 + ssh[h] + nh[h] * (ybh[h] - theta[h])^2)
               
               sig2[h] <- 1 / rgamma(n = 1, shape = a, rate = b)
          }
     }
     
     return(sig2)
}

sample_xi <- function(H, omega, theta, sig2, n, y) {
     lp <- outer(
          X   = y,
          Y   = 1:H,
          FUN = function(y_i, h) {
               log(omega[h]) +
                    dnorm(
                         x    = y_i,
                         mean = theta[h],
                         sd   = sqrt(sig2[h]),
                         log  = TRUE
                    )
          }
     )
     
     xi <- apply(
          X      = lp,
          MARGIN = 1,
          FUN    = function(row) {
               prob <- exp(row - max(row))
               prob <- prob / sum(prob)
               
               sample(x = 1:H, size = 1, prob = prob)
          }
     )
     
     return(xi)
}

sample_omega <- function(nh, alpha0) {
     omega <- c(gtools::rdirichlet(n = 1, alpha = alpha0 + nh))
     
     return(omega)
}

# Muestreador de Gibbs ---------------------------------------------------------

mcmc <- function(y, H, mu0, gam02, nu0, sig02, alpha0, n_sams, n_burn, n_skip, verbose = TRUE) {
     # Ajustes
     y  <- scale(y)
     yb <- attr(y, "scaled:center")
     sy <- attr(y, "scaled:scale")
     y  <- as.numeric(y)
     n  <- length(y)
     
     # Número de iteraciones
     B    <- n_burn + n_sams * n_skip
     ncat <- max(1, floor(0.1 * B))
     
     # Valores iniciales
     xi    <- kmeans(x = y, centers = H, nstart = 25)$cluster
     omega <- as.numeric(table(factor(x = xi, levels = 1:H))) / n
     theta <- rnorm(n = H, mean = mu0, sd = sqrt(gam02))
     sig2  <- 1 / rgamma(n = H, shape = nu0 / 2, rate = nu0 * sig02 / 2)
     
     # Almacenamiento
     THETA <- matrix(data = NA, nrow = n_sams, ncol = H)
     SIG2  <- matrix(data = NA, nrow = n_sams, ncol = H)
     OMEGA <- matrix(data = NA, nrow = n_sams, ncol = H)
     XI    <- matrix(data = NA, nrow = n_sams, ncol = n)
     LL    <- rep(NA, n_sams)
     
     # Cadena
     for (i in 1:B) {
          # Actualizar estadísticos suficientes
          nh  <- as.numeric(table(factor(x = xi, levels = 1:H)))
          ybh <- ssh <- rep(NA, H)
          
          for (h in 1:H) {
               if (nh[h] > 0) {
                    indexh <- xi == h
                    ybh[h] <- mean(y[indexh])
                    ssh[h] <- sum((y[indexh] - ybh[h])^2)
               }
          }
          
          # Actualizar parámetros
          sig2  <- sample_sig2(nh, ybh, ssh, H, nu0, sig02, theta, sig2)
          theta <- sample_theta(nh, ybh, H, mu0, gam02, theta, sig2)
          omega <- sample_omega(nh, alpha0)
          xi    <- sample_xi(H, omega, theta, sig2, n, y)
          
          # Almacenar y log-verosimilitud
          if (i > n_burn && (i - n_burn) %% n_skip == 0) {
               k <- (i - n_burn) / n_skip
               
               THETA[k, ] <- sy * theta + yb
               SIG2[k, ]  <- sy^2 * sig2
               OMEGA[k, ] <- omega
               XI[k, ]    <- xi
               LL[k]      <- sum(dnorm(y, mean = theta[xi], sd = sqrt(sig2[xi]), log = TRUE))
          }
          
          # Progreso
          if (verbose && i %% ncat == 0) {
               cat(sprintf("%.1f%% completado\n", 100 * i / B))
          }
     }
     
     # Salida
     return(
          list(
               THETA  = THETA,
               SIG2   = SIG2,
               OMEGA  = OMEGA,
               XI     = XI,
               LL     = LL,
               center = yb,
               scale  = sy
          )
     )
}

# Ajuste del modelo ------------------------------------------------------------

# Número de grupos
H <- 5

# Hiperparámetros
alpha0 <- rep(1 / H, H)
mu0    <- 0
gam02  <- 1
nu0    <- 1
sig02  <- 1

# Número de iteraciones
n_sams <- 10000
n_burn <- 10000
n_skip <- 10

# Muestreador de Gibbs
set.seed(123)

chain <- mcmc(
     y       = y,
     H       = H,
     mu0     = mu0,
     gam02   = gam02,
     nu0     = nu0,
     sig02   = sig02,
     alpha0  = alpha0,
     n_sams  = n_sams,
     n_burn  = n_burn,
     n_skip  = n_skip,
     verbose = TRUE
)

# Guardar muestras
save(chain, file = "simulacion_modelo_mezcla_finita_muestras_mcmc.RData")

# Cargar muestras --------------------------------------------------------------

load("simulacion_modelo_mezcla_finita_muestras_mcmc.RData")

# Función para graficar cadenas ------------------------------------------------

plot_cadena <- function(x, ylab, col, cex, file) {
     pdf(
          file      = file,
          width     = 6,
          height    = 4,
          pointsize = 15
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

# Cadena de la log-verosimilitud ----------------------------------------------

plot_cadena(
     x    = chain$LL,
     ylab = "Log-verosimilitud",
     col  = 1,
     cex  = 0.5,
     file = "simulacion_modelo_mezcla_finita_cadena_logverosimilitud.pdf"
)

# Número posterior de componentes ocupados -------------------------------------

n_sams <- nrow(chain$XI)
H      <- ncol(chain$OMEGA)

H_post <- apply(
     X      = chain$XI,
     MARGIN = 1,
     FUN    = function(x) length(unique(x))
)

prob_H <- table(
     factor(
          x      = H_post,
          levels = 1:H
     )
) / n_sams

H_hat <- as.integer(names(which.max(table(H_post))))

prob_H

H_hat

# Gráfico ----------------------------------------------------------------------

pdf(
     file      = "simulacion_modelo_mezcla_finita_numero_componentes_ocupados.pdf",
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
     x    = 1:H,
     y    = as.numeric(prob_H),
     type = "h",
     lwd  = 5,
     xlim = c(0.5, H + 0.5),
     ylim = c(0, max(prob_H) * 1.1),
     xaxt = "n",
     yaxt = "n",
     main = "",
     xlab = "Número de componentes ocupados",
     ylab = "Probabilidad posterior"
)

axis(
     side   = 1,
     at     = 1:H,
     labels = 1:H
)

axis(
     side   = 2,
     at     = seq(0, max(prob_H) * 1.1, length.out = 6),
     labels = round(seq(0, max(prob_H) * 1.1, length.out = 6), 2)
)

dev.off()

# Densidad posterior predictiva ------------------------------------------------

# Grilla donde se evalúa la densidad predictiva posterior
M <- 1000

x0 <- seq(
     from       = min(y) - 2,
     to         = max(y) + 2,
     length.out = M
)

# Número de muestras posteriores y número máximo de componentes
B <- nrow(chain$OMEGA)
H <- ncol(chain$OMEGA)

# Densidad verdadera evaluada en la grilla
y0 <- rep(NA, M)

# Matriz para almacenar la densidad predictiva en cada iteración
FE <- matrix(data = NA, nrow = B, ncol = M)

for (i in 1:M) {
     # Evaluar la densidad verdadera en x0[i]
     y0[i] <- f_true(x0[i])
     
     for (b in 1:B) {
          # Evaluar la densidad de mezcla en x0[i] para la iteración b
          FE[b, i] <- sum(
               chain$OMEGA[b, ] *
                    dnorm(
                         x    = x0[i],
                         mean = chain$THETA[b, ],
                         sd   = sqrt(chain$SIG2[b, ])
                    )
          )
     }
}

# Resúmenes posteriores de la densidad ----------------------------------------

# Media posterior de la densidad predictiva
f_hat <- colMeans(FE)

# Intervalo de credibilidad puntual del 95%
f_inf <- apply(
     X      = FE,
     MARGIN = 2,
     FUN    = quantile,
     probs  = 0.025
)

f_sup <- apply(
     X      = FE,
     MARGIN = 2,
     FUN    = quantile,
     probs  = 0.975
)

# Visualización de muestras posteriores ---------------------------------------

pdf(
     file      = "simulacion_modelo_mezcla_finita_densidades_posteriores.pdf",
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
     x        = x0,
     y        = f_hat,
     type     = "n",
     xlim     = range(x0),
     ylim     = c(0, max(FE) * 1.01),
     cex.axis = 0.85,
     cex.lab  = 0.95,
     main     = "",
     xlab     = "y",
     ylab     = "Densidad"
)

set.seed(123)

index_plot <- sample(
     x       = 1:nrow(FE),
     size    = min(100, nrow(FE)),
     replace = FALSE
)

for (b in index_plot) {
     lines(
          x   = x0,
          y   = FE[b, ],
          lwd = 1,
          col = adjustcolor("black", alpha.f = 0.1)
     )
}

lines(
     x   = x0,
     y   = y0,
     lwd = 3,
     col = 2
)

dev.off()

# Media posterior e intervalo de credibilidad al 95% --------------------------

pdf(
     file      = "simulacion_modelo_mezcla_finita_densidad_media_ic95.pdf",
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
     x        = x0,
     y        = f_hat,
     type     = "n",
     xlim     = range(x0),
     ylim     = c(0, max(FE) * 1.01),
     cex.axis = 0.85,
     cex.lab  = 0.95,
     main     = "",
     xlab     = "y",
     ylab     = "Densidad"
)

polygon(
     x      = c(x0, rev(x0)),
     y      = c(f_inf, rev(f_sup)),
     col    = adjustcolor("gray70", alpha.f = 0.35),
     border = NA
)

lines(
     x   = x0,
     y   = y0,
     lwd = 2,
     col = 2
)

lines(
     x   = x0,
     y   = f_hat,
     lwd = 3
)

legend(
     "topright",
     legend = c("Media", "IC 95%", "Verdad"),
     lwd    = 2,
     lty    = 1,
     col    = c(1, "gray70", 2),
     bty    = "n",
     cex    = 0.85
)

dev.off()

# Matriz de probabilidades posteriores de coagrupamiento -----------------------

B <- nrow(chain$XI)
n <- ncol(chain$XI)

A <- matrix(data = 0, nrow = n, ncol = n)

for (b in 1:B) {
     # Indicadoras de coagrupamiento para la iteración b
     A <- A + outer(
          X   = chain$XI[b, ],
          Y   = chain$XI[b, ],
          FUN = "=="
     ) / B
}

diag(A) <- 1

# Ordenar las observaciones de acuerdo con la partición verdadera
indices <- order(xi_true)
A       <- A[indices, indices]

# Visualización de la matriz de coagrupamiento --------------------------------

pdf(
     file      = "simulacion_modelo_mezcla_finita_matriz_probabilidades_coagrupamiento.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(2.75, 2.75, 0.5, 0.5),
     mgp   = c(1.7, 0.7, 0)
)

corrplot::corrplot(
     corr        = A,
     is.corr     = FALSE,
     method      = "color",
     addgrid.col = NA,
     tl.pos      = "n",
     cl.cex      = 0.8
)

dev.off()

# Partición estimada por clustering jerárquico ---------------------------------

D_A <- as.dist(1 - A)

hc_A <- hclust(
     d      = D_A,
     method = "average"
)

particion <- cutree(
     tree = hc_A,
     k    = H_hat
)

# Resumen de la partición
print(table(particion))

# Matriz de grupos
n <- length(particion)

AA <- matrix(data = 0, nrow = n, ncol = n)

diag(AA) <- 1

for (i in 1:(n - 1)) {
     for (j in (i + 1):n) {
          if (particion[i] == particion[j]) {
               AA[i, j] <- 1
               AA[j, i] <- 1
          }
     }
}

# Ordenar la matriz de grupos
order_indices <- order(particion)
AA            <- AA[order_indices, order_indices]

# Visualización de la matriz de partición -------------------------------------

pdf(
     file      = "simulacion_modelo_mezcla_finita_matriz_particion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(2.75, 2.75, 0.5, 0.5),
     mgp   = c(1.7, 0.7, 0)
)

corrplot::corrplot(
     corr        = AA,
     is.corr     = FALSE,
     method      = "color",
     addgrid.col = NA,
     tl.pos      = "n",
     cl.pos      = "n"
)

dev.off()

# Adjusted Rand Index ----------------------------------------------------------

B   <- nrow(chain$XI)
ari <- rep(NA, B)

for (b in 1:B) {
     ari[b] <- aricode::ARI(c1 = chain$XI[b, ], c2 = as.numeric(xi_true))
}

# Media posterior
round(mean(ari), 3)

# Intervalo de credibilidad al 95%
round(quantile(x = ari, probs = c(0.025, 0.975)), 3)

# Distribución posterior del ARI
pdf(
     file      = "simulacion_modelo_mezcla_finita_ari.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(2.75, 2.75, 0.5, 0.5),
     mgp   = c(1.7, 0.7, 0)
)

hist(
     x        = ari,
     breaks   = 20,
     col      = "gray90",
     border   = "white",
     prob     = TRUE,
     main     = "",
     xlab     = "ARI",
     ylab     = "Densidad",
     cex.axis = 0.85,
     cex.lab  = 0.95
)

# Resumen posterior
percentiles <- quantile(x = ari, probs = c(0.025, 0.975))
media_ari   <- mean(ari)

abline(v = percentiles[1], col = 4, lty = 2, lwd = 2)
abline(v = percentiles[2], col = 4, lty = 2, lwd = 2)
abline(v = media_ari, col = 2, lty = 2, lwd = 2)

legend(
     "topleft",
     legend = c("Media", "IC 95%"),
     col    = c(2, 4),
     fill   = c(2, 4),
     border = c(2, 4),
     bty    = "n"
)

dev.off()

# Distribución posterior de H --------------------------------------------------

B <- nrow(chain$XI)

H_post <- apply(
     X      = chain$XI,
     MARGIN = 1,
     FUN    = function(x) length(unique(x))
)

H_max <- max(H_post)

H_tab <- table(
     factor(
          x      = H_post,
          levels = 1:H_max,
          labels = 1:H_max
     )
)

H_hat <- which.max(H_tab)
H     <- H_hat

# Estimación de medias y varianzas de los componentes --------------------------

# Iteraciones con H_hat componentes ocupados
index_H <- which(H_post == H_hat)
B_H     <- length(index_H)
H       <- H_hat

THETA_H <- matrix(data = NA, nrow = B_H, ncol = H)
SIG2_H  <- matrix(data = NA, nrow = B_H, ncol = H)
OMEGA_H <- matrix(data = NA, nrow = B_H, ncol = H)

for (r in 1:B_H) {
     b <- index_H[r]
     
     active_clusters <- sort(unique(chain$XI[b, ]))
     THETA_H[r, ]    <- chain$THETA[b, active_clusters]
     SIG2_H[r, ]     <- chain$SIG2[b, active_clusters]
     OMEGA_H[r, ]    <- chain$OMEGA[b, active_clusters]
}

# Valor de referencia para corregir el cambio de etiquetas
theta_pos <- colMeans(THETA_H)

# Generar permutaciones
permu <- gtools::permutations(n = H, r = H)

# Promediar en los espacios permutados
THETA_corrected <- matrix(data = NA, nrow = B_H, ncol = H)
SIG2_corrected  <- matrix(data = NA, nrow = B_H, ncol = H)
OMEGA_corrected <- matrix(data = NA, nrow = B_H, ncol = H)

for (r in 1:B_H) {
     theta_current <- THETA_H[r, ]
     sig2_current  <- SIG2_H[r, ]
     omega_current <- OMEGA_H[r, ]
     
     # Se escoge la permutación más cercana al valor de referencia
     dist <- apply(
          X      = permu,
          MARGIN = 1,
          FUN    = function(p) {
               permuted_theta <- theta_current[p]
               sum((permuted_theta - theta_pos)^2)
          }
     )
     
     best_permu <- permu[which.min(dist), ]
     
     THETA_corrected[r, ] <- theta_current[best_permu]
     SIG2_corrected[r, ]  <- sig2_current[best_permu]
     OMEGA_corrected[r, ] <- omega_current[best_permu]
}

# Ordenar los componentes según la media posterior
order_theta     <- order(colMeans(THETA_corrected))
THETA_corrected <- THETA_corrected[, order_theta]
SIG2_corrected  <- SIG2_corrected[, order_theta]
OMEGA_corrected <- OMEGA_corrected[, order_theta]

# Resumen posterior de las medias
theta_media <- colMeans(THETA_corrected)
theta_cv    <- apply(X = THETA_corrected, MARGIN = 2, FUN = sd) / abs(theta_media)

tab_theta <- rbind(
     theta_media,
     theta_cv,
     apply(
          X      = THETA_corrected,
          MARGIN = 2,
          FUN    = quantile,
          probs  = c(0.025, 0.975)
     )
)

colnames(tab_theta) <- paste("Cluster", 1:H)
rownames(tab_theta) <- c("Media", "CV", "2.5%", "97.5%")

round(t(tab_theta), 3)

# Resumen posterior de las varianzas
sig2_media <- colMeans(SIG2_corrected)
sig2_cv    <- apply(X = SIG2_corrected, MARGIN = 2, FUN = sd) / abs(sig2_media)

tab_sig2 <- rbind(
     sig2_media,
     sig2_cv,
     apply(
          X      = SIG2_corrected,
          MARGIN = 2,
          FUN    = quantile,
          probs  = c(0.025, 0.975)
     )
)

colnames(tab_sig2) <- paste("Cluster", 1:H)
rownames(tab_sig2) <- c("Media", "CV", "2.5%", "97.5%")

round(t(tab_sig2), 3)

# Resumen posterior de los pesos de mezcla
omega_media <- colMeans(OMEGA_corrected)
omega_cv    <- apply(X = OMEGA_corrected, MARGIN = 2, FUN = sd) / abs(omega_media)

tab_omega <- rbind(
     omega_media,
     omega_cv,
     apply(
          X      = OMEGA_corrected,
          MARGIN = 2,
          FUN    = quantile,
          probs  = c(0.025, 0.975)
     )
)

colnames(tab_omega) <- paste("Cluster", 1:H)
rownames(tab_omega) <- c("Media", "CV", "2.5%", "97.5%")

round(t(tab_omega), 3)

# Fin --------------------------------------------------------------------------