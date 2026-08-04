# Configuración ----------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book")

# Librerías
library(readxl)
library(dplyr)
library(stringr)
library(stringi)
library(mvtnorm)
library(coda)

# Datos ------------------------------------------------------------------------

# Importar archivos
ideal_points <- read_excel("parapolitica_ideal_points.xlsx")
respuesta    <- read_excel("parapolitica_respuesta.xlsx")

# Estandarizar nombre de la primera columna
names(ideal_points)[1] <- "legislador"
names(respuesta)[1]    <- "legislador"

# Función para construir llave de emparejamiento
normalizar_nombre <- function(x) {
     x |>
          str_replace("\\s*:.*$", "") |>                 # elimina partido después de :
          str_squish() |>                                # elimina espacios extra
          str_to_upper() |>                              # convierte a mayúsculas
          stringi::stri_trans_general("Latin-ASCII") |>  # elimina tildes y ñ
          str_replace_all("[^A-Z ]", "") |>              # elimina signos raros
          str_squish()
}

# Crear llaves
ideal_points <- ideal_points |>
     mutate(id_legislador = normalizar_nombre(legislador))

respuesta <- respuesta |>
     mutate(id_legislador = normalizar_nombre(legislador))

# Revisar nombres que no hacen match
sin_match <- ideal_points |>
     anti_join(respuesta, by = "id_legislador") |>
     select(legislador, id_legislador)

sin_match

# Unir bases
parapolitica_modelo <- ideal_points |>
     left_join(
          respuesta |>
               select(-legislador),
          by = "id_legislador"
     ) |>
     select(
          legislador,
          Mean,
          parapolitica0,
          parapolitica1
     )

# Ajustar los nombres de las variables
parapolitica_modelo <- parapolitica_modelo |>
     rename(
          ideal   = Mean,
          parapol = parapolitica0
     ) |>
     select(
          legislador,
          ideal,
          parapol
     )

# Asegurar formato de la variable respuesta
parapolitica_modelo <- parapolitica_modelo |>
     mutate(
          parapol = as.integer(parapol)
     )

# Revisar resultado
glimpse(parapolitica_modelo)

# Ajustar modelo logit frecuentista --------------------------------------------

mod_logit_freq <- glm(
     parapol ~ ideal,
     data   = parapolitica_modelo,
     family = binomial(link = "logit")
)

# Resumen del modelo
summary(mod_logit_freq)

# Coeficientes estimados
coef(mod_logit_freq)

# Razones de odds
exp(coef(mod_logit_freq))

# Intervalos de confianza del 95% para los coeficientes
confint(mod_logit_freq)

# Intervalos de confianza del 95% para las razones de odds
exp(confint(mod_logit_freq))

# Datos para ajustar el modelo -------------------------------------------------

# Variable respuesta
y <- parapolitica_modelo$parapol

# Matriz de diseño: intercepto y punto ideal
X <- cbind(
     1,
     ideal = parapolitica_modelo$ideal
)

# Rango de los puntos ideales
round(range(parapolitica_modelo$ideal), 3)

# Tamaño de la muestra y número de predictores
n <- nrow(X)
p <- ncol(X)

# Análisis exploratorio --------------------------------------------------------

# Distribución de la variable respuesta
tabla_y <- table(y)

prop_y <- prop.table(tabla_y)

resumen_y <- data.frame(
     parapolitica = names(tabla_y),
     frecuencia   = as.numeric(tabla_y),
     proporcion   = round(as.numeric(prop_y), 3)
)

resumen_y

# Resumen global del punto ideal
resumen_ideal <- c(
     media   = mean(parapolitica_modelo$ideal),
     mediana = median(parapolitica_modelo$ideal),
     sd      = sd(parapolitica_modelo$ideal),
     minimo  = min(parapolitica_modelo$ideal),
     maximo  = max(parapolitica_modelo$ideal)
)

round(resumen_ideal, 3)

# Resumen del punto ideal según parapolítica
resumen_por_grupo <- aggregate(
     ideal ~ parapol,
     data = parapolitica_modelo,
     FUN  = function(x) {
          c(
               n       = length(x),
               media   = mean(x),
               mediana = median(x),
               sd      = sd(x),
               minimo  = min(x),
               maximo  = max(x)
          )
     }
)

resumen_por_grupo <- do.call(
     data.frame,
     resumen_por_grupo
)

names(resumen_por_grupo) <- c(
     "parapolitica",
     "n",
     "media",
     "mediana",
     "sd",
     "minimo",
     "maximo"
)

round(resumen_por_grupo, 3)

# Histograma de puntos ideales
pdf(
     file      = "parapolitica_eda_histograma_puntos_ideales.pdf",
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
     x      = X[, "ideal"],
     freq   = FALSE,
     col    = "gray85",
     border = "white",
     xlab   = "Punto ideal",
     ylab   = "Densidad",
     main   = ""
)

lines(
     density(X[, "ideal"]),
     lwd = 2
)

rug(
     X[, "ideal"],
     col = adjustcolor("black", 0.4)
)

abline(
     v   = mean(X[, "ideal"]),
     col = 2,
     lwd = 2,
     lty = 2
)

box()

dev.off()

# Diagrama de cajas de puntos ideales por grupo
pdf(
     file      = "parapolitica_eda_boxplot_puntos_ideales.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

grupo_parapol <- factor(
     parapolitica_modelo$parapol,
     levels = c(0, 1),
     labels = c("No involucrado", "Involucrado")
)

ylim_box <- range(parapolitica_modelo$ideal)
ylim_box <- ylim_box + c(-0.08, 0.08) * diff(ylim_box)

boxplot(
     parapolitica_modelo$ideal ~ grupo_parapol,
     outline  = FALSE,
     col      = "gray90",
     border   = "gray30",
     boxwex   = 0.45,
     lwd      = 1.4,
     medcol   = 1,
     medlwd   = 2,
     whisklty = 1,
     staplelty = 1,
     cex.axis = 0.85,
     xlab     = "",
     ylab     = "Punto ideal",
     ylim     = ylim_box,
     main     = ""
)

stripchart(
     parapolitica_modelo$ideal ~ grupo_parapol,
     vertical = TRUE,
     method   = "jitter",
     jitter   = 0.12,
     pch      = 16,
     cex      = 0.55,
     col      = adjustcolor("black", 0.45),
     add      = TRUE
)

points(
     x   = 1:2,
     y   = tapply(parapolitica_modelo$ideal, grupo_parapol, mean),
     pch = 18,
     cex = 1.25,
     col = 1
)

abline(
     h   = 0,
     lty = 3,
     col = "gray60"
)

box()

dev.off()

# MCMC modelo logit ------------------------------------------------------------

log1pexp <- function(x) {
     ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

log_likelihood_logit <- function(beta, y, X) {
     eta <- as.numeric(X %*% beta)
     sum(y * eta - log1pexp(eta))
}

log_posterior_logit <- function(beta, y, X, beta0, Sigma0) {
     log_lik <- log_likelihood_logit(beta, y, X)
     
     log_prior <- mvtnorm::dmvnorm(
          x     = beta,
          mean  = beta0,
          sigma = Sigma0,
          log   = TRUE
     )
     
     log_lik + log_prior
}

mcmc_logit <- function(
          y,
          X,
          beta0,
          Sigma0,
          B,
          burn_in,
          thin,
          seed = 123,
          progress = TRUE
) {
     # Verificaciones básicas
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n <- nrow(X)
     p <- ncol(X)
     
     beta_ini <- rep(0, p)
     
     # Calibración de la propuesta
     X_cov <- X[, -1, drop = FALSE]
     
     mod_logit_freq <- glm(
          y ~ X_cov,
          family = binomial(link = "logit")
     )
     
     Delta <- (2.38^2 / p) * vcov(mod_logit_freq)
     
     # Almacenamiento
     B_post <- burn_in + B * thin
     
     BETA_post <- matrix(NA, nrow = B, ncol = p)
     LL_post   <- numeric(B)
     
     colnames(BETA_post) <- paste0("beta", seq_len(p) - 1)
     
     ncat <- max(1, floor(B_post / 10))
     
     # Cadena
     beta <- beta_ini
     acr  <- 0
     s    <- 0
     
     set.seed(seed)
     
     for (b in seq_len(B_post)) {
          # 1. Propuesta
          beta_p <- c(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = beta,
                    sigma = Delta
               )
          )
          
          # 2. Razón de aceptación en escala logarítmica
          log_post_star <- log_posterior_logit(beta_p, y, X, beta0, Sigma0)
          log_post_curr <- log_posterior_logit(beta, y, X, beta0, Sigma0)
          
          log_r <- log_post_star - log_post_curr
          
          # 3. Aceptar o rechazar
          if (log(runif(1)) <= min(0, log_r)) {
               beta <- beta_p
               acr  <- acr + 1
          }
          
          # 4. Almacenar solo después del calentamiento y cada thin iteraciones
          if (b > burn_in && (b - burn_in) %% thin == 0) {
               s <- s + 1
               
               BETA_post[s, ] <- beta
               LL_post[s]     <- log_likelihood_logit(beta, y, X)
          }
          
          # 5. Progreso
          if (progress && b %% ncat == 0) {
               cat(round(100 * b / B_post, 1), "% completado ...\n", sep = "")
          }
     }
     
     # Tasa empírica de aceptación
     tasa_aceptacion <- acr / B_post
     
     # Retornar resultados
     list(
          BETA            = BETA_post,
          LL              = LL_post,
          tasa_aceptacion = tasa_aceptacion,
          Delta           = Delta
     )
}

# MCMC modelo probit aumentado -------------------------------------------------

log_likelihood_probit <- function(beta, y, X) {
     eta <- as.numeric(X %*% beta)
     
     sum(
          ifelse(
               y == 1,
               pnorm(eta, log.p = TRUE),
               pnorm(eta, lower.tail = FALSE, log.p = TRUE)
          )
     )
}

rnorm_trunc_probit <- function(mu, y) {
     n <- length(mu)
     z <- numeric(n)
     
     # y = 1: Normal truncada en (0, Inf)
     id1 <- which(y == 1)
     
     if (length(id1) > 0) {
          z[id1] <- truncnorm::rtruncnorm(
               n    = length(id1),
               a    = 0,
               b    = Inf,
               mean = mu[id1],
               sd   = 1
          )
     }
     
     # y = 0: Normal truncada en (-Inf, 0]
     id0 <- which(y == 0)
     
     if (length(id0) > 0) {
          z[id0] <- truncnorm::rtruncnorm(
               n    = length(id0),
               a    = -Inf,
               b    = 0,
               mean = mu[id0],
               sd   = 1
          )
     }
     
     z
}

mcmc_probit <- function(
          y,
          X,
          beta0,
          Sigma0,
          B,
          burn_in,
          thin,
          seed = 123,
          progress = TRUE
) {
     # Verificaciones básicas
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n <- nrow(X)
     p <- ncol(X)
     
     beta_ini <- rep(0, p)
     
     # Matriz de precisión previa
     P0 <- solve(Sigma0)
     
     # Matriz de covarianzas condicional de beta
     V_beta <- solve(P0 + t(X) %*% X)
     
     # Almacenamiento
     B_post <- burn_in + B * thin
     
     BETA_post <- matrix(NA, nrow = B, ncol = p)
     LL_post   <- numeric(B)
     
     colnames(BETA_post) <- paste0("beta", seq_len(p) - 1)
     
     ncat <- max(1, floor(B_post / 10))
     
     # Cadena
     beta <- beta_ini
     z    <- ifelse(y == 1, 1, -1)
     s    <- 0
     
     set.seed(seed)
     
     for (b in seq_len(B_post)) {
          # 1. Simular variables latentes z
          eta <- as.numeric(X %*% beta)
          z   <- rnorm_trunc_probit(mu = eta, y = y)
          
          # 2. Simular beta de su dcc Normal multivariada
          m_beta <- V_beta %*% (P0 %*% beta0 + t(X) %*% z)
          
          beta <- c(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = c(m_beta),
                    sigma = V_beta
               )
          )
          
          # 3. Almacenar solo después del calentamiento y cada thin iteraciones
          if (b > burn_in && (b - burn_in) %% thin == 0) {
               s <- s + 1
               
               BETA_post[s, ] <- beta
               LL_post[s]     <- log_likelihood_probit(beta, y, X)
          }
          
          # 4. Progreso
          if (progress && b %% ncat == 0) {
               cat(round(100 * b / B_post, 1), "% completado ...\n", sep = "")
          }
     }
     
     # Retornar resultados
     list(
          BETA   = BETA_post,
          LL     = LL_post,
          V_beta = V_beta
     )
}

# Ajuste del modelo logit ------------------------------------------------------

# Hiperparámetros
beta0  <- rep(0, p)
Sigma0 <- 10 * diag(1, p)

# Configuración del algoritmo
B       <- 10000
burn_in <- 10000
thin    <- 10

# Ajuste del modelo
muestras_logit <- mcmc_logit(
     y        = y,
     X        = X,
     beta0    = beta0,
     Sigma0   = Sigma0,
     B        = B,
     burn_in  = burn_in,
     thin     = thin,
     seed     = 123,
     progress = TRUE
)

# Tasa de aceptación
round(muestras_logit$tasa_aceptacion, 3)

# Extraer resultados
BETA_logit <- muestras_logit$BETA
LL_logit   <- muestras_logit$LL
acr_logit  <- muestras_logit$tasa_aceptacion

# Ajuste del modelo probit -----------------------------------------------------

# Hiperparámetros
beta0  <- rep(0, p)
Sigma0 <- 10 * diag(1, p)

# Configuración del algoritmo
B       <- 10000
burn_in <- 10000
thin    <- 10

# Ajuste del modelo
muestras_probit <- mcmc_probit(
     y        = y,
     X        = X,
     beta0    = beta0,
     Sigma0   = Sigma0,
     B        = B,
     burn_in  = burn_in,
     thin     = thin,
     seed     = 123,
     progress = TRUE
)

# Extraer resultados
BETA_probit <- muestras_probit$BETA
LL_probit   <- muestras_probit$LL

# Convergencia y resumen posterior modelo logit --------------------------------

# Matriz de cantidades posteriores: coeficientes y log-verosimilitud
MCMC_logit <- cbind(
     BETA_logit,
     logLik = LL_logit
)

# Tamaños efectivos de muestra
ESS_logit <- effectiveSize(MCMC_logit)

# Errores estándar de Monte Carlo
MCSE_logit <- apply(
     MCMC_logit,
     2,
     function(x) {
          sd(x) / sqrt(effectiveSize(x))
     }
)

# Resúmenes posteriores
resumen_posterior_logit <- data.frame(
     ess   = ESS_logit,
     mcse  = MCSE_logit,
     media = colMeans(MCMC_logit),
     sd    = apply(MCMC_logit, 2, sd),
     q025  = apply(MCMC_logit, 2, quantile, probs = 0.025),
     q975  = apply(MCMC_logit, 2, quantile, probs = 0.975)
)

round(resumen_posterior_logit, 3)

# Convergencia y resumen posterior modelo probit -------------------------------

# Matriz de cantidades posteriores: coeficientes y log-verosimilitud
MCMC_probit <- cbind(
     BETA_probit,
     logLik = LL_probit
)

# Tamaños efectivos de muestra
ESS_probit <- effectiveSize(MCMC_probit)

# Errores estándar de Monte Carlo
MCSE_probit <- apply(
     MCMC_probit,
     2,
     function(x) {
          sd(x) / sqrt(effectiveSize(x))
     }
)

# Resúmenes posteriores
resumen_posterior_probit <- data.frame(
     ess   = ESS_probit,
     mcse  = MCSE_probit,
     media = colMeans(MCMC_probit),
     sd    = apply(MCMC_probit, 2, sd),
     q025  = apply(MCMC_probit, 2, quantile, probs = 0.025),
     q975  = apply(MCMC_probit, 2, quantile, probs = 0.975)
)

round(resumen_posterior_probit, 3)

# Densidades posteriores de la log-verosimilitud: logit vs probit --------------

pdf(
     file      = paste0("parapolitica_logit_probit_log_verosimilitud.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

dens_LL_logit  <- density(LL_logit)
dens_LL_probit <- density(LL_probit)

ylim_LL <- range(dens_LL_logit$y, dens_LL_probit$y)

plot(
     dens_LL_logit,
     type = "l",
     lwd  = 3,
     col  = 4,
     xlab = "Log-verosimilitud",
     ylab = "Densidad",
     ylim = ylim_LL,
     main = ""
)

lines(
     dens_LL_probit,
     lwd = 3,
     col = 2
)

legend(
     "topleft",
     legend = c("Logit", "Probit"),
     col    = c(4, 2),
     lwd    = 3,
     bty    = "n",
     cex    = 1.5
)

dev.off()

# Densidades posteriores de los coeficientes: logit vs probit ------------------

parametros <- colnames(BETA_logit)

for (j in seq_along(parametros)) {
     pdf(
          file = paste0(
               "parapolitica_logit_probit_posterior_",
               parametros[j],
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     dens_logit  <- density(BETA_logit[, j])
     dens_probit <- density(BETA_probit[, j])
     
     xlim_j <- range(dens_logit$x, dens_probit$x)
     ylim_j <- range(dens_logit$y, dens_probit$y)
     
     plot(
          dens_logit,
          type = "l",
          lwd  = 3,
          col  = 4,
          xlim = xlim_j,
          ylim = ylim_j,
          main = "",
          xlab = parametros[j],
          ylab = "Densidad"
     )
     
     lines(
          dens_probit,
          lwd = 3,
          col = 2
     )
     
     box()
     
     dev.off()
}

# Curvas posteriores de probabilidad: modelo logit -----------------------------

# Grilla de puntos ideales
ideal_grid <- seq(
     from       = -3,
     to         = 3,
     length.out = 1000
)

# Matriz de diseño para predicción
X_grid <- cbind(
     1,
     ideal = ideal_grid
)

# Probabilidades posteriores para todas las iteraciones
PROB_logit <- sapply(
     seq_len(nrow(BETA_logit)),
     function(b) {
          plogis(X_grid %*% BETA_logit[b, ])
     }
)

# Media posterior de la curva
PROB_logit_media <- rowMeans(PROB_logit)

# Selección aleatoria de 100 curvas posteriores
set.seed(123)

ind_muestras <- sample(
     x    = seq_len(nrow(BETA_logit)),
     size = 100
)

PROB_logit_muestras <- PROB_logit[, ind_muestras]

# Gráfico
pdf(
     file      = "parapolitica_logit_curva_probabilidad_muestras.pdf",
     width     = 7,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = ideal_grid,
     y    = PROB_logit_media,
     type = "n",
     ylim = c(0, 1),
     xlab = "Punto ideal",
     ylab = "Probabilidad de parapolítica",
     main = ""
)

# Curvas posteriores simuladas
for (j in seq_len(ncol(PROB_logit_muestras))) {
     lines(
          x   = ideal_grid,
          y   = PROB_logit_muestras[, j],
          col = adjustcolor("gray40", 0.18),
          lwd = 1
     )
}

# Media posterior
lines(
     x   = ideal_grid,
     y   = PROB_logit_media,
     col = 4,
     lwd = 3
)

# Datos observados
points(
     x   = parapolitica_modelo$ideal,
     y   = jitter(parapolitica_modelo$parapol, amount = 0.025),
     pch = 16,
     cex = 0.6,
     col = adjustcolor("black", 0.5)
)

box()

dev.off()

# Curvas posteriores de probabilidad: logit vs probit --------------------------

# Grilla de puntos ideales
ideal_grid <- seq(
     from       = -3,
     to         = 3,
     length.out = 1000
)

# Matriz de diseño para predicción
X_grid <- cbind(
     1,
     ideal = ideal_grid
)

# Probabilidades posteriores para todas las iteraciones
PROB_logit <- sapply(
     seq_len(nrow(BETA_logit)),
     function(b) {
          plogis(X_grid %*% BETA_logit[b, ])
     }
)

PROB_probit <- sapply(
     seq_len(nrow(BETA_probit)),
     function(b) {
          pnorm(X_grid %*% BETA_probit[b, ])
     }
)

# Resúmenes posteriores de las curvas
PROB_logit_media <- rowMeans(PROB_logit)

PROB_logit_IC <- apply(
     PROB_logit,
     1,
     quantile,
     probs = c(0.025, 0.975)
)

PROB_probit_media <- rowMeans(PROB_probit)

PROB_probit_IC <- apply(
     PROB_probit,
     1,
     quantile,
     probs = c(0.025, 0.975)
)

# Gráfico
pdf(
     file      = "parapolitica_logit_probit_curva_probabilidad.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = ideal_grid,
     y    = PROB_logit_media,
     type = "n",
     ylim = c(0, 1),
     xlab = "Punto ideal",
     ylab = "Probabilidad de parapolítica"
)

# Bandas posteriores del 95%: logit
lines(
     x   = ideal_grid,
     y   = PROB_logit_IC[1, ],
     col = 4,
     lwd = 1,
     lty = 2
)

lines(
     x   = ideal_grid,
     y   = PROB_logit_IC[2, ],
     col = 4,
     lwd = 1,
     lty = 2
)

# Bandas posteriores del 95%: probit
lines(
     x   = ideal_grid,
     y   = PROB_probit_IC[1, ],
     col = 2,
     lwd = 1,
     lty = 2
)

lines(
     x   = ideal_grid,
     y   = PROB_probit_IC[2, ],
     col = 2,
     lwd = 1,
     lty = 2
)

# Medias posteriores
lines(
     x   = ideal_grid,
     y   = PROB_logit_media,
     col = 4,
     lwd = 3
)

lines(
     x   = ideal_grid,
     y   = PROB_probit_media,
     col = 2,
     lwd = 3
)

# Datos observados
points(
     x   = parapolitica_modelo$ideal,
     y   = jitter(parapolitica_modelo$parapol, amount = 0.025),
     pch = 16,
     cex = 0.6,
     col = adjustcolor("black", 0.5)
)

legend(
     "topleft",
     legend = c("Logit", "Probit"),
     col    = c(4, 2),
     lwd    = 3,
     lty    = 1,
     cex    = 1.5,
     bty    = "n"
)

box()

dev.off()

# AUC para cada modelo: logit vs probit ----------------------------------------

# Función para calcular AUC usando la media posterior de beta
auc_media_posterior <- function(BETA, y, X, link = c("logit", "probit")) {
     require(pROC)
     
     link <- match.arg(link)
     
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     # Media posterior de los coeficientes
     beta_media <- colMeans(BETA)
     
     # Predictor lineal
     eta <- as.numeric(X %*% beta_media)
     
     # Probabilidades estimadas de éxito
     if (link == "logit") {
          prob_hat <- plogis(eta)
     }
     
     if (link == "probit") {
          prob_hat <- pnorm(eta)
     }
     
     # AUC
     AUC <- as.numeric(
          pROC::auc(
               pROC::roc(
                    response  = y,
                    predictor = prob_hat,
                    levels    = c(0, 1),
                    direction = "<",
                    quiet     = TRUE
               )
          )
     )
     
     list(
          beta_media = beta_media,
          prob_hat   = prob_hat,
          AUC        = AUC
     )
}

AUC_logit_media <- auc_media_posterior(
     BETA = BETA_logit,
     y    = y,
     X    = X,
     link = "logit"
)

AUC_probit_media <- auc_media_posterior(
     BETA = BETA_probit,
     y    = y,
     X    = X,
     link = "probit"
)

# Resumen
resumen_AUC_media <- data.frame(
     AUC = c(
          AUC_logit_media$AUC,
          AUC_probit_media$AUC
     )
)

rownames(resumen_AUC_media) <- c("Logit", "Probit")

round(resumen_AUC_media, 3)

# DIC y WAIC para modelos logit y probit ---------------------------------------

# Función estable para log(1 + exp(x))
log1pexp <- function(x) {
     ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

# Matriz de log-verosimilitudes punto a punto
# Retorna una matriz de dimensión n x B, donde cada columna corresponde
# a una iteración posterior
log_lik_matrix_binary <- function(BETA, y, X, link = c("logit", "probit")) {
     link <- match.arg(link)
     
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     B <- nrow(BETA)
     n <- length(y)
     
     LL_mat <- matrix(NA, nrow = n, ncol = B)
     
     for (b in seq_len(B)) {
          eta_b <- as.numeric(X %*% BETA[b, ])
          
          if (link == "logit") {
               LL_mat[, b] <- y * eta_b - log1pexp(eta_b)
          }
          
          if (link == "probit") {
               LL_mat[, b] <- ifelse(
                    y == 1,
                    pnorm(eta_b, log.p = TRUE),
                    pnorm(eta_b, lower.tail = FALSE, log.p = TRUE)
               )
          }
     }
     
     LL_mat
}

# Función log-verosimilitud en la media posterior
log_likelihood_binary <- function(beta, y, X, link = c("logit", "probit")) {
     link <- match.arg(link)
     
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     eta <- as.numeric(X %*% beta)
     
     if (link == "logit") {
          ll <- sum(y * eta - log1pexp(eta))
     }
     
     if (link == "probit") {
          ll <- sum(
               ifelse(
                    y == 1,
                    pnorm(eta, log.p = TRUE),
                    pnorm(eta, lower.tail = FALSE, log.p = TRUE)
               )
          )
     }
     
     ll
}

# Función auxiliar log-mean-exp
log_mean_exp <- function(x) {
     m <- max(x)
     m + log(mean(exp(x - m)))
}

# Calcular DIC y WAIC
dic_waic_binary <- function(BETA, y, X, link = c("logit", "probit")) {
     link <- match.arg(link)
     
     # Log-verosimilitudes punto a punto
     LL_mat <- log_lik_matrix_binary(
          BETA = BETA,
          y    = y,
          X    = X,
          link = link
     )
     
     # Log-verosimilitud total por iteración
     LL_total <- colSums(LL_mat)
     
     # DIC
     D_bar <- mean(-2 * LL_total)
     
     beta_bar <- colMeans(BETA)
     
     D_hat <- -2 * log_likelihood_binary(
          beta = beta_bar,
          y    = y,
          X    = X,
          link = link
     )
     
     p_D <- D_bar - D_hat
     DIC <- D_bar + p_D
     
     # WAIC
     lppd <- sum(
          apply(
               LL_mat,
               1,
               log_mean_exp
          )
     )
     
     p_WAIC <- sum(
          apply(
               LL_mat,
               1,
               var
          )
     )
     
     WAIC <- -2 * (lppd - p_WAIC)
     
     data.frame(
          D_bar  = D_bar,
          D_hat  = D_hat,
          p_D    = p_D,
          DIC    = DIC,
          lppd   = lppd,
          p_WAIC = p_WAIC,
          WAIC   = WAIC
     )
}

# Aplicar a cada modelo
criterios_logit <- dic_waic_binary(
     BETA = BETA_logit,
     y    = y,
     X    = X,
     link = "logit"
)

criterios_probit <- dic_waic_binary(
     BETA = BETA_probit,
     y    = y,
     X    = X,
     link = "probit"
)

# Tabla comparativa
criterios_modelos <- rbind(
     Logit  = criterios_logit,
     Probit = criterios_probit
)

round(criterios_modelos, 3)

# Fin --------------------------------------------------------------------------