# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos analizados en
#    Dunn, Peter K., and Gordon K. Smyth. Generalized linear models with
#    examples in R. Vol. 53. New York: Springer, 2018.

# Datos ------------------------------------------------------------------------

library(GLMsData)

data(lungcap)

dat <- lungcap

dat$Gender <- ifelse(dat$Gender == "F", 0, 1)

dat$Smoke_f <- factor(
     dat$Smoke,
     levels = c(0, 1),
     labels = c("No fumador", "Fumador")
)

dat$Gender_f <- factor(
     dat$Gender,
     levels = c(0, 1),
     labels = c("Femenino", "Masculino")
)

# Análisis exploratorio --------------------------------------------------------

# Variable respuesta en escala logarítmica
dat$logFEV <- log(dat$FEV)

# Revisión básica de la base
summary(dat[, c("FEV", "logFEV", "Age", "Ht", "Gender", "Smoke")])

# Revisión de valores faltantes
colSums(is.na(dat[, c("FEV", "logFEV", "Age", "Ht", "Gender", "Smoke")]))

# Revisión de variables categóricas
round(prop.table(table(dat$Smoke_f)), 3)
round(prop.table(table(dat$Gender_f)), 3)

# Revisión de la variable respuesta original
summary(dat$FEV)
round(sd(dat$FEV), 3)

# Revisión de la variable respuesta transformada
summary(dat$logFEV)
round(sd(dat$logFEV), 3)

# Histograma de log(FEV) -------------------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_histograma_fev.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

logfev <- na.omit(dat$logFEV)

hist(
     logfev,
     freq   = FALSE,
     breaks = 25,
     main   = "",
     xlab   = "log(FEV)",
     ylab   = "Densidad",
     col    = "gray85",
     border = "white"
)

lines(density(logfev), lwd = 2)

curve(
     expr = dnorm(x, mean = mean(logfev), sd = sd(logfev)),
     lwd  = 2,
     lty  = 2,
     add  = TRUE
)

legend(
     "topright",
     legend = c("Kernel", "Normal"),
     lwd    = c(2, 2),
     lty    = c(1, 2),
     bty    = "n"
)

dev.off()

# Boxplot de log(FEV) por tabaquismo -------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_boxplot_fev_smoke.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

boxplot(
     logFEV ~ Smoke_f,
     data    = dat,
     main    = "",
     xlab    = "Estado de tabaquismo",
     ylab    = "log(FEV)",
     col     = "gray85",
     border  = "gray30",
     lwd     = 1.5,
     boxwex  = 0.45,
     outline = FALSE
)

stripchart(
     logFEV ~ Smoke_f,
     data     = dat,
     vertical = TRUE,
     method   = "jitter",
     pch      = 16,
     cex      = 0.6,
     col      = adjustcolor("black", alpha.f = 0.25),
     add      = TRUE
)

points(
     x   = 1:2,
     y   = tapply(dat$logFEV, dat$Smoke_f, mean, na.rm = TRUE),
     pch = 19,
     cex = 1.2
)

dev.off()

# Boxplot de log(FEV) por género -----------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_boxplot_fev_gender.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

boxplot(
     logFEV ~ Gender_f,
     data    = dat,
     main    = "",
     xlab    = "Género",
     ylab    = "log(FEV)",
     col     = "gray85",
     border  = "gray30",
     lwd     = 1.5,
     boxwex  = 0.45,
     outline = FALSE
)

stripchart(
     logFEV ~ Gender_f,
     data     = dat,
     vertical = TRUE,
     method   = "jitter",
     pch      = 16,
     cex      = 0.6,
     col      = adjustcolor("black", alpha.f = 0.25),
     add      = TRUE
)

points(
     x   = 1:2,
     y   = tapply(dat$logFEV, dat$Gender_f, mean, na.rm = TRUE),
     pch = 19,
     cex = 1.2
)

dev.off()

# Diagrama de dispersión de log(FEV) y edad ------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_scatter_fev_age.pdf",
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
     dat$Age,
     dat$logFEV,
     pch  = ifelse(dat$Smoke == 1, 17, 16),
     col  = ifelse(dat$Smoke == 1, "gray20", "gray85"),
     xlab = "Edad en años",
     ylab = "log(FEV)",
     main = ""
)

legend(
     "topleft",
     legend = c("No fumador", "Fumador"),
     pch    = c(16, 17),
     col    = c("gray85", "gray20"),
     bty    = "n"
)

dev.off()

# Diagrama de dispersión de log(FEV) y estatura --------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_scatter_fev_ht.pdf",
     width     = 5,
     height    = 5,
     pointsize = 13
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     dat$Ht,
     dat$logFEV,
     pch  = ifelse(dat$Smoke == 1, 17, 16),
     col  = ifelse(dat$Smoke == 1, "gray20", "gray85"),
     xlab = "Estatura en pulgadas",
     ylab = "log(FEV)",
     main = ""
)

legend(
     "topleft",
     legend = c("No fumador", "Fumador"),
     pch    = c(16, 17),
     col    = c("gray85", "gray20"),
     bty    = "n"
)

dev.off()

# Respuesta y matriz de diseño -------------------------------------------------

# Variable respuesta en escala logarítmica
dat$logFEV <- log(dat$FEV)

y <- as.numeric(dat$logFEV)

# Matriz de diseño
X <- cbind(
     intercepto = 1,
     Age        = dat$Age,
     Ht         = dat$Ht,
     Gender     = dat$Gender,
     Smoke      = dat$Smoke
)

X <- as.matrix(X)

# Dimensiones
n <- length(y)
p <- ncol(X)

# Mínimos cuadrados ordinarios -------------------------------------------------

# Estimación de los coeficientes beta por OLS
XtX <- crossprod(X)
Xty <- crossprod(X, y)

beta_ols <- solve(XtX, Xty)

rownames(beta_ols) <- colnames(X)
colnames(beta_ols) <- "Estimación"

print(round(beta_ols, 3))

# Estimación de la varianza residual por OLS
residuos <- as.numeric(y - X %*% beta_ols)

sig2_ols <- sum(residuos^2) / (n - p)

print(round(sig2_ols, 3))

# Muestreador de Gibbs ---------------------------------------------------------

mcmc <- function(
          y, X, beta0, Sigma0, nu0, sigma20,
          n_sams, n_burn, n_skip, verbose = TRUE
) {
     # Ajustes
     y <- as.numeric(y)
     X <- as.matrix(X)
     n <- nrow(X)
     p <- ncol(X)
     
     beta0  <- as.numeric(beta0)
     Sigma0 <- as.matrix(Sigma0)
     
     # Número de iteraciones
     B    <- n_burn + n_sams * n_skip
     ncat <- max(1, floor(0.1 * B))
     
     # Cantidades fijas
     XtX     <- crossprod(X)
     Xty     <- crossprod(X, y)
     iSigma0 <- solve(Sigma0)
     
     # Valores iniciales
     beta   <- as.numeric(qr.solve(X, y))
     resid  <- as.numeric(y - X %*% beta)
     sigma2 <- sum(resid^2) / (n - p)
     
     # Almacenamiento
     BETA   <- matrix(data = NA, nrow = n_sams, ncol = p)
     SIGMA2 <- rep(NA, n_sams)
     LL     <- rep(NA, n_sams)
     
     colnames(BETA) <- colnames(X)
     
     # Cadena
     for (i in 1:B) {
          # Actualizar beta
          V_beta <- solve(iSigma0 + XtX / sigma2)
          m_beta <- V_beta %*% (iSigma0 %*% beta0 + Xty / sigma2)
          
          beta <- as.numeric(
               MASS::mvrnorm(
                    n     = 1,
                    mu    = as.numeric(m_beta),
                    Sigma = V_beta
               )
          )
          
          # Actualizar sigma2
          resid   <- as.numeric(y - X %*% beta)
          rss     <- sum(resid^2)
          a_sigma <- (nu0 + n) / 2
          b_sigma <- (nu0 * sigma20 + rss) / 2
          
          sigma2 <- 1 / rgamma(n = 1, shape = a_sigma, rate = b_sigma)
          
          # Almacenar y log-verosimilitud
          if (i > n_burn && (i - n_burn) %% n_skip == 0) {
               k <- (i - n_burn) / n_skip
               
               ll <- sum(
                    dnorm(
                         x    = y,
                         mean = as.numeric(X %*% beta),
                         sd   = sqrt(sigma2),
                         log  = TRUE
                    )
               )
               
               BETA[k, ] <- beta
               SIGMA2[k] <- sigma2
               LL[k]     <- ll
          }
          
          # Progreso
          if (verbose && i %% ncat == 0) {
               cat(sprintf("%.1f%% completado\n", 100 * i / B))
          }
     }
     
     # Salida
     return(
          list(
               BETA   = BETA,
               SIGMA2 = SIGMA2,
               LL     = LL
          )
     )
}

# Ajuste del modelo lineal Normal ----------------------------------------------

# Hiperparámetros
beta0   <- rep(0, p)
Sigma0  <- diag(100, p)
nu0     <- 1
sigma20 <- sig2_ols

# Número de iteraciones
n_sams <- 10000
n_burn <- 10000
n_skip <- 10

# Ajuste del modelo Bayesiano
set.seed(123)

chain <- mcmc(
     y       = y,
     X       = X,
     beta0   = beta0,
     Sigma0  = Sigma0,
     nu0     = nu0,
     sigma20 = sigma20,
     n_sams  = n_sams,
     n_burn  = n_burn,
     n_skip  = n_skip,
     verbose = TRUE
)

# Traza de la log-verosimilitud
plot(
     chain$LL,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = ""
)

dev.off()

# Trazas de los parámetros del modelo
par(
     mfrow = c(2, 3),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

for (j in 1:ncol(chain$BETA)) {
     plot(
          chain$BETA[, j],
          type = "p",
          cex  = 0.3,
          col  = adjustcolor(1, 0.5),
          main = colnames(chain$BETA)[j],
          xlab = "Iteración",
          ylab = ""
     )
}

plot(
     chain$SIGMA2,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     main = expression(sigma^2),
     xlab = "Iteración",
     ylab = ""
)

dev.off()

# Inferencia -------------------------------------------------------------------

# Muestras posteriores conjuntas de beta y sigma2
POST <- cbind(chain$BETA, sigma2 = chain$SIGMA2)

# Tabla resumen para beta y sigma2
post_media <- apply(POST, 2, mean)
post_sd    <- apply(POST, 2, sd)
post_q025  <- apply(POST, 2, quantile, probs = 0.025)
post_q975  <- apply(POST, 2, quantile, probs = 0.975)

post_summary <- data.frame(
     media = post_media,
     sd    = post_sd,
     cv    = post_sd / abs(post_media),
     q025  = post_q025,
     q975  = post_q975
)

print(round(post_summary, 3))

# Probabilidades posteriores útiles para la interpretación
prob_summary <- data.frame(
     prob_mayor_0 = apply(chain$BETA, 2, function(x) mean(x > 0)),
     prob_menor_0 = apply(chain$BETA, 2, function(x) mean(x < 0))
)

print(round(prob_summary, 3))

# Gráfico de intervalos de credibilidad para los coeficientes
beta_mean <- apply(chain$BETA[, -1], 2, mean)
beta_q025 <- apply(chain$BETA[, -1], 2, quantile, probs = 0.025)
beta_q975 <- apply(chain$BETA[, -1], 2, quantile, probs = 0.975)

pdf(
     file      = "lungcap_regresion_lineal_normal_intervalos_beta.pdf",
     width     = 6,
     height    = 5,
     pointsize = 14
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 4, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

y_pos <- seq_along(beta_mean)

plot(
     beta_mean,
     y_pos,
     type = "n",
     xlim = range(beta_q025, beta_q975),
     ylim = range(y_pos),
     yaxt = "n",
     xlab = "Coeficiente",
     ylab = "",
     main = ""
)

abline(
     v   = 0,
     lty = 2,
     col = "gray60",
     lwd = 2
)

segments(
     x0  = beta_q025,
     y0  = y_pos,
     x1  = beta_q975,
     y1  = y_pos,
     lwd = 2.5,
     col = "gray35"
)

points(
     beta_mean,
     y_pos,
     pch = 19,
     cex = 1.2,
     col = "black"
)

axis(
     side   = 2,
     at     = y_pos,
     labels = names(beta_mean),
     las    = 1
)

dev.off()

# Evaluación del ajuste del modelo ---------------------------------------------

# Métricas de ajuste por iteración posterior
B <- nrow(chain$BETA)

R2   <- rep(NA, B)
RMSE <- rep(NA, B)
MAE  <- rep(NA, B)

y_bar <- mean(y)
tss   <- sum((y - y_bar)^2)

for (b in 1:B) {
     yhat_b <- as.numeric(X %*% chain$BETA[b, ])
     e_b    <- y - yhat_b
     
     RMSE[b] <- sqrt(mean(e_b^2))
     MAE[b]  <- mean(abs(e_b))
     R2[b]   <- 1 - sum(e_b^2) / tss
}

# Resumen posterior
FIT <- cbind(
     R2   = R2,
     RMSE = RMSE,
     MAE  = MAE
)

fit_summary <- data.frame(
     media = apply(FIT, 2, mean),
     cv    = apply(FIT, 2, sd) / abs(apply(FIT, 2, mean)),
     q025  = apply(FIT, 2, quantile, probs = 0.025),
     q975  = apply(FIT, 2, quantile, probs = 0.975)
)

print(round(fit_summary, 3))

# Diagnósticos del modelo ajustado a log(FEV) ----------------------------------

# Cantidades posteriores básicas
B <- nrow(chain$BETA)

beta_post_mean   <- colMeans(chain$BETA)
sigma2_post_mean <- mean(chain$SIGMA2)

yhat <- as.numeric(X %*% beta_post_mean)
res  <- as.numeric(y - yhat)

res_std <- res / sqrt(sigma2_post_mean)

# Matriz sombrero y medidas de influencia
H_mat <- X %*% solve(crossprod(X)) %*% t(X)
h_ii  <- diag(H_mat)

cook <- (res^2 / (p * sigma2_post_mean)) * (h_ii / (1 - h_ii)^2)

# Resumen numérico de residuos en escala logarítmica
resumen_residuos <- c(
     media = mean(res),
     sd    = sd(res),
     q025  = unname(quantile(res, 0.025)),
     q500  = unname(quantile(res, 0.500)),
     q975  = unname(quantile(res, 0.975)),
     min   = min(res),
     max   = max(res)
)

print(round(resumen_residuos, 3))

# Métricas simples para revisar patrones residuales
metricas_supuestos <- c(
     cor_res_yhat             = cor(res, yhat),
     cor_absres_yhat          = cor(abs(res), yhat),
     cor_absres_age           = cor(abs(res), dat$Age),
     cor_absres_ht            = cor(abs(res), dat$Ht),
     max_leverage             = max(h_ii),
     max_cook                 = max(cook),
     prop_abs_resstd_mayor_2  = mean(abs(res_std) > 2),
     prop_abs_resstd_mayor_3  = mean(abs(res_std) > 3)
)

print(round(metricas_supuestos, 3))

# Residuos vs valores ajustados ------------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_residuos_ajustados.pdf",
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
     yhat,
     res,
     type = "n",
     xlab = "Valores ajustados de log(FEV)",
     ylab = "Residuos",
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     col = "gray60",
     lwd = 1.5
)

points(
     yhat,
     res,
     pch = 16,
     cex = 0.65,
     col = adjustcolor("black", alpha.f = 0.45)
)

lines(
     lowess(yhat, res),
     lwd = 2,
     col = "gray20"
)

dev.off()

# Q-Q plot de residuos estandarizados con bandas de confianza ------------------

res_std_ord <- sort(res_std)
n_res       <- length(res_std_ord)

p_q    <- ppoints(n_res)
q_theo <- qnorm(p_q)

q_y <- quantile(res_std, probs = c(0.25, 0.75), na.rm = TRUE)
q_x <- qnorm(c(0.25, 0.75))

slope_qq <- diff(q_y) / diff(q_x)
int_qq   <- q_y[1] - slope_qq * q_x[1]

pdf(
     file      = "lungcap_regresion_lineal_normal_qqplot_residuos.pdf",
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
     q_theo,
     res_std_ord,
     type = "n",
     xlab = "Cuantiles teóricos",
     ylab = "Cuantiles muestrales",
     main = ""
)

abline(
     a   = int_qq,
     b   = slope_qq,
     lwd = 2,
     col = "gray35"
)

points(
     q_theo,
     res_std_ord,
     pch = 16,
     cex = 0.65,
     col = adjustcolor("black", alpha.f = 0.55)
)

dev.off()

# Escala-localización -----------------------------------------------------------

sqrt_abs_res_std <- sqrt(abs(res_std))

pdf(
     file      = "lungcap_regresion_lineal_normal_escala_localizacion.pdf",
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
     yhat,
     sqrt_abs_res_std,
     type = "n",
     xlab = "Valores ajustados de log(FEV)",
     ylab = expression(sqrt("|residuo estandarizado|")),
     main = ""
)

points(
     yhat,
     sqrt_abs_res_std,
     pch = 16,
     cex = 0.65,
     col = adjustcolor("black", alpha.f = 0.45)
)

lines(
     lowess(yhat, sqrt_abs_res_std),
     lwd = 2,
     col = "gray20"
)

dev.off()

# Distancia de Cook ------------------------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_distancia_cook.pdf",
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
     seq_along(cook),
     cook,
     type = "n",
     xlab = "Observación",
     ylab = "Distancia de Cook",
     main = ""
)

segments(
     x0  = seq_along(cook),
     y0  = 0,
     x1  = seq_along(cook),
     y1  = cook,
     lwd = 1.2,
     col = adjustcolor("black", alpha.f = 0.55)
)

abline(
     h   = 4 / n,
     lty = 2,
     col = "gray60",
     lwd = 1.5
)

dev.off()

# Residuos frente a edad -------------------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_residuos_age.pdf",
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
     dat$Age,
     res,
     type = "n",
     xlab = "Edad en años",
     ylab = "Residuos de log(FEV)",
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     col = "gray60",
     lwd = 1.5
)

points(
     dat$Age,
     res,
     pch = 16,
     cex = 0.65,
     col = adjustcolor("black", alpha.f = 0.45)
)

lines(
     lowess(dat$Age, res),
     lwd = 2,
     col = "gray20"
)

dev.off()

# Residuos frente a estatura ---------------------------------------------------

pdf(
     file      = "lungcap_regresion_lineal_normal_residuos_ht.pdf",
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
     dat$Ht,
     res,
     type = "n",
     xlab = "Estatura en pulgadas",
     ylab = "Residuos de log(FEV)",
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     col = "gray60",
     lwd = 1.5
)

points(
     dat$Ht,
     res,
     pch = 16,
     cex = 0.65,
     col = adjustcolor("black", alpha.f = 0.45)
)

lines(
     lowess(dat$Ht, res),
     lwd = 2,
     col = "gray20"
)

dev.off()

# Verificación predictiva posterior con todas las iteraciones ------------------

B <- nrow(chain$BETA)

T_obs <- c(
     media   = mean(y),
     mediana = median(y),
     sd      = sd(y),
     riq     = IQR(y),
     min     = min(y),
     max     = max(y)
)

T_rep <- matrix(
     data = NA,
     nrow = B,
     ncol = length(T_obs)
)

colnames(T_rep) <- names(T_obs)

set.seed(123)

for (b in seq_len(B)) {
     beta_b   <- chain$BETA[b, ]
     sigma2_b <- chain$SIGMA2[b]
     
     yhat_b <- as.numeric(X %*% beta_b)
     
     yrep_b <- rnorm(
          n    = n,
          mean = yhat_b,
          sd   = sqrt(sigma2_b)
     )
     
     T_rep[b, ] <- c(
          media   = mean(yrep_b),
          mediana = median(yrep_b),
          sd      = sd(yrep_b),
          riq     = IQR(yrep_b),
          min     = min(yrep_b),
          max     = max(yrep_b)
     )
}

ppp <- colMeans(
     sweep(T_rep, 2, T_obs, FUN = "<=")
)

ppc_summary <- data.frame(
     observado          = as.numeric(T_obs),
     media_rep          = colMeans(T_rep),
     q025_rep           = apply(T_rep, 2, quantile, probs = 0.025),
     q975_rep           = apply(T_rep, 2, quantile, probs = 0.975),
     valor_p_predictivo = ppp
)

print(round(ppc_summary, 3))

# Gráfico de verificación predictiva posterior --------------------------------

# Histogramas individuales para la verificación predictiva posterior -----------

nombres_archivos <- c(
     media   = "lungcap_regresion_lineal_normal_ppc_media.pdf",
     mediana = "lungcap_regresion_lineal_normal_ppc_mediana.pdf",
     sd      = "lungcap_regresion_lineal_normal_ppc_sd.pdf",
     riq     = "lungcap_regresion_lineal_normal_ppc_riq.pdf",
     min     = "lungcap_regresion_lineal_normal_ppc_min.pdf",
     max     = "lungcap_regresion_lineal_normal_ppc_max.pdf"
)

etiquetas_x <- c(
     media   = "Media de log(FEV)",
     mediana = "Mediana de log(FEV)",
     sd      = "Desviación estándar de log(FEV)",
     riq     = "Rango intercuartílico de log(FEV)",
     min     = "Mínimo de log(FEV)",
     max     = "Máximo de log(FEV)"
)

for (j in colnames(T_rep)) {
     pdf(
          file      = nombres_archivos[j],
          width     = 5,
          height    = 5,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     hist_j <- hist(
          T_rep[, j],
          breaks = 25,
          plot   = FALSE
     )
     
     dens_j <- density(T_rep[, j])
     
     ylim_j <- c(
          0,
          1.05 * max(hist_j$density, dens_j$y, na.rm = TRUE)
     )
     
     xlim_j <- range(
          hist_j$breaks,
          T_obs[j],
          na.rm = TRUE
     )
     
     plot(
          NA,
          xlim = xlim_j,
          ylim = ylim_j,
          xlab = etiquetas_x[j],
          ylab = "Densidad",
          main = "",
          axes = TRUE
     )
     
     rect(
          xleft   = hist_j$breaks[-length(hist_j$breaks)],
          ybottom = 0,
          xright  = hist_j$breaks[-1],
          ytop    = hist_j$density,
          col     = "gray85",
          border  = "white"
     )
     
     abline(
          v   = T_obs[j],
          lwd = 2,
          lty = 2,
          col = "gray40"
     )
     
     usr <- par("usr")
     
     text(
          x      = usr[2] - 0.02 * diff(usr[1:2]),
          y      = usr[4] - 0.06 * diff(usr[3:4]),
          labels = paste0(
               "p = ",
               formatC(
                    unname(ppp[j]),
                    format = "f",
                    digits = 3
               )
          ),
          adj = c(1, 1),
          cex = 1.1
     )
     
     box()
     
     dev.off()
}

# Fin --------------------------------------------------------------------------