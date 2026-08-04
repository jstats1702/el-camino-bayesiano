# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos simulados --------------------------------------------------------------

set.seed(123)

n       <- 80
p       <- 50
n_train <- round(0.75 * n)

rho <- 0.75

Sigma_x <- outer(
     seq_len(p),
     seq_len(p),
     function(i, j) rho^abs(i - j)
)

Z <- matrix(
     rnorm(n * p),
     nrow = n,
     ncol = p
)

X_raw <- Z %*% chol(Sigma_x)

colnames(X_raw) <- paste0("x", seq_len(p))

id_train <- sample(
     seq_len(n),
     size    = n_train,
     replace = FALSE
)

id_test <- setdiff(
     seq_len(n),
     id_train
)

# Escalamiento usando solo el conjunto de entrenamiento
x_mean <- colMeans(X_raw[id_train, , drop = FALSE])
x_sd   <- apply(X_raw[id_train, , drop = FALSE], 2, sd)

X <- sweep(X_raw, 2, x_mean, "-")
X <- sweep(X, 2, x_sd, "/")

beta_true <- rep(0, p)

beta_true[c(1, 2, 5, 12, 25, 40)] <- c(
     3.0,
     -2.5,
     2.0,
     -1.8,
     1.5,
     -1.2
)

names(beta_true) <- colnames(X)

sigma_true <- 2.5

y <- as.numeric(
     X %*% beta_true +
          rnorm(n, mean = 0, sd = sigma_true)
)

X_train <- X[id_train, , drop = FALSE]
y_train <- y[id_train]

X_test <- X[id_test, , drop = FALSE]
y_test <- y[id_test]

# Muestreador de Gibbs previa difusa -------------------------------------------

gibbs_difusa <- function(
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

# Muestreador de Gibbs regresión Ridge -----------------------------------------

gibbs_ridge <- function(
          y,
          X,
          B,
          n_burn,
          n_skip,
          a_sigma,
          b_sigma,
          a_lambda,
          b_lambda,
          verbose = TRUE
) {
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n          <- length(y)
     p          <- ncol(X)
     total_iter <- n_burn + B * n_skip
     
     XtX <- crossprod(X)
     Xty <- crossprod(X, y)
     Ip  <- diag(p)
     
     beta_samples <- matrix(
          NA_real_,
          nrow = B,
          ncol = p
     )
     
     sigma2_samples  <- numeric(B)
     lambda2_samples <- numeric(B)
     loglik_samples  <- numeric(B)
     
     beta    <- rep(0, p)
     sigma2  <- var(y)
     lambda2 <- 1
     
     progress_step <- max(1, floor(total_iter / 10))
     
     for (b in seq_len(total_iter)) {
          # Actualizar beta usando la precisión posterior
          precision_beta <- XtX + lambda2 * Ip
          chol_precision  <- chol(precision_beta)
          
          V_beta <- chol2inv(chol_precision)
          m_beta <- V_beta %*% Xty
          
          beta <- as.numeric(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = as.numeric(m_beta),
                    sigma = sigma2 * V_beta
               )
          )
          
          # Actualizar sigma^2
          resid <- y - drop(X %*% beta)
          
          sigma2 <- 1 / rgamma(
               n     = 1,
               shape = a_sigma + 0.5 * (n + p),
               rate  = b_sigma + 0.5 * (sum(resid^2) + lambda2 * sum(beta^2))
          )
          
          # Actualizar lambda^2
          lambda2 <- rgamma(
               n     = 1,
               shape = a_lambda + 0.5 * p,
               rate  = b_lambda + 0.5 * sum(beta^2) / sigma2
          )
          
          if (b > n_burn && (b - n_burn) %% n_skip == 0) {
               i <- (b - n_burn) / n_skip
               
               beta_samples[i, ]  <- beta
               sigma2_samples[i]  <- sigma2
               lambda2_samples[i] <- lambda2
               
               loglik_samples[i] <- sum(
                    dnorm(
                         x    = y,
                         mean = drop(X %*% beta),
                         sd   = sqrt(sigma2),
                         log  = TRUE
                    )
               )
          }
          
          if (verbose && b %% progress_step == 0) {
               cat(
                    "Iteración", b, "de", total_iter,
                    paste0("(", round(100 * b / total_iter), "%)\n")
               )
          }
     }
     
     beta_names <- colnames(X)
     
     if (is.null(beta_names)) {
          beta_names <- paste0("beta", seq_len(p))
     }
     
     colnames(beta_samples) <- beta_names
     
     list(
          beta    = beta_samples,
          sigma2  = sigma2_samples,
          lambda2 = lambda2_samples,
          loglik  = loglik_samples
     )
}

# Muestreador de Gibbs regresión Lasso -----------------------------------------

# Función para simular una Gaussiana inversa
# Parametrización: media = mu, parámetro de forma = lambda
rinvgauss_msh <- function(n, mu, lambda) {
     z <- rnorm(n)^2
     
     x <- mu +
          (mu^2 * z) / (2 * lambda) -
          (mu / (2 * lambda)) *
          sqrt(4 * mu * lambda * z + mu^2 * z^2)
     
     u <- runif(n)
     
     ifelse(
          u <= mu / (mu + x),
          x,
          mu^2 / x
     )
}

gibbs_lasso <- function(
          y,
          X,
          B,
          n_burn,
          n_skip,
          a_sigma,
          b_sigma,
          a_lambda,
          b_lambda,
          verbose = TRUE
) {
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n          <- length(y)
     p          <- ncol(X)
     total_iter <- n_burn + B * n_skip
     
     XtX <- crossprod(X)
     Xty <- crossprod(X, y)
     
     beta_samples <- matrix(
          NA_real_,
          nrow = B,
          ncol = p
     )
     
     sigma2_samples  <- numeric(B)
     lambda2_samples <- numeric(B)
     loglik_samples  <- numeric(B)
     
     beta    <- rep(0, p)
     sigma2  <- var(y)
     lambda2 <- 1
     tau2    <- rep(1, p)
     
     progress_step <- max(1, floor(total_iter / 10))
     
     for (b in seq_len(total_iter)) {
          # Actualizar beta usando la precisión posterior
          Dtau_inv       <- diag(1 / tau2, nrow = p, ncol = p)
          precision_beta <- XtX + Dtau_inv
          chol_precision <- chol(precision_beta)
          
          V_beta <- chol2inv(chol_precision)
          m_beta <- V_beta %*% Xty
          
          beta <- as.numeric(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = as.numeric(m_beta),
                    sigma = sigma2 * V_beta
               )
          )
          
          # Actualizar las varianzas locales tau_j^2
          mu_inv_tau2 <- sqrt(
               lambda2 * sigma2 / pmax(beta^2, .Machine$double.eps)
          )
          
          inv_tau2 <- rinvgauss_msh(
               n      = p,
               mu     = mu_inv_tau2,
               lambda = lambda2
          )
          
          tau2 <- 1 / inv_tau2
          
          # Actualizar sigma^2
          resid <- y - drop(X %*% beta)
          
          sigma2 <- 1 / rgamma(
               n     = 1,
               shape = a_sigma + 0.5 * (n + p),
               rate  = b_sigma + 0.5 * (sum(resid^2) + sum(beta^2 / tau2))
          )
          
          # Actualizar lambda^2
          lambda2 <- rgamma(
               n     = 1,
               shape = a_lambda + p,
               rate  = b_lambda + 0.5 * sum(tau2)
          )
          
          if (b > n_burn && (b - n_burn) %% n_skip == 0) {
               i <- (b - n_burn) / n_skip
               
               beta_samples[i, ]  <- beta
               sigma2_samples[i]  <- sigma2
               lambda2_samples[i] <- lambda2
               
               loglik_samples[i] <- sum(
                    dnorm(
                         x    = y,
                         mean = drop(X %*% beta),
                         sd   = sqrt(sigma2),
                         log  = TRUE
                    )
               )
          }
          
          if (verbose && b %% progress_step == 0) {
               cat(
                    "Iteración", b, "de", total_iter,
                    paste0("(", round(100 * b / total_iter), "%)\n")
               )
          }
     }
     
     beta_names <- colnames(X)
     
     if (is.null(beta_names)) {
          beta_names <- paste0("beta", seq_len(p))
     }
     
     colnames(beta_samples) <- beta_names
     
     list(
          beta    = beta_samples,
          sigma2  = sigma2_samples,
          lambda2 = lambda2_samples,
          loglik  = loglik_samples
     )
}

# Ajuste del modelo previa difusa ----------------------------------------------

# OLS
XtX      <- crossprod(X_train)
Xty      <- crossprod(X_train, y_train)
beta_ols <- solve(XtX, Xty)
residuos <- as.numeric(y_train - X_train %*% beta_ols)
sig2_ols <- sum(residuos^2) / (n_train - p)

print(round(sig2_ols, 3))

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

muestras_difusa <- gibbs_difusa(
     y       = y_train,
     X       = X_train,
     beta0   = beta0,
     Sigma0  = Sigma0,
     nu0     = nu0,
     sigma20 = sigma20,
     n_sams  = n_sams,
     n_burn  = n_burn,
     n_skip  = n_skip,
     verbose = TRUE
)

# Ajuste del modelo Ridge ------------------------------------------------------

# Hiperparámetros
a_sigma <- 2
b_sigma <- sig2_ols

a_lambda <- 1
b_lambda <- 1

# Configuración del muestreador
B      <- 10000
n_burn <- 10000
n_skip <- 10

# Ajuste del modelo
set.seed(123)

muestras_ridge <- gibbs_ridge(
     y        = y_train,
     X        = X_train,
     B        = B,
     n_burn   = n_burn,
     n_skip   = n_skip,
     a_sigma  = a_sigma,
     b_sigma  = b_sigma,
     a_lambda = a_lambda,
     b_lambda = b_lambda,
     verbose  = TRUE
)

# Ajuste del modelo Lasso ------------------------------------------------------

# Hiperparámetros
a_sigma <- 2
b_sigma <- sig2_ols

a_lambda <- 1
b_lambda <- 1

# Configuración del muestreador
B      <- 10000
n_burn <- 10000
n_skip <- 10

# Ajuste del modelo
set.seed(123)

muestras_lasso <- gibbs_lasso(
     y        = y_train,
     X        = X_train,
     B        = B,
     n_burn   = n_burn,
     n_skip   = n_skip,
     a_sigma  = a_sigma,
     b_sigma  = b_sigma,
     a_lambda = a_lambda,
     b_lambda = b_lambda,
     verbose  = TRUE
)

# Convergencia -----------------------------------------------------------------

# Cadena log-verosimilitud previa difusa
plot(
     muestras_difusa$LL,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = ""
)

# Cadena log-verosimilitud Ridge
plot(
     muestras_ridge$loglik,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = ""
)

# Cadena log-verosimilitud Lasso
plot(
     muestras_lasso$loglik,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = ""
)

dev.off()

# Comparación de beta verdadero con beta estimado ------------------------------

beta_true <- as.numeric(beta_true)

beta_difusa <- as.matrix(muestras_difusa$BETA)
beta_ridge  <- as.matrix(muestras_ridge$beta)
beta_lasso  <- as.matrix(muestras_lasso$beta)

muestras_beta <- list(
     Difusa = beta_difusa,
     Ridge  = beta_ridge,
     Lasso  = beta_lasso
)

resumir_estimacion_beta <- function(beta_samples, beta_true, modelo) {
     beta_hat <- colMeans(beta_samples)
     errores  <- beta_hat - beta_true
     
     q025 <- apply(
          beta_samples,
          2,
          quantile,
          probs = 0.025
     )
     
     q975 <- apply(
          beta_samples,
          2,
          quantile,
          probs = 0.975
     )
     
     longitud_ic <- q975 - q025
     
     post_sd <- apply(
          beta_samples,
          2,
          sd
     )
     
     data.frame(
          Modelo           = modelo,
          RMSE             = sqrt(mean(errores^2)),
          MAE              = mean(abs(errores)),
          L2               = sqrt(sum(errores^2)),
          MaxAE            = max(abs(errores)),
          Sesgo            = abs(mean(errores)),
          Longitud_IC95    = mean(longitud_ic),
          SD_promedio      = mean(post_sd),
          Cobertura_IC95   = mean(q025 <= beta_true & beta_true <= q975),
          Correlacion_beta = cor(beta_hat, beta_true),
          row.names        = NULL
     )
}

resumen_beta_media <- do.call(
     rbind,
     lapply(
          names(muestras_beta),
          function(modelo) {
               resumir_estimacion_beta(
                    beta_samples = muestras_beta[[modelo]],
                    beta_true    = beta_true,
                    modelo       = modelo
               )
          }
     )
)

resumen_beta_media[, -1] <- round(
     resumen_beta_media[, -1],
     3
)

print(resumen_beta_media)

# Predicción en el conjunto de prueba -----------------------------------------

X_test <- as.matrix(X_test)
y_test <- as.numeric(y_test)

beta_hat_difusa <- colMeans(as.matrix(muestras_difusa$BETA))
beta_hat_ridge  <- colMeans(as.matrix(muestras_ridge$beta))
beta_hat_lasso  <- colMeans(as.matrix(muestras_lasso$beta))

y_hat_difusa <- as.numeric(X_test %*% beta_hat_difusa)
y_hat_ridge  <- as.numeric(X_test %*% beta_hat_ridge)
y_hat_lasso  <- as.numeric(X_test %*% beta_hat_lasso)

# Métricas predictivas
calcular_metricas_pred <- function(y_obs, y_hat, modelo) {
     error <- y_obs - y_hat
     
     data.frame(
          Modelo  = modelo,
          RMSE    = sqrt(mean(error^2)),
          MAE     = mean(abs(error)),
          Sesgo   = mean(error),
          R2_pred = 1 - sum(error^2) / sum((y_obs - mean(y_obs))^2),
          Corr    = cor(y_obs, y_hat)
     )
}

metricas_pred <- rbind(
     calcular_metricas_pred(
          y_obs  = y_test,
          y_hat  = y_hat_difusa,
          modelo = "Difusa"
     ),
     calcular_metricas_pred(
          y_obs  = y_test,
          y_hat  = y_hat_ridge,
          modelo = "Ridge"
     ),
     calcular_metricas_pred(
          y_obs  = y_test,
          y_hat  = y_hat_lasso,
          modelo = "Lasso"
     )
)

metricas_pred[, -1] <- round(metricas_pred[, -1], 4)

print(metricas_pred)

# Comparación de selección de coeficientes no nulos ----------------------------

beta_true <- as.numeric(beta_true)

epsilon  <- 0.10
gamma    <- 0.95
nivel_ic <- 0.95

muestras_beta <- list(
     Difusa = as.matrix(muestras_difusa$BETA),
     Ridge  = as.matrix(muestras_ridge$beta),
     Lasso  = as.matrix(muestras_lasso$beta)
)

stopifnot(
     all(
          sapply(muestras_beta, ncol) == length(beta_true)
     )
)

calcular_auc <- function(score, verdadero) {
     verdadero <- as.logical(verdadero)
     n_pos     <- sum(verdadero)
     n_neg     <- sum(!verdadero)
     
     if (n_pos == 0 || n_neg == 0) {
          return(NA_real_)
     }
     
     ranks <- rank(
          score,
          ties.method = "average"
     )
     
     (
          sum(ranks[verdadero]) -
               n_pos * (n_pos + 1) / 2
     ) / (n_pos * n_neg)
}

resumen_seleccion_modelo <- function(
          beta_samples,
          beta_true,
          modelo,
          epsilon,
          gamma,
          nivel_ic
) {
     alfa               <- 1 - nivel_ic
     beta_true          <- as.numeric(beta_true)
     verdadero_no_nulo  <- beta_true != 0
     beta_names         <- colnames(beta_samples)
     
     if (is.null(beta_names)) {
          beta_names <- paste0("beta", seq_along(beta_true))
     }
     
     media   <- colMeans(beta_samples)
     mediana <- apply(beta_samples, 2, median)
     
     q025 <- apply(
          beta_samples,
          2,
          quantile,
          probs = alfa / 2
     )
     
     q975 <- apply(
          beta_samples,
          2,
          quantile,
          probs = 1 - alfa / 2
     )
     
     prob_pos        <- colMeans(beta_samples > 0)
     prob_neg        <- colMeans(beta_samples < 0)
     prob_signo      <- pmax(prob_pos, prob_neg)
     prob_relevancia <- colMeans(abs(beta_samples) > epsilon)
     
     seleccion_ic               <- (q025 > 0) | (q975 < 0)
     seleccion_relevancia       <- prob_relevancia > gamma
     seleccion_signo_relevancia <- (prob_signo > gamma) &
          (prob_relevancia > gamma)
     
     data.frame(
          Modelo                     = modelo,
          Coeficiente                = beta_names,
          beta_true                  = beta_true,
          verdadero_no_nulo          = verdadero_no_nulo,
          media                      = media,
          mediana                    = mediana,
          q025                       = q025,
          q975                       = q975,
          prob_pos                   = prob_pos,
          prob_neg                   = prob_neg,
          prob_signo                 = prob_signo,
          prob_relevancia            = prob_relevancia,
          seleccion_ic               = seleccion_ic,
          seleccion_relevancia       = seleccion_relevancia,
          seleccion_signo_relevancia = seleccion_signo_relevancia,
          row.names                  = NULL
     )
}

safe_div <- function(a, b) {
     ifelse(
          b == 0,
          NA_real_,
          a / b
     )
}

evaluar_regla <- function(tabla, columna_seleccion, nombre_regla) {
     verdadero    <- tabla$verdadero_no_nulo
     seleccionado <- tabla[[columna_seleccion]]
     
     TP <- sum(seleccionado & verdadero)
     FP <- sum(seleccionado & !verdadero)
     TN <- sum(!seleccionado & !verdadero)
     FN <- sum(!seleccionado & verdadero)
     
     data.frame(
          Modelo       = tabla$Modelo[1],
          Regla        = nombre_regla,
          TP           = TP,
          FP           = FP,
          TN           = TN,
          FN           = FN,
          Sensibilidad = safe_div(TP, TP + FN),
          Especificidad = safe_div(TN, TN + FP),
          Precision    = safe_div(TP, TP + FP),
          F1           = safe_div(2 * TP, 2 * TP + FP + FN),
          Exactitud    = safe_div(TP + TN, TP + FP + TN + FN),
          AUC          = calcular_auc(
               score     = tabla$prob_relevancia,
               verdadero = verdadero
          )
     )
}

resumen_coeficientes <- do.call(
     rbind,
     lapply(
          names(muestras_beta),
          function(modelo) {
               resumen_seleccion_modelo(
                    beta_samples = muestras_beta[[modelo]],
                    beta_true    = beta_true,
                    modelo       = modelo,
                    epsilon      = epsilon,
                    gamma        = gamma,
                    nivel_ic     = nivel_ic
               )
          }
     )
)

metricas_seleccion <- do.call(
     rbind,
     lapply(
          split(resumen_coeficientes, resumen_coeficientes$Modelo),
          function(tabla) {
               rbind(
                    evaluar_regla(
                         tabla             = tabla,
                         columna_seleccion = "seleccion_ic",
                         nombre_regla      = "IC excluye cero"
                    ),
                    evaluar_regla(
                         tabla             = tabla,
                         columna_seleccion = "seleccion_relevancia",
                         nombre_regla      = paste0(
                              "P(|beta| > ", epsilon, ") > ", gamma
                         )
                    ),
                    evaluar_regla(
                         tabla             = tabla,
                         columna_seleccion = "seleccion_signo_relevancia",
                         nombre_regla      = paste0(
                              "Signo y relevancia > ", gamma
                         )
                    )
               )
          }
     )
)

metricas_seleccion[, c(
     "Sensibilidad",
     "Especificidad",
     "Precision",
     "F1",
     "Exactitud",
     "AUC"
)] <- round(
     metricas_seleccion[, c(
          "Sensibilidad",
          "Especificidad",
          "Precision",
          "F1",
          "Exactitud",
          "AUC"
     )],
     4
)

resumen_coeficientes[, c(
     "media",
     "mediana",
     "q025",
     "q975",
     "prob_pos",
     "prob_neg",
     "prob_signo",
     "prob_relevancia"
)] <- round(
     resumen_coeficientes[, c(
          "media",
          "mediana",
          "q025",
          "q975",
          "prob_pos",
          "prob_neg",
          "prob_signo",
          "prob_relevancia"
     )],
     4
)

print(metricas_seleccion)

print(
     resumen_coeficientes[
          order(
               resumen_coeficientes$Modelo,
               -abs(resumen_coeficientes$beta_true),
               resumen_coeficientes$Coeficiente
          ),
     ]
)

# Intervalos de credibilidad posteriores por modelo ----------------------------

# Muestras posteriores
beta_difusa <- as.matrix(muestras_difusa$BETA)
beta_ridge  <- as.matrix(muestras_ridge$beta)
beta_lasso  <- as.matrix(muestras_lasso$beta)

beta_true <- as.numeric(beta_true)

p     <- length(beta_true)
x_pos <- seq_len(p)

# Estadísticas posteriores
media_difusa <- colMeans(beta_difusa)

ic95_difusa <- apply(
     beta_difusa,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

media_ridge <- colMeans(beta_ridge)

ic95_ridge <- apply(
     beta_ridge,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

media_lasso <- colMeans(beta_lasso)

ic95_lasso <- apply(
     beta_lasso,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

ylim_all <- range(
     c(
          ic95_difusa,
          ic95_ridge,
          ic95_lasso,
          beta_true
     ),
     na.rm = TRUE
)

# Función auxiliar para evitar repetir código
graficar_intervalos_beta <- function(media, ic95, beta_true, titulo_modelo) {
     plot(
          NA,
          NA,
          xlab = "Indice del coeficiente",
          ylab = expression(beta[j]),
          xlim = c(0.5, length(beta_true) + 0.5),
          ylim = ylim_all,
          main = ""
     )
     
     abline(
          h   = 0,
          col = "gray70",
          lwd = 2
     )
     
     for (j in seq_along(beta_true)) {
          segments(
               x0  = j,
               y0  = ic95[1, j],
               x1  = j,
               y1  = ic95[2, j],
               lwd = 1
          )
          
          points(
               x   = j,
               y   = media[j],
               pch = 16,
               cex = 0.6
          )
          
          points(
               x   = j,
               y   = beta_true[j],
               pch = 15,
               col = "red3",
               cex = 0.8
          )
     }
     
     legend(
          "topright",
          legend = titulo_modelo,
          bty    = "n",
          cex    = 1.5
     )
}

# Gráfico 1: modelo semiconjugado con previa difusa
pdf(
     file      = "simulacion_regularizacion_difusa.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

graficar_intervalos_beta(
     media         = media_difusa,
     ic95          = ic95_difusa,
     beta_true     = beta_true,
     titulo_modelo = "Previa difusa"
)

dev.off()

# Gráfico 2: Ridge Bayesiano
pdf(
     file      = "simulacion_regularizacion_ridge.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

graficar_intervalos_beta(
     media         = media_ridge,
     ic95          = ic95_ridge,
     beta_true     = beta_true,
     titulo_modelo = "Ridge Bayesiano"
)

dev.off()

# Gráfico 3: Lasso Bayesiano
pdf(
     file      = "simulacion_regularizacion_lasso.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

graficar_intervalos_beta(
     media         = media_lasso,
     ic95          = ic95_lasso,
     beta_true     = beta_true,
     titulo_modelo = "Lasso Bayesiano"
)

dev.off()

# Resultados -------------------------------------------------------------------

print(resumen_beta_media)

print(metricas_pred)

print(metricas_seleccion)

# Fin --------------------------------------------------------------------------