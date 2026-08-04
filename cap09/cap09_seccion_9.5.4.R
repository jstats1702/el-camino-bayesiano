# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

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

# Muestreador de Gibbs selección de variables ----------------------------------

lpy_X <- function(y, X, g, nu0, s20) {
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n <- length(y)
     p <- ncol(X)
     
     if (p == 0) {
          ssr_g <- sum(y^2)
          s20   <- mean(y^2)
     } else {
          XtX_inv <- solve(crossprod(X))
          H_g     <- (g / (g + 1)) * X %*% XtX_inv %*% t(X)
          ssr_g   <- as.numeric(crossprod(y, (diag(n) - H_g) %*% y))
     }
     
     log_marg_lik <-
          -0.5 * n * log(2 * pi) +
          lgamma(0.5 * (nu0 + n)) -
          lgamma(0.5 * nu0) -
          0.5 * p * log(1 + g) +
          0.5 * nu0 * log(0.5 * nu0 * s20) -
          0.5 * (nu0 + n) * log(0.5 * (nu0 * s20 + ssr_g))
     
     as.numeric(log_marg_lik)
}

lm_gprior <- function(y, X, g, nu0, s20, B) {
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n <- length(y)
     p <- ncol(X)
     
     if (p == 0) {
          ssr_g <- sum(y^2)
          
          s2 <- 1 / rgamma(
               B,
               shape = (nu0 + n) / 2,
               rate  = (nu0 * s20 + ssr_g) / 2
          )
          
          return(
               list(
                    beta = matrix(numeric(0), nrow = B, ncol = 0),
                    s2   = s2
               )
          )
     }
     
     # Cantidades posteriores bajo el prior g
     XtX_inv <- solve(crossprod(X))
     H_g     <- (g / (g + 1)) * X %*% XtX_inv %*% t(X)
     ssr_g   <- as.numeric(crossprod(y, (diag(n) - H_g) %*% y))
     
     V_beta <- (g / (g + 1)) * XtX_inv
     E_beta <- as.numeric(V_beta %*% crossprod(X, y))
     
     # Simulación posterior de sigma^2
     s2 <- 1 / rgamma(
          B,
          shape = (nu0 + n) / 2,
          rate  = (nu0 * s20 + ssr_g) / 2
     )
     
     # Simulación posterior de beta
     Z    <- matrix(rnorm(B * p), nrow = B, ncol = p)
     beta <- sweep(Z %*% chol(V_beta), 1, sqrt(s2), "*")
     beta <- sweep(beta, 2, E_beta, "+")
     
     list(
          beta = beta,
          s2   = s2
     )
}

loglik_normal <- function(y, X, beta, s2) {
     y    <- as.numeric(y)
     X    <- as.matrix(X)
     beta <- as.numeric(beta)
     
     n     <- length(y)
     resid <- y - as.numeric(X %*% beta)
     
     -0.5 * n * log(2 * pi * s2) - 0.5 * sum(resid^2) / s2
}

gibbs_seleccion_variables <- function(
          y,
          X,
          B,
          n_burn,
          n_skip,
          g = nrow(X),
          nu0 = 1,
          s20 = NULL,
          a_theta = 1,
          b_theta = 1,
          verbose = TRUE
) {
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n          <- length(y)
     p          <- ncol(X)
     total_iter <- n_burn + B * n_skip
     
     if (is.null(s20)) {
          s20 <- summary(lm(y ~ -1 + X))$sigma^2
     }
     
     beta_samples   <- matrix(NA_real_, nrow = B, ncol = p)
     z_samples      <- matrix(NA_real_, nrow = B, ncol = p)
     sigma2_samples <- numeric(B)
     theta_samples  <- numeric(B)
     loglik_samples <- numeric(B)
     
     # Estado inicial
     z     <- rep(1, p)
     theta <- a_theta / (a_theta + b_theta)
     
     progress_step <- max(1, floor(total_iter / 10))
     
     for (b in seq_len(total_iter)) {
          # Actualizar los indicadores de inclusión
          for (j in sample(seq_len(p))) {
               z_plus    <- z
               z_plus[j] <- 1
               
               z_minus    <- z
               z_minus[j] <- 0
               
               lpy_plus <- lpy_X(
                    y   = y,
                    X   = X[, z_plus == 1, drop = FALSE],
                    g   = g,
                    nu0 = nu0,
                    s20 = s20
               )
               
               lpy_minus <- lpy_X(
                    y   = y,
                    X   = X[, z_minus == 1, drop = FALSE],
                    g   = g,
                    nu0 = nu0,
                    s20 = s20
               )
               
               log_odds <-
                    lpy_plus -
                    lpy_minus +
                    log(theta) -
                    log(1 - theta)
               
               prob_inclusion <- plogis(log_odds)
               
               z[j] <- rbinom(
                    n    = 1,
                    size = 1,
                    prob = prob_inclusion
               )
          }
          
          # Actualizar theta
          theta <- rbeta(
               n      = 1,
               shape1 = a_theta + sum(z),
               shape2 = b_theta + p - sum(z)
          )
          
          # Simular beta y sigma^2 dado el modelo activo
          X_z <- X[, z == 1, drop = FALSE]
          
          fit <- lm_gprior(
               y   = y,
               X   = X_z,
               g   = g,
               nu0 = nu0,
               s20 = s20,
               B   = 1
          )
          
          beta <- rep(0, p)
          
          if (sum(z) > 0) {
               beta[z == 1] <- fit$beta[1, ]
          }
          
          sigma2 <- fit$s2
          
          # Almacenar las muestras posteriores
          if (b > n_burn && (b - n_burn) %% n_skip == 0) {
               i <- (b - n_burn) / n_skip
               
               beta_samples[i, ]   <- beta
               z_samples[i, ]      <- z
               sigma2_samples[i]   <- sigma2
               theta_samples[i]    <- theta
               
               loglik_samples[i] <- loglik_normal(
                    y    = y,
                    X    = X,
                    beta = beta,
                    s2   = sigma2
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
     colnames(z_samples)    <- beta_names
     
     list(
          beta   = beta_samples,
          z      = z_samples,
          sigma2 = sigma2_samples,
          theta  = theta_samples,
          loglik = loglik_samples
     )
}

# Iteraciones muestreadores ----------------------------------------------------

# Número de iteraciones
n_sams <- 10000
n_burn <- 10000
n_skip <- 10

B <- n_sams

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

# Ajuste del modelo de selección de variables ---------------------------------

# Hiperparámetros
g   <- nrow(X_train)
nu0 <- 1
s20 <- sig2_ols

a_theta <- 1
b_theta <- 1

# Ajuste del modelo
set.seed(123)

muestras_seleccion <- gibbs_seleccion_variables(
     y       = y_train,
     X       = X_train,
     B       = B,
     n_burn  = n_burn,
     n_skip  = n_skip,
     g       = g,
     nu0     = nu0,
     s20     = s20,
     a_theta = a_theta,
     b_theta = b_theta,
     verbose = TRUE
)

save(
     muestras_seleccion,
     file = "simualcion_seleccion_variables_muestras_seleccion_variables.RData"
)

# Convergencia -----------------------------------------------------------------

# Cadena log-verosimilitud selección de variables
plot(
     muestras_seleccion$loglik,
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

beta_difusa    <- as.matrix(muestras_difusa$BETA)
beta_ridge     <- as.matrix(muestras_ridge$beta)
beta_lasso     <- as.matrix(muestras_lasso$beta)
beta_seleccion <- as.matrix(muestras_seleccion$beta)

muestras_beta <- list(
     Difusa    = beta_difusa,
     Ridge     = beta_ridge,
     Lasso     = beta_lasso,
     Seleccion = beta_seleccion
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
          Sesgo            = mean(errores),
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

# Predicción en el conjunto de prueba ------------------------------------------

resumen_modelo_bayesiano <- function(muestras_beta, nombre_modelo) {
     # Convertir las muestras posteriores a matriz
     muestras_beta <- as.matrix(muestras_beta)
     
     # Media posterior de los coeficientes
     beta_hat <- colMeans(muestras_beta)
     
     # Predicción sobre el conjunto de prueba
     y_test_hat <- as.numeric(X_test %*% beta_hat)
     
     # Medidas de desempeño predictivo
     errores <- y_test - y_test_hat
     
     amse  <- mean(errores^2)
     rmse  <- sqrt(amse)
     mae   <- mean(abs(errores))
     corr  <- cor(y_test, y_test_hat)
     r2    <- 1 - sum(errores^2) / sum((y_test - mean(y_test))^2)
     sesgo <- abs(mean(y_test_hat - y_test))
     
     # Devolver resultados
     list(
          beta_hat  = beta_hat,
          y_test_hat = y_test_hat,
          metricas = c(
               RMSE  = rmse,
               MAE   = mae,
               Corr  = corr,
               R2    = r2,
               Sesgo = sesgo
          )
     )
}

resultados_modelos <- list(
     Difusa = resumen_modelo_bayesiano(
          muestras_beta = muestras_difusa$BETA,
          nombre_modelo = "difusa"
     ),
     Ridge = resumen_modelo_bayesiano(
          muestras_beta = muestras_ridge$beta,
          nombre_modelo = "ridge"
     ),
     Lasso = resumen_modelo_bayesiano(
          muestras_beta = muestras_lasso$beta,
          nombre_modelo = "lasso"
     ),
     Seleccion = resumen_modelo_bayesiano(
          muestras_beta = muestras_seleccion$beta,
          nombre_modelo = "seleccion_variables"
     )
)

tabla_metricas <- do.call(
     rbind,
     lapply(
          names(resultados_modelos),
          function(modelo) {
               data.frame(
                    Modelo = modelo,
                    t(resultados_modelos[[modelo]]$metricas),
                    row.names = NULL
               )
          }
     )
)

tabla_metricas[, -1] <- round(
     tabla_metricas[, -1],
     3
)

print(tabla_metricas)

# Comparación de selección de coeficientes no nulos ----------------------------

beta_true <- as.numeric(beta_true)

epsilon  <- 0.10
gamma    <- 0.95
nivel_ic <- 0.95

muestras_modelos <- list(
     Difusa = list(
          beta = as.matrix(muestras_difusa$BETA),
          z    = NULL
     ),
     Ridge = list(
          beta = as.matrix(muestras_ridge$beta),
          z    = NULL
     ),
     Lasso = list(
          beta = as.matrix(muestras_lasso$beta),
          z    = NULL
     ),
     Seleccion = list(
          beta = as.matrix(muestras_seleccion$beta),
          z    = as.matrix(muestras_seleccion$z)
     )
)

calcular_auc <- function(score, verdadero) {
     verdadero <- as.logical(verdadero)
     
     idx       <- is.finite(score)
     score     <- score[idx]
     verdadero <- verdadero[idx]
     
     n_pos <- sum(verdadero)
     n_neg <- sum(!verdadero)
     
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
          nivel_ic,
          z_samples = NULL
) {
     alfa              <- 1 - nivel_ic
     beta_true         <- as.numeric(beta_true)
     verdadero_no_nulo <- beta_true != 0
     beta_names        <- colnames(beta_samples)
     
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
     
     if (is.null(z_samples)) {
          prob_inclusion      <- rep(NA_real_, length(beta_true))
          seleccion_inclusion <- rep(NA, length(beta_true))
     } else {
          prob_inclusion      <- colMeans(z_samples)
          seleccion_inclusion <- prob_inclusion > gamma
     }
     
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
          prob_inclusion             = prob_inclusion,
          seleccion_ic               = seleccion_ic,
          seleccion_relevancia       = seleccion_relevancia,
          seleccion_signo_relevancia = seleccion_signo_relevancia,
          seleccion_inclusion        = seleccion_inclusion,
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

evaluar_regla <- function(
          tabla,
          columna_seleccion,
          nombre_regla,
          columna_score = "prob_relevancia"
) {
     verdadero    <- tabla$verdadero_no_nulo
     seleccionado <- tabla[[columna_seleccion]]
     
     if (all(is.na(seleccionado))) {
          return(NULL)
     }
     
     seleccionado <- as.logical(seleccionado)
     
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
               score     = tabla[[columna_score]],
               verdadero = verdadero
          )
     )
}

resumen_coeficientes <- do.call(
     rbind,
     lapply(
          names(muestras_modelos),
          function(modelo) {
               resumen_seleccion_modelo(
                    beta_samples = muestras_modelos[[modelo]]$beta,
                    beta_true    = beta_true,
                    modelo       = modelo,
                    epsilon      = epsilon,
                    gamma        = gamma,
                    nivel_ic     = nivel_ic,
                    z_samples    = muestras_modelos[[modelo]]$z
               )
          }
     )
)

metricas_seleccion <- do.call(
     rbind,
     lapply(
          split(resumen_coeficientes, resumen_coeficientes$Modelo),
          function(tabla) {
               metricas <- list(
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
                    ),
                    evaluar_regla(
                         tabla             = tabla,
                         columna_seleccion = "seleccion_inclusion",
                         nombre_regla      = paste0("P(z = 1) > ", gamma),
                         columna_score     = "prob_inclusion"
                    )
               )
               
               do.call(rbind, metricas[!sapply(metricas, is.null)])
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
     "prob_relevancia",
     "prob_inclusion"
)] <- round(
     resumen_coeficientes[, c(
          "media",
          "mediana",
          "q025",
          "q975",
          "prob_pos",
          "prob_neg",
          "prob_signo",
          "prob_relevancia",
          "prob_inclusion"
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

# Muestras posteriores del modelo de selección de variables
beta_seleccion <- as.matrix(muestras_seleccion$beta)
beta_true      <- as.numeric(beta_true)

p     <- length(beta_true)
x_pos <- seq_len(p)

# Estadísticas posteriores
media_seleccion <- colMeans(beta_seleccion)

ic95_seleccion <- apply(
     beta_seleccion,
     2,
     quantile,
     probs = c(0.025, 0.975)
)

ylim_all <- range(
     c(
          ic95_seleccion,
          beta_true
     ),
     na.rm = TRUE
)

# Gráfico: modelo de selección Bayesiana de variables
pdf(
     file      = "simulacion_seleccion_variables_intervalos.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     NA,
     NA,
     xlab = "Índice del coeficiente",
     ylab = expression(beta[j]),
     xlim = c(0.5, p + 0.5),
     ylim = ylim_all,
     main = ""
)

abline(
     h   = 0,
     col = "gray70",
     lwd = 2
)

for (j in seq_len(p)) {
     segments(
          x0  = j,
          y0  = ic95_seleccion[1, j],
          x1  = j,
          y1  = ic95_seleccion[2, j],
          lwd = 1
     )
     
     points(
          x   = j,
          y   = media_seleccion[j],
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
     legend = "Selección",
     bty    = "n",
     cex    = 1.5
)

dev.off()

# Resumen ----------------------------------------------------------------------

print(resumen_beta_media)

print(tabla_metricas)

print(metricas_seleccion)

# Inferencia posterior para selección Bayesiana de variables -------------------

beta_true <- as.numeric(beta_true)

Z_seleccion    <- as.matrix(muestras_seleccion$z)
BETA_seleccion <- as.matrix(muestras_seleccion$beta)

p <- length(beta_true)

beta_names <- colnames(Z_seleccion)

colnames(Z_seleccion)    <- beta_names
colnames(BETA_seleccion) <- beta_names

# Probabilidades posteriores de inclusión
inclusion_probs <- colMeans(Z_seleccion)

# Verdadero modelo activo
z_true   <- as.integer(beta_true != 0)
p_z_true <- sum(z_true)

# Media posterior de beta
beta_hat <- colMeans(BETA_seleccion)

# Reglas de selección
gamma_mediana <- 0.50
gamma_alta    <- 0.95

z_hat_mediana <- as.integer(inclusion_probs > gamma_mediana)
z_hat_alta    <- as.integer(inclusion_probs > gamma_alta)

# Gráfico: probabilidades posteriores de inclusión
pdf(
     file      = "simulacion_seleccion_variables_prob_inclusion_beta_true.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.2, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

col_barras <- ifelse(
     z_true == 0,
     "red3",
     "darkgreen"
)

plot(
     inclusion_probs,
     type = "h",
     lwd  = 4,
     col  = col_barras,
     ylim = c(0, 1.08),
     xlab = "Índice del regresor",
     ylab = expression(
          paste(
               "Pr(",
               italic(z[j] == 1),
               " | ",
               bold(y),
               ", ",
               bold(X),
               ")"
          )
     )
)

abline(
     h   = 0.50,
     lty = 2,
     lwd = 2.5,
     col = "gray25"
)

abline(
     h   = 0.95,
     lty = 3,
     lwd = 2.5,
     col = "gray25"
)

dev.off()

# Inferencia sobre theta
theta_samples <- as.numeric(muestras_seleccion$theta)

ic95_theta <- quantile(
     theta_samples,
     probs = c(0.025, 0.975)
)

resumen_theta <- c(
     media             = mean(theta_samples),
     mediana           = median(theta_samples),
     q025              = ic95_theta[1],
     q975              = ic95_theta[2],
     theta_true_aprox  = p_z_true / p
)

print(round(resumen_theta, 3))

pdf(
     file      = "simulacion_seleccion_variables_theta.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.1, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     theta_samples,
     breaks = 30,
     freq   = FALSE,
     xlab   = expression(theta),
     ylab   = "Densidad",
     main   = "",
     col    = "gray85",
     border = "white"
)

abline(
     v   = mean(theta_samples),
     lwd = 2
)

abline(
     v   = ic95_theta,
     col = "blue3",
     lwd = 2,
     lty = 3
)

abline(
     v   = p_z_true / p,
     col = "red3",
     lwd = 2,
     lty = 2
)

legend(
     "topright",
     legend = c(
          "Media",
          "IC 95%",
          expression(p[z]^{true} / p)
     ),
     lty = c(1, 3, 2),
     col = c("black", "blue3", "red3"),
     lwd = 2,
     bty = "n"
)

box()

dev.off()

# Inferencia sobre el número de variables activas

p_z_samples <- rowSums(Z_seleccion)

resumen_p_z <- c(
     media    = mean(p_z_samples),
     mediana  = median(p_z_samples),
     q025     = quantile(p_z_samples, 0.025),
     q975     = quantile(p_z_samples, 0.975),
     p_z_true = p_z_true
)

print(round(resumen_p_z, 3))

pdf(
     file      = "simulacion_seleccion_variables_tamano_modelo.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.1, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

tab_p_z <- table(
     factor(
          p_z_samples,
          levels = 0:p
     )
)

prob_p_z <- as.numeric(tab_p_z) / sum(tab_p_z)

plot(
     x    = 0:p,
     y    = prob_p_z,
     type = "h",
     lwd  = 4,
     col  = "gray50",
     xlim = c(-0.5, p + 0.5),
     ylim = c(0, max(prob_p_z) * 1.15),
     xlab = expression(p[bold(z)]),
     ylab = "Probabilidad posterior"
)

abline(
     v   = p_z_true,
     col = "red3",
     lwd = 2.5,
     lty = 2
)

legend(
     "topright",
     legend = expression(p[bold(z)]^{true}),
     lty    = 2,
     col    = "red3",
     lwd    = 2.5,
     bty    = "n"
)

dev.off()

# Fin --------------------------------------------------------------------------