# Configuración ----------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

# Cargar datos
nitro <- read.table(
     file   = "nitro.txt",
     header = FALSE,
     col.names = c(
          "nitro",
          "size",
          "farm"
     )
)

# Organizar variables
nitro$farm <- factor(
     nitro$farm,
     levels = sort(unique(nitro$farm))
)

y    <- as.matrix(nitro$size)
x    <- as.matrix(nitro$nitro)
farm <- nitro$farm

# Dimensiones del estudio
n_total <- nrow(nitro)
m       <- length(levels(farm))
n       <- as.integer(n_total / m)

# Matriz de diseño
X <- cbind(1, x)
p <- ncol(X)

# Grupos
grupo <- nitro$farm

colnames(X) <- c(
     "Intercepto",
     "Nitrogeno"
)

# Análisis exploratorio --------------------------------------------------------

# Histograma y densidad del tamaño de planta
densidad_y <- density(
     x      = as.numeric(y),
     adjust = 1.5,
     na.rm  = TRUE
)

pdf(
     file      = "nitro_histograma_size.pdf",
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
     x      = as.numeric(y),
     breaks = "FD",
     freq   = FALSE,
     col    = "gray92",
     border = "white",
     xlab   = "Tamaño de planta",
     ylab   = "Densidad",
     main   = "",
     ylim   = c(0, max(densidad_y$y, na.rm = TRUE) * 1.10)
)

lines(
     densidad_y,
     lwd = 2,
     col = "gray30"
)

box()

dev.off()

# Diagrama de dispersión con recta OLS
reg_ols <- lm(
     formula = size ~ nitro,
     data    = nitro
)

pdf(
     file      = "nitro_dispersion_tamano_nitro.pdf",
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
     x    = nitro$nitro,
     y    = nitro$size,
     pch  = 20,
     cex  = 1.2,
     col  = adjustcolor("gray30", 0.7),
     xlab = "Concentración de nitrógeno",
     ylab = "Tamaño de planta",
     main = ""
)

abline(
     reg_ols,
     col = "gray30",
     lwd = 2
)

box()

dev.off()

# Diagrama de dispersión por finca
col_farm <- grDevices::hcl.colors(
     n       = m,
     palette = "Dark 3"
)

pdf(
     file      = "nitro_dispersion_fincas.pdf",
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
     x    = nitro$nitro,
     y    = nitro$size,
     type = "n",
     xlab = "Concentración de nitrógeno",
     ylab = "Tamaño de planta",
     main = ""
)

for (j in seq_len(m)) {
     ind_j <- nitro$farm == levels(nitro$farm)[j]
     
     points(
          x   = nitro$nitro[ind_j],
          y   = nitro$size[ind_j],
          col = adjustcolor(col_farm[j], 0.7),
          cex = 1,
          pch = 16
     )
}

box()

dev.off()

# Diagramas de caja del tamaño de planta por finca
nitro <- as.data.frame(nitro)

nitro$farm <- factor(
     nitro$farm,
     levels = sort(unique(nitro$farm))
)

pdf(
     file      = "nitro_boxplot_size_farm.pdf",
     width     = 8,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

boxplot(
     size ~ farm,
     data     = nitro,
     col      = "gray92",
     border   = "gray30",
     boxwex   = 0.5,
     outline  = TRUE,
     outpch   = 16,
     outcex   = 0.5,
     cex.axis = 0.7,
     las      = 2,
     xlab     = "Finca",
     ylab     = "Tamaño de planta",
     main     = ""
)

box()

dev.off()

# Mínimos cuadrados ordinarios -------------------------------------------------

fit_ols <- lm(y ~ x)

summary(fit_ols)

coef_ols   <- coef(reg_ols)
sigma2_ols <- sigma(reg_ols)^2

round(coef_ols, 3)
round(sigma2_ols, 3)
round(sqrt(sigma2_ols), 3)

# Funciones auxiliares ---------------------------------------------------------

rigamma <- function(n, shape, scale) {
     1 / rgamma(
          n     = n,
          shape = shape,
          rate  = scale
     )
}

rdirichlet1 <- function(alpha) {
     x <- rgamma(
          n     = length(alpha),
          shape = alpha,
          rate  = 1
     )
     
     x / sum(x)
}

rmvnorm1 <- function(mu, Sigma) {
     mu <- as.numeric(mu)
     R  <- chol(Sigma)
     
     as.numeric(mu + t(R) %*% rnorm(length(mu)))
}

rinvwishart1 <- function(df, S) {
     W <- stats::rWishart(
          n     = 1,
          df    = df,
          Sigma = solve(S)
     )[, , 1]
     
     solve(W)
}

# Distribuciones condicionales completas --------------------------------------

rcond_gamma_chlrm <- function(omega, beta_k, sigma2_k, Xj, Yj) {
     m <- length(Yj)
     K <- length(omega)
     
     gamma <- integer(m)
     
     for (j in seq_len(m)) {
          log_q <- numeric(K)
          
          for (k in seq_len(K)) {
               log_q[k] <- log(omega[k]) +
                    sum(
                         dnorm(
                              x    = Yj[[j]],
                              mean = as.numeric(Xj[[j]] %*% beta_k[k, ]),
                              sd   = sqrt(sigma2_k[k]),
                              log  = TRUE
                         )
                    )
          }
          
          q <- exp(log_q - max(log_q))
          q <- q / sum(q)
          
          gamma[j] <- sample(
               x    = seq_len(K),
               size = 1,
               prob = q
          )
     }
     
     gamma
}

rcond_omega_chlrm <- function(gamma, alpha0, K) {
     m_k <- tabulate(
          gamma,
          nbins = K
     )
     
     list(
          omega  = rdirichlet1(alpha0 + m_k),
          m_k    = m_k,
          activo = m_k > 0,
          K_star = sum(m_k > 0)
     )
}

rcond_beta_k_chlrm <- function(
          k,
          activo,
          gamma,
          indices_grupo,
          X,
          y,
          beta,
          Sigma,
          Sigma_inv,
          Sigma_inv_beta,
          sigma2_k
) {
     if (activo[k]) {
          ind_k <- unlist(
               indices_grupo[gamma == k],
               use.names = FALSE
          )
          
          X_k <- X[ind_k, , drop = FALSE]
          y_k <- y[ind_k]
          
          V_beta_k <- solve(
               Sigma_inv + crossprod(X_k) / sigma2_k[k]
          )
          
          m_beta_k <- as.numeric(
               V_beta_k %*%
                    (
                         Sigma_inv_beta +
                              crossprod(X_k, y_k) / sigma2_k[k]
                    )
          )
          
          beta_k <- rmvnorm1(
               mu    = m_beta_k,
               Sigma = V_beta_k
          )
     } else {
          beta_k <- rmvnorm1(
               mu    = beta,
               Sigma = Sigma
          )
     }
     
     beta_k
}

rcond_beta_chlrm <- function(
          beta_k,
          activo,
          Lambda0_inv,
          Lambda0_inv_mu0,
          Sigma_inv,
          K_star
) {
     sum_beta_k <- colSums(
          beta_k[activo, , drop = FALSE]
     )
     
     V_beta <- solve(
          Lambda0_inv + K_star * Sigma_inv
     )
     
     m_beta <- as.numeric(
          V_beta %*%
               (
                    Lambda0_inv_mu0 +
                         Sigma_inv %*% sum_beta_k
               )
     )
     
     rmvnorm1(
          mu    = m_beta,
          Sigma = V_beta
     )
}

rcond_Sigma_chlrm <- function(
          beta_k,
          beta,
          activo,
          S0,
          n0,
          K_star
) {
     S_Sigma <- S0
     
     for (k in which(activo)) {
          dif_k <- matrix(
               beta_k[k, ] - beta,
               ncol = 1
          )
          
          S_Sigma <- S_Sigma + dif_k %*% t(dif_k)
     }
     
     rinvwishart1(
          df = n0 + K_star,
          S  = S_Sigma
     )
}

rcond_sigma2_k_chlrm <- function(
          k,
          activo,
          gamma,
          indices_grupo,
          X,
          y,
          beta_k,
          xi2,
          nu0
) {
     if (activo[k]) {
          ind_k <- unlist(
               indices_grupo[gamma == k],
               use.names = FALSE
          )
          
          X_k <- X[ind_k, , drop = FALSE]
          y_k <- y[ind_k]
          
          resid_k <- as.numeric(
               y_k - X_k %*% beta_k[k, ]
          )
          
          a_sigma_k <- 0.5 * (nu0 + length(y_k))
          b_sigma_k <- 0.5 * (nu0 * xi2 + sum(resid_k^2))
     } else {
          a_sigma_k <- 0.5 * nu0
          b_sigma_k <- 0.5 * nu0 * xi2
     }
     
     rigamma(
          n     = 1,
          shape = a_sigma_k,
          scale = b_sigma_k
     )
}

rcond_xi2_chlrm <- function(
          sigma2_k,
          activo,
          nu0,
          a0,
          b0,
          K_star
) {
     rgamma(
          n     = 1,
          shape = a0 + 0.5 * K_star * nu0,
          rate  = b0 + 0.5 * nu0 * sum(1 / sigma2_k[activo])
     )
}

loglik_chlrm <- function(gamma, beta_k, sigma2_k, Xj, Yj) {
     m  <- length(Yj)
     ll <- 0
     
     for (j in seq_len(m)) {
          k_j <- gamma[j]
          
          ll <- ll +
               sum(
                    dnorm(
                         x    = Yj[[j]],
                         mean = as.numeric(Xj[[j]] %*% beta_k[k_j, ]),
                         sd   = sqrt(sigma2_k[k_j]),
                         log  = TRUE
                    )
               )
     }
     
     ll
}

# Muestreador de Gibbs para el modelo CHLRM -----------------------------------

gibbs_chlrm <- function(
          y,
          X,
          grupo,
          K,
          g,
          alpha0,
          mu0,
          Lambda0,
          n0,
          S0,
          nu0,
          a0,
          b0,
          n_sams,
          n_burn,
          n_skip,
          seed = 123,
          verbose = TRUE
) {
     # Ajustes
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n_total <- length(y)
     p       <- ncol(X)
     
     grupo <- factor(
          grupo,
          levels = sort(unique(grupo))
     )
     
     grupos <- levels(grupo)
     m      <- length(grupos)
     
     indices_grupo <- split(
          x = seq_len(n_total),
          f = grupo
     )
     
     n_j <- as.numeric(table(grupo))
     
     beta_names <- colnames(X)
     
     if (is.null(beta_names) || any(beta_names == "")) {
          beta_names <- paste0(
               "beta_",
               seq_len(p)
          )
     }
     
     colnames(X) <- beta_names
     
     Xj <- lapply(
          indices_grupo,
          function(ind) {
               X[ind, , drop = FALSE]
          }
     )
     
     Yj <- lapply(
          indices_grupo,
          function(ind) {
               y[ind]
          }
     )
     
     # Cantidades OLS para inicialización y referencia
     beta_ols <- as.numeric(qr.solve(X, y))
     
     resid_ols <- as.numeric(y - X %*% beta_ols)
     
     sigma2_ols <- sum(resid_ols^2) / (n_total - p)
     
     # Hiperparámetros
     mu0     <- as.numeric(mu0)
     Lambda0 <- as.matrix(Lambda0)
     S0      <- as.matrix(S0)
     alpha0  <- as.numeric(alpha0)
     
     Lambda0_inv <- solve(Lambda0)
     
     Lambda0_inv_mu0 <- as.numeric(Lambda0_inv %*% mu0)
     
     # Número total de iteraciones
     B_total <- n_burn + n_sams * n_skip
     ncat    <- max(1, floor(0.1 * B_total))
     
     # Valores iniciales
     beta  <- mu0
     Sigma <- S0
     
     sigma2_k <- rep(
          max(sigma2_ols, .Machine$double.eps),
          K
     )
     
     xi2 <- max(sigma2_ols, .Machine$double.eps)
     
     beta_k <- matrix(
          rep(beta, each = K),
          nrow = K,
          ncol = p
     )
     
     gamma <- sample(
          x       = seq_len(K),
          size    = m,
          replace = TRUE
     )
     
     omega <- rdirichlet1(alpha0)
     
     Sigma_inv <- solve(Sigma)
     
     Sigma_inv_beta <- as.numeric(Sigma_inv %*% beta)
     
     # Almacenamiento
     BETA_K <- array(
          NA_real_,
          dim = c(n_sams, K, p),
          dimnames = list(
               NULL,
               paste0("cluster_", seq_len(K)),
               beta_names
          )
     )
     
     SIGMA2_K <- matrix(
          NA_real_,
          nrow = n_sams,
          ncol = K
     )
     
     BETA <- matrix(
          NA_real_,
          nrow = n_sams,
          ncol = p
     )
     
     SIGMA <- array(
          NA_real_,
          dim = c(n_sams, p, p),
          dimnames = list(
               NULL,
               beta_names,
               beta_names
          )
     )
     
     XI2 <- numeric(n_sams)
     
     OMEGA <- matrix(
          NA_real_,
          nrow = n_sams,
          ncol = K
     )
     
     GAMMA <- matrix(
          NA_integer_,
          nrow = n_sams,
          ncol = m
     )
     
     KSTAR <- integer(n_sams)
     LL    <- numeric(n_sams)
     
     colnames(SIGMA2_K) <- paste0("sigma2_", seq_len(K))
     colnames(BETA)     <- beta_names
     colnames(OMEGA)    <- paste0("omega_", seq_len(K))
     colnames(GAMMA)    <- paste0("gamma_", grupos)
     
     # Cadena
     for (b in seq_len(B_total)) {
          # Actualizar gamma_j
          gamma <- rcond_gamma_chlrm(
               omega    = omega,
               beta_k   = beta_k,
               sigma2_k = sigma2_k,
               Xj       = Xj,
               Yj       = Yj
          )
          
          # Actualizar omega
          omega_out <- rcond_omega_chlrm(
               gamma  = gamma,
               alpha0 = alpha0,
               K      = K
          )
          
          omega  <- omega_out$omega
          m_k    <- omega_out$m_k
          activo <- omega_out$activo
          K_star <- omega_out$K_star
          
          # Actualizar beta_k
          for (k in seq_len(K)) {
               beta_k[k, ] <- rcond_beta_k_chlrm(
                    k              = k,
                    activo         = activo,
                    gamma          = gamma,
                    indices_grupo  = indices_grupo,
                    X              = X,
                    y              = y,
                    beta           = beta,
                    Sigma          = Sigma,
                    Sigma_inv      = Sigma_inv,
                    Sigma_inv_beta = Sigma_inv_beta,
                    sigma2_k       = sigma2_k
               )
          }
          
          # Actualizar beta
          beta <- rcond_beta_chlrm(
               beta_k          = beta_k,
               activo          = activo,
               Lambda0_inv     = Lambda0_inv,
               Lambda0_inv_mu0 = Lambda0_inv_mu0,
               Sigma_inv       = Sigma_inv,
               K_star          = K_star
          )
          
          # Actualizar Sigma
          Sigma <- rcond_Sigma_chlrm(
               beta_k = beta_k,
               beta   = beta,
               activo = activo,
               S0     = S0,
               n0     = n0,
               K_star = K_star
          )
          
          Sigma_inv <- solve(Sigma)
          
          Sigma_inv_beta <- as.numeric(Sigma_inv %*% beta)
          
          # Actualizar sigma_k^2
          for (k in seq_len(K)) {
               sigma2_k[k] <- rcond_sigma2_k_chlrm(
                    k             = k,
                    activo        = activo,
                    gamma         = gamma,
                    indices_grupo = indices_grupo,
                    X             = X,
                    y             = y,
                    beta_k        = beta_k,
                    xi2           = xi2,
                    nu0           = nu0
               )
          }
          
          # Actualizar xi^2
          xi2 <- rcond_xi2_chlrm(
               sigma2_k = sigma2_k,
               activo   = activo,
               nu0      = nu0,
               a0       = a0,
               b0       = b0,
               K_star   = K_star
          )
          
          # Almacenar muestra posterior
          if (b > n_burn && (b - n_burn) %% n_skip == 0) {
               l <- (b - n_burn) / n_skip
               
               BETA_K[l, , ] <- beta_k
               SIGMA2_K[l, ] <- sigma2_k
               BETA[l, ]     <- beta
               SIGMA[l, , ]  <- Sigma
               XI2[l]        <- xi2
               OMEGA[l, ]    <- omega
               GAMMA[l, ]    <- gamma
               KSTAR[l]      <- K_star
               
               LL[l] <- loglik_chlrm(
                    gamma    = gamma,
                    beta_k   = beta_k,
                    sigma2_k = sigma2_k,
                    Xj       = Xj,
                    Yj       = Yj
               )
          }
          
          # Progreso
          if (verbose && b %% ncat == 0) {
               cat(
                    sprintf(
                         "%.1f%% completado\n",
                         100 * b / B_total
                    )
               )
          }
     }
     
     # Salida
     list(
          BETA_K   = BETA_K,
          SIGMA2_K = SIGMA2_K,
          BETA     = BETA,
          SIGMA    = SIGMA,
          XI2      = XI2,
          OMEGA    = OMEGA,
          GAMMA    = GAMMA,
          KSTAR    = KSTAR,
          LL       = LL,
          hyper    = list(
               g          = g,
               alpha0     = alpha0,
               mu0        = mu0,
               Lambda0    = Lambda0,
               n0         = n0,
               S0         = S0,
               nu0        = nu0,
               a0         = a0,
               b0         = b0,
               sigma2_ols = sigma2_ols,
               beta_ols   = beta_ols
          ),
          info = list(
               n_total    = n_total,
               m          = m,
               n_j        = n_j,
               p          = p,
               K          = K,
               grupos     = grupos,
               beta_names = beta_names
          )
     )
}

# Ajuste del modelo ------------------------------------------------------------

XtX        <- crossprod(X)
XtX_inv    <- solve(XtX)
beta_ols   <- as.numeric(qr.solve(X, y))
resid_ols  <- as.numeric(y - X %*% beta_ols)
sigma2_ols <- sum(resid_ols^2) / (n_total - p)

# Hiperparámetros
K       <- m
g       <- n_total
mu0     <- beta_ols
Lambda0 <- g * sigma2_ols * XtX_inv
n0      <- p + 2
S0      <- Lambda0
nu0     <- 1
a0      <- 1
b0      <- 1 / sigma2_ols
alpha0  <- rep(1 / K, K)

muestras_chlrm <- gibbs_chlrm(
     y,
     X,
     grupo,
     K,
     g,
     alpha0,
     mu0,
     Lambda0,
     n0,
     S0,
     nu0,
     a0,
     b0,
     n_sams = 10000,
     n_burn = 10000,
     n_skip = 50
)

save(
     muestras_chlrm,
     file = "nitro_muestras_distribucion_posterior.RData"
)

load(
     file = "nitro_muestras_distribucion_posterior.RData"
)

# Número de iteraciones
B <- length(muestras_chlrm$LL)

# Convergencia -----------------------------------------------------------------

# Cadena log-verosimilitud Ridge
plot(
     muestras_chlrm$LL,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = ""
)

dev.off()

# Distribución posterior del número de conglomerados ocupados ------------------

GAMMA <- muestras_chlrm$GAMMA

B <- nrow(GAMMA)

K_post <- apply(
     GAMMA,
     1,
     function(x) length(unique(x))
)

tabla_K <- table(K_post) / B

K_vals <- as.integer(
     names(tabla_K)
)

K_map <- K_vals[
     which.max(tabla_K)
]

pdf(
     file      = "nitro_numero_clusters_posterior.pdf",
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
     x    = K_vals,
     y    = as.numeric(tabla_K),
     type = "h",
     lwd  = 5,
     lend = 1,
     col  = "gray30",
     xlab = "Número de conglomerados",
     ylab = "Probabilidad posterior",
     main = "",
     xaxt = "n",
     ylim = c(0, max(tabla_K) * 1.10)
)

axis(
     side   = 1,
     at     = seq(from = min(K_vals), to = max(K_vals), by = 1),
     labels = seq(from = min(K_vals), to = max(K_vals), by = 1)
)

box()

dev.off()

round(tabla_K, 3)

K_map

# Matriz de incidencia posterior ----------------------------------------------

matriz_incidencia <- function(GAMMA, verbose = TRUE) {
     GAMMA <- as.matrix(GAMMA)
     
     B <- nrow(GAMMA)
     m <- ncol(GAMMA)
     
     A <- matrix(
          0,
          nrow = m,
          ncol = m
     )
     
     ncat <- max(1, floor(0.1 * B))
     
     for (b in seq_len(B)) {
          A <- A + outer(
               GAMMA[b, ],
               GAMMA[b, ],
               FUN = "=="
          )
          
          if (verbose && b %% ncat == 0) {
               cat(
                    sprintf(
                         "%.1f%% completado\n",
                         100 * b / B
                    )
               )
          }
     }
     
     A <- A / B
     
     diag(A) <- 1
     
     colnames(A) <- colnames(GAMMA)
     rownames(A) <- colnames(GAMMA)
     
     A
}

A <- matriz_incidencia(
     GAMMA   = muestras_chlrm$GAMMA,
     verbose = TRUE
)

# Partición estimada -----------------------------------------------------------

estimar_particion_hclust <- function(
          A,
          KSTAR = NULL,
          K_hat = NULL,
          metodo = "average"
) {
     A <- as.matrix(A)
     
     m <- nrow(A)
     
     if (is.null(K_hat)) {
          if (is.null(KSTAR)) {
               stop("Debe especificar K_hat o proporcionar KSTAR.")
          }
          
          K_hat <- as.integer(
               names(
                    which.max(
                         table(KSTAR)
                    )
               )
          )
     }
     
     distancia <- as.dist(1 - A)
     
     hc <- hclust(
          d      = distancia,
          method = metodo
     )
     
     cluster <- cutree(
          tree = hc,
          k    = K_hat
     )
     
     tamanos <- sort(
          table(cluster),
          decreasing = TRUE
     )
     
     niveles_cluster <- as.integer(
          names(tamanos)
     )
     
     cluster_relab <- cluster
     
     for (k in seq_along(niveles_cluster)) {
          cluster_relab[cluster == niveles_cluster[k]] <- k
     }
     
     orden <- order(cluster_relab)
     
     cluster_ordenado <- cluster_relab[orden]
     
     cortes <- which(
          cluster_ordenado[-length(cluster_ordenado)] !=
               cluster_ordenado[-1]
     )
     
     list(
          A                = A,
          A_ordenada       = A[orden, orden],
          cluster          = cluster_relab,
          cluster_ordenado = cluster_ordenado,
          orden            = orden,
          cortes           = cortes,
          K_hat            = K_hat,
          hclust           = hc
     )
}

# Estimar partición usando clustering jerárquico -------------------------------

particion_hc <- estimar_particion_hclust(
     A      = A,
     KSTAR  = muestras_chlrm$KSTAR,
     metodo = "average"
)

A_ordenada  <- particion_hc$A_ordenada
cluster_hat <- particion_hc$cluster
indices     <- particion_hc$orden
cortes      <- particion_hc$cortes
K_hat       <- particion_hc$K_hat

# Graficar matriz de incidencia posterior -------------------------------------

# Escala de colores
col_incidencia <- colorRampPalette(
     colors = c(
          "white",
          "yellow",
          "red"
     )
)(100)

pdf(
     file      = "nitro_matriz_incidencia_hclust.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

# Matriz de incidencia
par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

image(
     x    = seq_len(nrow(A_ordenada)),
     y    = seq_len(ncol(A_ordenada)),
     z    = A_ordenada,
     col  = col_incidencia,
     zlim = c(0, 1),
     axes = FALSE,
     xlab = "Finca",
     ylab = "Finca",
     main = ""
)

axis(
     side     = 1,
     at       = seq_len(length(indices)),
     labels   = indices,
     las      = 2,
     cex.axis = 0.75
)

axis(
     side     = 2,
     at       = seq_len(length(indices)),
     labels   = indices,
     las      = 1,
     cex.axis = 0.75
)

box()

if (length(cortes) > 0) {
     abline(
          h   = cortes + 0.5,
          v   = cortes + 0.5,
          col = "gray20",
          lwd = 1
     )
}

dev.off()

# Resumen de la partición estimada --------------------------------------------

split(
     x = indices,
     f = cluster_hat[indices]
)

table(cluster_hat)

# Matriz binaria de la partición estimada --------------------------------------

P_hat <- outer(
     cluster_hat,
     cluster_hat,
     FUN = "=="
)

P_hat <- 1 * P_hat

rownames(P_hat) <- names(cluster_hat)
colnames(P_hat) <- names(cluster_hat)

P_hat_ordenada <- P_hat[
     indices,
     indices
]

# Gráfico de la partición estimada --------------------------------------------

col_particion <- c(
     "white",
     "red"
)

pdf(
     file      = "nitro_particion_estimada_hclust.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

image(
     x      = seq_len(nrow(P_hat_ordenada)),
     y      = seq_len(ncol(P_hat_ordenada)),
     z      = P_hat_ordenada,
     col    = col_particion,
     breaks = c(-0.5, 0.5, 1.5),
     axes   = FALSE,
     xlab   = "Finca",
     ylab   = "Finca",
     main   = ""
)

axis(
     side     = 1,
     at       = seq_len(length(indices)),
     labels   = indices,
     las      = 2,
     cex.axis = 0.75
)

axis(
     side     = 2,
     at       = seq_len(length(indices)),
     labels   = indices,
     las      = 1,
     cex.axis = 0.75
)

box()

if (length(cortes) > 0) {
     abline(
          h   = cortes + 0.5,
          v   = cortes + 0.5,
          col = "gray20",
          lwd = 1
     )
}

dev.off()

# Estimación de coeficientes por conglomerado estimado -------------------------

GAMMA <- as.matrix(
     muestras_chlrm$GAMMA
)

BETA_K <- muestras_chlrm$BETA_K

SIGMA2_K <- as.matrix(
     muestras_chlrm$SIGMA2_K
)

B <- nrow(GAMMA)
p <- dim(BETA_K)[3]

Xi_hat <- as.integer(
     cluster_hat
)

K_hat <- max(
     Xi_hat
)

cluster_ref <- vector(
     mode   = "list",
     length = K_hat
)

for (k in seq_len(K_hat)) {
     cluster_ref[[k]] <- sort(
          which(Xi_hat == k)
     )
}

K_post <- apply(
     GAMMA,
     1,
     function(x) length(unique(x))
)

BETAk_post  <- NULL
SIGMAk_post <- NULL

for (b in seq_len(B)) {
     gamma_b <- as.integer(
          GAMMA[b, ]
     )
     
     gamma_unique <- sort(
          unique(gamma_b)
     )
     
     if (K_post[b] == K_hat) {
          cluster_b <- vector(
               mode   = "list",
               length = K_hat
          )
          
          for (l in seq_len(K_hat)) {
               cluster_b[[l]] <- sort(
                    which(gamma_b == gamma_unique[l])
               )
          }
          
          labels_b <- rep(
               NA_integer_,
               K_hat
          )
          
          for (k in seq_len(K_hat)) {
               for (l in seq_len(K_hat)) {
                    if (
                         length(cluster_ref[[k]]) == length(cluster_b[[l]]) &&
                         all(cluster_ref[[k]] == cluster_b[[l]])
                    ) {
                         labels_b[k] <- gamma_unique[l]
                    }
               }
          }
          
          if (!any(is.na(labels_b))) {
               beta_b <- matrix(
                    BETA_K[b, labels_b, , drop = FALSE],
                    nrow = K_hat,
                    ncol = p
               )
               
               sigma_b <- sqrt(
                    SIGMA2_K[
                         b,
                         labels_b
                    ]
               )
               
               BETAk_post <- rbind(
                    BETAk_post,
                    as.vector(beta_b)
               )
               
               SIGMAk_post <- rbind(
                    SIGMAk_post,
                    as.vector(sigma_b)
               )
          }
     }
}

if (is.null(BETAk_post)) {
     stop("No se encontraron iteraciones compatibles con la partición estimada.")
}

betak <- matrix(
     colMeans(BETAk_post),
     nrow = K_hat,
     ncol = p
)

betak_sd <- matrix(
     apply(BETAk_post, 2, sd),
     nrow = K_hat,
     ncol = p
)

betak_q025 <- matrix(
     apply(BETAk_post, 2, quantile, probs = 0.025),
     nrow = K_hat,
     ncol = p
)

betak_q975 <- matrix(
     apply(BETAk_post, 2, quantile, probs = 0.975),
     nrow = K_hat,
     ncol = p
)

sigmak <- colMeans(
     SIGMAk_post
)

sigmak_sd <- apply(
     SIGMAk_post,
     2,
     sd
)

sigmak_q025 <- apply(
     SIGMAk_post,
     2,
     quantile,
     probs = 0.025
)

sigmak_q975 <- apply(
     SIGMAk_post,
     2,
     quantile,
     probs = 0.975
)

colnames(betak)      <- dimnames(BETA_K)[[3]]
colnames(betak_sd)   <- dimnames(BETA_K)[[3]]
colnames(betak_q025) <- dimnames(BETA_K)[[3]]
colnames(betak_q975) <- dimnames(BETA_K)[[3]]

rownames(betak)      <- paste0("cluster_", seq_len(K_hat))
rownames(betak_sd)   <- paste0("cluster_", seq_len(K_hat))
rownames(betak_q025) <- paste0("cluster_", seq_len(K_hat))
rownames(betak_q975) <- paste0("cluster_", seq_len(K_hat))

names(sigmak)      <- paste0("cluster_", seq_len(K_hat))
names(sigmak_sd)   <- paste0("cluster_", seq_len(K_hat))
names(sigmak_q025) <- paste0("cluster_", seq_len(K_hat))
names(sigmak_q975) <- paste0("cluster_", seq_len(K_hat))

resumen_parametros_conglomerados <- data.frame(
     conglomerado = paste0("cluster_", seq_len(K_hat)),
     
     beta0_media = betak[, 1],
     beta0_sd    = betak_sd[, 1],
     beta0_q025  = betak_q025[, 1],
     beta0_q975  = betak_q975[, 1],
     
     beta1_media = betak[, 2],
     beta1_sd    = betak_sd[, 2],
     beta1_q025  = betak_q025[, 2],
     beta1_q975  = betak_q975[, 2],
     
     sigma_media = sigmak,
     sigma_sd    = sigmak_sd,
     sigma_q025  = sigmak_q025,
     sigma_q975  = sigmak_q975,
     
     row.names = NULL
)

resumen_parametros_conglomerados[, -1] <- round(
     resumen_parametros_conglomerados[, -1],
     3
)

resumen_parametros_conglomerados

# Gráficos por conglomerado estimado ------------------------------------------

x <- as.numeric(
     X[, 2]
)

y <- as.numeric(
     y
)

farm <- factor(
     farm,
     levels = sort(unique(farm))
)

m <- length(
     levels(farm)
)

col_cluster <- grDevices::hcl.colors(
     n       = K_hat,
     palette = "Dark 3"
)

for (k in seq_len(K_hat)) {
     pdf(
          file = paste0(
               "nitro_dispersion_cluster_",
               k,
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
     
     plot(
          x    = x,
          y    = y,
          type = "n",
          xlab = "Concentración de nitrógeno",
          ylab = "Tamaño de planta",
          main = ""
     )
     
     points(
          x   = x,
          y   = y,
          pch = 20,
          col = "gray85",
          cex = 1
     )
     
     abline(
          a   = betak[k, 1],
          b   = betak[k, 2],
          lwd = 2,
          col = "gray30"
     )
     
     for (j in seq_len(m)) {
          if (Xi_hat[j] == k) {
               ind_j <- farm == levels(farm)[j]
               
               text(
                    x      = x[ind_j],
                    y      = y[ind_j],
                    labels = j,
                    col    = col_cluster[k],
                    cex    = 1.1
               )
          }
     }
     
     box()
     
     dev.off()
}

# Chequeo predictivo posterior -------------------------------------------------

estadisticos_ppp <- function(y, X, grupo) {
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     grupo <- factor(
          grupo,
          levels = sort(unique(grupo))
     )
     
     media      <- mean(y)
     desviacion <- sd(y)
     iqr        <- IQR(y)
     
     asimetria <- mean((y - mean(y))^3) / sd(y)^3
     curtosis  <- mean((y - mean(y))^4) / sd(y)^4
     
     media_grupo <- tapply(
          y,
          grupo,
          mean
     )
     
     var_entre_grupos <- var(
          as.numeric(media_grupo)
     )
     
     c(
          media            = media,
          desviacion       = desviacion,
          iqr              = iqr,
          asimetria        = asimetria,
          curtosis         = curtosis,
          var_entre_grupos = var_entre_grupos
     )
}

# Simulación predictiva posterior ---------------------------------------------

GAMMA <- as.matrix(
     muestras_chlrm$GAMMA
)

BETA_K <- muestras_chlrm$BETA_K

SIGMA2_K <- as.matrix(
     muestras_chlrm$SIGMA2_K
)

B <- nrow(GAMMA)

y <- as.numeric(
     nitro$size
)

X <- cbind(
     1,
     nitro$nitro
)

grupo <- factor(
     nitro$farm,
     levels = sort(unique(nitro$farm))
)

m <- length(
     levels(grupo)
)

indices_grupo <- split(
     x = seq_along(y),
     f = grupo
)

T_obs <- estadisticos_ppp(
     y     = y,
     X     = X,
     grupo = grupo
)

T_rep <- matrix(
     NA_real_,
     nrow = B,
     ncol = length(T_obs)
)

colnames(T_rep) <- names(T_obs)

set.seed(123)

for (b in seq_len(B)) {
     y_rep <- numeric(
          length(y)
     )
     
     for (j in seq_len(m)) {
          ind_j <- indices_grupo[[j]]
          k_j   <- GAMMA[b, j]
          
          media_j <- as.numeric(
               X[ind_j, , drop = FALSE] %*% BETA_K[b, k_j, ]
          )
          
          y_rep[ind_j] <- rnorm(
               n    = length(ind_j),
               mean = media_j,
               sd   = sqrt(SIGMA2_K[b, k_j])
          )
     }
     
     T_rep[b, ] <- estadisticos_ppp(
          y     = y_rep,
          X     = X,
          grupo = grupo
     )
}

# Valores p predictivos posteriores -------------------------------------------

ppp <- colMeans(
     sweep(
          T_rep,
          2,
          T_obs,
          FUN = "<"
     )
)

resumen_ppp <- data.frame(
     estadistico      = names(T_obs),
     observado        = as.numeric(T_obs),
     media_predictiva = colMeans(T_rep),
     q025_predictivo  = apply(T_rep, 2, quantile, probs = 0.025),
     q500_predictivo  = apply(T_rep, 2, quantile, probs = 0.500),
     q975_predictivo  = apply(T_rep, 2, quantile, probs = 0.975),
     ppp              = as.numeric(ppp),
     row.names        = NULL
)

resumen_ppp

# Histogramas del chequeo predictivo posterior --------------------------------

estadisticos <- colnames(T_rep)

labels_x <- c(
     media            = "Media",
     desviacion       = "Desviación estándar",
     iqr              = "Rango intercuartílico",
     asimetria        = "Asimetría",
     curtosis         = "Curtosis",
     var_entre_grupos = "Varianza entre fincas"
)

nombres_archivo <- c(
     media            = "media",
     desviacion       = "desviacion_estandar",
     iqr              = "rango_intercuartilico",
     asimetria        = "asimetria",
     curtosis         = "curtosis",
     var_entre_grupos = "varianza_entre_fincas"
)

for (s in estadisticos) {
     xlab_s <- labels_x[s]
     
     if (is.na(xlab_s)) {
          xlab_s <- s
     }
     
     archivo_s <- nombres_archivo[s]
     
     if (is.na(archivo_s)) {
          archivo_s <- s
     }
     
     ppp_s <- ppp[s]
     
     pdf(
          file = paste0(
               "nitro_histograma_ppp_",
               archivo_s,
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
     
     hist(
          x      = T_rep[, s],
          breaks = 25,
          freq   = FALSE,
          col    = "gray90",
          border = "white",
          xlab   = xlab_s,
          ylab   = "Densidad",
          main   = "",
          xlim   = range(
               c(
                    T_rep[, s],
                    T_obs[s]
               ),
               na.rm = TRUE
          )
     )
     
     abline(
          v   = T_obs[s],
          lwd = 2,
          lty = 2,
          col = "gray20"
     )
     
     legend(
          "topright",
          legend = paste0(
               "ppp = ",
               format(
                    round(ppp_s, 3),
                    nsmall = 3
               )
          ),
          bty = "n",
          cex = 1.1
     )
     
     box()
     
     dev.off()
}

# Chequeo residuales -----------------------------------------------------------

GAMMA <- as.matrix(
     muestras_chlrm$GAMMA
)

BETA_K <- muestras_chlrm$BETA_K

y_obs <- as.numeric(
     nitro$size
)

X_res <- cbind(
     1,
     nitro$nitro
)

colnames(X_res) <- c(
     "Intercepto",
     "Nitrogeno"
)

farm <- factor(
     nitro$farm,
     levels = sort(unique(nitro$farm))
)

indices_grupo <- split(
     x = seq_along(y_obs),
     f = farm
)

B       <- nrow(GAMMA)
m       <- length(indices_grupo)
n_total <- length(y_obs)

mu_sum <- numeric(n_total)

for (b in seq_len(B)) {
     mu_b <- numeric(n_total)
     
     for (j in seq_len(m)) {
          ind_j <- indices_grupo[[j]]
          k_j   <- GAMMA[b, j]
          
          mu_b[ind_j] <- as.numeric(
               X_res[ind_j, , drop = FALSE] %*%
                    BETA_K[b, k_j, ]
          )
     }
     
     mu_sum <- mu_sum + mu_b
}

mu_hat <- mu_sum / B

resid     <- y_obs - mu_hat
resid_std <- resid / sd(resid)

resumen_residuales <- data.frame(
     finca                  = farm,
     y                      = y_obs,
     ajustado               = mu_hat,
     residual               = resid,
     residual_estandarizado = resid_std
)

summary(
     resumen_residuales$residual_estandarizado
)

# Gráficos diagnósticos de residuales -----------------------------------------

pdf(
     file      = "nitro_residuales_ajustados.pdf",
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
     x    = mu_hat,
     y    = resid,
     pch  = 20,
     col  = "gray30",
     xlab = "Valor ajustado posterior",
     ylab = "Residual",
     ylim = 6 * c(-1, 1),
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2,
     col = "gray20"
)

box()

dev.off()

# Residuales estandarizados contra valores ajustados

pdf(
     file      = "nitro_residuales_estandarizados_ajustados.pdf",
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
     x    = mu_hat,
     y    = resid_std,
     pch  = 20,
     col  = "gray30",
     ylim = 4 * c(-1, 1),
     xlab = "Valor ajustado posterior",
     ylab = "Residual estandarizado",
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2,
     col = "gray20"
)

abline(
     h   = c(-2, 2),
     lty = 3,
     lwd = 1.5,
     col = "gray50"
)

box()

dev.off()

# Histograma de residuales estandarizados

pdf(
     file      = "nitro_residuales_estandarizados_histograma.pdf",
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
     x      = resid_std,
     breaks = "FD",
     freq   = FALSE,
     col    = "gray90",
     border = "white",
     xlab   = "Residual estandarizado",
     ylab   = "Densidad",
     xlim   = 4 * c(-1, 1),
     main   = ""
)

curve(
     dnorm(x),
     add = TRUE,
     n   = 1000,
     lwd = 2,
     col = "gray20"
)

box()

dev.off()

# Gráfico Q-Q Normal de residuales estandarizados

pdf(
     file      = "nitro_residuales_estandarizados_qqplot.pdf",
     width     = 5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

qqnorm(
     y    = resid_std,
     pch  = 20,
     col  = "gray30",
     xlab = "Cuantiles teóricos normales",
     ylab = "Cuantiles muestrales",
     xlim = 3 * c(-1, 1),
     ylim = 3 * c(-1, 1),
     main = ""
)

qqline(
     y   = resid_std,
     lty = 2,
     lwd = 2,
     col = "gray20"
)

box()

dev.off()

# Residuales estandarizados por finca

pdf(
     file      = "nitro_residuales_por_finca.pdf",
     width     = 7,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

boxplot(
     resid_std ~ farm,
     col      = "gray90",
     border   = "gray30",
     boxwex   = 0.55,
     outline  = TRUE,
     outpch   = 16,
     outcex   = 0.5,
     cex.axis = 0.7,
     las      = 2,
     xlab     = "Finca",
     ylab     = "Residual estandarizado",
     main     = ""
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2,
     col = "gray20"
)

abline(
     h   = c(-2, 2),
     lty = 3,
     lwd = 1.5,
     col = "gray50"
)

box()

dev.off()

# Fin --------------------------------------------------------------------------