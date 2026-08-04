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

# Muestreador de Gibbs (previa difusa y unitaria) ------------------------------

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

# Monte Carlo (previa g) -------------------------------------------------------

mc_previag <- function(
          y, X, g, nu0, sigma20,
          n_sams, verbose = TRUE
) {
     # Ajustes
     y <- as.numeric(y)
     X <- as.matrix(X)
     n <- nrow(X)
     p <- ncol(X)
     
     # Número de simulaciones
     B    <- n_sams
     ncat <- max(1, floor(0.1 * B))
     
     # Cantidades fijas
     XtX  <- crossprod(X)
     Xty  <- crossprod(X, y)
     iXtX <- solve(XtX)
     
     H_g <- (g / (g + 1)) * X %*% iXtX %*% t(X)
     
     SSR_g <- as.numeric(t(y) %*% (diag(n) - H_g) %*% y)
     
     V_beta <- (g / (g + 1)) * iXtX
     E_beta <- as.numeric(V_beta %*% Xty)
     
     a_sigma <- 0.5 * (nu0 + n)
     b_sigma <- 0.5 * (nu0 * sigma20 + SSR_g)
     
     # Almacenamiento
     BETA   <- matrix(data = NA, nrow = B, ncol = p)
     SIGMA2 <- rep(NA, B)
     LL     <- rep(NA, B)
     
     colnames(BETA) <- colnames(X)
     
     # Simulación Monte Carlo directa
     for (b in 1:B) {
          sigma2 <- 1 / rgamma(
               n     = 1,
               shape = a_sigma,
               rate  = b_sigma
          )
          
          beta <- as.numeric(
               MASS::mvrnorm(
                    n     = 1,
                    mu    = E_beta,
                    Sigma = sigma2 * V_beta
               )
          )
          
          ll <- sum(
               dnorm(
                    x    = y,
                    mean = as.numeric(X %*% beta),
                    sd   = sqrt(sigma2),
                    log  = TRUE
               )
          )
          
          BETA[b, ] <- beta
          SIGMA2[b] <- sigma2
          LL[b]     <- ll
          
          if (verbose && b %% ncat == 0) {
               cat(sprintf("%.1f%% completado\n", 100 * b / B))
          }
     }
     
     # Salida
     return(
          list(
               BETA   = BETA,
               SIGMA2 = SIGMA2,
               LL     = LL,
               g      = g,
               SSR_g  = SSR_g
          )
     )
}

# Monte Carlo (previa impropia) ------------------------------------------------

mc_impropia <- function(y, X, n_sams, verbose = TRUE) {
     # Ajustes
     y <- as.numeric(y)
     X <- as.matrix(X)
     n <- nrow(X)
     p <- ncol(X)
     
     # Número de simulaciones
     B    <- n_sams
     ncat <- max(1, floor(0.1 * B))
     
     # Cantidades fijas
     XtX  <- crossprod(X)
     Xty  <- crossprod(X, y)
     iXtX <- solve(XtX)
     
     beta_ols <- as.numeric(iXtX %*% Xty)
     resid    <- as.numeric(y - X %*% beta_ols)
     RSS_ols  <- sum(resid^2)
     
     V_beta <- iXtX
     
     a_sigma <- 0.5 * (n - p)
     b_sigma <- 0.5 * RSS_ols
     
     # Almacenamiento
     BETA   <- matrix(data = NA, nrow = B, ncol = p)
     SIGMA2 <- rep(NA, B)
     LL     <- rep(NA, B)
     
     colnames(BETA) <- colnames(X)
     
     # Simulación de Monte Carlo directa
     for (b in 1:B) {
          sigma2 <- 1 / rgamma(
               n     = 1,
               shape = a_sigma,
               rate  = b_sigma
          )
          
          beta <- as.numeric(
               MASS::mvrnorm(
                    n     = 1,
                    mu    = beta_ols,
                    Sigma = sigma2 * V_beta
               )
          )
          
          ll <- sum(
               dnorm(
                    x    = y,
                    mean = as.numeric(X %*% beta),
                    sd   = sqrt(sigma2),
                    log  = TRUE
               )
          )
          
          BETA[b, ] <- beta
          SIGMA2[b] <- sigma2
          LL[b]     <- ll
          
          if (verbose && b %% ncat == 0) {
               cat(sprintf("%.1f%% completado\n", 100 * b / B))
          }
     }
     
     # Salida
     return(
          list(
               BETA     = BETA,
               SIGMA2   = SIGMA2,
               LL       = LL,
               RSS_ols  = RSS_ols,
               beta_ols = beta_ols
          )
     )
}

# Ajuste del modelo lineal Normal previa difusa --------------------------------

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

chain_difusa <- mcmc(
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

# Ajuste del modelo lineal Normal previa unitaria ------------------------------

# Hiperparámetros
beta0   <- as.numeric(beta_ols)
Sigma0  <- n * sig2_ols * solve(XtX)
nu0     <- 1
sigma20 <- sig2_ols

# Número de iteraciones
n_sams <- 10000
n_burn <- 10000
n_skip <- 10

# Ajuste del modelo Bayesiano
set.seed(123)

chain_unitaria <- mcmc(
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

# Ajuste del modelo lineal Normal previa g -------------------------------------

# Hiperparámetros
g       <- n
nu0     <- 1
sigma20 <- sig2_ols

# Número de simulaciones
n_sams <- 10000

# Ajuste del modelo Bayesiano con previa g
set.seed(123)

chain_previag <- mc_previag(
     y       = y,
     X       = X,
     g       = g,
     nu0     = nu0,
     sigma20 = sigma20,
     n_sams  = n_sams,
     verbose = TRUE
)

# Ajuste del modelo lineal Normal previa impropia ------------------------------

# Hiperparámetros
g       <- n
nu0     <- 1
sigma20 <- sig2_ols

# Número de simulaciones
n_sams <- 10000

# Ajuste del modelo Bayesiano con previa g
set.seed(123)

chain_impropia <- mc_impropia(
     y       = y,
     X       = X,
     n_sams  = n_sams,
     verbose = TRUE
)

# Cadenas para el análisis de sensibilidad ------------------------------------

chains <- list(
     "Difusa"   = chain_difusa,
     "Unitaria" = chain_unitaria,
     "Previa g" = chain_previag,
     "Impropia" = chain_impropia
)

cols <- c(
     "black",
     "red3",
     "blue3",
     "darkgreen"
)

# Parámetros asociados con covariables
parametros_beta <- colnames(chain_difusa$BETA)

# Gráficos de densidad posterior por covariable --------------------------------

for (param in parametros_beta) {
     densidades <- lapply(
          chains,
          function(chain) density(chain$BETA[, param])
     )
     
     xlim_param <- range(
          unlist(
               lapply(densidades, function(d) d$x)
          )
     )
     
     ylim_param <- c(
          0,
          max(
               unlist(
                    lapply(densidades, function(d) d$y)
               )
          )
     )
     
     nombre_archivo <- paste0(
          "lungcap_regresion_lineal_normal_sensibilidad_densidad_beta_",
          param,
          ".pdf"
     )
     
     pdf(
          file      = nombre_archivo,
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
          densidades[[1]],
          type = "n",
          xlim = xlim_param,
          ylim = ylim_param,
          xlab = paste0("Coeficiente de ", param),
          ylab = "Densidad",
          main = ""
     )
     
     for (j in seq_along(densidades)) {
          lines(
               densidades[[j]],
               col = cols[j],
               lwd = 2
          )
     }
     
     dev.off()
}

# Densidades posteriores de sigma2 --------------------------------------------

densidades_sigma2 <- lapply(
     chains,
     function(chain) density(chain$SIGMA2)
)

xlim_sigma2 <- range(
     unlist(
          lapply(densidades_sigma2, function(d) d$x)
     )
)

ylim_sigma2 <- c(
     0,
     max(
          unlist(
               lapply(densidades_sigma2, function(d) d$y)
          )
     )
)

pdf(
     file      = "lungcap_regresion_lineal_normal_sensibilidad_densidad_sigma2.pdf",
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
     densidades_sigma2[[1]],
     type = "n",
     xlim = xlim_sigma2,
     ylim = ylim_sigma2,
     xlab = expression(sigma^2),
     ylab = "Densidad",
     main = ""
)

for (j in seq_along(densidades_sigma2)) {
     lines(
          densidades_sigma2[[j]],
          col = cols[j],
          lwd = 2
     )
}

legend(
     "topright",
     legend = names(chains),
     col    = cols,
     fill   = cols,
     border = cols,
     bty    = "n",
     cex    = 0.9
)

dev.off()

# Fin --------------------------------------------------------------------------