# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(boot)))

# Cargar datos
dat <- read.csv(
     file         = "~/Dropbox/UN/bayes_book/Ocupados.CSV",
     sep          = ";",
     comment.char = "#"
)

# Seleccionar ingresos laborales de individuos en Bogotá D.C.
y <- dat$INGLABO[dat$DPTO == 11]

rm(dat)

# Remover individuos que no reportaron ingresos
y <- y[!is.na(y)]

# Remover individuos con ingresos iguales a cero
y <- y[y > 0]

# Expresar ingresos en millones de pesos
y <- y / 1e6

# Análisis descriptivo ---------------------------------------------------------

summary(y)

n      <- length(y)
mean_y <- mean(y)
medi_y <- median(y)
sd_y   <- sd(y)
cv_y   <- 100 * sd_y / mean_y

n
round(mean_y, 3)
round(medi_y, 3)
round(sd_y, 3)
round(cv_y, 3)

# Visualización de los datos ---------------------------------------------------

pdf(file = "GEIH_hist_datos.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

hist(
     y,
     freq   = FALSE,
     xlim   = c(0, 10),
     breaks = 30,
     xlab   = "Ingresos",
     ylab   = "Densidad",
     main   = "",
     col    = "gray90",
     border = "white"
)

abline(v = mean_y, lty = 1, col = 2, lwd = 2)
abline(v = medi_y, lty = 1, col = 4, lwd = 2)

legend(
     "topright",
     bty    = "n",
     lty    = 1,
     lwd    = 2,
     col    = c(2, 4),
     legend = c("Media", "Mediana")
)

dev.off()

# Teorema del límite central ---------------------------------------------------

# Error estándar de la media muestral
se_y <- sd_y / sqrt(n)

# Intervalo de confianza aproximado al 95% basado en el TLC
ci_theta_TLC <- c(
     mean_y - qnorm(0.975) * se_y,
     mean_y + qnorm(0.975) * se_y
)

round(ci_theta_TLC, 3)

# Modelo Gamma -----------------------------------------------------------------

# Ecuación de verosimilitud para el modelo Gamma
score_beta_G <- function(beta) {
     n * log(beta) - n * digamma(beta * mean_y) + sum(log(y))
}

# Estimar beta
beta_mle_G <- uniroot(
     f         = score_beta_G,
     interval  = c(1e-8, 100),
     extendInt = "yes"
)$root

# Estimar alpha
alpha_mle_G <- beta_mle_G * mean_y

round(c(alpha_mle_G, beta_mle_G), 3)

# Log-verosimilitud del modelo Gamma
ll_G <- function(alpha, beta, y) {
     sum(dgamma(y, shape = alpha, rate = beta, log = TRUE))
}

# Log-verosimilitud negativa
m_ll_G <- function(pars, data) {
     alpha <- pars[1]
     beta  <- pars[2]
     
     -ll_G(alpha = alpha, beta = beta, y = data)
}

# Optimización numérica
fit_G <- optim(
     par     = c(alpha_mle_G, beta_mle_G),
     fn      = m_ll_G,
     data    = y,
     method  = "L-BFGS-B",
     lower   = c(1e-8, 1e-8),
     hessian = TRUE
)

# Verificar convergencia
if (fit_G$convergence != 0) {
     warning("El algoritmo de optimización no convergió para el modelo Gamma.")
}

# Estimar parámetros
alpha_mle_G <- fit_G$par[1]
beta_mle_G  <- fit_G$par[2]

round(c(alpha_mle_G, beta_mle_G), 3)

# Información observada
J_G <- fit_G$hessian

round(J_G, 3)

# Matriz Hessiana
H_G <- -J_G

H_G[1, 1]
det(H_G)

# Matriz de varianzas y covarianzas
SI_G <- solve(J_G)

round(SI_G, 6)

# MLE de la media bajo el modelo Gamma
theta_mle_G <- alpha_mle_G / beta_mle_G

round(theta_mle_G, 3)

# Gradiente de theta = alpha / beta
delta_G <- as.matrix(c(1 / beta_mle_G, -alpha_mle_G / beta_mle_G^2))

# Error estándar por método Delta
se_theta_mle_G <- sqrt(t(delta_G) %*% SI_G %*% delta_G)

round(se_theta_mle_G, 3)

# Intervalo de confianza aproximado al 95%
ci_theta_G <- c(
     theta_mle_G - qnorm(0.975) * se_theta_mle_G,
     theta_mle_G + qnorm(0.975) * se_theta_mle_G
)

round(ci_theta_G, 3)

# Modelo Log-normal ------------------------------------------------------------

# Ingresos en escala logarítmica
log_y <- log(y)

# Estimaciones MLE bajo el modelo Log-normal
mu_mle_LN     <- mean(log_y)
sigma2_mle_LN <- mean((log_y - mu_mle_LN)^2)
sigma_mle_LN  <- sqrt(sigma2_mle_LN)

round(mu_mle_LN, 3)
round(sigma2_mle_LN, 3)

# MLE de la media bajo el modelo Log-normal
theta_mle_LN <- exp(mu_mle_LN + sigma2_mle_LN / 2)

round(theta_mle_LN, 3)

# Matriz de varianzas y covarianzas
SI_LN <- rbind(
     c(sigma2_mle_LN / n, 0),
     c(0, sigma2_mle_LN / (2 * n))
)

# Gradiente de theta = exp(mu + sigma^2 / 2)
delta_LN <- as.matrix(c(theta_mle_LN, sigma_mle_LN * theta_mle_LN))

# Error estándar por método Delta
se_theta_mle_LN <- sqrt(t(delta_LN) %*% SI_LN %*% delta_LN)

round(se_theta_mle_LN, 3)

# Intervalo de confianza aproximado al 95%
ci_theta_LN <- c(
     theta_mle_LN - qnorm(0.975) * se_theta_mle_LN,
     theta_mle_LN + qnorm(0.975) * se_theta_mle_LN
)

round(ci_theta_LN, 3)

# Modelo Gaussiano inverso -----------------------------------------------------

# Estimar mu
mu_mle_GI <- mean_y

# Estimar lambda
lambda_mle_GI <- n / sum((y - mu_mle_GI)^2 / (mu_mle_GI^2 * y))

round(c(mu_mle_GI, lambda_mle_GI), 3)

# Log-verosimilitud del modelo Gaussiano inverso
ll_GI <- function(mu, lambda, y) {
     sum(
          0.5 * log(lambda) -
               0.5 * log(2 * pi) -
               1.5 * log(y) -
               lambda * (y - mu)^2 / (2 * mu^2 * y)
     )
}

# Log-verosimilitud negativa
m_ll_GI <- function(pars, data) {
     mu     <- pars[1]
     lambda <- pars[2]
     
     -ll_GI(mu = mu, lambda = lambda, y = data)
}

# Optimización numérica
fit_GI <- optim(
     par     = c(mu_mle_GI, lambda_mle_GI),
     fn      = m_ll_GI,
     data    = y,
     method  = "L-BFGS-B",
     lower   = c(1e-8, 1e-8),
     hessian = TRUE
)

# Verificar convergencia
if (fit_GI$convergence != 0) {
     warning("El algoritmo de optimización no convergió para el modelo Gaussiano inverso.")
}

# Estimar parámetros
mu_mle_GI     <- fit_GI$par[1]
lambda_mle_GI <- fit_GI$par[2]

round(c(mu_mle_GI, lambda_mle_GI), 3)

# Información observada
J_GI <- fit_GI$hessian

round(J_GI, 3)

# Matriz de varianzas y covarianzas
SI_GI <- solve(J_GI)

round(SI_GI, 6)

# MLE de la media bajo el modelo Gaussiano inverso
theta_mle_GI <- mu_mle_GI

round(theta_mle_GI, 3)

# Gradiente de theta = mu
delta_GI <- as.matrix(c(1, 0))

# Error estándar por método Delta
se_theta_mle_GI <- sqrt(t(delta_GI) %*% SI_GI %*% delta_GI)

round(se_theta_mle_GI, 3)

# Intervalo de confianza aproximado al 95%
ci_theta_GI <- c(
     theta_mle_GI - qnorm(0.975) * se_theta_mle_GI,
     theta_mle_GI + qnorm(0.975) * se_theta_mle_GI
)

round(ci_theta_GI, 3)

# Modelo Weibull ---------------------------------------------------------------

# Valores iniciales bajo el modelo Weibull
shape_init_W <- 1
scale_init_W <- mean_y

# Log-verosimilitud del modelo Weibull
ll_W <- function(shape, scale, y) {
     sum(dweibull(y, shape = shape, scale = scale, log = TRUE))
}

# Log-verosimilitud negativa
m_ll_W <- function(pars, data) {
     shape <- pars[1]
     scale <- pars[2]
     
     -ll_W(shape = shape, scale = scale, y = data)
}

# Optimización numérica
fit_W <- optim(
     par     = c(shape_init_W, scale_init_W),
     fn      = m_ll_W,
     data    = y,
     method  = "L-BFGS-B",
     lower   = c(1e-8, 1e-8),
     hessian = TRUE
)

# Verificar convergencia
if (fit_W$convergence != 0) {
     warning("El algoritmo de optimización no convergió para el modelo Weibull.")
}

# Estimar parámetros
shape_mle_W <- fit_W$par[1]
scale_mle_W <- fit_W$par[2]

round(c(shape_mle_W, scale_mle_W), 3)

# Información observada
J_W <- fit_W$hessian

round(J_W, 3)

# Matriz de varianzas y covarianzas
SI_W <- solve(J_W)

round(SI_W, 6)

# MLE de la media bajo el modelo Weibull
theta_mle_W <- scale_mle_W * gamma(1 + 1 / shape_mle_W)

round(theta_mle_W, 3)

# Gradiente de theta = scale * Gamma(1 + 1 / shape)
a_W <- 1 + 1 / shape_mle_W

delta_W <- as.matrix(c(
     -scale_mle_W * gamma(a_W) * digamma(a_W) / shape_mle_W^2,
     gamma(a_W)
))

# Error estándar por método Delta
se_theta_mle_W <- sqrt(t(delta_W) %*% SI_W %*% delta_W)

round(se_theta_mle_W, 3)

# Intervalo de confianza aproximado al 95%
ci_theta_W <- c(
     theta_mle_W - qnorm(0.975) * se_theta_mle_W,
     theta_mle_W + qnorm(0.975) * se_theta_mle_W
)

round(ci_theta_W, 3)

# Modelo Dagum -----------------------------------------------------------------

# Valores iniciales bajo el modelo Dagum
a_init_D <- 2
b_init_D <- medi_y
p_init_D <- 1

# Log-verosimilitud del modelo Dagum
ll_D <- function(a, b, p, y) {
     sum(
          log(a) + log(p) - log(y) +
               a * p * log(y / b) -
               (p + 1) * log1p((y / b)^a)
     )
}

# Log-verosimilitud negativa
m_ll_D <- function(pars, data) {
     a <- pars[1]
     b <- pars[2]
     p <- pars[3]
     
     -ll_D(a = a, b = b, p = p, y = data)
}

# Optimización numérica
fit_D <- optim(
     par     = c(a_init_D, b_init_D, p_init_D),
     fn      = m_ll_D,
     data    = y,
     method  = "L-BFGS-B",
     lower   = c(1 + 1e-8, 1e-8, 1e-8),
     hessian = TRUE
)

# Verificar convergencia
if (fit_D$convergence != 0) {
     warning("El algoritmo de optimización no convergió para el modelo Dagum.")
}

# Estimar parámetros
a_mle_D <- fit_D$par[1]
b_mle_D <- fit_D$par[2]
p_mle_D <- fit_D$par[3]

round(c(a_mle_D, b_mle_D, p_mle_D), 3)

# Información observada
J_D <- fit_D$hessian

round(J_D, 3)

# Matriz de varianzas y covarianzas
SI_D <- solve(J_D)

round(SI_D, 6)

# MLE de la media bajo el modelo Dagum
theta_mle_D <- b_mle_D *
     gamma(p_mle_D + 1 / a_mle_D) *
     gamma(1 - 1 / a_mle_D) /
     gamma(p_mle_D)

round(theta_mle_D, 3)

# Gradiente de theta respecto a a, b y p
delta_D <- as.matrix(c(
     theta_mle_D / a_mle_D^2 *
          (digamma(1 - 1 / a_mle_D) - digamma(p_mle_D + 1 / a_mle_D)),
     theta_mle_D / b_mle_D,
     theta_mle_D *
          (digamma(p_mle_D + 1 / a_mle_D) - digamma(p_mle_D))
))

# Error estándar por método Delta
se_theta_mle_D <- sqrt(t(delta_D) %*% SI_D %*% delta_D)

round(se_theta_mle_D, 3)

# Intervalo de confianza aproximado al 95%
ci_theta_D <- c(
     theta_mle_D - qnorm(0.975) * se_theta_mle_D,
     theta_mle_D + qnorm(0.975) * se_theta_mle_D
)

round(ci_theta_D, 3)

# Bootstrap --------------------------------------------------------------------

# Número de remuestras bootstrap
M <- 10000

# Almacenar medias bootstrap
bs_mean_vec <- rep(NA_real_, M)

# Bootstrap no paramétrico para la media
set.seed(42)

for (i in seq_len(M)) {
     bs_y          <- sample(y, size = n, replace = TRUE)
     bs_mean_vec[i] <- mean(bs_y)
}

# Media bootstrap
theta_boot <- mean(bs_mean_vec)

round(theta_boot, 3)

# Error estándar bootstrap
se_boot <- sd(bs_mean_vec)

round(se_boot, 3)

# Intervalo de confianza percentil al 95%
ci_boot <- quantile(
     bs_mean_vec,
     probs = c(0.025, 0.975),
     names = FALSE
)

round(ci_boot, 3)

# AIC y BIC --------------------------------------------------------------------

# Funciones para AIC y BIC
get_AIC <- function(loglik, k) {
     -2 * loglik + 2 * k
}

get_BIC <- function(loglik, k, n) {
     -2 * loglik + k * log(n)
}

# AIC y BIC para el modelo Gamma
loglik_G <- ll_G(alpha_mle_G, beta_mle_G, y)

k_G <- 2

AIC_G <- get_AIC(loglik_G, k_G)
BIC_G <- get_BIC(loglik_G, k_G, n)

round(c(loglik_G = loglik_G, AIC_G = AIC_G, BIC_G = BIC_G), 3)

# AIC y BIC para el modelo Log-normal
loglik_LN <- sum(dlnorm(y, meanlog = mu_mle_LN, sdlog = sigma_mle_LN, log = TRUE))

k_LN <- 2

AIC_LN <- get_AIC(loglik_LN, k_LN)
BIC_LN <- get_BIC(loglik_LN, k_LN, n)

round(c(loglik_LN = loglik_LN, AIC_LN = AIC_LN, BIC_LN = BIC_LN), 3)

# AIC y BIC para el modelo Gaussiano inverso
loglik_GI <- ll_GI(mu_mle_GI, lambda_mle_GI, y)

k_GI <- 2

AIC_GI <- get_AIC(loglik_GI, k_GI)
BIC_GI <- get_BIC(loglik_GI, k_GI, n)

round(c(loglik_GI = loglik_GI, AIC_GI = AIC_GI, BIC_GI = BIC_GI), 3)

# AIC y BIC para el modelo Weibull
loglik_W <- ll_W(shape_mle_W, scale_mle_W, y)

k_W <- 2

AIC_W <- get_AIC(loglik_W, k_W)
BIC_W <- get_BIC(loglik_W, k_W, n)

round(c(loglik_W = loglik_W, AIC_W = AIC_W, BIC_W = BIC_W), 3)

# AIC y BIC para el modelo Dagum
loglik_D <- ll_D(a_mle_D, b_mle_D, p_mle_D, y)

k_D <- 3

AIC_D <- get_AIC(loglik_D, k_D)
BIC_D <- get_BIC(loglik_D, k_D, n)

round(c(loglik_D = loglik_D, AIC_D = AIC_D, BIC_D = BIC_D), 3)

# Tabla de resumen -------------------------------------------------------------

# Tabla resumen de estimaciones, intervalos y criterios
results_table <- data.frame(
     model = c(
          "TLC", "Bootstrap", "Gamma", "Log-normal",
          "Gaussiano inverso", "Weibull", "Dagum"
     ),
     theta_hat = c(
          mean_y,
          as.numeric(theta_boot),
          as.numeric(theta_mle_G),
          as.numeric(theta_mle_LN),
          as.numeric(theta_mle_GI),
          as.numeric(theta_mle_W),
          as.numeric(theta_mle_D)
     ),
     se = c(
          as.numeric(se_y),
          as.numeric(se_boot),
          as.numeric(se_theta_mle_G),
          as.numeric(se_theta_mle_LN),
          as.numeric(se_theta_mle_GI),
          as.numeric(se_theta_mle_W),
          as.numeric(se_theta_mle_D)
     ),
     lower = c(
          ci_theta_TLC[1],
          ci_boot[1],
          ci_theta_G[1],
          ci_theta_LN[1],
          ci_theta_GI[1],
          ci_theta_W[1],
          ci_theta_D[1]
     ),
     upper = c(
          ci_theta_TLC[2],
          ci_boot[2],
          ci_theta_G[2],
          ci_theta_LN[2],
          ci_theta_GI[2],
          ci_theta_W[2],
          ci_theta_D[2]
     ),
     AIC = c(NA, NA, AIC_G, AIC_LN, AIC_GI, AIC_W, AIC_D),
     BIC = c(NA, NA, BIC_G, BIC_LN, BIC_GI, BIC_W, BIC_D)
)

# Longitud del intervalo
results_table$length <- results_table$upper - results_table$lower

# Reordenar columnas
results_table <- results_table[, c(
     "model", "theta_hat", "se", "lower", "upper", "length", "AIC", "BIC"
)]

# Redondear resultados
results_table_round       <- results_table
results_table_round[, -1] <- round(results_table_round[, -1], 3)

results_table_round

# Visualización de los modelos -------------------------------------------------

# Densidad del modelo Gaussiano inverso
dginv <- function(x, mu, lambda) {
     dens <- rep(0, length(x))
     idx  <- x > 0
     
     dens[idx] <- sqrt(lambda / (2 * pi * x[idx]^3)) *
          exp(-lambda * (x[idx] - mu)^2 / (2 * mu^2 * x[idx]))
     
     dens
}

# Densidad del modelo Dagum
ddagum <- function(x, a, b, p) {
     dens <- rep(0, length(x))
     idx  <- x > 0
     
     dens[idx] <- (a * p / x[idx]) *
          (x[idx] / b)^(a * p) /
          (1 + (x[idx] / b)^a)^(p + 1)
     
     dens
}

# Grilla para superponer densidades
x_grid <- seq(from = 0, to = 10, length.out = 10000)

# Histograma con curvas ajustadas
pdf(file = "GEIH_hist_modelos.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

hist(
     y,
     freq   = FALSE,
     xlim   = c(0, 10),
     breaks = 30,
     xlab   = "Ingresos",
     ylab   = "Densidad",
     main   = "",
     col    = "gray90",
     border = "white"
)

lines(
     x_grid,
     dgamma(x_grid, shape = alpha_mle_G, rate = beta_mle_G),
     col = 2,
     lwd = 2
)

lines(
     x_grid,
     dlnorm(x_grid, meanlog = mu_mle_LN, sdlog = sigma_mle_LN),
     col = 3,
     lwd = 2
)

lines(
     x_grid,
     dginv(x_grid, mu = mu_mle_GI, lambda = lambda_mle_GI),
     col = 4,
     lwd = 2
)

lines(
     x_grid,
     dweibull(x_grid, shape = shape_mle_W, scale = scale_mle_W),
     col = 5,
     lwd = 2
)

lines(
     x_grid,
     ddagum(x_grid, a = a_mle_D, b = b_mle_D, p = p_mle_D),
     col = 6,
     lwd = 2
)

legend(
     "topright",
     bty    = "n",
     lty    = 1,
     lwd    = 2,
     col    = 2:6,
     legend = c("Gamma", "L-Nor", "G-Inv", "Weibull", "Dagum")
)

dev.off()

# Fin --------------------------------------------------------------------------