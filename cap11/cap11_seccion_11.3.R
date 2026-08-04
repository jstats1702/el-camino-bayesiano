# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Paquetes
library(maps)
library(sf)
library(spdep)
library(Matrix)
library(ggplot2)


# Datos ------------------------------------------------------------------------

flu <- read.csv("flu_data.csv", header = TRUE)    # Base completa
flu <- flu[485:496, ]                             # Primeras 12 semanas de 2013
flu <- flu[, -c(1, 2, 12)]                        # Sin fecha, Alaska ni Hawái

states.names <- colnames(flu)

data <- as.matrix(flu)
p    <- nrow(data)                                # T: número de semanas
m    <- ncol(data)                                # I: número de estados


# Regiones del Departamento de Salud y Servicios Humanos de Estados Unidos ------

regions <- read.csv("flu_regions.csv", header = TRUE)     # Base completa
regions <- regions[-c(2, 12), ]                           # Sin Alaska ni Hawái
regions <- as.matrix(regions)

e <- length(unique(regions))                              # Número de regiones


# Distribución por raza --------------------------------------------------------

whites <- read.csv("flu_race.csv", header = TRUE)    # Base completa
whites <- whites[, 2]                                # Porcentaje de población blanca
whites <- whites[-c(2, 12)]                          # Sin Alaska ni Hawái
whites <- as.matrix(whites)


# Distribución por edad --------------------------------------------------------

old <- read.csv("flu_age.csv", header = TRUE)    # Base completa
old <- old[, 4]                                  # Porcentaje de población mayor de 65 años
old <- old[-c(2, 12)]                            # Sin Alaska ni Hawái
old <- as.matrix(old)


# Respuesta binaria ------------------------------------------------------------

y <- matrix(0, p, m)

for (i in 1:p) {
     for (j in 1:m) {
          if (data[i, j] > 7500) {
               y[i, j] <- 1
          }
     }
}

y <- as.matrix(c(y))


# Vecindad entre estados -------------------------------------------------------

# Usar geometría planar en lugar de s2
sf::sf_use_s2(FALSE)

# Cargar polígonos de estados
usa.state <- map(
     database = "state",
     fill     = TRUE,
     plot     = FALSE
)

# Extraer identificadores de estados
state.ID <- sapply(
     strsplit(usa.state$names, ":"),
     function(x) x[1]
)

# Convertir el mapa a objeto sf
usa.sf <- sf::st_as_sf(
     usa.state,
     IDs = state.ID
)

# Corregir geometrías inválidas
usa.sf <- sf::st_make_valid(usa.sf)

# Proyectar a coordenadas planas
usa.sf <- sf::st_transform(
     usa.sf,
     crs = 5070
)

# Construir lista de vecinos
usa.nb <- spdep::poly2nb(
     usa.sf,
     queen = TRUE
)

# Convertir vecindades a matriz binaria
usa.adj.mat <- A <- spdep::nb2mat(
     usa.nb,
     style = "B"
)   # Matriz de vecindad


# Vecindad entre regiones ------------------------------------------------------

A2 <- c(
     0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
     1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
     0, 0, 0, 1, 0, 0, 1, 1, 1, 0,
     0, 0, 0, 1, 1, 1, 0, 1, 0, 0,
     0, 0, 0, 0, 1, 1, 1, 0, 1, 1,
     0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
     0, 0, 0, 0, 0, 0, 0, 1, 1, 0
)

A2 <- matrix(
     A2,
     nrow  = e,
     ncol  = e,
     byrow = TRUE
)


# Gráfico de la matriz de vecindad entre estados -------------------------------

pdf(
     file      = "flu_matriz_vecindad_estados.pdf",
     height    = 5,
     width     = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(1.4, 1.4, 0.4, 0.4),
     mgp   = c(1.75, 0.75, 0)
)

colorscale <- c("white", rev(heat.colors(100)))

image(
     x    = 1:m,
     y    = 1:m,
     z    = A,
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "",
     col  = colorscale[seq(floor(100 * min(A)), floor(100 * max(A)))]
)

axis(
     side     = 1,
     at       = 1:m,
     labels   = NA,
     las      = 2,
     cex.axis = 0.4
)

axis(
     side     = 2,
     at       = 1:m,
     labels   = NA,
     las      = 2,
     cex.axis = 0.4
)

box()

dev.off()


# Gráfico de la matriz de vecindad entre regiones ------------------------------

pdf(
     file      = "flu_matriz_vecindad_regiones.pdf",
     height    = 5,
     width     = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(1.4, 1.4, 0.4, 0.4),
     mgp   = c(1.75, 0.75, 0)
)

colorscale <- c("white", rev(heat.colors(100)))

image(
     x    = 1:e,
     y    = 1:e,
     z    = A2,
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "",
     col  = colorscale[seq(floor(100 * min(A2)), floor(100 * max(A2)))]
)

axis(
     side     = 1,
     at       = 1:e,
     labels   = 1:e,
     las      = 2,
     cex.axis = 0.75
)

axis(
     side     = 2,
     at       = 1:e,
     labels   = 1:e,
     las      = 2,
     cex.axis = 0.75
)

box()

dev.off()


# Matriz de precisión de estados ----------------------------------------------

W <- -A

diag(W) <- rowSums(A)

r     <- rankMatrix(W)[1]
n.nei <- rowSums(A)


# Matriz de precisión de regiones ---------------------------------------------

W2 <- -A2

diag(W2) <- rowSums(A2)

r2     <- rankMatrix(W2)[1]
n.nei2 <- rowSums(A2)


# Mapa de proporción de población blanca ---------------------------------------

# Mapa base de estados
all_states <- map_data("state")

# Unir covariable con mapa
whites <- read.csv("flu_race.csv", header = TRUE)    # Base completa
whites <- whites[, 2]                                # Porcentaje de población blanca
whites <- whites[-c(2, 12)]                          # Sin Alaska ni Hawái
whites <- as.data.frame(whites)
whites <- cbind(rownames(A), whites)

colnames(whites) <- c("region", "percent")

whites <- merge(
     all_states,
     whites,
     by = "region"
)

head(whites)

# Visualización
q <- ggplot() +
     geom_polygon(
          data   = whites,
          aes(
               x     = long,
               y     = lat,
               group = group,
               fill  = percent
          ),
          colour = "black"
     ) +
     scale_fill_continuous(
          low   = "thistle2",
          high  = "darkred",
          guide = "colorbar"
     ) +
     labs(
          fill  = "%",
          title = "",
          x     = "",
          y     = ""
     ) +
     theme(
          panel.background = element_rect(colour = "black"),
          panel.grid.major = element_line(size = 1),
          plot.margin      = margin(0, 0, 0, 0)
     ) +
     theme_minimal(base_size = 30)

q

ggsave(
     filename = "flu_mapa_proporcion_blancos.pdf",
     plot     = q,
     width    = 10,
     height   = 7
)

dev.off()


# Mapa de proporción de población mayor de 65 años -----------------------------

# Mapa base de estados
all_states <- map_data("state")

# Unir covariable con mapa
old <- read.csv("flu_age.csv", header = TRUE)    # Base completa
old <- old[, 4]                                  # Porcentaje de población mayor de 65 años
old <- old[-c(2, 12)]                            # Sin Alaska ni Hawái
old <- as.data.frame(old)
old <- cbind(rownames(A), old)

colnames(old) <- c("region", "percent")

old <- merge(
     all_states,
     old,
     by = "region"
)

head(old)

# Visualización
q <- ggplot() +
     geom_polygon(
          data   = old,
          aes(
               x     = long,
               y     = lat,
               group = group,
               fill  = percent
          ),
          colour = "black"
     ) +
     scale_fill_continuous(
          low   = "yellow1",
          high  = "yellow4",
          guide = "colorbar"
     ) +
     labs(
          fill  = "%",
          title = "",
          x     = "",
          y     = ""
     ) +
     theme(
          panel.background = element_rect(colour = "black"),
          panel.grid.major = element_line(size = 1),
          plot.margin      = margin(0, 0, 0, 0)
     ) +
     theme_minimal(base_size = 30)

q

ggsave(
     filename = "flu_mapa_proporcion_mayores_65.pdf",
     plot     = q,
     width    = 10,
     height   = 7
)

dev.off()


# Mapa de la respuesta por semana ----------------------------------------------

# Mapa base de estados
all_states <- map_data("state")

# Respuesta binaria
y <- matrix(0, p, m)

for (i in 1:p) {
     for (j in 1:m) {
          if (data[i, j] > 7500) {
               y[i, j] <- 1
          }
     }
}

y <- t(y)

# Visualización por semana
for (tt in 1:p) {
     
     # Unir respuesta con mapa
     yt <- as.data.frame(as.factor(y[, tt]))
     yt <- cbind(rownames(A), yt)
     
     colnames(yt) <- c("region", paste("week_", tt, sep = ""))
     
     yt <- merge(
          all_states,
          yt,
          by = "region"
     )
     
     q <- ggplot() +
          geom_polygon(
               data   = yt,
               aes(
                    x     = long,
                    y     = lat,
                    group = group,
                    fill  = yt[, 7]
               ),
               colour = "black"
          ) +
          scale_fill_manual(
               values = c("white", "red"),
               guide  = "none"
          ) +
          labs(
               title = "",
               x     = "",
               y     = ""
          ) +
          theme_minimal(base_size = 30) +
          theme(
               legend.position  = "none",
               axis.text        = element_blank(),
               axis.ticks       = element_blank(),
               panel.grid       = element_blank(),
               panel.background = element_rect(colour = "black"),
               plot.margin      = margin(0, 30, 0, 0)
          )
     
     q
     
     ggsave(
          filename = paste("flu_mapa_respuesta_semana_", tt, ".pdf", sep = ""),
          plot     = q,
          width    = 10,
          height   = 7
     )
}

# Fin --------------------------------------------------------------------------

# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Paquetes
library(maps)
library(sf)
library(spdep)
library(Matrix)
library(ggplot2)


# Datos ------------------------------------------------------------------------

flu <- read.csv("flu_data.csv", header = TRUE)    # Base completa
flu <- flu[485:496, ]                             # Primeras 12 semanas de 2013
flu <- flu[, -c(1, 2, 12)]                        # Sin fecha, Alaska ni Hawái

states.names <- colnames(flu)

data <- as.matrix(flu)
p    <- nrow(data)                                # T: número de semanas
m    <- ncol(data)                                # I: número de estados


# Regiones del Departamento de Salud y Servicios Humanos de Estados Unidos ------

regions <- read.csv("flu_regions.csv", header = TRUE)     # Base completa
regions <- regions[-c(2, 12), ]                           # Sin Alaska ni Hawái
regions <- as.matrix(regions)

e <- length(unique(regions))                               # Número de regiones


# Distribución por raza --------------------------------------------------------

whites <- read.csv("flu_race.csv", header = TRUE)    # Base completa
whites <- whites[, 2]                                # Porcentaje de población blanca
whites <- whites[-c(2, 12)]                          # Sin Alaska ni Hawái
whites <- as.matrix(whites)


# Distribución por edad --------------------------------------------------------

old <- read.csv("flu_age.csv", header = TRUE)    # Base completa
old <- old[, 4]                                  # Porcentaje de población mayor de 65 años
old <- old[-c(2, 12)]                            # Sin Alaska ni Hawái
old <- as.matrix(old)


# Respuesta binaria ------------------------------------------------------------

y <- matrix(0, p, m)

for (i in 1:p) {
     for (j in 1:m) {
          if (data[i, j] > 7500) {
               y[i, j] <- 1
          }
     }
}

y <- as.matrix(c(y))


# Matriz de vecindades estados -------------------------------------------------

# Mapa base de estados
usa.state <- maps::map(
     database = "state",
     fill     = TRUE,
     plot     = FALSE
)

# Convertir a objeto sf
usa.sf <- sf::st_as_sf(usa.state)

# Identificador de estado
usa.sf$region <- usa.sf$ID

# Corregir geometrías
usa.sf <- sf::st_make_valid(usa.sf)

# Desactivar temporalmente el motor s2
s2.original <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)

# Vecindad entre estados
usa.nb <- spdep::poly2nb(
     usa.sf,
     row.names = usa.sf$region,
     queen     = TRUE
)

# Matriz de adyacencia
usa.adj.mat <- A <- spdep::nb2mat(
     usa.nb,
     style       = "B",
     zero.policy = TRUE
)

rownames(A) <- usa.sf$region
colnames(A) <- usa.sf$region

# Restaurar la configuración original
sf::sf_use_s2(s2.original)


# Matriz de vecindades regiones -----------------------------------------------

A2 <- matrix(
     c(
          0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
          0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
          0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
          0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
          0, 0, 0, 1, 0, 0, 1, 1, 1, 0,
          0, 0, 0, 1, 1, 1, 0, 1, 0, 0,
          0, 0, 0, 0, 1, 1, 1, 0, 1, 1,
          0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
          0, 0, 0, 0, 0, 0, 0, 1, 1, 0
     ),
     nrow  = e,
     ncol  = e,
     byrow = TRUE
)


# Matriz de precisión W para estados ------------------------------------------

W <- -A
diag(W) <- rowSums(A)

r     <- rankMatrix(W)[1]
n.nei <- rowSums(A)


# Matriz de precisión W2 para regiones ----------------------------------------

W2 <- -A2
diag(W2) <- rowSums(A2)

r2     <- rankMatrix(W2)[1]
n.nei2 <- rowSums(A2)


# Matriz de diseño X para efectos fijos ----------------------------------------

k <- 4

X <- matrix(0, m * p, k)

X[, 1] <- 1                                           # Intercepto
X[, 2] <- whites %x% as.matrix(rep(1, p))             # Porcentaje de población blanca
X[, 3] <- old %x% as.matrix(rep(1, p))                # Porcentaje de población mayor
X[, 4] <- as.matrix(rep(1, m)) %x% as.matrix(1:p)     # Semana

M <- dim(X)[1]
k <- dim(X)[2]


# Matriz de diseño Z para efectos espaciales -----------------------------------

k2 <- e

Z <- matrix(0, m * p, k2)

for (i in 1:m) {
     Z[((i - 1) * p + 1):(i * p), regions[i]] <- 1
}

rm(i)

M  <- dim(Z)[1]
k2 <- dim(Z)[2]

head(Z, 2 * p)


# Muestreador de Gibbs ---------------------------------------------------------

# Nota sobre la implementación
# Los datos se almacenan por estado, dejando que las semanas varíen dentro de cada
# estado; por esta razón, las precisiones temporales se organizan como rep(ka.w, m).
# Además, aunque xi se presenta teóricamente como un parámetro separado, en el
# código se incorpora como la cuarta columna de X y se actualiza junto con beta.

sample.w <- function(X, Z, beta, gama, theta, phi, ka.w, y)
{
     a <- b <- rep(0, M)
     
     a[y == 0] <- -Inf
     b[y == 1] <- Inf
     
     mu <- X %*% beta + Z %*% gama + theta + phi
     w  <- matrix(NA, M, 1)
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          w[idx] <- truncnorm::rtruncnorm(
               length(idx),
               a[idx],
               b[idx],
               mu[idx],
               sqrt(1 / ka.w[tt])
          )
     }
     
     as.numeric(w)
}


sample.beta <- function(X, Z, sig2.0, w, gama, theta, phi, ka.w)
{
     V.inv <- diag(rep(ka.w, m))
     
     SIGMA <- chol2inv(
          chol(
               diag(k) / sig2.0 +
                    t(X) %*% V.inv %*% X
          )
     )
     
     MEAN <- SIGMA %*% t(X) %*% V.inv %*% (
          w - Z %*% gama - theta - phi
     )
     
     as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
}


sample.gama <- function(X, Z, W2, sig2.0, w, beta, theta, phi, ka.w)
{
     V.inv <- diag(rep(ka.w, m))
     
     SIGMA <- chol2inv(
          chol(
               W2 / sig2.0 +
                    t(Z) %*% V.inv %*% Z
          )
     )
     
     MEAN <- SIGMA %*% t(Z) %*% V.inv %*% (
          w - X %*% beta - theta - phi
     )
     
     gama <- as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
     gama <- gama - mean(gama)
     
     gama
}


sample.theta <- function(X, Z, beta, gama, w, phi, ka.h, ka.w)
{
     theta <- matrix(NA, M, 1)
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          SIGMA <- diag(m) / (ka.h[tt] + ka.w[tt])
          
          MEAN <- SIGMA %*% (
               ka.w[tt] * (
                    w[idx] -
                         X[idx, ] %*% beta -
                         Z[idx, ] %*% gama -
                         phi[idx]
               )
          )
          
          theta[idx] <- as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
          # theta[idx] <- theta[idx] - mean(theta[idx])
     }
     
     as.numeric(theta)
}


sample.phi <- function(X, Z, W, beta, gama, w, theta, ka.c, ka.w)
{
     phi <- matrix(NA, M, 1)
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          SIGMA <- chol2inv(
               chol(
                    ka.c[tt] * W +
                         ka.w[tt] * diag(m)
               )
          )
          
          MEAN <- SIGMA %*% (
               ka.w[tt] * (
                    w[idx] -
                         X[idx, ] %*% beta -
                         Z[idx, ] %*% gama -
                         theta[idx]
               )
          )
          
          phi[idx] <- as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
          phi[idx] <- phi[idx] - mean(phi[idx])
     }
     
     as.numeric(phi)
}


sample.ka.w <- function(nu.0, X, Z, beta, gama, w, theta, phi)
{
     ka.w  <- matrix(NA, p, 1)
     shape <- nu.0 / 2 + m / 2
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          rate <- nu.0 / 2 +
               0.5 * sum(
                    (
                         w[idx] -
                              X[idx, ] %*% beta -
                              Z[idx, ] %*% gama -
                              theta[idx] -
                              phi[idx]
                    )^2
               )
          
          ka.w[tt] <- rgamma(1, shape, rate)
     }
     
     as.numeric(ka.w)
}


sample.ka.h <- function(a.h, b.h, theta)
{
     ka.h  <- matrix(NA, p, 1)
     shape <- a.h + m / 2
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          rate <- b.h + 0.5 * sum(theta[idx]^2)
          
          ka.h[tt] <- rgamma(1, shape, rate)
     }
     
     as.numeric(ka.h)
}


sample.ka.c <- function(a.c, b.c, W, phi)
{
     ka.c  <- matrix(NA, p, 1)
     shape <- a.c + r / 2
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          rate <- b.c + 0.5 * as.numeric(
               t(phi[idx]) %*% W %*% phi[idx]
          )
          
          ka.c[tt] <- rgamma(1, shape, rate)
     }
     
     as.numeric(ka.c)
}

MCMC <- function(y, X, Z, W, W2, sig2.0, nu.0, a.h, b.h, a.c, b.c,
                 n.sams, n.burn, n.skip)
{
     require(Matrix)
     library(mvtnorm)
     
     # Dimensiones
     m  <- dim(W)[1]          # Número de estados
     M  <- dim(X)[1]          # Número de observaciones
     k  <- dim(X)[2]          # Número de covariables fijas
     k2 <- dim(Z)[2]          # Número de covariables espaciales
     p  <- M / m              # Número de tiempos
     r  <- rankMatrix(W)[1]   # Rango de W, estados
     r2 <- rankMatrix(W2)[1]  # Rango de W2, regiones
     
     # Número total de iteraciones
     n.total <- n.burn + n.skip * n.sams
     
     # Almacenamiento de muestras posteriores
     beta.chain   <- matrix(NA, n.sams, k)
     gama.chain   <- matrix(NA, n.sams, k2)
     w.chain      <- matrix(NA, n.sams, M)
     theta.chain  <- matrix(NA, n.sams, M)
     phi.chain    <- matrix(NA, n.sams, M)
     ka.w.chain   <- matrix(NA, n.sams, p)
     ka.h.chain   <- matrix(NA, n.sams, p)
     ka.c.chain   <- matrix(NA, n.sams, p)
     loglik.chain <- rep(NA, n.sams)
     
     # Inicialización de parámetros
     beta  <- rep(0, k)
     gama  <- rep(0, k2)
     w     <- rep(0, M)
     theta <- rep(0, M)
     phi   <- rep(0, M)
     ka.w  <- rep(1, p)
     ka.h  <- rep(1, p)
     ka.c  <- rep(1, p)
     
     # Puntos de progreso
     progreso <- unique(round(seq(0.10, 1, by = 0.10) * n.total))
     
     # MCMC
     s <- 0
     
     for (l in 1:n.total) {
          
          # Muestreo
          w     <- sample.w(X, Z, beta, gama, theta, phi, ka.w, y)
          beta  <- sample.beta(X, Z, sig2.0, w, gama, theta, phi, ka.w)
          gama  <- sample.gama(X, Z, W2, sig2.0, w, beta, theta, phi, ka.w)
          theta <- sample.theta(X, Z, beta, gama, w, phi, ka.h, ka.w)
          phi   <- sample.phi(X, Z, W, beta, gama, w, theta, ka.c, ka.w)
          ka.w  <- sample.ka.w(nu.0, X, Z, beta, gama, w, theta, phi)
          ka.h  <- sample.ka.h(a.h, b.h, theta)
          ka.c  <- sample.ka.c(a.c, b.c, W, phi)
          
          # Almacenamiento posterior
          if (l > n.burn && (l - n.burn) %% n.skip == 0) {
               s <- s + 1
               
               beta.chain[s, ]  <- beta
               gama.chain[s, ]  <- gama
               w.chain[s, ]     <- w
               theta.chain[s, ] <- theta
               phi.chain[s, ]   <- phi
               ka.w.chain[s, ]  <- ka.w
               ka.h.chain[s, ]  <- ka.h
               ka.c.chain[s, ]  <- ka.c
               
               eta <- X %*% beta + Z %*% gama + theta + phi
               pi  <- pnorm(sqrt(rep(ka.w, m)) * eta)
               
               loglik.chain[s] <- sum(dbinom(y, size = 1, prob = pi, log = TRUE))
          }
          
          # Progreso
          if (l %in% progreso || l == n.total) {
               cat(
                    "Progreso: ",
                    round(100 * l / n.total, 1),
                    "%\n",
                    sep = ""
               )
          }
     }
     
     # Salida
     list(
          beta.chain   = beta.chain,
          gama.chain   = gama.chain,
          w.chain      = w.chain,
          theta.chain  = theta.chain,
          phi.chain    = phi.chain,
          ka.w.chain   = ka.w.chain,
          ka.h.chain   = ka.h.chain,
          ka.c.chain   = ka.c.chain,
          loglik.chain = loglik.chain,
          n.eff        = n.sams
     )
}


# Hiperparámetros --------------------------------------------------------------

# Especificación previa
sig2.0 <- 100  # psi0 = 1/sig2.0
nu.0   <- 3

# Escala previa objetivo para los efectos aleatorios en la escala latente
s0 <- 2

# Precisión de theta_it
a.h <- 2
b.h <- (a.h - 1) * s0^2

# Precisión de phi_it
a.c <- 2
b.c <- (0.7^2) * mean(rowSums(A)) * b.h  # BCG, p. 156


# Ajuste del modelo ------------------------------------------------------------

n.sams <- 10000
n.burn <- 10000
n.skip <- 10

# Ajuste del modelo

set.seed(123)

fit <- MCMC(
     y, X, Z, W, W2,
     sig2.0, nu.0,
     a.h, b.h,
     a.c, b.c,
     n.sams, n.burn, n.skip
)

save(fit, file = "flu_muestras_distribucion_posterior.RData")

load(file = "flu_muestras_distribucion_posterior.RData")


# Traza de la log-verosimilitud ------------------------------------------------

plot(
     x    = fit$loglik.chain,
     type = "p",
     pch  = 16,
     cex  = 0.3,
     xlab = "Iteración",
     ylab = "Log-verosimilitud",
     main = "Log-verosimilitud"
)

dev.off()


# Inferencia posterior de los efectos fijos y xi -------------------------------

resumen_posterior <- function(x)
{
     c(
          media = mean(x),
          sd    = sd(x),
          li95  = unname(quantile(x, 0.025)),
          ls95  = unname(quantile(x, 0.975))
     )
}

tabla.efectos <- t(
     apply(
          fit$beta.chain,
          2,
          resumen_posterior
     )
)

rownames(tabla.efectos) <- c(
     "beta_1",
     "beta_2",
     "beta_3",
     "xi"
)

tabla.efectos <- round(tabla.efectos, 3)

tabla.efectos


# Inferencia posterior gamma ---------------------------------------------------

# Resumen posterior de los efectos aleatorios
gama.media <- apply(fit$gama.chain, 2, mean)
gama.li95  <- apply(fit$gama.chain, 2, quantile, probs = 0.025)
gama.ls95  <- apply(fit$gama.chain, 2, quantile, probs = 0.975)

# Colores según posición del intervalo
col.intervalo <- rep("black", length(gama.media))
col.intervalo[gama.li95 > 0] <- "darkgreen"
col.intervalo[gama.ls95 < 0] <- "red"

pdf(
     file      = "flu_inferencia_efectos_aleatorios.pdf",
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
     x    = 1:length(gama.media),
     y    = gama.media,
     type = "p",
     pch  = 16,
     col  = col.intervalo,
     ylim = range(c(gama.li95, gama.ls95)),
     xaxt = "n",
     xlab = "Efecto aleatorio",
     ylab = "Media posterior",
     main = ""
)

axis(
     side   = 1,
     at     = 1:length(gama.media),
     labels = 1:length(gama.media),
     las    = 1
)

arrows(
     x0     = 1:length(gama.media),
     y0     = gama.li95,
     x1     = 1:length(gama.media),
     y1     = gama.ls95,
     angle  = 90,
     code   = 3,
     length = 0,
     col    = col.intervalo
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2
)

box()

dev.off()


# Inferencia posterior efectos espacio-temporales ------------------------------

# Efecto espacio-temporal total
efecto.et <- fit$theta.chain + fit$phi.chain

# Media posterior por estado y semana
media.et <- matrix(NA, nrow = m, ncol = p)

for (i in 1:m) {
     idx <- ((i - 1) * p + 1):(i * p)
     media.et[i, ] <- apply(efecto.et[, idx], 2, mean)
}

# Región de cada estado
region.estado <- as.numeric(
     Z[seq(1, m * p, by = p), ] %*% (1:ncol(Z))
)

# Regiones
regiones <- sort(unique(region.estado))

# Paneles por región
for (rr in regiones) {
     
     idx.region <- which(region.estado == rr)
     
     pdf(
          file      = paste0("flu_inferencia_efectos_espaciotemporales_region_", rr, ".pdf"),
          width     = 7,
          height    = 5,
          pointsize = 25
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = 1:p,
          y    = media.et[idx.region[1], ],
          type = "n",
          ylim = range(media.et),
          xlab = "Semana",
          ylab = "Efecto",
          main = ""
     )
     
     for (i in idx.region) {
          lines(
               x   = 1:p,
               y   = media.et[i, ],
               col = 2,
               lwd = 1.5
          )
     }
     
     abline(
          h   = 0,
          lty = 2,
          lwd = 1.2
     )
     
     box()
     
     dev.off()
}


# Inferencia posterior componentes de precisión --------------------------------

# Medias posteriores de los componentes de precisión
ka.w.media <- apply(fit$ka.w.chain, 2, mean)
ka.h.media <- apply(fit$ka.h.chain, 2, mean)
ka.c.media <- apply(fit$ka.c.chain, 2, mean)

# Archivo de salida
pdf(
     file      = "flu_inferencia_componentes_precision.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Gráfico
plot(
     x    = 1:p,
     y    = ka.w.media,
     type = "n",
     ylim = range(c(ka.w.media, ka.h.media, ka.c.media)),
     xaxt = "n",
     xlab = "Semana",
     ylab = "Precisión",
     main = ""
)

axis(
     side   = 1,
     at     = 1:p,
     labels = 1:p
)

lines(
     x   = 1:p,
     y   = ka.w.media,
     col = "black",
     lwd = 2
)

lines(
     x   = 1:p,
     y   = ka.h.media,
     col = "red",
     lwd = 2
)

lines(
     x   = 1:p,
     y   = ka.c.media,
     col = "blue",
     lwd = 2
)

legend(
     "topleft",
     legend = c(
          expression(kappa[t]),
          expression(tau[t]),
          expression(lambda[t])
     ),
     col = c("black", "red", "blue"),
     lty = 1,
     lwd = 2,
     bty = "n"
)

box()

dev.off()


# Predicción -------------------------------------------------------------------

post.pred <- function(i, tt, fit)
{
     n.sam    <- fit$n.eff
     beta.sam <- fit$beta.chain
     gama.sam <- fit$gama.chain
     w.it.sam <- rep(NA, n.sam)
     
     x.it    <- rep(0, k)
     x.it[1] <- 1
     x.it[2] <- whites[i]  # Porcentaje de población blanca
     x.it[3] <- old[i]     # Porcentaje de población mayor
     x.it[4] <- tt         # Semana
     
     z.it <- rep(0, k2)
     z.it[regions[i]] <- 1
     
     idx <- (i - 1) * p + tt
     
     for (b in 1:n.sam) {
          kw.t     <- fit$ka.w.chain[b, tt]
          theta.it <- fit$theta.chain[b, idx]
          phi.it   <- fit$phi.chain[b, idx]
          
          mu <- sum(x.it * beta.sam[b, ]) +
               sum(z.it * gama.sam[b, ]) +
               theta.it +
               phi.it
          
          w.it.sam[b] <- rnorm(1, mu, sqrt(1 / kw.t))
     }
     
     y.it <- rep(0, n.sam)
     y.it[w.it.sam > 0] <- 1
     
     list(pr = mean(y.it), sam = y.it)
}

# Mapa base de estados
all_states <- map_data("state")

# Mapas de probabilidad predictiva posterior
for (tt in 1:p) {
     
     post.prob <- rep(NA, m)
     
     for (i in 1:m) {
          post.prob[i] <- post.pred(i, tt, fit)$pr
     }
     
     # Unir probabilidad predictiva con mapa
     post.prob <- as.data.frame(post.prob)
     post.prob <- cbind(rownames(A), post.prob)
     
     colnames(post.prob) <- c("region", "prob")
     
     post.prob <- merge(
          all_states,
          post.prob,
          by = "region"
     )
     
     q <- ggplot() +
          geom_polygon(
               data = post.prob,
               aes(
                    x     = long,
                    y     = lat,
                    group = group,
                    fill  = prob
               ),
               colour = "black"
          ) +
          scale_fill_continuous(
               low    = "white",
               high   = "red",
               limits = c(0, 1),
               guide  = "none"
          ) +
          labs(
               title = "",
               x     = "",
               y     = ""
          ) +
          theme_minimal(base_size = 30) +
          theme(
               legend.position  = "none",
               axis.text        = element_blank(),
               axis.ticks       = element_blank(),
               panel.grid       = element_blank(),
               panel.background = element_rect(colour = "black"),
               plot.margin      = margin(0, 30, 0, 0)
          )
     
     q
     
     ggsave(
          filename = paste("flu_inferencia_prediccion_semana_", tt, ".pdf", sep = ""),
          plot     = q,
          width    = 10,
          height   = 7
     )
     
     rm(post.prob)
}

# Fin --------------------------------------------------------------------------