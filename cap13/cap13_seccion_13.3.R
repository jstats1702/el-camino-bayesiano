# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(MASS)))

# Datos ------------------------------------------------------------------------

data("galaxies", package = "MASS")

dat <- data.frame(
     id       = seq_along(galaxies),
     velocity = as.numeric(galaxies) / 1000
)

# Ordenar únicamente para visualización
dat <- dat[
     order(dat$velocity),
     ,
     drop = FALSE
]

rownames(dat) <- NULL

y <- dat$velocity
n <- length(y)

# Análisis exploratorio --------------------------------------------------------

# Resumen numérico
descriptive_summary <- c(
     n       = n,
     minimum = min(y),
     q025    = quantile(y, 0.025, names = FALSE),
     q25     = quantile(y, 0.250, names = FALSE),
     median  = median(y),
     mean    = mean(y),
     q75     = quantile(y, 0.750, names = FALSE),
     q975    = quantile(y, 0.975, names = FALSE),
     maximum = max(y),
     sd      = sd(y),
     iqr     = IQR(y)
)

print(
     round(
          descriptive_summary,
          digits = 3
     )
)

# Histograma y estimación mediante núcleos
pdf(
     file = "dpm_galaxias_exploratorio_histograma.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

histogram <- hist(
     y,
     breaks = 25,
     probability = TRUE,
     col = "gray90",
     border = "white",
     main = "",
     xlab = "Velocidad",
     ylab = "Densidad"
)

rug(y)

dev.off()

# Velocidades ordenadas
pdf(
     file = "dpm_galaxias_exploratorio_orden.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     seq_len(n),
     sort(y),
     pch = 16,
     cex = 0.75,
     xlab = "Índice de orden",
     ylab = "Velocidad",
     main = ""
)

lines(
     seq_len(n),
     sort(y),
     col = "gray65"
)

dev.off()

# Funciones auxiliares ---------------------------------------------------------

# Reetiqueta los componentes con enteros consecutivos.
relabel_partition <- function(xi) {
     as.integer(
          match(
               xi,
               unique(xi)
          )
     )
}

# Construye una partición inicial a partir de los datos ordenados.
initialize_partition <- function(y, K_init = 8L) {
     n <- length(y)
     K <- min(as.integer(K_init), n)
     
     if (K < 1L) {
          stop("K_init debe ser un entero positivo.")
     }
     
     if (K == 1L) {
          return(rep.int(1L, n))
     }
     
     order_y <- order(y)
     
     groups <- as.integer(
          floor(
               (seq_len(n) - 1L) * K / n
          ) + 1L
     )
     
     xi <- integer(n)
     xi[order_y] <- groups
     
     xi
}

# Calcula tamaños, sumas y sumas de cuadrados por componente.
partition_statistics <- function(y, xi) {
     xi <- relabel_partition(xi)
     K <- max(xi)
     
     counts <- tabulate(
          xi,
          nbins = K
     )
     
     sums <- vapply(
          X = seq_len(K),
          FUN = function(ell) {
               sum(y[xi == ell])
          },
          FUN.VALUE = numeric(1)
     )
     
     sums_sq <- vapply(
          X = seq_len(K),
          FUN = function(ell) {
               sum(y[xi == ell]^2)
          },
          FUN.VALUE = numeric(1)
     )
     
     list(
          n      = as.integer(counts),
          sum    = sums,
          sum_sq = sums_sq
     )
}

# Calcula los parámetros posteriores Normal--Gamma Inversa.
nig_posterior_stats <- function(
          n_l,
          sum_l,
          sum_sq_l,
          mu,
          kappa,
          beta,
          cc
) {
     if (n_l == 0L) {
          return(
               list(
                    kappa_tilde = kappa,
                    mu_tilde    = mu,
                    c_tilde     = cc,
                    beta_tilde  = beta
               )
          )
     }
     
     mean_l <- sum_l / n_l
     
     sse_l <- max(sum_sq_l - sum_l^2 / n_l, 0)
     
     kappa_tilde <- kappa + n_l
     
     mu_tilde <- (kappa * mu + sum_l) / kappa_tilde
     
     c_tilde <- cc + n_l / 2
     
     beta_tilde <- beta +
          sse_l / 2 +
          (kappa * n_l / (2 * kappa_tilde)) * (mean_l - mu)^2
     
     list(
          kappa_tilde = kappa_tilde,
          mu_tilde    = mu_tilde,
          c_tilde     = c_tilde,
          beta_tilde  = beta_tilde
     )
}

# Evalúa la log-densidad predictiva t de un componente.
log_predictive_nig_stats <- function(
          y_new,
          n_l,
          sum_l,
          sum_sq_l,
          mu,
          kappa,
          beta,
          cc
) {
     posterior <- nig_posterior_stats(
          n_l      = n_l,
          sum_l    = sum_l,
          sum_sq_l = sum_sq_l,
          mu       = mu,
          kappa    = kappa,
          beta     = beta,
          cc       = cc
     )
     
     df <- 2 * posterior$c_tilde
     
     scale2 <- (
          posterior$beta_tilde *
               (posterior$kappa_tilde + 1)
     ) / (
          posterior$c_tilde *
               posterior$kappa_tilde
     )
     
     dt(
          x = (y_new - posterior$mu_tilde) / sqrt(scale2),
          df = df,
          log = TRUE
     ) -
          0.5 * log(scale2)
}

# Evalúa la log-verosimilitud marginal de un componente.
log_marginal_nig_stats <- function(
          n_l,
          sum_l,
          sum_sq_l,
          mu,
          kappa,
          beta,
          cc
) {
     if (n_l == 0L) {
          return(0)
     }
     
     posterior <- nig_posterior_stats(
          n_l      = n_l,
          sum_l    = sum_l,
          sum_sq_l = sum_sq_l,
          mu       = mu,
          kappa    = kappa,
          beta     = beta,
          cc       = cc
     )
     
     -0.5 * n_l * log(2 * pi) +
          0.5 * (
               log(kappa) -
                    log(posterior$kappa_tilde)
          ) +
          cc * log(beta) -
          posterior$c_tilde * log(posterior$beta_tilde) +
          lgamma(posterior$c_tilde) -
          lgamma(cc)
}

# Suma las log-verosimilitudes marginales de la partición.
partition_log_marginal <- function(
          statistics,
          mu,
          kappa,
          beta,
          cc
) {
     K <- length(statistics$n)
     
     if (K == 0L) {
          return(0)
     }
     
     sum(
          vapply(
               X = seq_len(K),
               FUN = function(ell) {
                    log_marginal_nig_stats(
                         n_l      = statistics$n[ell],
                         sum_l    = statistics$sum[ell],
                         sum_sq_l = statistics$sum_sq[ell],
                         mu       = mu,
                         kappa    = kappa,
                         beta     = beta,
                         cc       = cc
                    )
               },
               FUN.VALUE = numeric(1)
          )
     )
}

# Muestrea una categoría a partir de log-pesos no normalizados.
sample_log_weights <- function(log_weights) {
     infinite_weights <- which(log_weights == Inf)
     
     if (length(infinite_weights) > 0L) {
          return(
               sample(
                    infinite_weights,
                    size = 1L
               )
          )
     }
     
     maximum <- max(log_weights)
     
     weights <- exp(
          log_weights -
               maximum
     )
     
     weights[!is.finite(weights)] <- 0
     
     total <- sum(weights)
     
     sample.int(
          n = length(weights),
          size = 1L,
          prob = weights
     )
}

# Actualiza secuencialmente las asignaciones del Gibbs colapsado.
update_assignments <- function(
          y,
          xi,
          statistics,
          alpha,
          mu,
          kappa,
          beta,
          cc,
          random_scan = TRUE
) {
     n <- length(y)
     
     counts  <- statistics$n
     sums    <- statistics$sum
     sums_sq <- statistics$sum_sq
     
     scan_order <- if (random_scan) {
          sample.int(n)
     } else {
          seq_len(n)
     }
     
     for (i in scan_order) {
          old_label <- xi[i]
          
          # Retirar la observación
          counts[old_label] <- counts[old_label] - 1L
          sums[old_label] <- sums[old_label] - y[i]
          sums_sq[old_label] <- sums_sq[old_label] - y[i]^2
          
          # Eliminar componentes vacíos
          if (counts[old_label] == 0L) {
               counts <- counts[-old_label]
               sums <- sums[-old_label]
               sums_sq <- sums_sq[-old_label]
               
               xi[xi > old_label] <- xi[xi > old_label] - 1L
          }
          
          K_minus <- length(counts)
          
          log_weights <- numeric(K_minus + 1L)
          
          # Componentes ocupados
          if (K_minus > 0L) {
               for (ell in seq_len(K_minus)) {
                    log_weights[ell] <- log(counts[ell]) +
                         log_predictive_nig_stats(
                              y_new     = y[i],
                              n_l       = counts[ell],
                              sum_l     = sums[ell],
                              sum_sq_l  = sums_sq[ell],
                              mu        = mu,
                              kappa     = kappa,
                              beta      = beta,
                              cc        = cc
                         )
               }
          }
          
          # Componente nuevo
          log_weights[K_minus + 1L] <- log(alpha) +
               log_predictive_nig_stats(
                    y_new     = y[i],
                    n_l       = 0L,
                    sum_l     = 0,
                    sum_sq_l  = 0,
                    mu        = mu,
                    kappa     = kappa,
                    beta      = beta,
                    cc        = cc
               )
          
          new_label <- sample_log_weights(log_weights)
          
          xi[i] <- new_label
          
          # Incorporar la observación
          if (new_label == K_minus + 1L) {
               counts <- c(
                    counts,
                    1L
               )
               
               sums <- c(
                    sums,
                    y[i]
               )
               
               sums_sq <- c(
                    sums_sq,
                    y[i]^2
               )
          } else {
               counts[new_label] <- counts[new_label] + 1L
               sums[new_label] <- sums[new_label] + y[i]
               sums_sq[new_label] <- sums_sq[new_label] + y[i]^2
          }
     }
     
     list(
          xi = xi,
          statistics = list(
               n      = as.integer(counts),
               sum    = sums,
               sum_sq = sums_sq
          )
     )
}

# Actualiza alpha mediante el esquema de Escobar--West.
update_alpha_escobar_west <- function(
          alpha,
          K,
          n,
          a_alpha,
          b_alpha
) {
     auxiliary <- rbeta(
          n = 1L,
          shape1 = alpha + 1,
          shape2 = n
     )
     
     auxiliary <- max(
          auxiliary,
          .Machine$double.xmin
     )
     
     rate <- b_alpha -
          log(auxiliary)
     
     numerator <- a_alpha +
          K -
          1
     
     mixture_probability <- numerator / (
          numerator +
               n * rate
     )
     
     shape <- if (
          runif(1) <
          mixture_probability
     ) {
          a_alpha + K
     } else {
          a_alpha + K - 1
     }
     
     rgamma(
          n = 1L,
          shape = shape,
          rate = rate
     )
}

# Evalúa la log-posterior de mu, kappa y beta.
log_hyperparameter_posterior <- function(
          mu,
          kappa,
          beta,
          statistics,
          hyperparameters
) {
     log_prior <- dnorm(
          x = mu,
          mean = hyperparameters$a_mu,
          sd = sqrt(hyperparameters$b_mu),
          log = TRUE
     ) +
          dgamma(
               x = kappa,
               shape = hyperparameters$a_kappa,
               rate = hyperparameters$b_kappa,
               log = TRUE
          ) +
          dgamma(
               x = beta,
               shape = hyperparameters$a_beta,
               rate = hyperparameters$b_beta,
               log = TRUE
          )
     
     log_prior +
          partition_log_marginal(
               statistics = statistics,
               mu = mu,
               kappa = kappa,
               beta = beta,
               cc = hyperparameters$cc
          )
}

# Decide la aceptación de una propuesta de Metropolis--Hastings.
accept_metropolis <- function(log_ratio) {
     if (is.na(log_ratio)) {
          return(FALSE)
     }
     
     if (log_ratio >= 0) {
          return(TRUE)
     }
     
     log(runif(1)) < log_ratio
}

# Actualiza mu, kappa y beta mediante Metropolis--Hastings.
update_hyperparameters <- function(
          mu,
          kappa,
          beta,
          statistics,
          hyperparameters,
          proposal_sd
) {
     accepted <- c(
          mu        = FALSE,
          log_kappa = FALSE,
          log_beta  = FALSE
     )
     
     current_target <- log_hyperparameter_posterior(
          mu = mu,
          kappa = kappa,
          beta = beta,
          statistics = statistics,
          hyperparameters = hyperparameters
     )
     
     # Actualizar mu
     mu_proposal <- rnorm(
          n = 1L,
          mean = mu,
          sd = proposal_sd["mu"]
     )
     
     proposal_target <- log_hyperparameter_posterior(
          mu = mu_proposal,
          kappa = kappa,
          beta = beta,
          statistics = statistics,
          hyperparameters = hyperparameters
     )
     
     if (
          accept_metropolis(
               proposal_target -
               current_target
          )
     ) {
          mu <- mu_proposal
          current_target <- proposal_target
          accepted["mu"] <- TRUE
     }
     
     # Actualizar kappa en escala logarítmica
     log_kappa <- log(kappa)
     
     log_kappa_proposal <- rnorm(
          n = 1L,
          mean = log_kappa,
          sd = proposal_sd["log_kappa"]
     )
     
     kappa_proposal <- exp(log_kappa_proposal)
     
     if (
          is.finite(kappa_proposal) &&
          kappa_proposal > 0
     ) {
          proposal_target <- log_hyperparameter_posterior(
               mu = mu,
               kappa = kappa_proposal,
               beta = beta,
               statistics = statistics,
               hyperparameters = hyperparameters
          )
          
          log_ratio <- proposal_target -
               current_target +
               log_kappa_proposal -
               log_kappa
          
          if (accept_metropolis(log_ratio)) {
               kappa <- kappa_proposal
               current_target <- proposal_target
               accepted["log_kappa"] <- TRUE
          }
     }
     
     # Actualizar beta en escala logarítmica
     log_beta <- log(beta)
     
     log_beta_proposal <- rnorm(
          n = 1L,
          mean = log_beta,
          sd = proposal_sd["log_beta"]
     )
     
     beta_proposal <- exp(log_beta_proposal)
     
     if (
          is.finite(beta_proposal) &&
          beta_proposal > 0
     ) {
          proposal_target <- log_hyperparameter_posterior(
               mu = mu,
               kappa = kappa,
               beta = beta_proposal,
               statistics = statistics,
               hyperparameters = hyperparameters
          )
          
          log_ratio <- proposal_target -
               current_target +
               log_beta_proposal -
               log_beta
          
          if (accept_metropolis(log_ratio)) {
               beta <- beta_proposal
               current_target <- proposal_target
               accepted["log_beta"] <- TRUE
          }
     }
     
     list(
          mu         = mu,
          kappa      = kappa,
          beta       = beta,
          accepted   = accepted,
          log_target = current_target
     )
}

# Genera los parámetros posteriores de los componentes ocupados.
sample_component_parameters <- function(
          y,
          xi,
          mu,
          kappa,
          beta,
          cc
) {
     xi <- relabel_partition(xi)
     
     statistics <- partition_statistics(
          y = y,
          xi = xi
     )
     
     K <- length(statistics$n)
     
     mu_star <- numeric(K)
     sigma2_star <- numeric(K)
     
     for (ell in seq_len(K)) {
          posterior <- nig_posterior_stats(
               n_l = statistics$n[ell],
               sum_l = statistics$sum[ell],
               sum_sq_l = statistics$sum_sq[ell],
               mu = mu,
               kappa = kappa,
               beta = beta,
               cc = cc
          )
          
          sigma2_star[ell] <- 1 / rgamma(
               n = 1L,
               shape = posterior$c_tilde,
               rate = posterior$beta_tilde
          )
          
          mu_star[ell] <- rnorm(
               n = 1L,
               mean = posterior$mu_tilde,
               sd = sqrt(
                    sigma2_star[ell] /
                         posterior$kappa_tilde
               )
          )
     }
     
     data.frame(
          component   = seq_len(K),
          n           = statistics$n,
          mu_star     = mu_star,
          sigma2_star = sigma2_star
     )
}

# Muestreador completamente colapsado ------------------------------------------

# Ajusta una mezcla Normal de localización y escala mediante un Gibbs completamente colapsado.
dpm_ls_collapsed <- function(
          y,
          B = 110000L,
          burn = 10000L,
          thin = 10L,
          K_init = 8L,
          hyperparameters,
          initial = NULL,
          proposal_sd = NULL,
          adapt = TRUE,
          adapt_interval = 50L,
          target_accept = 0.44,
          random_scan = TRUE,
          seed = NULL,
          verbose = TRUE,
          progress_every = 10000L
) {
     y <- as.numeric(y)
     n <- length(y)
     
     B <- as.integer(B)
     burn <- as.integer(burn)
     thin <- as.integer(thin)
     adapt_interval <- as.integer(adapt_interval)
     progress_every <- as.integer(progress_every)
     
     required_hyperparameters <- c(
          "a_alpha",
          "b_alpha",
          "a_mu",
          "b_mu",
          "a_kappa",
          "b_kappa",
          "a_beta",
          "b_beta",
          "cc"
     )
     
     missing_hyperparameters <- setdiff(
          required_hyperparameters,
          names(hyperparameters)
     )
     
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     if (is.null(initial)) {
          initial <- list()
     }
     
     # Inicializar la partición y los hiperparámetros
     xi <- if (is.null(initial$xi)) {
          initialize_partition(
               y = y,
               K_init = K_init
          )
     } else {
          relabel_partition(
               as.integer(initial$xi)
          )
     }
     
     alpha <- if (is.null(initial$alpha)) {
          hyperparameters$a_alpha / hyperparameters$b_alpha
     } else {
          as.numeric(initial$alpha)
     }
     
     mu <- if (is.null(initial$mu)) {
          hyperparameters$a_mu
     } else {
          as.numeric(initial$mu)
     }
     
     kappa <- if (is.null(initial$kappa)) {
          hyperparameters$a_kappa / hyperparameters$b_kappa
     } else {
          as.numeric(initial$kappa)
     }
     
     beta <- if (is.null(initial$beta)) {
          hyperparameters$a_beta / hyperparameters$b_beta
     } else {
          as.numeric(initial$beta)
     }
     
     # Definir las escalas de las propuestas
     if (is.null(proposal_sd)) {
          mu_proposal_sd <- 0.10 * sd(y)
          
          if (!is.finite(mu_proposal_sd) || mu_proposal_sd <= 0) {
               mu_proposal_sd <- 0.10
          }
          
          proposal_sd <- c(
               mu        = mu_proposal_sd,
               log_kappa = 0.20,
               log_beta  = 0.20
          )
     }
     
     if (is.null(names(proposal_sd)) && length(proposal_sd) == 3L) {
          names(proposal_sd) <- c(
               "mu",
               "log_kappa",
               "log_beta"
          )
     }
     
     proposal_sd <- proposal_sd[
          c(
               "mu",
               "log_kappa",
               "log_beta"
          )
     ]
     
     statistics <- partition_statistics(
          y = y,
          xi = xi
     )
     
     # Preparar el almacenamiento
     keep_iterations <- seq.int(
          from = burn + thin,
          to = B,
          by = thin
     )
     
     n_keep <- length(keep_iterations)
     
     draws <- data.frame(
          iteration    = keep_iterations,
          alpha        = numeric(n_keep),
          mu           = numeric(n_keep),
          kappa        = numeric(n_keep),
          beta         = numeric(n_keep),
          K            = integer(n_keep),
          log_marginal = numeric(n_keep)
     )
     
     partitions <- matrix(
          data = NA_integer_,
          nrow = n_keep,
          ncol = n
     )
     
     colnames(partitions) <- paste0(
          "observation_",
          seq_len(n)
     )
     
     acceptance_post <- c(
          mu        = 0,
          log_kappa = 0,
          log_beta  = 0
     )
     
     batch_acceptance <- c(
          mu        = 0,
          log_kappa = 0,
          log_beta  = 0
     )
     
     batch_size <- 0L
     adaptation_index <- 0L
     keep_index <- 0L
     
     for (b in seq_len(B)) {
          # Actualizar la partición
          assignment_update <- update_assignments(
               y = y,
               xi = xi,
               statistics = statistics,
               alpha = alpha,
               mu = mu,
               kappa = kappa,
               beta = beta,
               cc = hyperparameters$cc,
               random_scan = random_scan
          )
          
          xi <- assignment_update$xi
          statistics <- assignment_update$statistics
          K <- length(statistics$n)
          
          # Actualizar el parámetro de concentración
          alpha <- update_alpha_escobar_west(
               alpha = alpha,
               K = K,
               n = n,
               a_alpha = hyperparameters$a_alpha,
               b_alpha = hyperparameters$b_alpha
          )
          
          # Actualizar los hiperparámetros de la distribución base
          hyperparameter_update <- update_hyperparameters(
               mu = mu,
               kappa = kappa,
               beta = beta,
               statistics = statistics,
               hyperparameters = hyperparameters,
               proposal_sd = proposal_sd
          )
          
          mu <- hyperparameter_update$mu
          kappa <- hyperparameter_update$kappa
          beta <- hyperparameter_update$beta
          
          # Adaptar las propuestas durante el calentamiento
          if (adapt && b <= burn) {
               batch_acceptance <- batch_acceptance +
                    hyperparameter_update$accepted
               
               batch_size <- batch_size + 1L
               
               if (batch_size == adapt_interval || b == burn) {
                    adaptation_index <- adaptation_index + 1L
                    
                    batch_rate <- batch_acceptance / batch_size
                    
                    adaptation_gain <- min(
                         0.05,
                         1 / sqrt(adaptation_index)
                    )
                    
                    proposal_sd <- exp(
                         log(proposal_sd) +
                              adaptation_gain *
                              (batch_rate - target_accept)
                    )
                    
                    proposal_sd <- pmax(proposal_sd, 1e-6)
                    
                    names(proposal_sd) <- c(
                         "mu",
                         "log_kappa",
                         "log_beta"
                    )
                    
                    batch_acceptance[] <- 0
                    batch_size <- 0L
               }
          }
          
          # Registrar las aceptaciones posteriores al calentamiento
          if (b > burn) {
               acceptance_post <- acceptance_post +
                    hyperparameter_update$accepted
          }
          
          log_marginal <- partition_log_marginal(
               statistics = statistics,
               mu = mu,
               kappa = kappa,
               beta = beta,
               cc = hyperparameters$cc
          )
          
          # Almacenar las iteraciones seleccionadas
          if (b > burn && (b - burn) %% thin == 0L) {
               keep_index <- keep_index + 1L
               
               draws$alpha[keep_index] <- alpha
               draws$mu[keep_index] <- mu
               draws$kappa[keep_index] <- kappa
               draws$beta[keep_index] <- beta
               draws$K[keep_index] <- K
               draws$log_marginal[keep_index] <- log_marginal
               
               partitions[keep_index, ] <- xi
          }
          
          # Mostrar el avance de la cadena
          if (
               verbose &&
               (
                    b == 1L ||
                    b %% progress_every == 0L ||
                    b == B
               )
          ) {
               cat(
                    sprintf(
                         "Iteración %d de %d | K = %d | alpha = %.3f\n",
                         b,
                         B,
                         K,
                         alpha
                    )
               )
          }
     }
     
     acceptance_rate <- acceptance_post / (B - burn)
     
     result <- list(
          call = match.call(),
          y = y,
          hyperparameters = hyperparameters,
          settings = list(
               B = B,
               burn = burn,
               thin = thin,
               K_init = K_init,
               random_scan = random_scan,
               adapt = adapt,
               adapt_interval = adapt_interval,
               target_accept = target_accept,
               seed = seed
          ),
          draws = draws,
          partitions = partitions,
          acceptance = acceptance_rate,
          proposal_sd = proposal_sd,
          last_state = list(
               xi = xi,
               alpha = alpha,
               mu = mu,
               kappa = kappa,
               beta = beta,
               statistics = statistics
          )
     )
     
     result
}

# Hiperparámetros --------------------------------------------------------------

hyperparameters <- list(
     a_alpha = 1,
     b_alpha = 1,
     a_mu    = mean(y),
     b_mu    = 2 * var(y),
     a_kappa = 2,
     b_kappa = 1,
     a_beta  = 2,
     b_beta  = 1,
     cc      = 2
)

# Ajuste del modelo ------------------------------------------------------------

fit_dpm_galaxies <- dpm_ls_collapsed(
     y               = y,
     B               = 110000L,
     burn            = 10000L,
     thin            = 10L,
     K_init          = 8L,
     hyperparameters = hyperparameters,
     adapt           = TRUE,
     adapt_interval  = 50L,
     target_accept   = 0.44,
     random_scan     = TRUE,
     seed            = 1702,
     verbose         = TRUE,
     progress_every  = 10000L
)


# Organización de las cadenas --------------------------------------------------

xi_draws           <- fit_dpm_galaxies$partitions
alpha_draws        <- fit_dpm_galaxies$draws$alpha
K_draws            <- fit_dpm_galaxies$draws$K
log_marginal_draws <- fit_dpm_galaxies$draws$log_marginal

n_draws <- nrow(xi_draws)

# Diagnósticos de convergencia -------------------------------------------------

# Resume los diagnósticos de una única cadena MCMC.
diagnostic_summary <- function(draws) {
     draws <- as.numeric(draws)
     
     n_draws <- length(draws)
     draw_sd <- sd(draws)
     
     autocorrelation <- acf(
          draws,
          plot = FALSE,
          lag.max = min(n_draws - 1L, 1000L)
     )$acf[-1L]
     
     pair_count <- floor(length(autocorrelation) / 2L)
     
     pair_sums <- autocorrelation[2L * seq_len(pair_count) - 1L] +
          autocorrelation[2L * seq_len(pair_count)]
     
     first_nonpositive <- which(pair_sums <= 0)[1L]
     
     if (!is.na(first_nonpositive)) {
          pair_sums <- pair_sums[seq_len(first_nonpositive - 1L)]
     }
     
     effective_sample_size <- min(
          n_draws,
          n_draws / (1 + 2 * sum(pair_sums))
     )
     
     c(
          mean   = mean(draws),
          sd     = draw_sd,
          q025   = quantile(draws, 0.025, names = FALSE),
          median = median(draws),
          q975   = quantile(draws, 0.975, names = FALSE),
          acf1   = autocorrelation[1L],
          ESS    = effective_sample_size,
          MCSE   = draw_sd / sqrt(effective_sample_size)
     )
}

diagnostics <- rbind(
     alpha = diagnostic_summary(alpha_draws),
     K = diagnostic_summary(K_draws),
     log_marginal = diagnostic_summary(log_marginal_draws)
)

print(
     round(
          diagnostics,
          digits = 3
     )
)

# Cadena de la log-verosimilitud
pdf(
     file = "dpm_galaxias_logverosimilitud.pdf",
     width = 7,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     fit_dpm_galaxies$draws$iteration,
     log_marginal_draws,
     type = "p",
     pch  = 16,
     cex  = 0.5,
     col  = adjustcolor(1, 0.3),
     xlab = "Iteración",
     ylab = "Log-verosimilitud marginal",
     main = ""
)

dev.off()

# Matriz posterior de coagrupamiento ------------------------------------------

posterior_similarity <- matrix(
     0,
     nrow = n,
     ncol = n
)

for (b in seq_len(n_draws)) {
     xi <- xi_draws[b, ]
     
     posterior_similarity <- posterior_similarity +
          outer(
               xi,
               xi,
               FUN = "=="
          )
}

posterior_similarity <- posterior_similarity / n_draws


# Estimación de Dahl -----------------------------------------------------------

dahl_loss <- numeric(n_draws)
constant_term <- sum(posterior_similarity^2)

for (b in seq_len(n_draws)) {
     incidence <- outer(
          xi_draws[b, ],
          xi_draws[b, ],
          FUN = "=="
     )
     
     dahl_loss[b] <- sum(incidence) -
          2 * sum(incidence * posterior_similarity) +
          constant_term
}

dahl_index <- which.min(dahl_loss)
xi_hat <- xi_draws[dahl_index, ]
K_hat <- max(xi_hat)

# Región de credibilidad para la partición ------------------------------------

# Calcula la distancia de Binder entre dos particiones.
binder_distance <- function(xi_1, xi_2) {
     incidence_1 <- outer(
          xi_1,
          xi_1,
          FUN = "=="
     )
     
     incidence_2 <- outer(
          xi_2,
          xi_2,
          FUN = "=="
     )
     
     upper_indices <- upper.tri(incidence_1)
     
     sum(
          incidence_1[upper_indices] !=
               incidence_2[upper_indices]
     )
}

# Calcula la moda de una variable discreta.
discrete_mode <- function(x) {
     frequencies <- table(x)
     
     as.numeric(
          names(
               frequencies[which.max(frequencies)]
          )
     )
}

binder_distances <- vapply(
     X = seq_len(n_draws),
     FUN = function(b) {
          binder_distance(
               xi_draws[b, ],
               xi_hat
          )
     },
     FUN.VALUE = numeric(1)
)

binder_radius <- quantile(
     binder_distances,
     probs = 0.95,
     type = 1,
     names = FALSE
)

inside_credible_ball <- binder_distances <= binder_radius

credible_mass <- mean(inside_credible_ball)

credible_K_range <- range(
     K_draws[inside_credible_ball]
)

K_summary <- c(
     mean   = mean(K_draws),
     sd     = sd(K_draws),
     mode   = discrete_mode(K_draws),
     q025   = quantile(K_draws, 0.025, names = FALSE),
     median = median(K_draws),
     q975   = quantile(K_draws, 0.975, names = FALSE)
)

print(
     round(
          c(
               K_summary,
               K_hat         = K_hat,
               binder_radius = binder_radius,
               credible_mass = credible_mass,
               K_min_ball    = credible_K_range[1],
               K_max_ball    = credible_K_range[2]
          ),
          digits = 3
     )
)

# Gráfico de la matriz de coagrupamiento ---------------------------------------

partition_order <- order(
     ave(y, xi_hat, FUN = mean),
     y
)

P_ordered <- posterior_similarity[
     partition_order,
     partition_order
]

colorscale <- c(
     "white",
     rev(heat.colors(100))
)

pdf(
     file = "dpm_galaxias_particion.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mar = c(2.75, 2.75, 0.5, 0.5),
     mgp = c(1.7, 0.7, 0)
)

corrplot::corrplot(
     P_ordered,
     is.corr = FALSE,
     method = "color",
     col = colorscale,
     tl.pos = "n",
     addgrid.col = NA
)

dev.off()

# Distribución posterior del número de componentes -----------------------------

K_probabilities <- table(K_draws) / length(K_draws)

pdf(
     file = "dpm_galaxias_numero_componentes.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

barplot(
     K_probabilities,
     border = NA,
     col = "gray85",
     xlab = expression(K[n]),
     ylab = "Probabilidad posterior",
     ylim = c(0,0.13),
     main = ""
)

box()

dev.off()

# Incertidumbre de la partición ------------------------------------------------

pdf(
     file = "dpm_galaxias_distancia_binder.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

hist(
     binder_distances,
     breaks = 25,
     probability = TRUE,
     col = "gray85",
     border = "white",
     xlab = "Núm. pares dicordantes",
     ylab = "Densidad",
     main = ""
)

abline(
     v = binder_radius,
     lwd = 2,
     lty = 2
)

box()

dev.off()

# Extracciones de densidades y réplicas ---------------------------------------

# Genera un vector de probabilidades con distribución Dirichlet.
rdirichlet_one <- function(parameters) {
     gamma_draws <- rgamma(
          n = length(parameters),
          shape = parameters,
          rate = 1
     )
     
     gamma_draws / sum(gamma_draws)
}

# Genera pesos mediante una representación stick-breaking truncada.
stick_weights <- function(alpha, H) {
     breaks <- rbeta(
          n = H,
          shape1 = 1,
          shape2 = alpha
     )
     
     breaks[H] <- 1
     
     breaks * c(
          1,
          cumprod(1 - breaks[-H])
     )
}

# Posprocesamiento posterior
mu_draws    <- fit_dpm_galaxies$draws$mu
kappa_draws <- fit_dpm_galaxies$draws$kappa
beta_draws  <- fit_dpm_galaxies$draws$beta

set.seed(2468)

B_post <- min(
     5000L,
     n_draws
)

post_indices <- sort(
     sample.int(
          n = n_draws,
          size = B_post,
          replace = FALSE
     )
)

grid_margin <- 0.10 * diff(range(y))

density_grid <- seq(
     min(y) - grid_margin,
     max(y) + grid_margin,
     length.out = 500
)

density_draws <- matrix(
     NA_real_,
     nrow = B_post,
     ncol = length(density_grid)
)

predictive_density_draws <- matrix(
     NA_real_,
     nrow = B_post,
     ncol = length(density_grid)
)

y_new_draws <- numeric(B_post)

y_rep <- matrix(
     NA_real_,
     nrow = B_post,
     ncol = n
)

H_residual <- 30L

post_progress <- max(
     1L,
     floor(B_post / 10L)
)

for (s in seq_len(B_post)) {
     index <- post_indices[s]
     
     xi    <- xi_draws[index, ]
     alpha <- alpha_draws[index]
     mu    <- mu_draws[index]
     kappa <- kappa_draws[index]
     beta  <- beta_draws[index]
     
     statistics <- partition_statistics(
          y = y,
          xi = xi
     )
     
     K <- length(statistics$n)
     
     counts  <- statistics$n
     sums    <- statistics$sum
     sums_sq <- statistics$sum_sq
     
     mu_occupied     <- numeric(K)
     sigma2_occupied <- numeric(K)
     
     # Densidad predictiva asociada con un componente nuevo
     predictive_density <- alpha / (alpha + n) *
          exp(
               log_predictive_nig_stats(
                    y_new = density_grid,
                    n_l = 0L,
                    sum_l = 0,
                    sum_sq_l = 0,
                    mu = mu,
                    kappa = kappa,
                    beta = beta,
                    cc = hyperparameters$cc
               )
          )
     
     for (ell in seq_len(K)) {
          posterior <- nig_posterior_stats(
               n_l = counts[ell],
               sum_l = sums[ell],
               sum_sq_l = sums_sq[ell],
               mu = mu,
               kappa = kappa,
               beta = beta,
               cc = hyperparameters$cc
          )
          
          sigma2_occupied[ell] <- 1 / rgamma(
               n = 1L,
               shape = posterior$c_tilde,
               rate = posterior$beta_tilde
          )
          
          mu_occupied[ell] <- rnorm(
               n = 1L,
               mean = posterior$mu_tilde,
               sd = sqrt(
                    sigma2_occupied[ell] /
                         posterior$kappa_tilde
               )
          )
          
          # Densidad predictiva del componente ocupado
          predictive_density <- predictive_density +
               counts[ell] / (alpha + n) *
               exp(
                    log_predictive_nig_stats(
                         y_new = density_grid,
                         n_l = counts[ell],
                         sum_l = sums[ell],
                         sum_sq_l = sums_sq[ell],
                         mu = mu,
                         kappa = kappa,
                         beta = beta,
                         cc = hyperparameters$cc
                    )
               )
     }
     
     # Masas de los componentes ocupados y de la parte residual
     posterior_masses <- rdirichlet_one(
          c(
               alpha,
               counts
          )
     )
     
     # Aproximación truncada de la parte residual del DP
     residual_weights <- stick_weights(
          alpha = alpha,
          H = H_residual
     )
     
     sigma2_residual <- 1 / rgamma(
          n = H_residual,
          shape = hyperparameters$cc,
          rate = beta
     )
     
     mu_residual <- rnorm(
          n = H_residual,
          mean = mu,
          sd = sqrt(
               sigma2_residual / kappa
          )
     )
     
     density_draw <- numeric(
          length(density_grid)
     )
     
     for (ell in seq_len(K)) {
          density_draw <- density_draw +
               posterior_masses[ell + 1L] *
               dnorm(
                    density_grid,
                    mean = mu_occupied[ell],
                    sd = sqrt(sigma2_occupied[ell])
               )
     }
     
     residual_density <- numeric(
          length(density_grid)
     )
     
     for (h in seq_len(H_residual)) {
          residual_density <- residual_density +
               residual_weights[h] *
               dnorm(
                    density_grid,
                    mean = mu_residual[h],
                    sd = sqrt(sigma2_residual[h])
               )
     }
     
     density_draw <- density_draw +
          posterior_masses[1L] *
          residual_density
     
     mixture_weights <- c(
          posterior_masses[-1L],
          posterior_masses[1L] * residual_weights
     )
     
     mixture_means <- c(
          mu_occupied,
          mu_residual
     )
     
     mixture_variances <- c(
          sigma2_occupied,
          sigma2_residual
     )
     
     component_draws <- sample.int(
          n = length(mixture_weights),
          size = n + 1L,
          replace = TRUE,
          prob = mixture_weights
     )
     
     predictive_sample <- rnorm(
          n = n + 1L,
          mean = mixture_means[component_draws],
          sd = sqrt(
               mixture_variances[component_draws]
          )
     )
     
     density_draws[s, ] <- density_draw
     predictive_density_draws[s, ] <- predictive_density
     y_new_draws[s] <- predictive_sample[1L]
     y_rep[s, ] <- predictive_sample[-1L]
     
     if (s %% post_progress == 0L) {
          cat(
               "Posprocesamiento: ",
               round(100 * s / B_post),
               "% completado\n",
               sep = ""
          )
     }
}

# Estimación posterior de la densidad -----------------------------------------

posterior_density_mean <- colMeans(
     predictive_density_draws
)

posterior_density_lower <- apply(
     density_draws,
     2,
     quantile,
     probs = 0.025,
     names = FALSE
)

posterior_density_upper <- apply(
     density_draws,
     2,
     quantile,
     probs = 0.975,
     names = FALSE
)

density_histogram <- hist(
     y,
     breaks = 25,
     plot = FALSE
)

density_ylim <- c(
     0,
     max(
          density_histogram$density,
          posterior_density_upper
     )
)

pdf(
     file = "dpm_galaxias_densidad.pdf",
     width = 7,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     density_histogram,
     freq = FALSE,
     col = "gray92",
     border = "white",
     ylim = density_ylim,
     xlab = "Velocidad",
     ylab = "Densidad",
     main = ""
)

polygon(
     x = c(
          density_grid,
          rev(density_grid)
     ),
     y = c(
          posterior_density_lower,
          rev(posterior_density_upper)
     ),
     col = adjustcolor(
          "#6BAED6",
          alpha.f = 0.35
     ),
     border = NA
)

lines(
     density_grid,
     posterior_density_mean,
     col = "#08519C",
     lwd = 2.5
)

rug(y)

dev.off()

# Función de distribución empírica ---------------------------------------------

replicated_ecdf <- matrix(
     NA_real_,
     nrow = B_post,
     ncol = length(density_grid)
)

for (b in seq_len(B_post)) {
     replicated_ecdf[b, ] <- findInterval(
          density_grid,
          sort(y_rep[b, ]),
          rightmost.closed = TRUE
     ) / n
}

observed_ecdf <- findInterval(
     density_grid,
     sort(y),
     rightmost.closed = TRUE
) / n

ecdf_mean <- colMeans(
     replicated_ecdf
)

ecdf_lower <- apply(
     replicated_ecdf,
     2,
     quantile,
     probs = 0.025,
     names = FALSE
)

ecdf_upper <- apply(
     replicated_ecdf,
     2,
     quantile,
     probs = 0.975,
     names = FALSE
)

pdf(
     file = "dpm_galaxias_ppc_ecdf.pdf",
     width = 5,
     height = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     density_grid,
     observed_ecdf,
     type = "n",
     ylim = c(0, 1),
     xlab = "Velocidad",
     ylab = "Función de distribución",
     main = ""
)

polygon(
     x = c(
          density_grid,
          rev(density_grid)
     ),
     y = c(
          ecdf_lower,
          rev(ecdf_upper)
     ),
     col = adjustcolor(
          "#6BAED6",
          alpha.f = 0.35
     ),
     border = NA
)

lines(
     density_grid,
     ecdf_mean,
     col = "#08519C",
     lwd = 1.5,
     lty = 2
)

lines(
     density_grid,
     observed_ecdf,
     type = "s",
     lwd = 2
)

dev.off()

# Fin --------------------------------------------------------------------------