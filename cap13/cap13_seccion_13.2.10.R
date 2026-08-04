# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Datos ------------------------------------------------------------------------

n          <- 300
omega_true <- c(0.25, 0.5, 0.25)
mu_true    <- c(-3, 0, 3)
sigma2     <- 0.36

set.seed(123)

z_true <- sample.int(
     n       = length(omega_true),
     size    = n,
     replace = TRUE,
     prob    = omega_true
)

set.seed(123)

y <- rnorm(
     n    = n,
     mean = mu_true[z_true],
     sd   = sqrt(sigma2)
)

# Visualización

pdf(
     file      = "dpm_simulacion_datos.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     y,
     breaks      = 35,
     probability = TRUE,
     main        = "",
     xlab        = expression(y),
     ylab        = "Densidad",
     col         = "grey92",
     border      = "white"
)

abline(
     v   = mu_true,
     lwd = 2,
     lty = 2,
     col = 2
)

box()

dev.off()

# Funciones auxiliares ---------------------------------------------------------

sample_discrete <- function(log_weight) {
     log_weight <- log_weight - max(log_weight)
     prob       <- exp(log_weight)
     prob       <- prob / sum(prob)
     
     sample.int(
          n    = length(prob),
          size = 1,
          prob = prob
     )
}

relabel_partition <- function(xi) {
     occupied <- sort(unique(xi[xi > 0]))
     out      <- integer(length(xi))
     
     if (length(occupied) > 0) {
          for (ell in seq_along(occupied)) {
               out[xi == occupied[ell]] <- ell
          }
     }
     
     out
}

initialize_partition <- function(y, K) {
     if (K == 1) {
          return(rep(1L, length(y)))
     }
     
     fit <- kmeans(
          x       = matrix(y, ncol = 1),
          centers = K,
          nstart  = 30
     )
     
     cluster_mean  <- tapply(y, fit$cluster, mean)
     order_cluster <- order(cluster_mean)
     map           <- integer(K)
     
     map[order_cluster] <- seq_len(K)
     
     as.integer(map[fit$cluster])
}

posterior_component <- function(y_group, sigma2, mu0, tau02) {
     n_group <- length(y_group)
     
     posterior_variance <- 1 / (
          1 / tau02 +
               n_group / sigma2
     )
     
     posterior_mean <- posterior_variance * (
          mu0 / tau02 +
               sum(y_group) / sigma2
     )
     
     list(
          mean     = posterior_mean,
          variance = posterior_variance
     )
}

log_predictive_existing <- function(
          y_i,
          y_group,
          sigma2,
          mu0,
          tau02
) {
     pars <- posterior_component(
          y_group = y_group,
          sigma2  = sigma2,
          mu0     = mu0,
          tau02   = tau02
     )
     
     dnorm(
          x    = y_i,
          mean = pars$mean,
          sd   = sqrt(sigma2 + pars$variance),
          log  = TRUE
     )
}

log_predictive_new <- function(y_i, sigma2, mu0, tau02) {
     dnorm(
          x    = y_i,
          mean = mu0,
          sd   = sqrt(sigma2 + tau02),
          log  = TRUE
     )
}

sample_component_parameters <- function(
          y,
          xi,
          sigma2,
          mu0,
          tau02
) {
     K     <- max(xi)
     theta <- numeric(K)
     
     for (ell in seq_len(K)) {
          pars <- posterior_component(
               y_group = y[xi == ell],
               sigma2  = sigma2,
               mu0     = mu0,
               tau02   = tau02
          )
          
          theta[ell] <- rnorm(
               n    = 1,
               mean = pars$mean,
               sd   = sqrt(pars$variance)
          )
     }
     
     theta
}

update_alpha_escobar_west <- function(
          alpha,
          K,
          n,
          a_alpha,
          b_alpha
) {
     u <- rbeta(
          n      = 1,
          shape1 = alpha + 1,
          shape2 = n
     )
     
     rate_alpha <- b_alpha - log(u)
     
     prob_first <- (a_alpha + K - 1) / (
          a_alpha + K - 1 +
               n * rate_alpha
     )
     
     first_component <- rbinom(
          n    = 1,
          size = 1,
          prob = prob_first
     )
     
     shape_alpha <- if (first_component == 1) {
          a_alpha + K
     } else {
          a_alpha + K - 1
     }
     
     rgamma(
          n     = 1,
          shape = shape_alpha,
          rate  = rate_alpha
     )
}

clamp_unit <- function(x) {
     pmin(
          pmax(x, .Machine$double.eps),
          1 - .Machine$double.eps
     )
}

stick_weights <- function(V) {
     H         <- length(V)
     omega     <- numeric(H)
     remaining <- 1
     
     for (ell in seq_len(H)) {
          omega[ell] <- V[ell] * remaining
          remaining  <- remaining * (1 - V[ell])
     }
     
     omega / sum(omega)
}

conditional_loglik <- function(y, xi, theta, sigma2) {
     sum(
          dnorm(
               x    = y,
               mean = theta[xi],
               sd   = sqrt(sigma2),
               log  = TRUE
          )
     )
}

# Muestreador completamente colapsado ------------------------------------------

mcmc_fully_collapsed <- function(
          y,
          sigma2,
          mu0,
          tau02,
          a_alpha,
          b_alpha,
          B,
          burn,
          thin,
          init_K,
          seed
) {
     set.seed(seed)
     
     # Inicialización de la cadena
     n     <- length(y)
     xi    <- initialize_partition(y, init_K)
     alpha <- rgamma(
          n     = 1,
          shape = a_alpha,
          rate  = b_alpha
     )
     
     # Iteraciones que se conservarán
     keep <- seq.int(
          from = burn + 1L,
          to   = B,
          by   = thin
     )
     
     keep_iteration       <- logical(B)
     keep_iteration[keep] <- TRUE
     
     n_keep <- length(keep)
     
     alpha_draw  <- numeric(n_keep)
     K_draw      <- integer(n_keep)
     loglik_draw <- numeric(n_keep)
     
     save_index <- 0L
     
     elapsed <- system.time({
          for (b in seq_len(B)) {
               # Actualización aleatoria de las asignaciones
               update_order <- sample.int(n)
               
               for (i in update_order) {
                    # Retirar temporalmente la observación
                    xi[i] <- 0L
                    xi    <- relabel_partition(xi)
                    
                    K_minus <- max(xi)
                    
                    log_weight <- numeric(K_minus + 1L)
                    
                    # Pesos de los componentes ocupados
                    if (K_minus > 0L) {
                         component_size <- tabulate(
                              xi,
                              nbins = K_minus
                         )
                         
                         for (ell in seq_len(K_minus)) {
                              log_weight[ell] <-
                                   log(component_size[ell]) +
                                   log_predictive_existing(
                                        y_i     = y[i],
                                        y_group = y[xi == ell],
                                        sigma2  = sigma2,
                                        mu0     = mu0,
                                        tau02   = tau02
                                   )
                         }
                    }
                    
                    # Peso correspondiente a un componente nuevo
                    log_weight[K_minus + 1L] <-
                         log(alpha) +
                         log_predictive_new(
                              y_i    = y[i],
                              sigma2 = sigma2,
                              mu0    = mu0,
                              tau02  = tau02
                         )
                    
                    # Generar la nueva asignación
                    xi[i] <- sample_discrete(log_weight)
               }
               
               K <- max(xi)
               
               # Actualizar el parámetro de concentración
               alpha <- update_alpha_escobar_west(
                    alpha   = alpha,
                    K       = K,
                    n       = n,
                    a_alpha = a_alpha,
                    b_alpha = b_alpha
               )
               
               # Guardar las extracciones seleccionadas
               if (keep_iteration[b]) {
                    save_index <- save_index + 1L
                    
                    # Recuperar los parámetros integrados
                    theta <- sample_component_parameters(
                         y      = y,
                         xi     = xi,
                         sigma2 = sigma2,
                         mu0    = mu0,
                         tau02  = tau02
                    )
                    
                    alpha_draw[save_index] <- alpha
                    K_draw[save_index]     <- K
                    
                    loglik_draw[save_index] <- conditional_loglik(
                         y      = y,
                         xi     = xi,
                         theta  = theta,
                         sigma2 = sigma2
                    )
               }
          }
     })["elapsed"]
     
     list(
          alpha   = alpha_draw,
          K       = K_draw,
          loglik  = loglik_draw,
          elapsed = as.numeric(elapsed)
     )
}

# Muestreador parcialmente colapsado -------------------------------------------

mcmc_partially_collapsed <- function(
          y,
          sigma2,
          mu0,
          tau02,
          a_alpha,
          b_alpha,
          B,
          burn,
          thin,
          init_K,
          seed
) {
     set.seed(seed)
     
     # Inicialización de la cadena
     n  <- length(y)
     xi <- initialize_partition(y, init_K)
     
     theta <- sample_component_parameters(
          y      = y,
          xi     = xi,
          sigma2 = sigma2,
          mu0    = mu0,
          tau02  = tau02
     )
     
     alpha <- rgamma(
          n     = 1,
          shape = a_alpha,
          rate  = b_alpha
     )
     
     # Iteraciones que se conservarán
     keep <- seq.int(
          from = burn + 1L,
          to   = B,
          by   = thin
     )
     
     keep_iteration       <- logical(B)
     keep_iteration[keep] <- TRUE
     
     n_keep <- length(keep)
     
     alpha_draw  <- numeric(n_keep)
     K_draw      <- integer(n_keep)
     loglik_draw <- numeric(n_keep)
     
     save_index <- 0L
     
     elapsed <- system.time({
          for (b in seq_len(B)) {
               # Actualización aleatoria de las asignaciones
               update_order <- sample.int(n)
               
               for (i in update_order) {
                    # Retirar temporalmente la observación
                    old_component <- xi[i]
                    xi[i]         <- 0L
                    
                    # Eliminar el componente si queda vacío
                    if (!any(xi == old_component)) {
                         theta <- theta[-old_component]
                         
                         xi[xi > old_component] <-
                              xi[xi > old_component] - 1L
                    }
                    
                    K_minus   <- length(theta)
                    log_weight <- numeric(K_minus + 1L)
                    
                    # Pesos de los componentes ocupados
                    if (K_minus > 0L) {
                         component_size <- tabulate(
                              xi,
                              nbins = K_minus
                         )
                         
                         for (ell in seq_len(K_minus)) {
                              log_weight[ell] <-
                                   log(component_size[ell]) +
                                   dnorm(
                                        x    = y[i],
                                        mean = theta[ell],
                                        sd   = sqrt(sigma2),
                                        log  = TRUE
                                   )
                         }
                    }
                    
                    # Peso correspondiente a un componente nuevo
                    log_weight[K_minus + 1L] <-
                         log(alpha) +
                         log_predictive_new(
                              y_i    = y[i],
                              sigma2 = sigma2,
                              mu0    = mu0,
                              tau02  = tau02
                         )
                    
                    # Generar la nueva asignación
                    new_component <- sample_discrete(log_weight)
                    xi[i]         <- new_component
                    
                    # Generar el parámetro del componente nuevo
                    if (new_component == K_minus + 1L) {
                         pars <- posterior_component(
                              y_group = y[i],
                              sigma2  = sigma2,
                              mu0     = mu0,
                              tau02   = tau02
                         )
                         
                         theta <- c(
                              theta,
                              rnorm(
                                   n    = 1,
                                   mean = pars$mean,
                                   sd   = sqrt(pars$variance)
                              )
                         )
                    }
               }
               
               K <- max(xi)
               
               # Actualizar los parámetros de los componentes
               theta <- sample_component_parameters(
                    y      = y,
                    xi     = xi,
                    sigma2 = sigma2,
                    mu0    = mu0,
                    tau02  = tau02
               )
               
               # Actualizar el parámetro de concentración
               alpha <- update_alpha_escobar_west(
                    alpha   = alpha,
                    K       = K,
                    n       = n,
                    a_alpha = a_alpha,
                    b_alpha = b_alpha
               )
               
               # Guardar las extracciones seleccionadas
               if (keep_iteration[b]) {
                    save_index <- save_index + 1L
                    
                    alpha_draw[save_index] <- alpha
                    K_draw[save_index]     <- K
                    
                    loglik_draw[save_index] <- conditional_loglik(
                         y      = y,
                         xi     = xi,
                         theta  = theta,
                         sigma2 = sigma2
                    )
               }
          }
     })["elapsed"]
     
     list(
          alpha   = alpha_draw,
          K       = K_draw,
          loglik  = loglik_draw,
          elapsed = as.numeric(elapsed)
     )
}

# Muestreador no colapsado con truncamiento ------------------------------------

mcmc_truncated_blocked <- function(
          y,
          sigma2,
          mu0,
          tau02,
          a_alpha,
          b_alpha,
          H,
          B,
          burn,
          thin,
          init_K,
          seed
) {
     set.seed(seed)
     
     # Inicialización de la cadena
     n           <- length(y)
     sd_obs      <- sqrt(sigma2)
     sd_base     <- sqrt(tau02)
     stick_index <- seq_len(H - 1L)
     
     xi <- initialize_partition(
          y = y,
          K = min(init_K, H)
     )
     
     alpha <- rgamma(
          n     = 1,
          shape = a_alpha,
          rate  = b_alpha
     )
     
     # Inicialización de los átomos
     theta <- rnorm(
          n    = H,
          mean = mu0,
          sd   = sd_base
     )
     
     K_init <- max(xi)
     
     for (ell in seq_len(K_init)) {
          pars <- posterior_component(
               y_group = y[xi == ell],
               sigma2  = sigma2,
               mu0     = mu0,
               tau02   = tau02
          )
          
          theta[ell] <- rnorm(
               n    = 1,
               mean = pars$mean,
               sd   = sqrt(pars$variance)
          )
     }
     
     # Inicialización de las variables de partición y los pesos
     counts <- tabulate(
          xi,
          nbins = H
     )
     
     V <- numeric(H)
     
     for (ell in stick_index) {
          V[ell] <- clamp_unit(
               rbeta(
                    n      = 1,
                    shape1 = 1 + counts[ell],
                    shape2 = alpha + sum(counts[(ell + 1L):H])
               )
          )
     }
     
     V[H]  <- 1
     omega <- stick_weights(V)
     
     # Iteraciones que se conservarán
     keep <- seq.int(
          from = burn + 1L,
          to   = B,
          by   = thin
     )
     
     keep_iteration       <- logical(B)
     keep_iteration[keep] <- TRUE
     
     n_keep <- length(keep)
     
     alpha_draw       <- numeric(n_keep)
     K_draw           <- integer(n_keep)
     loglik_draw      <- numeric(n_keep)
     last_weight_draw <- numeric(n_keep)
     max_label_draw   <- integer(n_keep)
     
     occupancy_draw <- matrix(
          0L,
          nrow = n_keep,
          ncol = H
     )
     
     save_index <- 0L
     
     elapsed <- system.time({
          for (b in seq_len(B)) {
               # Actualizar las asignaciones
               for (i in seq_len(n)) {
                    log_weight <-
                         log(omega) +
                         dnorm(
                              x    = y[i],
                              mean = theta,
                              sd   = sd_obs,
                              log  = TRUE
                         )
                    
                    xi[i] <- sample_discrete(log_weight)
               }
               
               counts <- tabulate(
                    xi,
                    nbins = H
               )
               
               # Actualizar los átomos
               for (ell in seq_len(H)) {
                    if (counts[ell] > 0L) {
                         pars <- posterior_component(
                              y_group = y[xi == ell],
                              sigma2  = sigma2,
                              mu0     = mu0,
                              tau02   = tau02
                         )
                         
                         theta[ell] <- rnorm(
                              n    = 1,
                              mean = pars$mean,
                              sd   = sqrt(pars$variance)
                         )
                    } else {
                         theta[ell] <- rnorm(
                              n    = 1,
                              mean = mu0,
                              sd   = sd_base
                         )
                    }
               }
               
               # Actualizar las variables de partición
               for (ell in stick_index) {
                    V[ell] <- clamp_unit(
                         rbeta(
                              n      = 1,
                              shape1 = 1 + counts[ell],
                              shape2 = alpha + sum(counts[(ell + 1L):H])
                         )
                    )
               }
               
               V[H] <- 1
               
               # Reconstruir los pesos de la mezcla
               omega <- stick_weights(V)
               
               # Actualizar el parámetro de concentración
               alpha <- rgamma(
                    n     = 1,
                    shape = a_alpha + H - 1,
                    rate  = b_alpha - sum(log(1 - V[stick_index]))
               )
               
               # Guardar las extracciones seleccionadas
               if (keep_iteration[b]) {
                    save_index <- save_index + 1L
                    occupied   <- counts > 0L
                    
                    alpha_draw[save_index] <- alpha
                    K_draw[save_index]     <- sum(occupied)
                    
                    loglik_draw[save_index] <- conditional_loglik(
                         y      = y,
                         xi     = xi,
                         theta  = theta,
                         sigma2 = sigma2
                    )
                    
                    # Diagnósticos del truncamiento
                    last_weight_draw[save_index] <- omega[H]
                    max_label_draw[save_index]   <- max(which(occupied))
                    occupancy_draw[save_index, ] <- as.integer(occupied)
               }
          }
     })["elapsed"]
     
     list(
          alpha       = alpha_draw,
          K           = K_draw,
          loglik      = loglik_draw,
          elapsed     = as.numeric(elapsed),
          last_weight = last_weight_draw,
          max_label   = max_label_draw,
          occupancy   = occupancy_draw
     )
}

# Ejecución de las cadenas -----------------------------------------------------

# Hiperparámetros
mu0     <- 0
tau02   <- 9
a_alpha <- 1
b_alpha <- 1

# Configuración MCMC
B        <- 11000
burn     <- 1000
thin     <- 1
n_chains <- 3
init_K   <- c(1, 4, 8)
H        <- 25

# Semillas para reproducibilidad
seeds_fully     <- 123 + seq_len(n_chains)
seeds_partial   <- 123 + seq_len(n_chains)
seeds_truncated <- 123 + seq_len(n_chains)

# Almacenamiento de las cadenas
fully_collapsed     <- vector("list", length = n_chains)
partially_collapsed <- vector("list", length = n_chains)
truncated_blocked   <- vector("list", length = n_chains)

# Ajuste de los tres muestreadores
for (chain in seq_len(n_chains)) {
     cat(
          "\nCadena",
          chain,
          "de",
          n_chains,
          "\n"
     )
     
     # Muestreador completamente colapsado
     cat("  Ejecutando muestreador completamente colapsado...\n")
     
     fully_collapsed[[chain]] <- mcmc_fully_collapsed(
          y       = y,
          sigma2  = sigma2,
          mu0     = mu0,
          tau02   = tau02,
          a_alpha = a_alpha,
          b_alpha = b_alpha,
          B       = B,
          burn    = burn,
          thin    = thin,
          init_K  = init_K[chain],
          seed    = seeds_fully[chain]
     )
     
     # Muestreador parcialmente colapsado
     cat("  Ejecutando muestreador parcialmente colapsado...\n")
     
     partially_collapsed[[chain]] <- mcmc_partially_collapsed(
          y       = y,
          sigma2  = sigma2,
          mu0     = mu0,
          tau02   = tau02,
          a_alpha = a_alpha,
          b_alpha = b_alpha,
          B       = B,
          burn    = burn,
          thin    = thin,
          init_K  = init_K[chain],
          seed    = seeds_partial[chain]
     )
     
     # Muestreador no colapsado truncado
     cat("  Ejecutando muestreador no colapsado truncado...\n")
     
     truncated_blocked[[chain]] <- mcmc_truncated_blocked(
          y       = y,
          sigma2  = sigma2,
          mu0     = mu0,
          tau02   = tau02,
          a_alpha = a_alpha,
          b_alpha = b_alpha,
          H       = H,
          B       = B,
          burn    = burn,
          thin    = thin,
          init_K  = init_K[chain],
          seed    = seeds_truncated[chain]
     )
     
     cat("  Cadena completada.\n")
}

# Resultados organizados por método
samplers <- list(
     "Completamente colapsado" = fully_collapsed,
     "Parcialmente colapsado"  = partially_collapsed,
     "No colapsado truncado"   = truncated_blocked
)

# Organizar la configuración
settings <- list(
     n          = length(y),
     omega_true = omega_true,
     mu_true    = mu_true,
     sigma2     = sigma2,
     mu0        = mu0,
     tau02      = tau02,
     a_alpha    = a_alpha,
     b_alpha    = b_alpha,
     B          = B,
     burn       = burn,
     thin       = thin,
     n_chains   = n_chains,
     init_K     = init_K,
     H          = H
)

# Guardar datos, configuración y resultados
save(
     y,
     z_true,
     settings,
     samplers,
     file = "dpm_simulacion_mcmc_comparacion.RData"
)

# Diagnósticos -----------------------------------------------------------------

ess_one_chain <- function(
          x,
          max_lag = min(length(x) - 1L, 1000L)
) {
     n_draw          <- length(x)
     sample_variance <- var(x)
     
     # Casos sin información sobre autocorrelación
     if (
          n_draw < 3L ||
          !is.finite(sample_variance) ||
          sample_variance <= 0
     ) {
          return(as.numeric(n_draw))
     }
     
     max_lag <- min(
          as.integer(max_lag),
          n_draw - 1L
     )
     
     # Autocorrelaciones muestrales
     rho <- as.numeric(
          acf(
               x       = x,
               plot    = FALSE,
               lag.max = max_lag,
               demean  = TRUE
          )$acf[-1L]
     )
     
     # Suma inicial de pares positivos
     pair_sum     <- 0
     n_paired_lags <- length(rho) - length(rho) %% 2L
     
     if (n_paired_lags >= 2L) {
          for (j in seq.int(1L, n_paired_lags, by = 2L)) {
               current_pair <- rho[j] + rho[j + 1L]
               
               if (
                    !is.finite(current_pair) ||
                    current_pair <= 0
               ) {
                    break
               }
               
               pair_sum <- pair_sum + current_pair
          }
     }
     
     ess <- n_draw / (1 + 2 * pair_sum)
     
     # Evitar valores fuera del rango admisible
     max(
          1,
          min(n_draw, ess)
     )
}

rhat_basic <- function(chain_list) {
     n_chain <- length(chain_list)
     
     if (n_chain < 2L) {
          return(NA_real_)
     }
     
     n_draw <- min(
          vapply(
               chain_list,
               length,
               integer(1)
          )
     )
     
     if (n_draw < 2L) {
          return(NA_real_)
     }
     
     # Igualar la longitud de las cadenas
     draw_matrix <- do.call(
          cbind,
          lapply(
               chain_list,
               function(x) x[seq_len(n_draw)]
          )
     )
     
     chain_means <- colMeans(draw_matrix)
     
     chain_variances <- apply(
          draw_matrix,
          2,
          var
     )
     
     # Variabilidad dentro y entre cadenas
     W         <- mean(chain_variances)
     B_between <- n_draw * var(chain_means)
     
     if (
          !is.finite(W) ||
          W <= 0 ||
          !is.finite(B_between)
     ) {
          return(NA_real_)
     }
     
     variance_hat <-
          ((n_draw - 1) / n_draw) * W +
          B_between / n_draw
     
     sqrt(variance_hat / W)
}

summarize_variable <- function(
          chain_results,
          variable,
          elapsed_total
) {
     # Extraer las cadenas de la cantidad analizada
     chain_draws <- lapply(
          chain_results,
          function(result) result[[variable]]
     )
     
     combined_draws <- unlist(
          chain_draws,
          use.names = FALSE
     )
     
     # Sumar los tamaños efectivos de las cadenas
     ess_by_chain <- vapply(
          chain_draws,
          ess_one_chain,
          numeric(1)
     )
     
     ess_total    <- sum(ess_by_chain)
     posterior_sd <- sd(combined_draws)
     
     posterior_quantiles <- quantile(
          combined_draws,
          probs = c(0.025, 0.500, 0.975),
          names = FALSE
     )
     
     mcse <- posterior_sd / sqrt(ess_total)
     
     ess_second <- if (elapsed_total > 0) {
          ess_total / elapsed_total
     } else {
          NA_real_
     }
     
     c(
          mean       = mean(combined_draws),
          sd         = posterior_sd,
          q025       = posterior_quantiles[1],
          q500       = posterior_quantiles[2],
          q975       = posterior_quantiles[3],
          rhat       = rhat_basic(chain_draws),
          ess        = ess_total,
          mcse       = mcse,
          elapsed    = elapsed_total,
          ess_second = ess_second
     )
}

# Etiquetas de las cantidades analizadas
variable_labels <- c(
     alpha  = "alpha",
     K      = "K_n",
     loglik = "log L_c"
)

# Almacenamiento de los resúmenes
n_rows <- length(samplers) * length(variable_labels)

summary_rows <- vector(
     mode   = "list",
     length = n_rows
)

row_index <- 0L

# Resumir cada muestreador y cantidad
for (sampler_name in names(samplers)) {
     chain_results <- samplers[[sampler_name]]
     
     elapsed_total <- sum(
          vapply(
               chain_results,
               function(result) result$elapsed,
               numeric(1)
          )
     )
     
     for (variable in names(variable_labels)) {
          row_index <- row_index + 1L
          
          stats <- summarize_variable(
               chain_results = chain_results,
               variable      = variable,
               elapsed_total = elapsed_total
          )
          
          summary_rows[[row_index]] <- data.frame(
               Muestreador = sampler_name,
               Cantidad    = unname(variable_labels[variable]),
               Media       = unname(stats["mean"]),
               DE          = unname(stats["sd"]),
               Q025        = unname(stats["q025"]),
               Q500        = unname(stats["q500"]),
               Q975        = unname(stats["q975"]),
               Rhat        = unname(stats["rhat"]),
               ESS         = unname(stats["ess"]),
               MCSE        = unname(stats["mcse"]),
               Tiempo      = unname(stats["elapsed"]),
               ESSseg      = unname(stats["ess_second"]),
               stringsAsFactors = FALSE
          )
     }
}

# Tabla final de resultados
summary_table <- do.call(
     rbind,
     summary_rows
)

rownames(summary_table) <- NULL

print(summary_table)

# Combinar los diagnósticos de todas las cadenas
occupancy_all <- do.call(
     rbind,
     lapply(
          truncated_blocked,
          function(result) result$occupancy
     )
)

last_weight_all <- unlist(
     lapply(
          truncated_blocked,
          function(result) result$last_weight
     ),
     use.names = FALSE
)

max_label_all <- unlist(
     lapply(
          truncated_blocked,
          function(result) result$max_label
     ),
     use.names = FALSE
)

# Probabilidad posterior de ocupación de cada componente
occupancy_prob <- colMeans(occupancy_all)

# Últimos tres componentes de la representación truncada
n_tail <- min(3L, H)

tail_index <- seq.int(
     from = H - n_tail + 1L,
     to   = H
)

tail_occupied <- rowSums(
     occupancy_all[, tail_index, drop = FALSE]
) > 0L

# Resumen de los diagnósticos del truncamiento
truncation_table <- data.frame(
     H                  = H,
     Media_peso_final   = mean(last_weight_all),
     Q95_peso_final     = unname(
          quantile(
               last_weight_all,
               probs = 0.95
          )
     ),
     Prob_ultimo_ocupado = occupancy_prob[H],
     Prob_algun_ultimo_tres = mean(tail_occupied),
     Q95_indice_maximo = unname(
          quantile(
               max_label_all,
               probs = 0.95
          )
     )
)

print(truncation_table)

# Figuras ----------------------------------------------------------------------

# Combinar las extracciones de todas las cadenas
alpha_box <- lapply(
     samplers,
     function(chains) {
          unlist(
               lapply(
                    chains,
                    function(chain) chain$alpha
               ),
               use.names = FALSE
          )
     }
)

K_box <- lapply(
     samplers,
     function(chains) {
          unlist(
               lapply(
                    chains,
                    function(chain) chain$K
               ),
               use.names = FALSE
          )
     }
)

# Distribución posterior de alpha
pdf(
     file      = "dpm_simulacion_muestrador_alpha_comparacion.pdf",
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
     alpha_box,
     names = c(
          "Completo",
          "Parcial",
          "Truncado"
     ),
     las        = 1,
     ylab       = expression(alpha),
     main       = "",
     col        = "grey90",
     border     = "grey30",
     boxwex     = 0.55,
     staplewex  = 0.6,
     medlwd     = 2
)

box()

dev.off()

# Distribución posterior de K_n
pdf(
     file      = "dpm_simulacion_muestreador_K_comparacion.pdf",
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
     K_box,
     names = c(
          "Completo",
          "Parcial",
          "Truncado"
     ),
     las       = 1,
     ylab      = expression(K[n]),
     main      = "",
     col       = "grey90",
     border    = "grey30",
     boxwex    = 0.55,
     staplewex = 0.6,
     medlwd    = 2
)

box()

dev.off()

# Probabilidad posterior de ocupación
pdf(
     file      = "dpm_simulacion_truncacion_diagnostico.pdf",
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
     x    = seq_len(H),
     y    = occupancy_prob,
     type = "h",
     lwd  = 4,
     xlab = "Índice del componente",
     ylab = "Probabilidad posterior",
     ylim = c(0, 1)
)

box()

dev.off()

# Fin --------------------------------------------------------------------------