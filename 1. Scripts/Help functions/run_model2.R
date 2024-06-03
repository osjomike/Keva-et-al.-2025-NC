run_model2<-function (run, mix, source, discr, model_filename, alpha.prior = 1, 
          resid_err = NULL, process_err = NULL) 
{
  err_raw <- read.table(model_filename, comment.char = "", 
                        sep = ":", skip = 7, nrows = 1, colClasses = "character")
  if (err_raw[1, 2] == " Residual only") 
    err <- "resid"
  if (err_raw[1, 2] == " Process only (MixSIR, for N = 1)") 
    err <- "process"
  if (err_raw[1, 2] == " Residual * Process") 
    err <- "mult"
  if (length(alpha.prior) == 1) {
    if (alpha.prior == 1) 
      alpha.prior = rep(1, source$n.sources)
  }
  if (!is.numeric(alpha.prior)) {
    stop(paste("*** Error: Your prior is not a numeric vector of length(n.sources).\n        Try again or choose the uninformative prior option. For example,\n        c(1,1,1,1) is a valid (uninformative) prior for 4 sources. ***", 
               sep = ""))
  }
  if (length(alpha.prior) != source$n.sources) {
    stop(paste("*** Error: Length of your prior does not match the\n        number of sources (", 
               source$n.sources, "). Try again. ***", sep = ""))
  }
  if (length(which(alpha.prior == 0)) != 0) {
    stop(paste("*** Error: You cannot set any alpha = 0.\n      Instead, set = 0.01.***", 
               sep = ""))
  }
  if (is.numeric(alpha.prior) == F) 
    alpha.prior = 1
  if (length(alpha.prior) == 1) 
    alpha = rep(alpha.prior, source$n.sources)
  if (length(alpha.prior) > 1 & length(alpha.prior) != source$n.sources) 
    alpha = rep(1, source$n.sources)
  if (length(alpha.prior) > 1 & length(alpha.prior) == source$n.sources) 
    alpha = alpha.prior
  if (is.list(run)) {
    mcmc <- run
  }
  else {
    if (run == "test") 
      mcmc <- list(chainLength = 1000, burn = 500, thin = 1, 
                   chains = 3, calcDIC = TRUE)
    if (run == "very short") 
      mcmc <- list(chainLength = 10000, burn = 5000, thin = 5, 
                   chains = 3, calcDIC = TRUE)
    if (run == "short") 
      mcmc <- list(chainLength = 50000, burn = 25000, thin = 25, 
                   chains = 3, calcDIC = TRUE)
    if (run == "normal") 
      mcmc <- list(chainLength = 1e+05, burn = 50000, thin = 50, 
                   chains = 3, calcDIC = TRUE)
    if (run == "long") 
      mcmc <- list(chainLength = 3e+05, burn = 2e+05, thin = 100, 
                   chains = 3, calcDIC = TRUE)
    if (run == "very long") 
      mcmc <- list(chainLength = 1e+06, burn = 5e+05, thin = 500, 
                   chains = 3, calcDIC = TRUE)
    if (run == "extreme") 
      mcmc <- list(chainLength = 3e+06, burn = 1500000, 
                   thin = 500, chains = 3, calcDIC = TRUE)
  }
  n.sources <- source$n.sources
  N <- mix$N
  e <- matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, 
              ncol = (n.sources - 1))
  for (i in 1:(n.sources - 1)) {
    e[, i] <- exp(c(rep(sqrt(1/(i * (i + 1))), i), -sqrt(i/(i + 
                                                              1)), rep(0, n.sources - i - 1)))
    e[, i] <- e[, i]/sum(e[, i])
  }
  cross <- array(data = NA, dim = c(N, n.sources, n.sources - 
                                      1))
  tmp.p <- array(data = NA, dim = c(N, n.sources))
  jags.params <- c("p.global", "loglik")
  f.data <- character(0)
  if (mix$n.effects > 0 & !mix$fere) {
    factor1_levels <- mix$FAC[[1]]$levels
    Factor.1 <- mix$FAC[[1]]$values
    cross.fac1 <- array(data = NA, dim = c(factor1_levels, 
                                           n.sources, n.sources - 1))
    tmp.p.fac1 <- array(data = NA, dim = c(factor1_levels, 
                                           n.sources))
    if (mix$FAC[[1]]$re) 
      jags.params <- c(jags.params, "p.fac1", "ilr.fac1", 
                       "fac1.sig")
    else jags.params <- c(jags.params, "p.fac1", "ilr.fac1")
    f.data <- c(f.data, "factor1_levels", "Factor.1", "cross.fac1", 
                "tmp.p.fac1")
  }
  if (mix$n.effects > 1 & !mix$fere) {
    factor2_levels <- mix$FAC[[2]]$levels
    Factor.2 <- mix$FAC[[2]]$values
    if (mix$fac_nested[1]) {
      factor2_lookup <- mix$FAC[[1]]$lookup
      f.data <- c(f.data, "factor2_lookup")
    }
    if (mix$fac_nested[2]) {
      factor1_lookup <- mix$FAC[[2]]$lookup
      f.data <- c(f.data, "factor1_lookup")
    }
    cross.fac2 <- array(data = NA, dim = c(factor2_levels, 
                                           n.sources, n.sources - 1))
    tmp.p.fac2 <- array(data = NA, dim = c(factor2_levels, 
                                           n.sources))
    if (mix$FAC[[2]]$re) 
      jags.params <- c(jags.params, "p.fac2", "ilr.fac2", 
                       "fac2.sig")
    else jags.params <- c(jags.params, "p.fac2", "ilr.fac2")
    f.data <- c(f.data, "factor2_levels", "Factor.2", "cross.fac2", 
                "tmp.p.fac2")
  }
  if (mix$fere) {
    factor1_levels <- mix$FAC[[1]]$levels
    Factor.1 <- mix$FAC[[1]]$values
    if (mix$n.re == 1) {
      cross.fac1 <- array(data = NA, dim = c(factor1_levels, 
                                             n.sources, n.sources - 1))
      tmp.p.fac1 <- array(data = NA, dim = c(factor1_levels, 
                                             n.sources))
      jags.params <- c(jags.params, "p.fac1")
      f.data <- c(f.data, "cross.fac1", "tmp.p.fac1")
    }
    factor2_levels <- mix$FAC[[2]]$levels
    Factor.2 <- mix$FAC[[2]]$values
    if (mix$FAC[[2]]$re) 
      jags.params <- c(jags.params, "fac2.sig")
    jags.params <- c(jags.params, "ilr.global", "ilr.fac1", 
                     "ilr.fac2")
    f.data <- c(f.data, "factor1_levels", "Factor.1", "factor2_levels", 
                "Factor.2")
  }
  if (source$data_type == "raw") {
    SOURCE_array <- source$SOURCE_array
    n.rep <- source$n.rep
    s.data <- c("SOURCE_array", "n.rep")
  }
  else {
    MU_array <- source$MU_array
    SIG2_array <- source$SIG2_array
    n_array <- source$n_array
    s.data <- c("MU_array", "SIG2_array", "n_array")
  }
  if (!is.na(source$by_factor)) {
    source_factor_levels <- source$S_factor_levels
    s.data <- c(s.data, "source_factor_levels")
  }
  if (source$conc_dep) {
    conc <- source$conc
    s.data <- c(s.data, "conc")
  }
  c.data <- rep(NA, mix$n.ce)
  if (mix$n.ce > 0) {
    for (ce in 1:mix$n.ce) {
      name <- paste("Cont.", ce, sep = "")
      assign(name, as.vector(mix$CE[[ce]]))
      c.data[ce] <- paste("Cont.", ce, sep = "")
      jags.params <- c(jags.params, "ilr.global", paste("ilr.cont", 
                                                        ce, sep = ""), "p.ind")
    }
  }
  X_iso <- mix$data_iso
  n.iso <- mix$n.iso
  frac_mu <- discr$mu
  frac_sig2 <- discr$sig2
  all.data <- c("X_iso", "N", "n.sources", "n.iso", "alpha", 
                "frac_mu", "e", "cross", "tmp.p")
  jags.data <- c(all.data, f.data, s.data, c.data)
  I <- diag(n.iso)
  if (err == "resid" && mix$n.iso > 1) 
    jags.data <- c(jags.data, "I")
  if (err != "resid") 
    jags.data <- c(jags.data, "frac_sig2")
  if (err == "mult") 
    jags.params <- c(jags.params, "resid.prop")
  jags.inits <- function() {
    list(p.global = as.vector(MCMCpack::rdirichlet(1, alpha)))
  }
  for (j in 1:n.iso) {
    if (source$data_type == "raw") {
      if (!is.na(source$by_factor)) {
        mean.pool <- mean(c(X_iso[, j], as.vector(SOURCE_array[, 
                                                               j, , ])), na.rm = T)
        sd.pool <- sd(c(X_iso[, j], as.vector(SOURCE_array[, 
                                                           j, , ])), na.rm = T)
        SOURCE_array[, j, , ] <- (SOURCE_array[, j, , 
        ] - mean.pool)/sd.pool
      }
      else {
        mean.pool <- mean(c(X_iso[, j], as.vector(SOURCE_array[, 
                                                               j, ])), na.rm = T)
        sd.pool <- sd(c(X_iso[, j], as.vector(SOURCE_array[, 
                                                           j, ])), na.rm = T)
        SOURCE_array[, j, ] <- (SOURCE_array[, j, ] - 
                                  mean.pool)/sd.pool
      }
    }
    else {
      if (!is.na(source$by_factor)) {
        mean.pool <- (N * mean(X_iso[, j], na.rm = T) + 
                        as.vector(as.vector(n_array) %*% as.vector(MU_array[, 
                                                                            j, ])))/sum(c(as.vector(n_array), N))
        if (N > 1) 
          sd.pool <- sqrt((sum((as.vector(n_array) - 
                                  1) * as.vector(SIG2_array[, j, ])) + as.vector(as.vector(n_array) %*% 
                                                                                   as.vector(MU_array[, j, ])^2) + (N - 1) * 
                             stats::var(X_iso[, j], na.rm = T) + N * mean(X_iso[, 
                                                                                j], na.rm = T)^2 - sum(c(as.vector(n_array), 
                                                                                                         N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                                                                                                                                   N)) - 1))
        if (N == 1) 
          sd.pool <- sqrt((sum((as.vector(n_array) - 
                                  1) * as.vector(SIG2_array[, j, ])) + as.vector(as.vector(n_array) %*% 
                                                                                   as.vector(MU_array[, j, ])^2) + N * mean(X_iso[, 
                                                                                                                                  j], na.rm = T)^2 - sum(c(as.vector(n_array), 
                                                                                                                                                           N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                                                                                                                                                                                     N)) - 1))
        MU_array[, j, ] <- (MU_array[, j, ] - mean.pool)/sd.pool
        SIG2_array[, j, ] <- SIG2_array[, j, ]/sd.pool^2
      }
      else {
        mean.pool <- (N * mean(X_iso[, j], na.rm = T) + 
                        as.vector(as.vector(n_array) %*% as.vector(MU_array[, 
                                                                            j])))/sum(c(as.vector(n_array), N))
        if (N > 1) 
          sd.pool <- sqrt((sum((as.vector(n_array) - 
                                  1) * as.vector(SIG2_array[, j])) + as.vector(as.vector(n_array) %*% 
                                                                                 as.vector(MU_array[, j])^2) + (N - 1) * stats::var(X_iso[, 
                                                                                                                                          j], na.rm = T) + N * mean(X_iso[, j], na.rm = T)^2 - 
                             sum(c(as.vector(n_array), N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                                                                                 N)) - 1))
        if (N == 1) 
          sd.pool <- sqrt((sum((as.vector(n_array) - 
                                  1) * as.vector(SIG2_array[, j])) + as.vector(as.vector(n_array) %*% 
                                                                                 as.vector(MU_array[, j])^2) + N * mean(X_iso[, 
                                                                                                                              j], na.rm = T)^2 - sum(c(as.vector(n_array), 
                                                                                                                                                       N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                                                                                                                                                                                 N)) - 1))
        MU_array[, j] <- (MU_array[, j] - mean.pool)/sd.pool
        SIG2_array[, j] <- SIG2_array[, j]/sd.pool^2
      }
    }
    X_iso[, j] <- (X_iso[, j] - mean.pool)/sd.pool
    frac_mu[, j] <- frac_mu[, j]/sd.pool
    frac_sig2[, j] <- frac_sig2[, j]/sd.pool^2
  }
  
  jags.params <- c(jags.params, "ilr.cont1.r", "cont1.sig")

    jags.1 <- R2jags::jags(jags.data, inits = jags.inits, parameters.to.save = jags.params, 
        model.file = model_filename, n.chains = mcmc$chains, 
        n.burnin = mcmc$burn, n.thin = mcmc$thin, n.iter = mcmc$chainLength, 
        DIC = mcmc$calcDIC)
  return(jags.1)
}