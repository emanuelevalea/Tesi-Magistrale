################
### 1. Setup ###
################

### 1.1 Libraries ###
library(ReplicationSuccess) # Data
library(ggplot2) # Plots
library(pwr) # Small Telescopes Approach
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)
library(cmdstanr) # STAN
library(posterior)
library(stringr)

options(mc.cores = parallel::detectCores()) # Useful for the FB part

### 1.2 Upload the data ###

data("RProjects", package = "ReplicationSuccess")
data <- RProjects %>%
  mutate(
    zo = fiso / se_fiso,
    zr = fisr / se_fisr,
    c  = (se_fiso^2) / (se_fisr^2),
    d  = fisr / fiso
  ) %>%
  filter(!is.na(zo), !is.na(d), !is.na(c))
summary(data)
help(RProjects)

#################
### Chapter 2 ###
#################

### 2.1 Frequentist Methods ###
## 2.1.1 Statistical Significance Criterion ##
idx_sign <- data$po < 0.05 # Statistical Significant Studies
data$ssc[idx_sign] <- (sign(data$fiso[idx_sign]) == sign(data$fisr[idx_sign])) & (data$pr[idx_sign]<0.05) # SSC

## 2.1.2 Meta-Analysis Statistical Significance Criterion ## 
data$se_fismeta[idx_sign] <- sqrt(1 / ( (1 / data$se_fiso[idx_sign]^2) + (1 / data$se_fisr[idx_sign]^2) )) # Meta-Analitic variance
data$fismeta[idx_sign] <- (data$fiso[idx_sign]*(1/data$se_fiso[idx_sign]^2) + data$fisr[idx_sign]*(1/data$se_fisr[idx_sign]^2)) /
  ( (1/data$se_fiso[idx_sign]^2) + (1/data$se_fisr[idx_sign]^2) ) # Meta-Analitic Estimate
data$pmeta[idx_sign] <- 2 * (1 - pnorm(abs(data$fismeta[idx_sign]/data$se_fismeta[idx_sign]))) # Meta-Analitic p-value
data$ssc_meta[idx_sign] <- (sign(data$fiso[idx_sign]) == sign(data$fismeta[idx_sign])) & (data$pmeta[idx_sign]<0.05) # SSC using the Meta-Analitic effect 

## 2.1.3 Prediction Interval Criterion ##
data$pi_inf <- data$fiso - sqrt(1/(data$no-3)+1/(data$nr-3))
data$pi_sup <- data$fiso + sqrt(1/(data$no-3)+1/(data$nr-3))
data$pic <- (data$fisr>=data$pi_inf) & (data$fisr<=data$pi_sup)

## 2.1.4 Small Telescope Approach ##
d33 <- function(n) pwr.r.test(power=1/3,n=n, sig.level = 0.05)$r # Evaluate d33
data$r_sta <- sapply(data$no, d33)
data$fisd33 <- atanh(data$r_sta) # Trasformo in fisher z

data$z <- (data$fisr - data$fisd33)/data$se_fisr
data$psta <- pnorm(data$z)
data$sta <- (data$psta <= 0.05)

### 2.2 Bayesian Methods ###
## 2.2.1 Default Bayes Factor ##
BF10 <- function(theta_o, theta_r, se_r, r = sqrt(2)/2) {
  if (theta_o > 0){
    lh_H1 <- integrate(
      function(theta) {
        dnorm(theta_r, mean = theta, sd = se_r) * 2 * dcauchy(theta, location = 0, scale = r)
      },lower = 0, upper = Inf)$value
  } else{
    lh_H1 <- integrate(
      function(theta) {
        dnorm(theta_r, mean = theta, sd = se_r) * 2 * dcauchy(theta, location = 0, scale = r)
      },lower = -Inf, upper = 0)$value
  }
  BF <- lh_H1 / dnorm(theta_r, mean = 0, sd = se_r)
  return(BF)
}
data$BF_def <- rep(NA, nrow(RProjects))
for (i in 1:nrow(RProjects)){
  data$BF_def <- BF10(RProjects$ro[i], RProjects$fisr[i], RProjects$se_fisr[i])
}

## 2.2.2 Replication Bayes Factor ##
BF0A <- function(zo, d, c){ # Pawel-Held)
  BF <- sqrt(1+c)*exp((-(zo^2)/2)*((d^2)*c-((1-d)^2)/(1/c+1)))
  return(BF)
}
data$BFR <- rep(NA, nrow(RProjects))
for (i in 1:nrow(RProjects)){
  data$BFR[i] <- BF0A(data$zo[i], data$d[i], data$c[i])
}

## 2.2.3 Sceptical Bayes Factor ##
BF0S <- function(zo, g) {
  res <- sqrt(1 + g) * exp(-0.5 * (g/(1 + g)) * zo^2)
  return(res)
}

BFSA <- function(zo, d, c, g) {
  res <- sqrt((1/c + 1)/(1/c + g)) * exp(-0.5 * zo^2 * ((d^2)/(1/c + g) - ((d - 1)^2)/(1/c + 1)))
  return(res)
}

BFS <- function(zo, d, c, g_seq = seq(0.01, 20000, by = 0.1)) {
  
  # Evaluate BF0S, BFSA for each g
  BF0S_vals <- BF0S(zo, g_seq)
  BFSA_vals <- BFSA(zo, d, c, g_seq)
  
  # 3 possible cases #
  # (a): BFSA > BF0S always -> BFS undefined 
  if (all(BFSA_vals > BF0S_vals)) {
    return(list(g_gamma = NA, bfs = NA, case = "a"))
  }
  
  # (b): BFSA < BF0S always -> BFS = min{BF0S}
  if (all(BFSA_vals < BF0S_vals)) {
    g_gamma <- optimize(function(g) BF0S(zo, g), 
                        interval = c(0, max(g_seq)))$minimum
    BFS_value <- BF0S(zo, g_gamma)
    
    return(list(g_gamma = g_gamma, bfs = BFS_value, case = "b"))
  }
  
  # (c): intersection
  intersection <- function(g) {BF0S(zo, g) - BFSA(zo, d, c, g)}
  upper <- max(g_seq)  # Research Range
  lower <- 0
  while ((intersection(lower) * intersection(upper) > 0) & (upper > 1)) {
    upper <- upper - 100
  }
  
  # if f(lower)*f(upper)>0, then there is no x|f(x)=0 - > no solution
  if (intersection(lower) * intersection(upper) > 0) {
    return(list(g_gamma = NA, bfs = NA, case = "no_intersection"))
  } 
  else{
    g_gamma <- uniroot(intersection, lower = lower, upper = upper)$root
    BFS_value <- BF0S(zo, g_gamma)
    return(list(g_gamma = g_gamma, bfs = BFS_value, case = "c"))
  }
}

data$g_gamma <- rep(NA,nrow(RProjects))
data$BFS_value <- rep(NA,nrow(RProjects))
data$BFS_case <- rep(NA,nrow(RProjects))

for (i in 1:nrow(data)){
  res <- BFS(data$zo[i], data$d[i], data$c[i])
  data$g_gamma[i] <- res$g_gamma
  data$BFS_value[i]  <- res$bfs
  data$BFS_case[i]   <- res$case
}

## 2.2.4 Sceptical Mixture Bayes Factor ##
BF0SM <- function(zo, g, psi){
  if (is.na(psi)) return(NA)
  else{ 
    res <- (psi + (1 - psi) * (BF0S(zo, g))^-1)^-1 # Egidi, Consonni (2025)
    return (res)}
}

BFSMA <- function(zo, d, c, g, psi){
  if (is.na(psi)) return(NA)
  else{
    res <- psi * BF0A(zo, d, c) + (1 - psi) * BFSA(zo, d, c, g)  # Egidi, Consonni (2025)
    return(res)}
}

get_psi <- function(zo, h, alpha = 0.05) {
  
  psi <- (alpha - (1 - pchisq(zo^2 / (1 + h), 1))) / 
    ((1 - pchisq(zo^2, 1)) - (1 - pchisq(zo^2 / (1 + h), 1)))
  
  psi <- pmax(0, pmin(1, psi))
  
  psi <- round(psi,3)
  return(psi)
}

BFSM <- function(zo, d, c, h_seq = seq(0.01, 2000, by = 0.01), alpha = 0.05){
  
  h_star <- (zo^2)/(qchisq(1-alpha,1)) - 1
  if(h_star < 0) h_star <- 0
  h_seq_valid <- h_seq[h_seq>=h_star]
  
  psi_seq_valid <- get_psi(zo, h_seq_valid, alpha)
  
  bf_osm <- sapply(seq_along(h_seq_valid), function(i) BF0SM(zo, h_seq_valid[i], psi_seq_valid[i]))
  bf_sma <- sapply(seq_along(h_seq_valid), function(i) BFSMA(zo, d, c, h_seq_valid[i], psi_seq_valid[i]))
  
  # (a)
  if(all(bf_sma > bf_osm, na.rm = TRUE)){
    return(list(h_gamma = NA, bfsm = NA, psi = NA, case = "a"))
  }
  
  # (b)
  else if(all(bf_sma < bf_osm, na.rm = TRUE)){
    sol_min_bf0 <- function(h){
      psi <- get_psi(zo, h, alpha)
      return(BF0SM(zo, h, psi))
    } 
    opt <- optim(c(max(h_star,0.01)),  sol_min_bf0, method = "L-BFGS-B", lower = c(max(h_star,0.01)), upper = max(h_seq))
    h_gamma <- opt$par
    psi_star <- get_psi(zo, h_gamma, alpha)
    return(list(h_gamma = h_gamma, bfsm = opt$value, psi = psi_star, case = "b"))
  }
  else{
    # (c)
    intersection <- function(h){
      psi <- get_psi(zo, h, alpha)
      return(BF0SM(zo, h, psi) - BFSMA(zo, d, c, h, psi))
    }
    upper <- max(h_seq)
    while (intersection(h_star+1e-02)*intersection(upper)>0){
      upper <- upper-1
    }
    h_gamma <- uniroot(intersection, lower = h_star+1e-02, upper = upper)$root
    psi_star <- get_psi(zo, h_gamma, alpha)
    bfsm_val <- BF0SM(zo, h_gamma, psi_star)
    return(list(h_gamma = h_gamma, bfsm = bfsm_val, psi = psi_star, case = "c"))
  }
}

sceptical_mixture_p_value <-  function(zo, h_gamma, psi){
  psm <- psi*(1-pchisq(zo^2, 1))+(1-psi)*(1-pchisq((zo)^2/(1+h_gamma), 1))
  return(psm)
}

data$g_gamma_mix <- rep(NA, nrow(RProjects))
data$psi_gamma   <- rep(NA, nrow(RProjects))
data$BFSM_value  <- rep(NA, nrow(RProjects))
data$BFSM_case   <- rep(NA, nrow(RProjects))
data$p_sm <- rep(NA, nrow(RProjects))

for (i in 1:nrow(data)){
  res <- BFSM(data$zo[i], data$d[i], data$c[i]) # alpha = 0.05
  data$g_gamma_mix[i] <- res$h_gamma
  data$psi_gamma[i]   <- res$psi
  data$BFSM_value[i]  <- res$bfsm
  data$BFSM_case[i]   <- res$case
  data$psm[i] <- sceptical_mixture_p_value(data$zo[i], res$h_gamma, res$psi)
}

## 2.2.5 Sceptical p-value ##
sceptical_p_value <- function(zo, zr, c = 1){
  z2_h <- 2/(1/zo^2 + 1/zr^2)
  z2_a <- (zo^2 + zr^2)/2
  if (c==1){
    z2_s <- z2_h/2
  } else {
    z2_s <- (1/(c-1))*sqrt(z2_a + (c-1)*z2_h) - z2_a
  }
  z_s <- sqrt(z2_a)
  ps <- 1 - pnorm(z_s)
  ps_rec <- 1 - pnorm(z_s * sqrt((1+sqrt(5))/2))
  return(list(sceptical_p = ps, recalibrated_sceptical_p = ps_rec))
}

data$ps <- rep(NA, nrow(RProjects))
data$ps_rec <- rep(NA, nrow(RProjects))

for (i in 1:nrow(data)){
  res <- sceptical_p_value(data$zo[i], data$zr[i], data$c[i])
  data$ps[i] <- res$sceptical_p
  data$ps_rec[i] <- res$recalibrated_sceptical_p
}

### 2.3 Summary ###
data_results <- data.frame(
  Statistical_Significance      = data$ssc,
  SSC_Meta_Analysis              = data$ssc_meta,
  Prediction_Interval            = data$pic,
  Small_Telescope                = data$sta,
  BF_Default                     = data$BF_def,
  BF_Replication                 = data$BFR,
  BF_Sceptical                   = data$BFS_value,
  Sceptical_p_value               = data$ps,
  Recalibrated_Sceptical_p_value  = data$ps_rec,
  BF_Sceptical_Mixture            = data$BFSM_value,
  Sceptical_Mixture_p_value       = data$psm
)

#################
### Chapter 3 ###
#################

### 3.1 KL Calibration ###
KL_gamma <- function(zo, d, c, gamma) {
  0.5 * (log(gamma * c) + 1/(gamma * c) + zo^2 * (gamma - d)^2 / gamma - 1)
}
gamma_star_KL <- function(zo, d, c) {
  A  <- zo^2
  B  <- 1
  Cc <- -1/c - (d^2)*(zo^2)
  root <- (-B + sqrt(B^2 - 4*A*Cc)) / (2*A)
  pmin(1, pmax(0, root))
}
data$gamma_KL <- gamma_star_KL(data$zo, data$d, data$c)

### 3.2 Empirical Bayes ###
l_gamma <- function(zo, d, c, gamma) {
  -0.5 * log(1 + c * gamma) - 0.5 * c * zo^2 * (d - gamma)^2 / (1 + c * gamma)
}
gamma_star_EB <- function(zo, d, c) {
  A  <- zo^2 * c
  B  <- c + 2*zo^2
  Cc <- 1 - 2*zo^2*d - zo^2*c*d^2
  root <- (-B + sqrt(B^2 - 4*A*Cc)) / (2*A)
  pmin(1, pmax(0, root))
}
data$gamma_EB <- gamma_star_EB(data$zo, data$d, data$c)

### 3.3 Full Bayes ###

## 3.3.1 Model Fit ##
model <- cmdstan_model("Mod_Gamma.stan")
fit_prior <- function(a, b, label, data, model){
  stan_data <- list(
    N = nrow(data),
    theta_o = data$fiso,
    theta_r = data$fisr,
    var_o = data$se_fiso^2,
    var_r = data$se_fisr^2,
    a = a,
    b = b
  )
  fit <- model$sample(
    data = stan_data,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = .95,
    refresh = 0
  )
  summary_gamma <- fit$summary(
    variables = "gamma",
    mean, median, sd,
    q025 = ~quantile(.x, .025),
    q975 = ~quantile(.x, .975)
  ) %>%
    rename(
      Mean   = mean,
      Median = median,
      SD     = sd,
      Q025   = `2.5%`,
      Q975   = `97.5%`
    ) %>%
    mutate(
      Study = row_number(),
      Prior = label
    )
  list(fit = fit, summary = summary_gamma)
}

fit11 <- fit_prior(a = 1, b = 1, label = "Beta(1,1)", data = data, model = model)
fit22 <- fit_prior(a = 2, b = 2, label = "Beta(2,2)", data = data, model = model)
fit25 <- fit_prior(a = 2, b = 5, label = "Beta(2,5)", data = data, model = model)
fit52 <- fit_prior(a = 5, b = 2, label = "Beta(5,2)", data = data, model = model)

fits        <- list(fit11, fit22, fit25, fit52)
prior_names <- c("Beta(1,1)", "Beta(2,2)", "Beta(2,5)", "Beta(5,2)")

data$gamma_FB <- fits[[4]]$summary$Mean # Pick Beta(5,2)

## 3.3.2 Diagnostics ##
diagnostics <- purrr::map_dfr(seq_along(fits), function(i){
  s  <- fits[[i]]$fit$summary("gamma")
  ds <- fits[[i]]$fit$diagnostic_summary()
  tibble(
    Prior           = prior_names[i],
    n_divergent     = sum(ds$num_divergent),
    n_max_treedepth = sum(ds$num_max_treedepth),
    min_ebfmi       = min(ds$ebfmi),
    max_rhat        = max(s$rhat, na.rm = TRUE),
    min_ess_bulk    = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail    = min(s$ess_tail, na.rm = TRUE)
  )
})

### 3.3 Datasets for plots ###
lookup <- data %>%
  mutate(Study = row_number()) %>%
  select(Study, study, project, zo, d, c)

summary_long <- bind_rows(lapply(fits, function(x) x$summary)) %>%
  mutate(
    CI_width = Q975 - Q025,
    Prior = factor(Prior, levels = prior_names)
  ) %>%
  left_join(lookup, by = "Study")

projects <- sort(unique(data$project))

df_forest <- summary_long %>%
  filter(Prior == "Beta(1,1)") %>%
  group_by(project) %>%
  arrange(Median, .by_group = TRUE) %>%
  ungroup()

selected_names <- c(
  "Aviezer", "Janssen", "Kovacs",
  "Pyc", "Gneezy", "Hauser"
)
idx_selected <- sapply(selected_names, \(x)  which(str_detect(data$study,regex(x, ignore_case = TRUE)))[1])

names(idx_selected) <- selected_names

## 3.4 Threshold Calibration ##
tab_soglie <- tibble(
  Criterio = c("KL", "Empirical Bayes", "Full Bayes (Beta(5,2))"),
  Soglia_050 = c(
    mean(data$gamma_KL >= 0.5, na.rm = TRUE),
    mean(data$gamma_EB >= 0.5, na.rm = TRUE),
    mean(data$gamma_FB >= 0.5, na.rm = TRUE)
  ),
  Soglia_067 = c(
    mean(data$gamma_KL >= 2/3, na.rm = TRUE),
    mean(data$gamma_EB >= 2/3, na.rm = TRUE),
    mean(data$gamma_FB >= 2/3, na.rm = TRUE)
  )
)

#################
### Chapter 4 ###
#################

### 4.1 Simulating the priors and posteriors ###
simulate_prior_predictive <- function(theta_o, var_o, a = 5, b = 2, n_sim = 10000){
  gamma_sim <- rbeta(n_sim, a, b) # Draw n_sim values from the Prior
  rnorm(n_sim, mean = gamma_sim * theta_o, sd = sqrt(gamma_sim * var_o)) # Generate n_sim values for the possible replication
}

theta_r_rep_Beta52 <- fits[[4]]$fit$draws(variables = "theta_r_rep", format = "matrix") # Posterior Predictive Replications for Prior Beta(5,2)
theta_r_rep_U01 <- fits[[1]]$fit$draws(variables = "theta_r_rep", format = "matrix") # Posterior Predictive Replications for Prior U(0,1)

### 4.2 Density curves for each study ###
build_ppc_curves <- function(draws, study_label, K = 20){
  draws <- draws[sample(length(draws))]
  m <- floor(length(draws) / K)
  out <- data.frame()
  for (k in 1:K){
    sub <- draws[((k - 1) * m + 1):(k * m)]
    d <- density(sub)
    out <- rbind(out, data.frame(theta = d$x, density = d$y, rep_id = k, study = study_label))
  }
  out
}

### 4.4 Dataframes for plots ###
df_ppc_lines_prior_Beta52 <- data.frame()
df_ppc_lines_post_Beta52  <- data.frame()
df_ppc_lines_prior_U01 <- data.frame()
df_ppc_lines_post_U01  <- data.frame()
df_obs_curve       <- data.frame()

for (nm in names(idx_selected)){
  idx <- idx_selected[nm]
  
  sim_prior_52 <- simulate_prior_predictive(data$fiso[idx], data$se_fiso[idx]^2, a = 5, b = 2)
  df_ppc_lines_prior_Beta52 <- rbind(df_ppc_lines_prior_Beta52, build_ppc_curves(sim_prior_52, nm))
  
  sim_prior_01 <- simulate_prior_predictive(data$fiso[idx], data$se_fiso[idx]^2, a = 1, b = 1)
  df_ppc_lines_prior_U01 <- rbind(df_ppc_lines_prior_U01, build_ppc_curves(sim_prior_01, nm))
  
  sim_post_52 <- theta_r_rep_Beta52[, paste0("theta_r_rep[", idx, "]")]
  df_ppc_lines_post_Beta52 <- rbind(df_ppc_lines_post_Beta52, build_ppc_curves(sim_post_52, nm))
  
  sim_post_01 <- theta_r_rep_U01[, paste0("theta_r_rep[", idx, "]")]
  df_ppc_lines_post_U01 <- rbind(df_ppc_lines_post_U01, build_ppc_curves(sim_post_01, nm))
  
  theta_r_obs <- data$fisr[idx]
  se_r_obs    <- data$se_fisr[idx]
  grid <- seq(theta_r_obs - 4 * se_r_obs, theta_r_obs + 4 * se_r_obs, length.out = 300)
  df_obs_curve <- rbind(df_obs_curve, data.frame(
    theta = grid, density = dnorm(grid, theta_r_obs, se_r_obs), study = nm
  ))
}
