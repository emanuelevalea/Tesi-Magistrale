#############
### Plots ###
#############

## Set Up ##
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)
library(stringr)
library(tibble)

options(scipen=999)
par(mar = c(4.2, 4.2, 1.2, 1.2))

theme_gg <- theme_minimal(base_size = 13) +
  theme(
    text = element_text(family = "serif"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.title = element_text(face = "italic")
  )

theme_gg_facet <- theme_gg +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "italic")
  )

### Fig. 1 ###
x <- seq(-5, 5, length.out = 10000)
df <- data.frame(x = x, Normal = dnorm(x, 0, 1), Cauchy_1 = dcauchy(x, 0, 1),
                 Cauchy_scaled = dcauchy(x, 0, sqrt(2)/2))
df_long <- melt(df, id.vars = "x", variable.name = "Distribuzione", value.name = "Densita")

p1 <- ggplot(df_long, aes(x, Densita, color = Distribuzione)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(Normal = "#E6AB02", Cauchy_1 = "#1B9E77", Cauchy_scaled = "#377EB8"),
                     labels = c(expression(N(0,1)), expression(Cauchy(0,1)), expression(Cauchy(0,sqrt(2)/2)))) +
  labs(x = expression(theta), y = expression(f(theta)), color = "Distribuzione") +
  theme_gg

### Fig. 2 ###
df <- data.frame(x = x, Cauchy_scaled = dcauchy(x, 0, sqrt(2)/2),
                 Cauchy_plus  = ifelse(x > 0, 2*dcauchy(x, 0, sqrt(2)/2), NA),
                 Cauchy_minus = ifelse(x < 0, 2*dcauchy(x, 0, sqrt(2)/2), NA))
df_long <- melt(df, id.vars = "x", variable.name = "Distribuzione", value.name = "Densita")

p2 <- ggplot(df_long, aes(x, Densita, color = Distribuzione)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(Cauchy_scaled = "#377EB8", Cauchy_plus = "#1B9E77", Cauchy_minus = "#E6AB02"),
                     labels = c(expression(Cauchy(0,sqrt(2)/2)), expression(Cauchy^"+"), expression(Cauchy^"-"))) +
  labs(x = expression(theta), y = expression(f(theta)), color = "Distribuzione") +
  theme_gg

### Fig. 3 ###
theta_o <- 1; se_r <- 1; d_vals <- c(0.5, 1, 2); r_grid <- seq(0.001, 3, length.out = 1000)

dfBF10 <- do.call(rbind, lapply(d_vals, function(d) {
  data.frame(d = d, r = r_grid, BF10 = sapply(r_grid, function(rv) BF10(theta_o, d*theta_o, se_r, r = rv)))
}))
dfBF10$d <- factor(dfBF10$d, labels = c("0.5", "1", "2"))

p3 <- ggplot(dfBF10, aes(r, BF10, color = d)) +
  geom_line(linewidth = 1.2) +
  scale_y_log10() +
  scale_color_manual(values = c("#E6AB02", "#1B9E77", "#377EB8")) +
  geom_vline(xintercept = sqrt(2)/2, linetype = "dashed", color = "#7B4F9D") +
  labs(x = expression(r), y = expression(BF[1:0]), color = expression(d)) +
  scale_x_continuous(breaks = c(0, sqrt(2)/2, 1, 2, 3),
                     labels = c(0, expression(frac(sqrt(2), 2)), 1, 2, 3)) +
  theme_gg

### Fig. 4 ###
zo <- 2; d_vals <- c(0.5, 1, 2); c_grid <- seq(0.0001, 4, length.out = 10000)

dfBF0A <- do.call(rbind, lapply(d_vals, function(d)
  data.frame(d = d, c = c_grid, BF0A = BF0A(zo, d, c_grid))))
dfBF0A$d <- factor(dfBF0A$d, labels = c("0.5", "1", "2"))

p4 <- ggplot(dfBF0A, aes(c, BF0A, color = d)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#E6AB02", "#1B9E77", "#377EB8")) +
  labs(x = expression(c), y = expression(BF[R]), color = expression(d)) +
  geom_hline(yintercept = 1/3, linetype = "dashed", color = "#7B4F9D") +
  scale_y_continuous(breaks = c(0, 1/3, 1), labels = c(0, expression(frac(1,3)), 1)) +
  theme_gg

### Fig. 5 ###
zo_vals <- c(1, 2, 3); g_grid <- seq(0.0001, 25, length.out = 100000)

dfBF0S <- do.call(rbind, lapply(zo_vals, function(z)
  data.frame(zo = z, g = g_grid, BF0S = BF0S(z, g_grid))))
dfBF0S$zo <- factor(dfBF0S$zo, labels = c("1", "2", "3"))

p5 <- ggplot(dfBF0S, aes(g, BF0S, color = zo)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#E6AB02", "#1B9E77", "#377EB8")) +
  labs(x = expression(g[gamma]), y = expression(BF[0:S](hat(theta)[o], g[gamma])), color = expression(abs(z[0]))) +
  scale_y_log10(breaks = c(0, 1/3, 1), labels = c(0, expression(frac(1,3)), 1)) +
  geom_hline(yintercept = 1/3, linetype = "dashed", color = "#7B4F9D") +
  theme_gg

### Fig. 6 ###
d <- 1; c <- 1; zo_vals <- c(1, 2, 3)

dfBFSA <- do.call(rbind, lapply(zo_vals, function(z)
  data.frame(zo = z, g = g_grid, BFSA = BFSA(z, d, c, g_grid))))
dfBFSA$zo <- factor(dfBFSA$zo, labels = c("1", "2", "3"))

p6 <- ggplot(dfBFSA, aes(g, BFSA, color = zo)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#E6AB02", "#1B9E77", "#377EB8")) +
  labs(x = expression(g[gamma]), y = expression(BF[S:A](hat(theta)[r], g[gamma])), color = expression(abs(z[o]))) +
  scale_y_continuous(breaks = c(0, 1/3, 1), labels = c(0, expression(frac(1,3)), 1)) +
  geom_hline(yintercept = 1/3, linetype = "dashed", color = "#7B4F9D") +
  theme_gg

### Fig. 7 ###
zo_vals <- c(1, 2, 3); d_vals <- c(0.5, 1, 2); c <- 1
g_grid7 <- seq(0.01, 5, length.out = 500)

df_all <- expand.grid(g = g_grid7, d = d_vals, zo = zo_vals) %>%
  mutate(BF0S = BF0S(zo, g), BFSA = BFSA(zo, d, c, g))

df_long <- df_all %>%
  pivot_longer(cols = c(BF0S, BFSA), names_to = "tipo", values_to = "BF") %>%
  mutate(d = factor(d, labels = c("d==0.5", "d==1", "d==2")),
         zo = factor(zo, labels = c("z[o]==1", "z[o]==2", "z[o]==3")))

intersections <- expand.grid(d = d_vals, zo = zo_vals) %>%
  rowwise() %>%
  mutate(res = list(BFS(zo, d, c, g_grid7)), g_gamma = res$g_gamma, BFS_value = res$bfs, case = res$case) %>%
  ungroup() %>%
  filter(case == "c", !is.na(g_gamma)) %>%
  mutate(d = factor(d, labels = c("d==0.5", "d==1", "d==2")),
         zo = factor(zo, labels = c("z[o]==1", "z[o]==2", "z[o]==3")))

p7 <- ggplot(df_long, aes(g, BF, color = tipo)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"), labels = c(expression(BF[0:S]), expression(BF[S:A]))) +
  geom_point(data = intersections, aes(g_gamma, BFS_value), color = "black", shape = 4, size = 2, stroke = 1,
             inherit.aes = FALSE) +
  facet_grid(d ~ zo, labeller = label_parsed) +
  scale_y_continuous(breaks = c(0, 1/3, 1), labels = c(0, expression(1/3), 1)) +
  scale_x_continuous(breaks = 0:5, labels = 0:5) +
  labs(x = expression(g[gamma]), y = expression(BF), color = NULL) +
  theme_gg_facet

### Fig. 8 ###
zo <- 3; zr <- 2.5; d <- zr/zo; c <- 1
alpha8 <- 0.05  # FIX: alpha unico (l'originale mescolava 0.01 e 0.1)

h_star8 <- (zo^2) / qchisq(1 - alpha8, 1) - 1
g_grid8 <- seq(0.0001, 25, length.out = 100000)
psi_grid <- get_psi(zo, g_grid8, alpha8)
psi_grid[g_grid8 < h_star8] <- NA
psi_grid[psi_grid == 0] <- NA

data_plot8 <- data.frame(
  g = g_grid8, bf0s = BF0S(zo, g_grid8), bfsa = BFSA(zo, d, c, g_grid8),
  bf0sm = sapply(seq_along(g_grid8), function(i) BF0SM(zo, g_grid8[i], psi_grid[i])),
  bfsma = sapply(seq_along(g_grid8), function(i) BFSMA(zo, d, c, g_grid8[i], psi_grid[i]))
)
data_long8 <- pivot_longer(data_plot8, cols = -g, names_to = "curve", values_to = "value")

point_bfs  <- which.min(abs(data_plot8$bf0s - data_plot8$bfsa))
point_bfsm <- which.min(abs(data_plot8$bf0sm - data_plot8$bfsma))
points_df8 <- data.frame(g = c(data_plot8$g[point_bfs], data_plot8$g[point_bfsm]),
                         value = c(data_plot8$bf0s[point_bfs], data_plot8$bf0sm[point_bfsm]))

colors8 <- c(bf0s = "#E6AB02", bfsa = "#D95F02", bf0sm = "#7570B3", bfsma = "#00BFC4")

p8 <- ggplot(data_long8, aes(g, value, color = curve, linetype = curve)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = points_df8, aes(g, value), inherit.aes = FALSE,
             color = c("black", "darkblue"), shape = 4, stroke = 1.2, size = 3) +
  scale_color_manual(values = colors8,
                     labels = c(bf0s = expression(BF[0*S]), bfsa = expression(BF[S*A]),
                                bf0sm = expression(BF[0*SM]), bfsma = expression(BF[SM*A]))) +
  scale_linetype_manual(values = c(bf0s = "solid", bfsa = "solid", bf0sm = "dashed", bfsma = "dashed"),
                        labels = c(bf0s = expression(BF[0*S]), bfsa = expression(BF[S*A]),
                                   bf0sm = expression(BF[0*SM]), bfsma = expression(BF[SM*A]))) +
  scale_y_log10(breaks = c(0, points_df8$value[2], points_df8$value[1], 1),
                labels = c(0, expression(BF[SM]), expression(BF[S]), 1)) +
  labs(x = "Varianza relativa", y = "Bayes Factor", color = NULL, linetype = NULL) +
  theme_gg

### Fig. 9 ###
alpha_seq9 <- seq(0, 50, length.out = 500)
df9 <- data.frame(alpha = alpha_seq9, gamma = alpha_seq9 / (alpha_seq9 + 1))

p9 <- ggplot(df9, aes(alpha, gamma)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  labs(x = expression(alpha), y = expression(delta == alpha/(alpha+1))) +
  theme_gg


### Fig. 10 ###
theta_o10 <- 2; sigma_o10 <- 1
gamma_vals10 <- c(0.1, 0.3, 0.6, 1.0)
theta_grid10 <- seq(-2, 5, length.out = 500)

df10 <- do.call(rbind, lapply(gamma_vals10, function(g)
  data.frame(theta = theta_grid10, density = dnorm(theta_grid10, g*theta_o10, sqrt(g)*sigma_o10),
             gamma = factor(g, labels = paste0("delta == ", g)))))

p10 <- ggplot(df10, aes(theta, density, color = gamma)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = theta_o10, linetype = "dashed", color = "gray40") +
  annotate("text", x = theta_o10, y = 0, label = expression(hat(theta)[o]), vjust = -1, hjust = -0.3) +
  labs(x = expression(theta), y = expression(f(theta~"|"~H[delta])), color = NULL) +
  scale_color_discrete(labels = scales::parse_format()) +
  theme_gg

save_plot(p10, "Figura 10.jpg")

### Fig. 11-16 ###
make_objective_grid <- function(obj_fun, star_fun, ylab_expr,
                                d_vals = c(0.1, 0.5, 0.9), zo_vals = c(1, 2, 3),
                                c_fix = 1, gamma_seq = seq(0.01, 1, length.out = 500)) {
  df_curves <- expand.grid(gamma = gamma_seq, d = d_vals, zo = zo_vals) %>%
    mutate(y = obj_fun(zo, d, c_fix, gamma), d = factor(d), zo = factor(zo))
  df_points <- expand.grid(d = d_vals, zo = zo_vals) %>%
    mutate(gamma_star = star_fun(zo, d, c_fix), y_star = obj_fun(zo, d, c_fix, gamma_star),
           d = factor(d), zo = factor(zo))
  
  ggplot(df_curves, aes(gamma, y)) +
    geom_line(linewidth = 0.8, color = "steelblue") +
    geom_point(data = df_points, aes(gamma_star, y_star), color = "firebrick", size = 2.2) +
    geom_vline(data = df_points, aes(xintercept = gamma_star),
               linetype = "dashed", color = "firebrick", linewidth = 0.4) +
    facet_grid(d ~ zo, scales = "free_y",
               labeller = labeller(
                 d  = as_labeller(function(x) paste0("d == ", x), label_parsed),
                 zo = as_labeller(function(x) paste0("z[o] == ", x), label_parsed))) +
    labs(x = expression(delta), y = ylab_expr) +
    theme_gg_facet
}

make_gamma_vs_d <- function(star_fun, ylab_expr, d_seq = seq(0, 1.25, length.out = 300),
                            zo_vals = c(1, 2, 3, 4), c_fix = 1) {
  df <- expand.grid(d = d_seq, zo = zo_vals) %>%
    mutate(gamma_star = star_fun(zo, d, c_fix), zo = factor(zo, labels = paste0("z[o] == ", zo_vals)))
  ggplot(df, aes(d, gamma_star, color = zo)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray50") +
    labs(x = "d", y = ylab_expr, color = NULL) +
    scale_color_discrete(labels = scales::parse_format()) +
    theme_gg
}

make_gamma_vs_c <- function(star_fun, ylab_expr, c_seq = seq(0.1, 5, length.out = 300),
                            zo_fix = 1, d_vals = c(0.1, 0.5, 0.9, 1.2)) {
  df <- expand.grid(c = c_seq, d = d_vals) %>%
    mutate(gamma_star = star_fun(zo_fix, d, c), d = factor(d, labels = paste0("d = ", d_vals)))
  ggplot(df, aes(c, gamma_star, color = d)) +
    geom_line(linewidth = 1) +
    labs(x = "c", y = ylab_expr, color = NULL) +
    theme_gg
}

p11 <- make_objective_grid(KL_gamma, gamma_star_KL, expression(KL[delta]))

p12 <- make_gamma_vs_d(gamma_star_KL, expression(hat(delta)[KL]))

p13 <- make_gamma_vs_c(gamma_star_KL, expression(hat(delta)[KL]))

p14 <- make_objective_grid(l_gamma, gamma_star_EB, expression(l(delta)))

p15 <- make_gamma_vs_d(gamma_star_EB, expression(hat(delta)[EB]))

p16 <- make_gamma_vs_c(gamma_star_EB, expression(hat(delta)[EB]))

### Fig. 17 ###
gamma_grid17 <- seq(0.001, 1, length.out = 1000)
pars17 <- tibble(alpha = c(1, 2, 5, 2), beta = c(1, 2, 2, 5),
                 prior = c("Beta(1,1)", "Beta(2,2)", "Beta(5,2)", "Beta(2,5)"))

df17 <- pars17 %>% rowwise() %>%
  reframe(gamma = gamma_grid17, density = dbeta(gamma_grid17, alpha, beta), prior = prior)

p17 <- ggplot(df17, aes(gamma, density, color = prior)) +
  geom_line(linewidth = 1.1) +
  labs(x = expression(delta), y = expression(p(delta)), color = "Prior") +
  theme_gg

### Fig. 18 ###
p18 <- ggplot(data, aes(gamma_KL, gamma_EB)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "firebrick") +
  coord_equal() +
  labs(x = expression(hat(delta)[KL]), y = expression(hat(delta)[EB])) +
  theme_gg
save_plot(p18, "Figura 18.jpg")

### Fig.19 ###
df19 <- data %>% pivot_longer(cols = c(gamma_KL, gamma_EB), names_to = "Method", values_to = "Estimate")

p19 <- ggplot(df19, aes(x = Method, y = Estimate, fill = Method)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.4, width = 0.5) +
  scale_x_discrete(labels = c("Stima EB", "Stima KL")) +
  scale_fill_discrete(labels = c("Stima EB", "Stima KL")) +
  labs(x = NULL, y = expression(hat(delta)), fill = NULL) +
  theme_gg +
  theme(legend.position = "none")
save_plot(p19, "Figura 19.jpg")

### Fig. 20 ###
p20 <- ggplot(df_forest, aes(Median, fct_reorder(factor(Study), Median))) +
  geom_errorbarh(aes(xmin = Q025, xmax = Q975), height = 0, color = "steelblue", alpha = 0.7) +
  geom_point(size = 1.2, color = "steelblue") +
  geom_vline(xintercept = c(0, 1), linetype = "dotted", color = "gray50") +
  facet_wrap(~project, scales = "free_y", ncol = 2) +
  labs(x = expression(delta), y = NULL) +
  theme_gg +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_plot(p20, "Figura 20.jpg")

### Fig. 21a - 21b
make_ciwidth <- function(proj_name) {
  df_p <- summary_long %>% filter(project == proj_name)
  ggplot(df_p, aes(Prior, CI_width, fill = Prior)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
    labs(x = NULL, y = expression(delta), title = proj_name) +
    theme_gg +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
}

p21 <- (make_ciwidth(projects[1]) + make_ciwidth(projects[2])) /
  (make_ciwidth(projects[3]) + make_ciwidth(projects[4])) +
  plot_layout(axis_titles = "collect")

### Fig. 22 - 23 ###
gamma_draws_mat <- fits[[4]]$fit$draws(
  variables = "gamma",
  format = "matrix"
)

compose_theta_posterior <- function(i, n_draws = 10000) {
  vo <- data$se_fiso[i]^2
  vr <- data$se_fisr[i]^2
  g <- sample(gamma_draws_mat[, paste0("gamma[", i, "]")],n_draws,replace = TRUE)
  vp <- 1 / (1 / (g * vo) + 1 / vr)
  mp <- vp * (data$fiso[i] / vo + data$fisr[i] / vr)
  rnorm(
    n_draws,
    mean = mp,
    sd = sqrt(vp)
  )
}

build_study_panel <- function(i, label) {
  theta_o <- data$fiso[i]
  theta_r <- data$fisr[i]
  sigma_o <- data$se_fiso[i]
  sigma_r <- data$se_fisr[i]
  gamma_hat <- median(
    gamma_draws_mat[, paste0("gamma[", i, "]")]
  )
  grid <- seq(min(theta_o - 4 * sigma_o,theta_r - 4 * sigma_r), max(theta_o + 4 * sigma_o, theta_r + 4 * sigma_r), length.out = 500)
  curves <- bind_rows(
    tibble(theta = grid, density = dnorm(grid,mean = gamma_hat * theta_o, sd = sqrt(gamma_hat) * sigma_o), curva = "Prior ricalibrata"
    ),
    tibble(theta = grid, density = dnorm(grid, mean = theta_r,sd = sigma_r), curva = "Replica"
    )
  ) %>%
    mutate(study = label)
  posterior <- tibble(
    theta = compose_theta_posterior(i),
    study = label
  )
  list(
    curves = curves,
    posterior = posterior
  )
}

panels <- lapply(
  names(idx_selected),
  \(x) build_study_panel(idx_selected[x], x)
)

df_curves_all <- bind_rows(
  lapply(panels, `[[`, "curves")
)

df_posterior_all <- bind_rows(
  lapply(panels, `[[`, "posterior")
)
p22 <- ggplot() +
  geom_density(data = df_posterior_all, aes(x = theta),fill = "darkseagreen3", color = "darkgreen",
    alpha = 0.35, linewidth = 0.8,adjust = 1.2) +
  geom_line(data = df_curves_all, aes(x = theta,y = density, color = curva, linetype = curva), linewidth = 0.9) +
  facet_wrap(~ study, scales = "free",ncol = 3) +
  scale_color_manual(values = c("Prior ricalibrata" = "steelblue4", "Replica" = "firebrick3")) +
  scale_linetype_manual(values = c("Prior ricalibrata" = "dashed", "Replica" = "solid")) +
  labs(x = expression(theta), y = "Densità", color = NULL, linetype = NULL, title = "") +
  theme_gg_facet +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

study_idx <- 50
df50 <- build_study_panel(
  study_idx,
  data$study[study_idx]
)

p23 <- ggplot() +
  geom_density(data = df50$posterior, aes(x = theta), fill = "darkseagreen3", color = "darkgreen",
    alpha = 0.35, linewidth = 0.8, adjust = 1.2) +
  geom_line( data = df50$curves, aes(x = theta, y = density, color = curva, linetype = curva), linewidth = 1) +
  scale_color_manual(values = c("Prior ricalibrata" = "steelblue4","Replica" = "firebrick3")
  ) +
  scale_linetype_manual(values = c("Prior ricalibrata" = "dashed", "Replica" = "solid")) +
  labs(x = expression(theta),y = "Densità",color = NULL,linetype = NULL, title = ""
  ) +
  theme_gg +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

### Fig. 24 ###
gamma_FB_multi <- summary_long %>%
  filter(Prior %in% c("Beta(1,1)", "Beta(2,5)", "Beta(5,2)")) %>%
  select(Study, Prior, gamma_FB = Mean)

data_compare <- data %>%
  mutate(Study = row_number()) %>%
  select(Study, gamma_KL, gamma_EB) %>%
  left_join(gamma_FB_multi, by = "Study")  # join "many-to-many": una riga per Prior

df24 <- data_compare %>%
  pivot_longer(cols = c(gamma_KL, gamma_EB), names_to = "Metodo", values_to = "gamma_stima") %>%
  mutate(Metodo = recode(Metodo, gamma_KL = "hat(delta)[KL]", gamma_EB = "hat(delta)[EB]"))

p24<- ggplot(df24, aes(gamma_stima, gamma_FB)) +
  geom_point(alpha = 0.7, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "firebrick") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  facet_grid(Prior ~ Metodo, labeller = labeller(Metodo = label_parsed)) +
  labs(x = expression(hat(delta)), y = expression(E(delta ~ "|" ~ hat(theta)[r]))) +
  theme_gg +
  theme(strip.background = element_blank(), strip.text = element_text(size = 11),
        panel.spacing = unit(1.2, "lines"))

### Fig. 25 - 28 ###
plot_ppc <- function(df_lines, line_color){
  ggplot() +
    geom_line(data = df_lines, aes(theta, density, group = rep_id),
              color = line_color, alpha = 0.25, linewidth = 0.4) +
    geom_line(data = df_obs_curve, aes(theta, density), color = "firebrick", linewidth = 1.1) +
    facet_wrap(~study, scales = "free", ncol = 3) +
    labs(x = expression(theta), y = "Densita") +
    theme_gg_facet +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.ticks.y = element_blank(), axis.text.y = element_blank())
}

p25 <- plot_ppc(df_ppc_lines_prior_Beta52, "steelblue")

p26 <- plot_ppc(df_ppc_lines_prior_U01, "steelblue")

p27 <- plot_ppc(df_ppc_lines_post_Beta52, "darkgreen")

p28 <- plot_ppc(df_ppc_lines_post_U01, "darkgreen")
