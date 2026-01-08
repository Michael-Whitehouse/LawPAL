library(ggplot2)
library(patchwork)

# Theme with reduced vertical spacing between facets
tight_strip_theme <- theme(
  panel.spacing.y = unit(0.00, "lines"),   # very small vertical spacing
  strip.text.y = element_text(
    size = 50,                               # smaller 'n = ...'
    face = "bold",
    margin = margin(t = -1, b = -10)           # move strip closer to panel
  ),
  strip.background = element_blank(),
)


p_beta <- ggplot(beta_long, aes(x = T, y = beta)) +
  geom_boxplot(width = 0.85, linewidth = 0.5) +
  geom_hline(yintercept = 0.15, colour = "red", linewidth = 0.5) +
  facet_wrap(
    ~ N,
    ncol = 1,
    labeller = labeller(N = function(x) paste("n =", x))
  ) +
  coord_cartesian(ylim = c(0.125, 0.2)) +
  scale_x_discrete(expand = expansion(mult = 0.2)) +
  labs(x = NULL, y = expression(beta)) +
  theme_bw() +
  tight_strip_theme 

p_gamma <- ggplot(gamma_long, aes(x = T, y = gamma)) +
  geom_boxplot(width = 0.85, linewidth = 0.5) +
  geom_hline(yintercept = 0.1, colour = "red", linewidth = 0.5) +
  facet_wrap(
    ~ N,
    ncol = 1,
    labeller = labeller(N = function(x) paste("n =", x))
  ) +
  coord_cartesian(ylim = c(0.075, 0.175)) +
  scale_x_discrete(expand = expansion(mult = 0.2)) +
  labs(x = NULL, y = expression(gamma)) +
  theme_bw() +
  tight_strip_theme 

p_mu <- ggplot(mu_long, aes(x = T, y = mu)) +
  geom_boxplot(width = 0.85, linewidth = 0.5) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  facet_wrap(
    ~ N,
    ncol = 1,
    labeller = labeller(N = function(x) paste("n =", x))
  ) +
  coord_cartesian(ylim = c(0.4, 0.6)) +
  scale_x_discrete(expand = expansion(mult = 0.2)) +
  labs(x = paste("Time horizon ", expression(T)), y = expression(mu[q])) +
  theme_bw() +
  tight_strip_theme

p_sigma <- ggplot(sigma_long, aes(x = T, y = sigma)) +
  geom_boxplot(width = 0.85, linewidth = 0.5) +
  geom_hline(yintercept = 0.1, colour = "red", linewidth = 0.5) +
  facet_wrap(
    ~ N,
    ncol = 1,
    labeller = labeller(N = function(x) paste("n =", x))
  ) +
  coord_cartesian(ylim = c(0.05, 0.15)) +
  scale_x_discrete(expand = expansion(mult = 0.2)) +
  labs(x = "Time horizon T", y = expression(sigma[q]^2)) +
  theme_bw() +
  tight_strip_theme +
  theme(
    axis.title.y = element_text(
      margin = margin(r = -7)   # ← CONTROLS DISTANCE TO AXIS
    )
  )


(p_beta | p_gamma) / (p_mu | p_sigma)
