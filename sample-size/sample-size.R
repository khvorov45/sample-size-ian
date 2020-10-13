cat("sample size")

library(tidyverse)
library(furrr)

plan(multiprocess)

sample_size_dir <- "sample-size"

# Functions ===================================================================

# Returns OG scale titres
censor_logtitres <- function(logtitres) {
  cuts <- cut(exp(logtitres), c(-Inf, 5 * 2^(1:10), Inf)) %>% as.integer()
  5 * 2^(cuts - 1)
}

# Returns log titres
midpoint_titres <- function(titres) {
  if_else(titres == 5L, log(titres), log(titres) + log(2) / 2)
}

percent_diff_to_b_adjuvant <- function(percent_diff) {
  log((percent_diff / 100) + 1)
}

gen_one_pop <- function(n_per_group = 50,
                        b0 = 1.1, # Postvax for no adjuvant for base 1 (3-fold)
                        b_logbaseline = 0.8, # Add to logpostvax for 1 logbase
                        percent_diff = 10, # % diff b/w adjuvant and no adjuvant
                        b_adjuvant = percent_diff_to_b_adjuvant(percent_diff),
                        mean_logbaseline = 3.7, sd_logbaseline = 1,
                        sd_logpostvax = sd_logbaseline) {
  titre_vals <- 10 * 2^(0:10)
  n <- n_per_group * 2
  tibble(
    id = 1:n,
    adjuvant = extraDistr::rbern(n, 0.5),
    adjuvant_lbl = recode(
      adjuvant,
      "0" = "Vaccine", "1" = "Vaccine + adjuvant"
    ),
    logbaseline = rnorm(n, mean_logbaseline, sd_logbaseline),
    logpostvax = rnorm(
      n, b0 + b_logbaseline * logbaseline + b_adjuvant * adjuvant, sd_logpostvax
    ),
    baseline = censor_logtitres(logbaseline),
    postvax = censor_logtitres(logpostvax),
    logbaseline_mid = midpoint_titres(baseline),
    logpostvax_mid = midpoint_titres(postvax),
  )
}

fit_one <- function(...) {
  lm(logpostvax_mid ~ logbaseline_mid + adjuvant, gen_one_pop(...)) %>%
    broom::tidy()
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(sample_size_dir, paste0(name, ".pdf")), plot,
    units = "cm",
    ...
  )
}

# Script ======================================================================

relevant_diffs_percent <- seq(5, 30, 5)

# Plots of an example of a simulated sample at each of the
# relevant true differences
s <- function(percent_diff, ...) {
  gen_one_pop(percent_diff = percent_diff, ...) %>%
    mutate(
      percent_diff = percent_diff,
      percent_diff_lbl = glue::glue(
        "Diff {percent_diff}%"
      )
    )
}
samples <- map_dfr(relevant_diffs_percent, s) %>%
  mutate(percent_diff_lbl = reorder(percent_diff_lbl, percent_diff))

one_example_plot <- samples %>%
  select(id, percent_diff_lbl, adjuvant_lbl, baseline, postvax) %>%
  pivot_longer(
    c(baseline, postvax),
    names_to = "timepoint", values_to = "titre"
  ) %>%
  ggplot(aes(timepoint, titre)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.spacing = unit(0, "null")
  ) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:10)) +
  scale_x_discrete(expand = c(0.1, 0)) +
  facet_grid(percent_diff_lbl ~ adjuvant_lbl) +
  geom_line(aes(group = id), alpha = 0.3) +
  geom_point(alpha = 0.5, shape = 16)

save_plot(one_example_plot, "example-sample-titres", width = 10, height = 15)

one_example_diff_plot <- samples %>%
  mutate(diff = exp(logpostvax_mid - logbaseline_mid)) %>%
  ggplot(aes(adjuvant_lbl, diff)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.spacing = unit(0, "null")
  ) +
  scale_y_log10(
    "Fold-change postvax vs baseline",
    breaks = c(1, 2, 3, 5, 10, 20, 40)
  ) +
  xlab("Group") +
  facet_wrap(~percent_diff_lbl) +
  geom_boxplot()

save_plot(
  one_example_diff_plot, "example-sample-diffs",
  width = 10, height = 12
)
