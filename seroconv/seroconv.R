cat("sample size")

library(tidyverse)
library(furrr)

plan(multiprocess)

seroconv_dir <- "seroconv"

# Functions ===================================================================

expit <- function(x) 1 - 1 / (1 + exp(x))

gen_one_pop <- function(n_per_group = 50,
                        b0 = log(1 / 3), # Logodds of seroconv for no adjuvant
                        b_adjuvant = log(1.3)) { # log-OR for seroconv
  n <- n_per_group * 2
  pop <- tibble(
    id = 1:n,
    adjuvant = rep(c(0, 1), n_per_group),
    adjuvant_lbl = recode(
      adjuvant,
      "0" = "Vaccine", "1" = "Vaccine + adjuvant"
    ),
    seroconv = extraDistr::rbern(n, expit(b0 + b_adjuvant * adjuvant))
  )
  attr(pop, "true_vals") <- list(
    n_per_group = n_per_group,
    b0 = b0,
    b_adjuvant = b_adjuvant
  )
  pop
}

fit_one <- function(...) {
  pop <- gen_one_pop(...)
  res <- glm(seroconv ~ adjuvant, binomial, pop) %>%
    broom::tidy() %>%
    bind_cols(as_tibble(attr(pop, "true_vals")))
  res
}

fit_many <- function(nsim = 10, ...) {
  map_dfr(1:nsim, function(i) fit_one(...) %>% mutate(i = i))
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(seroconv_dir, paste0(name, ".pdf")), plot,
    units = "cm",
    ...
  )
}

save_data <- function(data, name) {
  write_csv(data, file.path(seroconv_dir, paste0(name, ".csv")))
}

# Script ======================================================================

# Plots of an example of a simulated sample at each of the
# relevant true differences
s <- function(b_adjuvant, ...) {
  gen_one_pop(b_adjuvant = b_adjuvant, ...) %>%
    mutate(
      OR = exp(b_adjuvant),
      seroconv_lbl = recode(
        seroconv,
        "0" = "Not seroconverted", "1" = "Seroconverted"
      ),
    )
}

samples <- map_dfr(log(seq(1.7, 2.2, 0.1)), s, n_per_group = 1e6)

one_example_table <- samples %>%
  count(OR, adjuvant, seroconv_lbl) %>%
  mutate(n = n / 1e6) %>%
  pivot_wider(names_from = "seroconv_lbl", values_from = n)

save_data(one_example_table, "example-table")

# Simulate studies

pars <- tribble(
  ~b_adjuvant, ~n_per_group,

  log(2.2), 125,
  log(2.1), 125,
  log(2), 125,

  log(1.9), 150,
  log(2), 150,
  log(2.1), 150,

  log(1.7), 250,
  log(1.8), 250,
  log(1.9), 250,
)

sim_res <- future_pmap_dfr(
  pars,
  function(b_adjuvant, n_per_group) {
    fit_many(1000, b_adjuvant = b_adjuvant, n_per_group = n_per_group)
  }
)

save_data(sim_res, "results")

summary <- sim_res %>%
  filter(term == "adjuvant") %>%
  group_by(b_adjuvant, n_per_group) %>%
  summarise(
    power = sum(estimate > 0 & p.value < 0.05) / n(),
    .groups = "drop"
  ) %>%
  mutate(OR = exp(b_adjuvant)) %>%
  select(-b_adjuvant) %>%
  arrange(n_per_group)

save_data(summary, "summary-seroconv")
