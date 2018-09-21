# libraries ----
library(ggplot2)
library(dplyr)
library(tidyr)

# load results ----
load("results/error_simulation/crc_zoom.Rdata")

# summaries ----
pop_N %>%
  group_by(type, r, i) %>%
  summarise(ll = quantile(pop_N, .025), ul = quantile(pop_N, .975)) %>%
  mutate(contains_truth = 1000 <= ul & 1000 >= ll) %>%
  group_by(type, r) %>%
  summarise(coverage = sum(contains_truth)/n()) %>%
  spread(type, coverage)

quantile(true_pop_N, c(.025, .975))


# plots ----
ggplot() +
  geom_density(aes(true_pop_N)) +
  geom_vline(aes(xintercept = 1000))

pop_N %>%
  group_by(type, r, i) %>%
  summarise(ll = quantile(pop_N, .025), ul = quantile(pop_N, .975)) %>%
  mutate(contains_truth = 1000 <= ul & 1000 >= ll) %>% 
  ggplot() +
  geom_segment(aes(x = i, xend = i, y = ll, yend = ul, colour = contains_truth)) +
  geom_hline(aes(yintercept = 1000)) +
  facet_grid(type~r, scales = "free_y") +
  scale_colour_manual(values = c("red", "black"))