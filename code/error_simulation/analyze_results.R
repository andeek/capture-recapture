# libraries ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# load results ----
load("results/error_simulation/crc_fixed.Rdata")
fixed_pop_N <- pop_N
load("results/error_simulation/crc.Rdata")
rl_pop_N <- pop_N
load("results/error_simulation/crc_zoom.Rdata")
rl_zoom_pop_N <- pop_N


# fixed plots ----
ggplot() +
  geom_density(aes(true_pop_N)) +
  geom_vline(aes(xintercept = 1000))

fixed_pop_N %>%
  filter(type == "add") %>%
  group_by(type, bucket_type, r, i) %>%
  summarise(ll = quantile(pop_N, .025), ul = quantile(pop_N, .975)) %>%
  mutate(contains_truth = 1000 <= ul & 1000 >= ll) %>% 
  ggplot() +
  geom_segment(aes(x = i, xend = i, y = ll, yend = ul, colour = contains_truth)) +
  geom_hline(aes(yintercept = 1000)) +
  facet_grid(type+bucket_type~r) +
  scale_colour_manual(values = c("red", "black")) -> p.add

fixed_pop_N %>%
  filter(type == "remove") %>%
  group_by(type, bucket_type, r, i) %>%
  summarise(ll = quantile(pop_N, .025), ul = quantile(pop_N, .975)) %>%
  mutate(contains_truth = 1000 <= ul & 1000 >= ll) %>% 
  ggplot() +
  geom_segment(aes(x = i, xend = i, y = ll, yend = ul, colour = contains_truth)) +
  geom_hline(aes(yintercept = 1000)) +
  facet_grid(type+bucket_type~r) +
  scale_colour_manual(values = c("red", "black")) -> p.remove

grid.arrange(p.add, p.remove)

# rl plots ----
rl_pop_N %>%
  group_by(type, r, i) %>%
  summarise(ll = quantile(pop_N, .025), ul = quantile(pop_N, .975)) %>%
  mutate(contains_truth = 1000 <= ul & 1000 >= ll) %>% 
  ggplot() +
  geom_segment(aes(x = i, xend = i, y = ll, yend = ul, colour = contains_truth)) +
  geom_hline(aes(yintercept = 1000)) +
  facet_grid(type~r) +
  scale_colour_manual(values = c("red", "black"))
