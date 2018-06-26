## libraries ----
library(ggplot2)
library(dplyr) 
library(tidyr)

## load results ----
load("results/geco_sim/crc_5dup_5dist.Rdata")
pop_N_bayes_5 <- pop_N_bayes
pop_N <- data.frame(dist = 5, truth = pop_N_truth, mpmms = pop_N_mpmms)
mpmms_lambda_5 <- mpmms_lambda
min_5 <- min(c(pop_N_truth, pop_N_mpmms, pop_N_bayes))
max_5 <- max(c(pop_N_truth, pop_N_mpmms, pop_N_bayes))

load("results/geco_sim/crc_5dup_10dist.Rdata")
pop_N_bayes_10 <- pop_N_bayes
pop_N <- rbind(pop_N, data.frame(dist = 10, truth = pop_N_truth, mpmms = pop_N_mpmms))
mpmms_lambda_10 <- mpmms_lambda
min_10 <- min(c(pop_N_truth, pop_N_mpmms, pop_N_bayes))
max_10 <- max(c(pop_N_truth, pop_N_mpmms, pop_N_bayes))

load("results/geco_sim/crc_5dup_15dist.Rdata")
pop_N_bayes_15 <- pop_N_bayes
pop_N <- rbind(pop_N, data.frame(dist = 15, truth = pop_N_truth, mpmms = pop_N_mpmms))
mpmms_lambda_15 <- mpmms_lambda
min_15 <- min(c(pop_N_truth, pop_N_mpmms, pop_N_bayes))
max_15 <- max(c(pop_N_truth, pop_N_mpmms, pop_N_bayes))


## plots ----
pop_N %>%
  gather(method, N, -dist) %>%
  ggplot() +
  geom_vline(aes(xintercept = 1000), lty = 2) +
  geom_density(aes(N, fill = method, colour = method), alpha = .5) +
  facet_grid(~dist) +
  xlim(c(0, 2500))

## summaries ---
pop_N %>%
  gather(method, N, -dist) %>%
  group_by(dist, method) %>%
  summarise(quantile_truth = sum(N < 1000)/n()) %>%
  spread(method, quantile_truth)


## posterior expected distribution ----
## bin, get empirical density
## this will be U()

post_exp_dsn <- function(pop_N_bayes, min, max, n) {
  p <- ncol(pop_N_bayes)
  from <- min(pop_N_bayes)
  to <- max(pop_N_bayes)
  emp_dens <- density(pop_N_bayes[, 1], from = min, to = max, n = n)
  U <- data.frame(x = emp_dens$x, y = emp_dens$y)
  
  for(i in 2:p) {
    temp_dens <- density(pop_N_bayes[, i], from = min, to = max, n = n)  
    U$x <- U$x + temp_dens$x
  }
  U$x <- U$x/p
  U
}

n_bin <- 5000
min_tot <- min(min_5, min_10)
max_tot <- max(max_5, max_10)

## truth
pop_N_density <- rbind(data.frame(dist = 5, data.frame(density(pop_N[pop_N$dist == 5,]$truth, from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])),
                       data.frame(dist = 10, data.frame(density(pop_N[pop_N$dist == 10,]$truth, from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])),
                       data.frame(dist = 15, data.frame(density(pop_N[pop_N$dist == 15,]$truth, from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])))
names(pop_N_density) <- c("dist", "x", "truth")
pop_N_density$x <- round(pop_N_density$x, 5)

## add mmpms
rbind(data.frame(dist = 5, data.frame(density(pop_N[pop_N$dist == 5,]$mpmms, from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])),
      data.frame(dist = 10, data.frame(density(pop_N[pop_N$dist == 10,]$mpmms, from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])),
      data.frame(dist = 15, data.frame(density(pop_N[pop_N$dist == 15,]$mpmms, from = min_tot, to = max_tot, n = n_bin)[c("x", "y")]))) %>%
  rename(mpmms = y) %>% 
  mutate(x = round(x, 5)) %>%
  left_join(pop_N_density, by = c("dist", "x")) -> pop_N_density

## add marginal posterior
# rbind(data.frame(dist = 5, post_exp_dsn(pop_N_bayes_5, min = min_tot, max = max_tot, n = n_bin)),
#       data.frame(dist = 10, post_exp_dsn(pop_N_bayes_10, min = min_tot, max = max_tot, n = n_bin))) %>%
#   rename(marginal_post = y) %>%
#   mutate(x = round(x, 5)) %>%
#   left_join(pop_N_density, by = c("dist", "x")) -> pop_N_density

## add joint posterior
rbind(data.frame(dist = 5, data.frame(density(as.numeric(pop_N_bayes_5), from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])),
      data.frame(dist = 10, data.frame(density(as.numeric(pop_N_bayes_10), from = min_tot, to = max_tot, n = n_bin)[c("x", "y")])),
      data.frame(dist = 15, data.frame(density(as.numeric(pop_N_bayes_15), from = min_tot, to = max_tot, n = n_bin)[c("x", "y")]))) %>%
  rename(joint_post = y) %>% 
  mutate(x = round(x, 5)) %>%
  left_join(pop_N_density, by = c("dist", "x")) -> pop_N_density

pop_N_density %>% 
  gather(method, y, -dist, -x) %>%
  mutate(xp = x*y) %>%
  group_by(method, dist) %>%
  summarise(mean_N = sum(xp)/sum(y)) -> pop_N_mean


pop_N_density %>%
  gather(method, y, -dist, -x) %>% 
  ggplot() +
  geom_line(aes(x, y, colour = method)) +
  geom_polygon(aes(x, y, fill = method), alpha = .2) +
  geom_vline(aes(xintercept = mean_N, colour = method), lty = 2, data = pop_N_mean) +
  geom_vline(aes(xintercept = 1000), lty = 1) +
  facet_wrap(~dist) +
  xlim(c(750, 1500))
              

## everything below is crackpot ----
## look at performance ----
load("results/geco_sim/eber_5dup_5dist.Rdata")
load("results/geco_sim/eber_diag_5dup_5dist.Rdata")
load("data/geco_sim/geco_5dup_5dist.Rdata")

identity_vec <- do.call(c, identity)
identity_pair <- eber::clust2pairs(eber::memb2clust(identity_vec))

conf <- eber::confusion.matrix(eber::clust2pairs(eber::memb2clust(mpmms_lambda_5)), identity_pair)
mpmms_recall = conf$TP/(conf$TP + conf$FN)
mpmms_precision = conf$TP/(conf$TP + conf$FP)

m <- ncol(lambda)
burnin <- 10000
thinning <- 40

diag_res$eval[seq_len(m) > burnin & seq_len(m) %% thinning == 0,] %>% 
  mutate(iter = 1:n()) %>%
  ggplot() +
  geom_line(aes(iter, precision)) +
  geom_hline(aes(yintercept = mpmms_precision), colour = "red")

diag_res$eval[seq_len(m) > burnin & seq_len(m) %% thinning == 0,] %>% 
  mutate(iter = 1:n()) %>%
  ggplot() +
  geom_line(aes(iter, recall)) +
  geom_hline(aes(yintercept = mpmms_recall), colour = "red")


## relationship with precision and recall ----
load("results/geco_sim/eber_5dup_5dist.Rdata")
load("results/geco_sim/eber_diag_5dup_5dist.Rdata")
load("data/geco_sim/geco_5dup_5dist.Rdata")

identity_vec <- do.call(c, identity)
identity_pair <- eber::clust2pairs(eber::memb2clust(identity_vec))

conf <- eber::confusion.matrix(eber::clust2pairs(eber::memb2clust(mpmms_lambda_5)), identity_pair)
mpmms_recall = conf$TP/(conf$TP + conf$FN)
mpmms_precision = conf$TP/(conf$TP + conf$FP)

m <- ncol(lambda)
burnin <- 10000
thinning <- 40


diag_res$eval[seq_len(m) > burnin & seq_len(m) %% thinning == 0, c("precision", "recall")] %>%
  mutate(MSE = colSums((pop_N_bayes - 1000)^2)/nrow(pop_N_bayes)) %>%
  mutate(method = "full_bayes") %>%
  bind_rows(data.frame(precision = mpmms_precision, recall = mpmms_recall, MSE = mean((pop_N$mpmms[pop_N$dist == 5] - 1000)^2), method = "mpmms")) %>%
  ggplot() +
  geom_point(aes(precision, recall, colour = MSE, shape = method)) +
  scale_color_continuous(low = "blue", high = "white")





