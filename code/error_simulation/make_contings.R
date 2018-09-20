# libraries + seed ----
set.seed(1234)

# global params ----
M <- 1000 # true pop size
D <- 5 # number of databases
strata_prop <- .75 # proportion of records in each strata
error_levels <- seq(.05, .5, by = .05)
rep <- 100

# inclusion table ----
# sample D "databases" from "population" with specified levels of overlap/inclusion
# 2 strata - 1: large and hard to find, 2: small and easy to find
inclusion <- data.frame(strata = rep(c(1, 2), each = D), prop = rep(c(strata_prop, 1 - strata_prop), each = D), db = seq_len(D), inclusion = c(rbeta(D, 3, 15), rbeta(D, 15, 3)))

# population ----
# ids + split pop into strata
pop <- data.frame(id = seq_len(M), strata = rep(c(1, 2), c(floor(strata_prop*M), M - floor(strata_prop*M))))

# make true contingency table ----
true_assignment <- pop # data frame showing how each record is assigned
for(i in seq_len(D)) {
  split_pop <- split(true_assignment, true_assignment$strata)
  pop_list <- lapply(split_pop, function(p) {
    stra <- unique(p$strata)
    inc <- as.logical(rbinom(nrow(p), 1, inclusion[inclusion$db == i & inclusion$strata == stra, "inclusion"]))
    true_assignment[true_assignment$strata == stra, paste0("db", i)] <- inc
    true_assignment[true_assignment$strata == stra,]
  })
  true_assignment <- do.call(rbind, pop_list)
  rownames(true_assignment) <- NULL
}
true_conting <- as.data.frame(table(true_assignment[, paste0("db", seq_len(D))])) # true contingengy table

# distort ----
# function for distorting a contingency table 
# takes a given percent of singletons and re assigns them to a given bucket
# or takes a given percent of non singletones and makes them singletons
# input: conting - a contingency table (df), r - a proportion to be distorted, error_type - "remove" or "add"
# output: distorted contingency table
distort_conting <- function(conting, r, error_type = "remove") {
  stopifnot(error_type %in% c("remove", "add"))
  stopifnot(is.numeric(r) & r <= 1 & r >= 0)
  
  tab <- conting # make copy
  tab[-ncol(tab)] <- apply(tab[-ncol(tab)], 2, function(x) as.logical(as.character(x))) # change factors to logical
  tab$type <- ifelse(rowSums(tab[, -ncol(tab)]) == 0, "zero", ifelse(rowSums(tab[, -ncol(tab)]) == 1, "single", "multi")) # classify buckets
  
  num_error <- ceiling(sum(tab[tab$type == "single", "Freq"])*r) # number of singletons to remove or add
  
  # choose where to take or add the singletons from
  where_single <- table(sample(which(tab$type == "single"), num_error, replace = TRUE))
  
  # choose where to move or take the error singletons from
  where_multi <- table(sample(which(tab$type == "multi"), num_error, replace = TRUE))
  
  if(error_type == "remove") {
    tab[names(where_single), "Freq"] <- tab[names(where_single), "Freq"] - where_single
    tab[names(where_multi), "Freq"] <- tab[names(where_multi), "Freq"] + where_multi
  } else if(error_type == "add") {
    tab[names(where_single), "Freq"] <- tab[names(where_single), "Freq"] + where_single
    tab[names(where_multi), "Freq"] <- tab[names(where_multi), "Freq"] - where_multi
  }
  return(tab[, -ncol(tab)])  
}

error_remove_conting <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_conting(true_conting, r, "remove")))
error_add_conting <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_conting(true_conting, r, "add")))
names(error_remove_conting) <- names(error_add_conting) <- error_levels

# save ----
save(true_assignment, true_conting, error_add_conting, error_remove_conting, M, D, strata_prop, error_levels, rep, inclusion,
     file = "results/error_simulation/contings.Rdata")



