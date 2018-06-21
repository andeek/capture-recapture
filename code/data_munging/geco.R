## command line args ---- 
## within db duplication levels 5; string distortion levels 1, 5, 10
# Rscript code/data_munging/jasa_pubs_10yr.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in the dup level and dist level.", call.=FALSE)
dup_level <- as.numeric(args[1])/100
dist_level <- as.numeric(args[2])/100

# 0. libraries ----
library(dplyr) # general manipulation
library(tidyr)
library(purrr)
set.seed(12345)

# 1. Get population
population <- read.csv("data/geco_sim/population.csv", stringsAsFactors = FALSE) 

## add birth dates
population[, "bdate"] <- sample(seq(as.Date('1980/01/01'), as.Date('2000/01/01'), by="day"), nrow(population))

# 2. Inclusion matrix
# sample D "databases" from "population" with specified levels of overlap/inclusion
D <- 4
strata_prop <- .75

# make the inclusion table
inclusion <- data.frame(strata = rep(c(1, 2), each = D),
                        prop = rep(c(strata_prop, 1 - strata_prop), each = D),
                        db = seq_len(D),
                        inclusion = c(rbeta(D, 1, 10), rbeta(D, 5, 1)))


# split population into strata
population$strata <- (seq_len(nrow(population)) > strata_prop*nrow(population)) + 1

# randomly select records from the population
record_db <- population
for(i in seq_len(D)) {
  split_pop <- split(record_db, record_db$strata)
  pop_list <- lapply(split_pop, function(pop) {
    stra <- unique(pop$strata)
    inc <- as.logical(rbinom(nrow(pop), 1, inclusion[inclusion$db == i & inclusion$strata == stra, "inclusion"]))
    record_db[record_db$strata == stra, paste0("db", i)] <- inc
    record_db[record_db$strata == stra,]
  })
  record_db <- do.call(rbind, pop_list)
}


# construct undublicated dbs
clean_db <- list(D)
truth_db <- list(D) # keep track of true inclusion by "id"
keep_cols <- c("fname", "lname", "bdate")
for(i in seq_len(D)){
  idx <- record_db[,paste0("db", i)]
  truth_db[[i]] <- data.frame(rec_id = population[idx, "rec_id"])
  clean_db[[i]] <- population[idx, c("rec_id", keep_cols)] # don't forget to remove id later
}

# 4. Add distortion and duplicates
# how many records will be duplicated 1x, 2x, 3x in each dataset?
# pull this from a multinomial distribution where p1 = .7, p2 = .2, p3 = .1
rep_prob <- c(.7, .2, .1)

dup_db <- lapply(clean_db, function(db) {
  num_dup <- ceiling(nrow(db)*dup_level)
  how_many <- rmultinom(1, size = num_dup, prob = rep_prob) 
  db$duplicate <- sample(c(rep(c(1, 2, 3), times = how_many), rep(0, times = nrow(db) - num_dup)))
  db[rep(row.names(db), db$duplicate + 1), -ncol(db)]
})

# create identity (population row_num) for each record
identity <- lapply(dup_db, function(db) {
  merge(data.frame(id = db$rec_id), # these are the duplicated ids
        data.frame(id = population$rec_id, rownum = seq_len(nrow(population))))$rownum # population row_num
})

# store truth of dups for later
truth_dup_db <- lapply(dup_db, function(db) {
  identity <- data.frame(record1 = numeric(0), record2 = numeric(0), rec_id = character(0))
  for(i in seq_len(nrow(db))) {
    for(j in seq_len(nrow(db))) {
      if(i < j) {
        if(db[i, "rec_id"] == db[j, "rec_id"]) {
          identity <- rbind(identity, cbind(i, j, rec_id = db[i, "rec_id"]))
        }
      }
    }
  }
  identity
})

# distort all of the duplicated records between and within dbs
# get idx of all duplicated records
idx_dup <- list()
for(i in seq_len(D)) {
  # which ids are in the duplicated db
  ids <- dup_db[[i]]$rec_id
  
  # flag those rows as duplicates if they appear in a "previous" db
  idx_dup[[i]] <- ids %in% as.character(do.call(rbind, truth_db[1:D < i])$rec_id)
  
  # now go find duplicates within each db
  idx_dup[[i]][unique(as.numeric(as.character(unlist(truth_dup_db[[i]][, 1:2]))))] <- TRUE
  
  idx_dup[[i]] <- which(idx_dup[[i]])
}

# remove rec_id from dup_db
dup_db <- lapply(dup_db, function(db) {
  db[, -which("rec_id" == names(db))]
})

# function to add distortion to a subset of the fields in each db
## TODO: this is where things can be changed for experimentation
add_noise <- function(db, idx, col_types) {
  p <- ncol(db)
  n <- nrow(db)
  
  for(i in seq_len(p)) {
    type <- col_types[i]
    if(type == "integer") {
      # if column in integer, add Normal noise then round to integer
      db[idx, i] <- round(db[idx, i] + rnorm(length(idx), mean = 0, sd = sqrt(var(db[, i], na.rm = TRUE))))
    } else if(type == "numeric") {
      # if column in numeric, add Normal noise
      db[idx, i] <- db[idx, i] + rnorm(length(idx), mean = 0, sd = sqrt(var(db[, i], na.rm = TRUE)))
    } else if(type == "string") {
      # if column is string, add string distortion
      strings <- db[idx, i]
      
      possible_replace <- c(letters, LETTERS, 
                            " ", "~", "`", "!", "@", "#", "$", "%", "^", "&", "*", "(", ")", ",", ".", "/", "?", ";", ":", "'", "\"")
      
      # replace randomly chosen characters in each chosen string at a given rate
      for(s in seq_len(length(strings))) {
        len <- nchar(strings[s])
        idx_rep <- sample(seq_len(len), ceiling(len*dist_level))
        for(id in idx_rep) {
          substr(strings[s], id, id) <- sample(possible_replace, 1)
        }
      }
      
      db[idx, i] <- strings
      
    } else if(type == "categorical") {
      # if column is categorial, use empirical distribution
      emp_dist <- prop.table(table(db[, i]))
      db[idx, i] <- sample(names(emp_dist), size = length(idx), prob = emp_dist)
    } else if(type == "date") {
      db[idx, i] <- as.Date(db[idx, i] + round(rnorm(length(idx), 0, 5)))
    } else {
      stop("Columns must be type string, integer, numeric, date, or categorical.")
    }
  }
  db
}

col_types <- c("string", "string", "date")
noisy_dup_db <- list()
for(i in seq_len(D)) {
  noisy_dup_db[[i]] <- add_noise(dup_db[[i]], idx = idx_dup[[i]], col_types = col_types) %>%
    separate(bdate, into = c("by", "bm", "bd"))
}

save(noisy_dup_db, identity, population, inclusion,
     file = paste0("data/geco_sim/geco_", dup_level*100, "dup_", dist_level*100, "dist.Rdata"))
