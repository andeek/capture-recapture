#TODO: this is where things can be changed for experimentation
# duplication level: start with 50, 30, 10 (can vary between dbs later)
dup_level <- .05

## TODO: if a string is chosen to be "distorted", how distorted do we want it?
# start lower and ramp up the string distortion .05, .1, .15
dist_level <- .15


# 0. libraries ----
library(tidyverse) # general manipulation
library(rcrossref) # get metadata from journals crossref
set.seed(12345)

# 1. Get population - JASA Articles from the past 10 years
# only care about some columns
cols <- c('ISSN', 'DOI', 'issue', 'issued', 'page', 'publisher', 'title', 'volume', 'type', 'subject', 'author')

# get JASA articles
jasa <- cr_journals(query = "Journal of the American Statistical Association")
jasa_10yr <- cr_journals(issn = jasa$data$issn, works = TRUE, filter = c(from_pub_date='2008-02-01'), select = cols, cursor = "*")

# authors are their own tibble, change to character column
# separate issued date into year, month, day
jasa_10yr$data %>%
  mutate(author = map_chr(author, ~ paste(.$given, .$family, collapse = " "))) %>%
  separate(issued, into = c("year", "month", "day")) %>%
  separate(page, into = c("page_begin", "page_end")) %>%
  mutate(year = as.integer(year), month = as.integer(month), issue = as.integer(issue), 
         page_begin = as.integer(page_begin), page_end = as.integer(page_end), volume = as.integer(volume)) %>%
  filter(author != "") -> population

# 2. Inclusion matrix
# sample D "databases" from "population" with specified levels of overlap/inclusion
D <- 4

# make the inclusion table
inclusion <- expand.grid(rep(list(c(FALSE, TRUE)), times = D)) %>% mutate(partition = seq_len(2^D))
names(inclusion)[seq_len(D)] <- paste0("db", seq_len(D))

# get incusion probabilities from generative model
## TODO: this is where things can be changed for experimentation
inclusion$p <- prop.table(runif(2^D))

# partition the unique records
splits <- ceiling(inclusion$p*nrow(population))
population$partition <- sample(rep(seq_len(2^D), times = splits))[seq_len(nrow(population))] 

# 3. Get clean databases
# join to the inclusion table to see which records are in which db
population %>% 
  left_join(inclusion) %>%
  select(-p, -partition) -> record_db

# construct undublicated dbs
clean_db <- list(D)
truth_db <- list(D) # keep track of true inclusion by "id" - doi
keep_cols <- c("issue", "year", "month", "page_begin", "page_end", "title", "volume", "author")
for(i in seq_len(D)){
  idx <- record_db[,paste0("db", i)]
  truth_db[[i]] <- data.frame(DOI = population[idx, "DOI"])
  clean_db[[i]] <- population[idx, c("DOI", keep_cols)] # don't forget to remove DOI later
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

# store truth of dups for later
truth_dup_db <- lapply(dup_db, function(db) {
  identity <- data.frame(record1 = numeric(0), record2 = numeric(0), doi = character(0))
  for(i in seq_len(nrow(db))) {
    for(j in seq_len(nrow(db))) {
      if(i < j) {
        if(db[i, "DOI"] == db[j, "DOI"]) {
          identity <- rbind(identity, cbind(i, j, DOI = db[i, "DOI"]))
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
  # which DOIs are in the duplicated db
  dois <- dup_db[[i]]$DOI
  
  # flag those rows as duplicates if they appear in a "previous" db
  idx_dup[[i]] <- dois %in% as.character(do.call(rbind, truth_db[1:D < i])$DOI)
  
  # now go find duplicates within each db
  idx_dup[[i]][unique(as.numeric(as.character(unlist(truth_dup_db[[i]][, 1:2]))))] <- TRUE
  
  idx_dup[[i]] <- which(idx_dup[[i]])
}

# remove DOI from dup_db
dup_db <- lapply(dup_db, function(db) {
  db[, -which("DOI" == names(db))]
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
    } else {
      stop("Columns must be type string, integer, numeric, or categorical.")
    }
  }
  db
}

col_types <- c("integer", "integer", "integer", "integer", "integer", "string", "integer", "string")
noisy_dup_db <- list()
for(i in seq_len(D)) {
  noisy_dup_db[[i]] <- add_noise(dup_db[[i]], idx = idx_dup[[i]], col_types = col_types)
}

save(noisy_dup_db, truth_db, truth_dup_db, population,
     file = paste0("../../data/jasa_sim_datasets/jasa_pubs_10yr_", dup_level*100, "dup_", dist_level*100, "dist.Rdata"))
