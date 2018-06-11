## command line args ---- 
## duplication levels 5, 10, 30, 50; string distortion levels 5, 10, 15
# Rscript code/data_munging/jasa_pubs_10yr.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in the dup level (5, 10, 30, 50) and dist level (5, 10, 15).", call.=FALSE)
if (!(args[1] %in% c(5, 10, 30, 50))) stop("Pass in the duplication level (5, 10, 30, 50)", call.=FALSE)
if (!(args[2] %in% c(5, 10, 15))) stop("Pass in the distortion level (5, 10, 15)", call.=FALSE)
dup_level <- as.numeric(args[1])/100
dist_level <- as.numeric(args[2])/100

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
  separate(issued, into = c("year", "month", "day"), fill = "right") %>%
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
  left_join(inclusion, by = "partition") %>%
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

# create identity (population row_num) for each record
identity <- lapply(dup_db, function(db) {
  merge(data.frame(id = db$DOI), # these are the duplicated DOIs
        data.frame(id = population$DOI, rownum = seq_len(nrow(population))))$rownum # population row_num
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

save(noisy_dup_db, identity, population, inclusion,
     file = paste0("data/jasa_sim/jasa_", dup_level*100, "dup_", dist_level*100, "dist.Rdata"))
