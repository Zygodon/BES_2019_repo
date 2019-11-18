---
title: "Notebook"
author: "John"
date: "09/11/2019"
output:
  # word_document: default
  # pdf_document: default
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---
Tidying up stuff for BES 2019 poster (Note: could poster have been made as a
presentation in a new .rmd document?). 

First candidate for tidying: nvc baplot.R

nvc barplot


```r
# Libraries
library("RMySQL")
```

```
## Loading required package: DBI
```

```r
library("RColorBrewer")
library(ggplot2)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

########################## MAIN ##############################

# GET DATA FROM DB
# Open connection to meadows data base to get the data
# My local DB
mydb = dbConnect(MySQL(), user='root', password='Mysql130641', dbname='meadows', host='127.0.0.1')
# Remote DB
#mydb = dbConnect(MySQL(), user='guest', dbname='meadows', port = #3306, host='sxouse.ddns.net')
rs1 = dbSendQuery(mydb, "select nvc from assemblies where nvc is not null and nvc != major_nvc_community;")
data <- fetch(rs1, n=-1)
dbDisconnectAll()
```

```
## Warning: Closing open result sets
```

```
## 1 connection(s) closed.
```

```r
# Tibble g: hand coded analysis groups - nvc with few samples grouped e.g. MG10.
# Also mires grouped - will be excluded from this analysis.
g <- tibble(grp= c(rep("M", 7), rep("MG10", 2), rep("MG1", 3), "MG5a", "MG5c", "MG6a", "MG6b", rep("MG7", 3), "MG9a"))
# grps: Distinct nvcs matched with their analysis groups.
grps <- bind_cols(data %>% distinct(nvc) %>% arrange(nvc), g)
# >left join attaches an analysis group to record for each sample.
d <- left_join (data, grps, by="nvc") %>% mutate_all(tolower)

# Bar plot. Note hand-coded palette corresponding to analysis groups 
# e.g. "M" and "MG9", white - excluded from analysis.
bp1 <- ggplot(d, aes(x=nvc, fill=grp)) + 
  geom_bar(colour="black", stat="count") +
  scale_fill_manual(values = 
  c("white", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#ffff33", "#a65628", "white")) +
  coord_flip() +
  labs(x = "assessed NVC", fill = "analysis group")
print(bp1)
```

![](Notebook_files/figure-html/unnamed-chunk-1-1.png)<!-- -->
Moving on to prox. This code gets sample and standard vectors in 22 dimensional species frequency space. The sample frequencies calculated from the raw DB data.
Results visualised as bubble plot, nearest standard frequency vs distance to MG5 for each sample; samples coloured according to assessed nvc (lover case, e.g. "mg5a").


```r
# Libraries
library("RMySQL")
library(tidyverse)
```

```
## -- Attaching packages ------------------------------------------ tidyverse 1.2.1 --
```

```
## v tibble  2.1.3     v purrr   0.3.2
## v tidyr   1.0.0     v stringr 1.4.0
## v readr   1.3.1     v forcats 0.4.0
```

```
## -- Conflicts --------------------------------------------- tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

edist <- function(s1,a1){
  x <- s1-a1
  sqrt(sum(x^2))
}

ppx <-  function(r)
{
  names(eds)[which(r[1:7] == min(r))]
}


########################## MAIN ##############################
# GET DATA FROM DB
# Open connection to meadows data base
con = dbConnect(MySQL(), user='root', password='Mysql130641', dbname='meadows', host='127.0.0.1')
# Remote DB with password - works Ok but table mg_standards6 is not available on PI. Should update.
# con <- dbConnect(MySQL(), 
#                  user  = "team",
#                  password    = rstudioapi::askForPassword("Database password"),
#                  port = 3306,
#                  host   = "sxouse.ddns.net")
# Get nvc standard frequencies
q <- sprintf('SELECT * FROM meadows.mg_standards6;')
rs1 = dbSendQuery(con, q)
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 2 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 3 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 4 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 5 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 6 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 7 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 8 imported as
## numeric
```

```r
stdFreqs <- fetch(rs1, n=-1) 
# Make string of species names to use in selecting species hits
s <- toString(stdFreqs %>% select(species))
s <- substring(s, 4, nchar(s)-2)

# Get trials for assemblies
# Start with MG5a,c, and MG6a,b
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG5a", "MG5c", "MG6a", "MG6b")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(con, q)
trials_mg5_6 <- as_tibble(fetch(rs1, n=-1))

# Separate query for MG7b, c and d (DON'T exclude the few MG7c records)
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG7b", "MG7c", "MG7d")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(con, q)
trials_mg7 <- as_tibble(fetch(rs1, n=-1))
trials_mg7 <- trials_mg7 %>% select(assemblies_id, quadrat_count) %>% mutate(nvc = "MG7")

# Separate query for amalgamated MG1
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG1b", "MG1c", "MG1e")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(con, q)
trials_mg1 <- as_tibble(fetch(rs1, n=-1))
trials_mg1 <- trials_mg1 %>% select(assemblies_id, quadrat_count) %>% mutate(nvc = "MG1")

# Separate query for amalgamated MG10
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG10", "MG10a", "MG10b")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(con, q)
trials_mg10 <- as_tibble(fetch(rs1, n=-1))
trials_mg10 <- trials_mg10 %>% select(assemblies_id, quadrat_count) %>% mutate(nvc = "MG10")

# Join them together with bind_rows
assembly_trials <- bind_rows(trials_mg5_6, trials_mg7, trials_mg1, trials_mg10)
colnames(assembly_trials) <- c("assemblies_id", "community", "trials")

# Get hits for each species, separate query required for each subgroup.
# Amalgamated by assembly
# Start with MG5a,c and MG6a,b
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG5a", "MG5c", "MG6a", "MG6b")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)
rs1 = dbSendQuery(con, q)
hits_mg5_6 <- as_tibble(fetch(rs1, n=-1))

# Then MG7
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG7b", "MG7c", "MG7d")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)
rs1 = dbSendQuery(con, q)
hits_mg7 <- as_tibble(fetch(rs1, n=-1))
hits_mg7 <- hits_mg7 %>% select(assemblies_id, nvc, species_name, 'count(species.species_id)') %>% mutate(nvc = "MG7")

# Then MG1
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG1b", "MG1c", "MG1e")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)

rs1 = dbSendQuery(con, q)
hits_mg1 <- as_tibble(fetch(rs1, n=-1))
hits_mg1 <- hits_mg1 %>% select(assemblies_id, nvc, species_name, 'count(species.species_id)') %>% mutate(nvc = "MG1")

# Then MG10
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG10", "MG10a", "MG10b")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)
rs1 = dbSendQuery(con, q)
hits_mg10 <- as_tibble(fetch(rs1, n=-1))
hits_mg10 <- hits_mg10 %>% select(assemblies_id, nvc, species_name, 'count(species.species_id)') %>% mutate(nvc = "MG10")

# Join them together with bind_rows
species_hits <- bind_rows(hits_mg5_6, hits_mg7, hits_mg1, hits_mg10)
colnames(species_hits) <- c("assembly", "assemblies_id", "community", "species", "hits")

assembly_hits <- species_hits %>% group_by(assemblies_id, community, species) %>% summarise(hits = sum(hits)) 

# Close the DB connection
dbDisconnectAll()
```

```
## Warning: Closing open result sets
```

```
## 1 connection(s) closed.
```

```r
# Ensure species not represented in a community get 0 hits
wide <- spread(assembly_hits, key=species, value=hits) # Coerces NAs
wide[is.na(wide)] <- 0 # Replaces NAs with 0
# And reconstitute assembly hits. Good general solution to the problem
# assembly_hits <- wide %>% gather(Anthoxanthum_odoratum:Trifolium_repens, key= species, value=hits)
assembly_hits <- wide %>% gather(colnames(wide[3:length(colnames(wide))])
, key= species, value=hits)

# Survey frequencies: join hits and trials and mutate to calculate frequencis
srvf <- (full_join(assembly_trials, assembly_hits, by = c("assemblies_id", "community")) 
         %>% mutate(freq = hits/trials)
         # And mutate to get 0.05, 0.5 and 0.95 quantiles (use 0.05 and 0.95 for CrI)
         %>%  mutate(CrI5 = qbeta(0.05,hits+1, 1+trials-hits))
         %>%  mutate(q50 = qbeta(0.5,hits+1, 1+trials-hits)) # For comparison with frequency as hits/trials
         %>%  mutate(CrI95 = qbeta(0.95,hits+1, 1+trials-hits))
         %>%  mutate(source="survey")) # Note pipe enclosed in () to ensure run onto new line

# Distance measures
output <- matrix(nrow=204, ncol=7) # The numer of analysis groups - we'll have 7 columns
t <-  rep(0,7)

for (i in 1:length(assembly_trials$assemblies_id)) {
  # a1: 22 values of q50 (mean)
  a1 <- srvf %>% filter(assemblies_id == assembly_trials$assemblies_id[i]) %>% select(q50)
  for (j in 1:length(t)) {  # number of communities times
    t1 <- stdFreqs[,j+2]
    t[j] <- edist(t1, a1$q50)
  }
  output[i,] <- t
}
eds <- as_tibble(output)
```

```
## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
## This warning is displayed once per session.
```

```r
colnames(eds) <- colnames(stdFreqs)[3:9]
eds$community <- tolower(assembly_trials$community)
eds$assemblies_id <- assembly_trials$assemblies_id

prox <- tibble(ass_id = assembly_trials$assemblies_id, d5 = 0, px="", dx=0)
prox$d5 <-  eds$MG5a                                  # d5: distance to MG5
prox$dx <- apply(eds[,1:7], 1, min)                   # distance to proximate nvc
prox$px <- apply(eds[1:7], 1, ppx)                    # proximate nvc

prox <- prox %>% mutate(source = "survey")
# # prox$community to lower case, to emphasize that these are assessed, not standard.
prox$community <- tolower(assembly_trials$community)
species_counts <- species_hits %>% group_by(assemblies_id) %>% tally()
species_counts <- rename(species_counts, ass_id = assemblies_id)
prox <- left_join (prox, species_counts, by = "ass_id") 
# Mean d5 for MG5a
avg <- prox %>% filter(px=="MG5a") %>% summarise(avg=mean(d5))
 
# Change the display order
prox1 <-  prox %>% mutate(px1=factor(px,levels(factor(px))[c(2,1,3,4,6,5,7)]))

g <- ggplot(data = prox1, aes(x = px1, y = d5, colour = community, size = exp(n))) +
  ylim(0, NA) +
  geom_jitter(stat = 'identity', shape=16, width=0.1, height=0) +  # alt shapes: 1, 16
  scale_shape_identity() +
  scale_size(range = c(2, 10)) +
  scale_colour_brewer(palette="Set1") +
  # ggtitle("Differences from MG5a standard") +
  xlab("Nearest standard community") +
  ylab("Distance from MG5a standard") +
  labs(color = "Assessed NVC") +
  labs(size = "exp(Species count) \n range 2:17") +
  # Increase size of legend symbols
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme_grey()
print(g)
```

![](Notebook_files/figure-html/unnamed-chunk-2-1.png)<!-- -->
Moving on to the confusion plot Figure 5 in poster.
This is essentially the same code as prox8, just copied and pasted over. Then followed
by code to make the base plot figure 3 and the final confusion plot.


```r
# Libraries
library("RMySQL")
library(tidyverse)

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

# Euclidean distance function
edist <- function(s1,a1){
  x <- s1-a1
  sqrt(sum(x^2))
}

# Returns the name of the nearest standard NVC given r, a row
# of distance values for one assembly. Wrapped as a function 
# for ease of using apply
ppx <-  function(r)
{
  names(eds)[which(r == min(r))]
}

########################## MAIN ##############################
# GET DATA FROM DB
# Open connection to meadows data base
mydb = dbConnect(MySQL(), user='root', password='Mysql130641', dbname='meadows', host='127.0.0.1')

# Get nvc standard frequencies
q <- sprintf('SELECT * FROM meadows.mg_standards6;')
rs1 = dbSendQuery(mydb, q)
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 2 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 3 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 4 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 5 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 6 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 7 imported as
## numeric
```

```
## Warning in .local(conn, statement, ...): Decimal MySQL column 8 imported as
## numeric
```

```r
stdFreqs <- fetch(rs1, n=-1) 
# Make string of species names to use in selecting species hits
s <- toString(stdFreqs %>% select(species))
s <- substring(s, 4, nchar(s)-2)

# Get trials for assemblies
# Start with MG5a,c, and MG6a,b
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG5a", "MG5c", "MG6a", "MG6b")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(mydb, q)
trials_mg5_6 <- as_tibble(fetch(rs1, n=-1))

# Separate query for MG7b, c, d (DON'T exclude the few MG7c records)
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG7b", "MG7c", "MG7d")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(mydb, q)
trials_mg7 <- as_tibble(fetch(rs1, n=-1))
trials_mg7 <- trials_mg7 %>% select(assemblies_id, quadrat_count) %>% mutate(nvc = "MG7")

# Separate query for amalgamated MG1
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG1b", "MG1c", "MG1e")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(mydb, q)
trials_mg1 <- as_tibble(fetch(rs1, n=-1))
trials_mg1 <- trials_mg1 %>% select(assemblies_id, quadrat_count) %>% mutate(nvc = "MG1")

# Separate query for amalgamated MG10
q <- sprintf('select assemblies_id, nvc, quadrat_count from assemblies
             where nvc in ("MG10", "MG10a", "MG10b")
             and quadrat_count > 0;') # Two assemblies have 0 quadrat count

rs1 = dbSendQuery(mydb, q)
trials_mg10 <- as_tibble(fetch(rs1, n=-1))
trials_mg10 <- trials_mg10 %>% select(assemblies_id, quadrat_count) %>% mutate(nvc = "MG10")

# Join them together with bind_rows
assembly_trials <- bind_rows(trials_mg5_6, trials_mg7, trials_mg1, trials_mg10)
colnames(assembly_trials) <- c("assemblies_id", "community", "trials")

# Get hits for each species, separate query required for each subgroup.
# Amalgamated by assembly
# Start with MG5a,c and MG6a,b
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG5a", "MG5c", "MG6a", "MG6b")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)
rs1 = dbSendQuery(mydb, q)
hits_mg5_6 <- as_tibble(fetch(rs1, n=-1))

# Then MG7
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG7b", "MG7c", "MG7d")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)
rs1 = dbSendQuery(mydb, q)
hits_mg7 <- as_tibble(fetch(rs1, n=-1))
hits_mg7 <- hits_mg7 %>% select(assemblies_id, nvc, species_name, 'count(species.species_id)') %>% mutate(nvc = "MG7")

# Then MG1
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG1b", "MG1c", "MG1e")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)

rs1 = dbSendQuery(mydb, q)
hits_mg1 <- as_tibble(fetch(rs1, n=-1))
hits_mg1 <- hits_mg1 %>% select(assemblies_id, nvc, species_name, 'count(species.species_id)') %>% mutate(nvc = "MG1")

# Then MG10
q <- sprintf('select assembly_name, assemblies_id, nvc, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where nvc in ("MG10", "MG10a", "MG10b")
             and species.species_name in (\"%s\")
             group by assemblies_id, species_name;', s)
rs1 = dbSendQuery(mydb, q)
hits_mg10 <- as_tibble(fetch(rs1, n=-1))
hits_mg10 <- hits_mg10 %>% select(assemblies_id, nvc, species_name, 'count(species.species_id)') %>% mutate(nvc = "MG10")

# Join them together with bind_rows
species_hits <- bind_rows(hits_mg5_6, hits_mg7, hits_mg1, hits_mg10)
colnames(species_hits) <- c("assembly", "assemblies_id", "community", "species", "hits")

assembly_hits <- species_hits %>% group_by(assemblies_id, community, species) %>% summarise(hits = sum(hits)) 

# Close the DB connection
dbDisconnectAll()
```

```
## Warning: Closing open result sets
```

```
## 1 connection(s) closed.
```

```r
# Ensure species not represented in a community get 0 hits
wide <- spread(assembly_hits, key=species, value=hits) # Coerces NAs
wide[is.na(wide)] <- 0 # Replaces NAs with 0
# And reconstitute assembly hits. Good general solution to the problem
# assembly_hits <- wide %>% gather(Anthoxanthum_odoratum:Trifolium_repens, key= species, value=hits)
assembly_hits <- wide %>% gather(colnames(wide[3:length(colnames(wide))])
, key= species, value=hits)

# Survey frequencies: join hits and trials and mutate to calculate frequencies
srvf <- (full_join(assembly_trials, assembly_hits, by = c("assemblies_id", "community")) 
         %>% mutate(freq = hits/trials)
         # And mutate to get 0.05, 0.5 and 0.95 quantiles (use 0.05 and 0.95 for CrI)
         %>%  mutate(CrI5 = qbeta(0.05,hits+1, 1+trials-hits))
         %>%  mutate(q50 = qbeta(0.5,hits+1, 1+trials-hits)) # For comparison with frequency as hits/trials
         %>%  mutate(CrI95 = qbeta(0.95,hits+1, 1+trials-hits))
         %>%  mutate(source="survey")) # Note pipe enclosed in () to ensure run onto new line

# Distance measures
output <- matrix(nrow=204, ncol=7) # The numer of analysis groups - we'll have 7 columns
t <-  rep(0,7)

for (i in 1:length(assembly_trials$assemblies_id)) {
  # a1: 22 values of q50 (mean)
  a1 <- srvf %>% filter(assemblies_id == assembly_trials$assemblies_id[i]) %>% select(q50)
  for (j in 1:length(t)) {  # number of communities times
    t1 <- stdFreqs[,j+2]
    t[j] <- edist(t1, a1$q50)
  }
  output[i,] <- t
}
# eds: Euclidean Distances to Survey data from standards  
eds <- as_tibble(output)
# Note use of actual column numbers here and below - should be replaced by lengths
colnames(eds) <- colnames(stdFreqs)[3:9] 

prox <- tibble(ass_id = assembly_trials$assemblies_id, d5 = 0, px="", dx=0)
prox$d5 <-  eds$MG5a                            # d5: distance to MG5
prox$dx <- apply(eds, 1, min)                   # distance to proximate nvc
prox$px <- apply(eds, 1, ppx)                   # proximate nvc

prox <- prox %>% mutate(source = "survey")
# prox$community to lower case, to emphasize that these are assessed, not standard.
prox$community <- tolower(assembly_trials$community)
species_counts <- species_hits %>% group_by(assemblies_id) %>% tally()
species_counts <- rename(species_counts, ass_id = assemblies_id)
prox <- left_join (prox, species_counts, by = "ass_id") 
# Mean d5 for MG5a
avg <- prox %>% filter(px=="MG5a") %>% summarise(avg=mean(d5))


# BASE DIAGRAM STANDARDS:
ref <- tibble(
  n = 0:6, 
  x = cos(n*0.897),
  y = sin(n*0.897)
)

# Determines order of nodes
ref <- ref %>%  mutate(lab=c("MG6b", "MG6a", "MG7", "MG10", "MG1", "MG5c", "MG5a"))

# Set which nodes are linked
start <- c("MG5a", "MG5a", "MG5a", "MG6b", "MG6b", "MG6a", "MG7", "MG6a")
end <- c("MG1", "MG5c", "MG6b", "MG10", "MG6a", "MG7", "MG6a", "MG10")
links <- tibble(start = start, end = end, x=0, y=0, xend=0, yend=0)
# links <- links %>% mutate(labs = c("neglect", "acid", "", "waterlogging", "", "seeding", "", "waterlogging"))

for (i in 1:length(start)){
  t <- ref %>% filter(lab==start[i])
  links$x[i] <- t$x
  links$y[i] <- t$y
  t <- ref %>% filter(lab==end[i])
  links$xend[i] <- t$x
  links$yend[i] <- t$y
}

# Base diagram plot
g3 <- ggplot(data=ref, aes(x=x, y=y)) + 
  ylim(-1.5, 1.5) + xlim(-1.5, 1.5) +
  coord_fixed() +
  geom_segment(data=links, aes(x = links$x, y = links$y, xend = links$xend, yend = links$yend), 
      size=6, colour="green2") +
  geom_point(shape = 16, colour = "grey50", size = 25) +
  # shapes16, 21
  geom_text(aes(label=lab),hjust=0.5, vjust=0, size=5, colour="white") +
  theme_void() #+ theme(legend.position="none")
print(g3)
```

![](Notebook_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
## Confusion matrix. Clever, eh?
crosstab1 <- (prox %>%
                group_by(px, community)%>%
                summarise(n=n())%>%
                spread(community, n))

crosstab2 <- (prox %>%
                group_by(community, px)%>%
                summarise(n=n())%>%
                spread(px, n))

ct2 <- crosstab2 %>% replace(is.na(.), 0)
ct1 <- crosstab1 %>% replace(is.na(.), 0)
m1 <- as.matrix(ct1[,2:8])
m2 <- as.matrix(ct2[,2:8])
gross_flow <- m2 + m1 # or m2 + t(m2)

gf <- as_tibble(gross_flow)
standard <- rep(colnames(gf),7)
gf <- gf %>% mutate(lab = colnames(gf))

# Make expected values
chi2 <- chisq.test(m2)
```

```
## Warning in chisq.test(m2): Chi-squared approximation may be incorrect
```

```r
chi1 <- chisq.test(m1)
```

```
## Warning in chisq.test(m1): Chi-squared approximation may be incorrect
```

```r
corr2 <- chi2$observed/chi2$expected
corr1 <- chi1$observed/chi1$expected
xflow <- corr2 + corr1 # expected flow/confusion. Use to select which links to draw
xf <- as_tibble(xflow)
xf <- xf %>% mutate(lab = colnames(xf))
# Make the xf matrix long and add community  names 
xf_long <-  xf %>% gather(MG1:MG7, key="lab", value="xratio") %>% mutate(standard = standard)

ref2 <- ref %>% rename(standard = lab)

# Make the gf matrix long and add community names
gf_long <-  gf %>% gather(MG1:MG7, key="lab", value="gross") %>% mutate(standard = standard)
# FULL join gf_long with ref - "lab" - gives slow changing "row" names
flows <- full_join(gf_long, ref, by = "lab") %>% select(-n)
# joining with ref2 "standard" gives fast changing "column" names
flows <-left_join(flows, ref2, by="standard") %>% select(-n)
flows <- flows %>% rename(x=x.x) %>% rename(y=y.x) %>% rename(xend=x.y) %>% rename(yend=y.y)
# Join in the obs/expected ratios, xratio.
flows <- left_join(flows, xf_long, by = c("lab", "standard"))

# Filter unwanted values: only include cross-tabs, reject xratio < 1. NOTE actually > 0.8 because
# Chisquare.test calculations of expected figures are not integer.
flows <- flows %>% filter(lab != standard) %>% filter(xratio>0.8)
comm_counts <- species_hits %>% group_by(community) %>% 
  summarise(count = n_distinct(assemblies_id)) %>%  ungroup() %>% rename(lab=community)
flows <- left_join(flows, comm_counts, by="lab") %>% rename(count_lab=count)

comm_counts <- rename(comm_counts, standard=lab)
flows <- left_join(flows, comm_counts, by="standard") %>% rename(count_std=count)
flows <-  flows %>% mutate(pool = count_lab + count_std) %>% mutate(gross_rel = 100*gross/pool)

# Confusion plot


g5 <- ggplot(data=flows, aes(x = x, y = y, xend = xend, yend = yend, size=gross_rel, colour=gross_rel)) +
  ylim(-1.5, 1.5) + xlim(-1.5, 1.5) +
  scale_colour_gradient(low = "grey99",high = "green2") +
  coord_fixed() +
  geom_segment() +
  geom_point(aes(fill=lab), shape = 21, size = 22, show.legend = F) +
  scale_fill_brewer(palette="Set1") +
  geom_text(aes(label=tolower(lab)),hjust=0.5, vjust=-0.6, size=5, colour="black") +
  geom_text(aes(label=count_lab),hjust=0.5, vjust=1, size=5, colour="black") +
  labs(colour = " exchanges\n% pool", size = "") +
  # Suppress the size legend, retaining only the colour legend
  scale_size(guide = 'none') +
  theme_void() + theme(legend.position="bottom", legend.box = "horizontal") 
print(g5)
```

![](Notebook_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

