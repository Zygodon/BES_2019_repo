# Libraries
library("RMySQL")
library("RColorBrewer")
library(ggplot2)
library(dplyr)

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

# barplot -----------------------------------------------------------------


# GET BAR PLOT DATA FROM DB
# Open connection to meadows data base to get the data
# Remote DB
mydb = dbConnect(MySQL(), user='guest', password = 'guest', dbname='meadows', port = 3306, host='sxouse.ddns.net')
rs1 = dbSendQuery(mydb, "select nvc from assemblies where nvc is not null and nvc != major_nvc_community;")
data <- fetch(rs1, n=-1)
dbDisconnectAll()

# Tibble g: hand coded analysis groups - nvc with few samples grouped e.g. MG10.
# Also mires grouped - will be excluded from this analysis.
g <- tibble(grp= c(rep("M", 7), rep("MG10", 2), rep("MG1", 3), "MG5a", "MG5c", "MG6a", "MG6b", rep("MG7", 3), "MG9a"))
# grps: Distinct nvcs matched with their analysis groups.
grps <- bind_cols(data %>% distinct(nvc) %>% arrange(nvc), g)
# >left join attaches an analysis group to record for each sample.
d <- left_join (data, grps, by="nvc") %>% mutate_all(tolower)

# Assessed NVC bar plot
# Note hand-coded palette corresponding to analysis groups 
# e.g. "M" and "MG9", white - excluded from analysis.
bp1 <- ggplot(d, aes(x=nvc, fill=grp)) + 
  geom_bar(colour="black", stat="count") +
  scale_fill_manual(values = 
  c("white", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#ffff33", "#a65628", "white")) +
  coord_flip() +
  labs(x = "assessed NVC", fill = "analysis group")
print(bp1)

# 



# Dissimilarities (euclidean distance) assessed nvc - standard nvc --------
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

# Ensure species not represented in a community get 0 hits
# NOTE: Hadley Wickham has something to say about this problem somewhere ...
wide <- spread(assembly_hits, key=species, value=hits) # Coerces NAs
wide[is.na(wide)] <- 0 # Replaces NAs with 0
# And reconstitute assembly hits. Good general solution to the problem
assembly_hits <- wide %>% gather(colnames(wide[3:length(colnames(wide))]), key= species, value=hits)

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


