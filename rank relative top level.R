# Not currently working. Needs updating getting species ids - table mg_constants
# doesn't exist any more. Use q <- sprintf('SELECT * FROM meadows.mg_standards6;').
# Probably other inconsistencies too.

# Libraries
library("RMySQL")
library(tidyverse)
library(ggrepel)
# library(factoextra) # may need later

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}


########################## MAIN ##############################

# GET DATA FROM DB
# Open connection to meadows data base
mydb = dbConnect(MySQL(), user='root', password='Mysql130641', dbname='meadows', host='127.0.0.1')

# Get nvc standard frequencies
q <- sprintf('SELECT * FROM meadows.mg_standards;')
rs1 = dbSendQuery(mydb, q)
stdFreqs <- fetch(rs1, n=-1)
# Useful species list
consts <- tibble(colnames(stdFreqs[3:18]))
colnames(consts) <- "species"
# Gather stdFreqs into more useful form and add columns for compatibility with srvf (survey freqs)
stdf <-( stdFreqs %>% gather(Anthoxanthum_odoratum:Trifolium_repens, key = species, value = "freq")
%>% select(-stdfreqs_id)
%>% mutate(CrI5 = 0, q50=0, CrI95=0, source="standard"))

# Get trials for the top_level nvc
q <- sprintf('select major_nvc_community, quadrat_count from assemblies
             where major_nvc_community in ("MG1", "MG5", "MG6", "MG7", "MG10")
             and major_nvc_community != nvc;')

rs1 = dbSendQuery(mydb, q)
quadrat_counts <- as_tibble(fetch(rs1, n=-1))
colnames(quadrat_counts) <- c("top_level", "q_count")
# Sum of quadrat counts for each top_level
top_level_trials <- quadrat_counts %>% group_by(top_level) %>% summarise(trials = sum(q_count))

# Get hits for each species and for each species find total hit count for each top level
# Grouped by major_nvc (top_level)
q <- sprintf('select major_nvc_community, species_name, count(species.species_id)
             from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where major_nvc_community in ("MG1", "MG5", "MG6", "MG7", "MG10")
             and major_nvc_community != nvc
             and species.species_id in (select species_id from mg_constants)
             group by major_nvc_community, species_name;')

rs1 = dbSendQuery(mydb, q)
species_hits <- as_tibble(fetch(rs1, n=-1))
colnames(species_hits) <- c("top_level", "species", "hits")
top_level_hits <- species_hits %>% group_by(top_level, species) %>% summarise(hits = sum(hits)) 

# Ensure species not represented in a community get 0 hits
wide <- spread(top_level_hits, key=species, value=hits) # Coerces NAs
wide[is.na(wide)] <- 0 # Replaces NAs with 0
# And reconstitute top_level hits. Good general solution to the problem
top_level_hits <- wide %>% gather(Anthoxanthum_odoratum:Trifolium_repens, key= species, value=hits)

# Survey frequencies: join hits and trials and mutate to calculate frequencis
srvf <- (full_join(top_level_trials, top_level_hits, by = "top_level") 
  %>% mutate(freq = hits/trials) 
  # And mutate to get 0.05, 0.5 and 0.95 quantiles (use 0.05 and 0.95 for CrI)
  %>%  mutate(CrI5 = qbeta(0.05,hits+1, 1+trials-hits))
  %>%  mutate(q50 = qbeta(0.5,hits+1, 1+trials-hits)) # For comparison with frequency as hits/trials
  %>%  mutate(CrI95 = qbeta(0.95,hits+1, 1+trials-hits))
  %>%  mutate(source="survey")) # Note pipe enclosed in () to ensure run onto new line

# Join tables preparatory to sorting on standard freqs for each top_level group
joined <- (left_join(stdf, srvf, by=c("top_level", "species"))
  %>% group_by(top_level) # Group on op_level so each top_level sorted correctly
  %>% arrange(freq.x, .by_group = TRUE)) # and then sort.

# Split the joined and sorted data preparatory to binding end-to end. Note not all columns
# included in the split. Just what are needed for plot.
split_srv <- select(joined, top_level, species, freq.y, CrI5.y, CrI95.y, source.y)
split_std <- select(joined, top_level, species, freq.x, CrI5.x, CrI95.x, source.x)
# Fix column names
colnames(split_srv) <- c("top_level", "species", "freq", "CrI5", "CrI95", "source")
colnames(split_std) <- c("top_level", "species", "freq", "CrI5", "CrI95", "source")

# Join the split data end-to end ready for plotting
df <- bind_rows(split_std, split_srv)

# Add necessary information for x axis (vertical) scale: std consts in red.
df <- mutate(df, std=ifelse((freq >= 0.7 & source == "standard"), "red", "black"))
# Fix some oddities - Plantago_lanceolata in MG1 standard should be red
df$std[121] <- "red"

# Cannot currently use facets with different scales.
# So instead of faceting must make five separate plots
tops <- df %>% distinct(top_level)
for(t in tops$top_level){
  dm <- filter(df, top_level==t)
  lim <- dm[[2]]
  p <- ggplot(dm, aes(x=species, y=freq, fill=source)) +
    scale_x_discrete(limits=lim[1:length(consts$species)]) + # Ensures species (x-axis) plotted in same order as df
    coord_flip() +
    geom_bar(stat="identity", position=position_dodge2(reverse=T)) +
    geom_errorbar(aes(ymin=CrI5, ymax=CrI95), width=.2,
                  position=position_dodge(-0.9)) +
    ggtitle(t)
  print(p + theme(axis.text.y = element_text(colour=dm$std)) + scale_fill_manual(values = c("tan2", "steelblue2")))
}




