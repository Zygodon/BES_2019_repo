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

########################## MAIN ##############################

# GET DATA FROM DB
# Open connection to meadows data base to get the data
# My local DB
#mydb = dbConnect(MySQL(), user='root', password='Mysql130641', dbname='meadows', host='127.0.0.1')
# Remote DB
mydb = dbConnect(MySQL(), user='guest', dbname='meadows', port = 3306, host='sxouse.ddns.net')
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
# 
