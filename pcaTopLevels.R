# Not currently working. Needs updating getting species ids - table mg_constants
# doesn't exist any more. Use q <- sprintf('SELECT * FROM meadows.mg_standards6;').
# Probably other inconsistencies too, especially: no longer using Reshape;
# other changes for Tidyverse.

# Libraries
library("RMySQL")
library(reshape2)
library(corrplot)
library(ggplot2)
library(deldir)
# library(colorspace) # May be needed later

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

ptSize <- function(z, flr, range, offset, scl){  # e.g. offset2, scl 10
  z1 <-(z-flr)/range
  (z1*scl) + offset
}

########################## MAIN ##############################

# GET DATA FROM DB
# Open connection to meadows data base
mydb = dbConnect(MySQL(), user='root', password='Mysql130641', dbname='meadows', host='127.0.0.1')

# Extract records of mg_constants with assessed NVC. Include quadrat count
q <- sprintf('select assemblies_id, major_nvc_community, nvc, species_name, count(species_name), 
             quadrat_count from assemblies
             join quadrats on quadrats.assembly_id = assemblies_id
             join records on records.quadrat_id = quadrats_id
             join species on species.species_id = records.species_id
             where major_nvc_community in ("MG1", "MG5", "MG6", "MG7", "MG10")
             and major_nvc_community != nvc
             and species.species_id in (select species_id from mg_constants)
             group by assemblies_id, species_name;')
rs1 = dbSendQuery(mydb, q)
dbData <- fetch(rs1, n=-1)

# Get nvc standard frequencies
q <- sprintf('SELECT * FROM meadows.mg_standards;')
rs1 = dbSendQuery(mydb, q)
stdFreqs <- fetch(rs1, n=-1)
dbDisconnectAll()

# Cosmetics
stdFreqs <- stdFreqs[,-1]

colnames(dbData)[2] <- "top_level"
colnames(dbData)[5] <- "records"
dbData$freq <- dbData$records/dbData$quadrat_count

d <- dbData[,c(2,4,7)]
d.melt <- melt(d, id = 1:2)
d.cast <- cast(d.melt, top_level ~ species_name, mean)
d.cast[is.na(d.cast)] <- 0
M <- cor(d.cast)
corrplot(M)
attributes(stdFreqs) <- attributes(d.cast)
M2 <- cor(stdFreqs)
corrplot(M2)

# PCA on DB
# pcaDb <- prcomp(d.cast, center = TRUE,scale. = TRUE)

# PCA on standards
pcaStd <- prcomp(stdFreqs, center = TRUE,scale. = TRUE)

# summary(pcaDb)
summary(pcaStd)
# screeplot(pcaDb,type="lines") # scree plot
screeplot(pcaStd,type="lines") # scree plot
# biplot(pcaDb)
biplot(pcaStd)

stdxyz <- as.data.frame(pcaStd$x[, 1:3])

# Project nvc classes assigned from data(DB) on to PCs of nvc standards
predDbFromStd <- predict(pcaStd, newdata=d.cast[, 2:17])
plot(predDbFromStd[,1:2], xlim=c(-5,5), ylim=c(-5,5))
text(x=predDbFromStd[,1], y=predDbFromStd[,2],labels=d.cast$top_level)
points(stdxyz[,1:2], col="red")
text(x=stdxyz[,1], y=stdxyz[,2],labels=stdFreqs$top_level, col="red")

surveyxyz <- as.data.frame(predDbFromStd[,1:3])

zVals <- c(stdxyz[,3], surveyxyz[,3])
zFloor <- min(zVals)
zRange <- max(zVals - min(zVals))


p1 <- ggplot(data = stdxyz[,1:2], aes(x=stdxyz[,1], y=stdxyz[,2], color=rownames(stdxyz))) +
  xlim(-5, 5) + ylim(-5, 5) +
  geom_segment(x=stdxyz[,1], xend=surveyxyz[,1], y=stdxyz[,2], yend=surveyxyz[,2],
               size=ptSize(stdxyz[,3], zFloor, zRange, 1, 1)) + 
  geom_point(x=surveyxyz[,1], y=surveyxyz[,2],color="red", size = ptSize(surveyxyz[,3], zFloor, zRange, 3, 12), shape=19) + #survey points bg
  geom_point(x=surveyxyz[,1], y=surveyxyz[,2],size = ptSize(surveyxyz[,3], zFloor, zRange, 3, 10), shape=19) + #survey points fill
  geom_point(color="grey40", size = ptSize(stdxyz[,3], zFloor, zRange, 3, 12), shape=19) + #std points bg
  geom_point(size = ptSize(stdxyz[,3], zFloor, zRange, 3, 10), shape=19)  #std points fill

print(p1) # Needed to ensure the plot is displayed when code run from "Source".

# ptScale <- floor(10*(max(predDbFromStd[,3]) - min(predDbFromStd[,3])))
# p1 <- ggplot(data = locsStd, aes(x = locsStd$PC1, y =locsStd$PC2)) +
#   xlim(-5, 5) + ylim(-5, 5) +
#   geom_segment(x=locsStd[,1], xend=predDbFromStd[,1], y=locsStd[,2], yend=predDbFromStd[,2]) +
#   # geom_point(color = "grey40", fill = "grey40", shape=21,
#   #            size = 2 + (predDbFromStd[,3] - min(predDbFromStd[,3]))* ptScale) +
#   geom_point(aes(color = rownames(locsStd)), fill = "grey40", shape=21,
#              size = 2 + (predDbFromStd[,3] - min(predDbFromStd[,3]))* ptScale) +
#   # geom_point(x=predDbFromStd[,1], y=predDbFromStd[,2], color = "red", fill="grey", shape=21,
#   #            size = 2 + (predDbFromStd[,3] - min(predDbFromStd[,3]))* ptScale)
#   geom_point(x=predDbFromStd[,1], y=predDbFromStd[,2], aes(color = rownames(locsStd)), fill="red", shape=21,
#              size = 2 + (predDbFromStd[,3] - min(predDbFromStd[,3]))* ptScale)
# print(p1) # Needed to ensure the plot is displayed when code run from "Source".


# Proceed to Voronoi on standards
# voronoi <- deldir(locsStd[,1], locsStd[,2])



# Correct approach is to Project the field results on to the standards
# ggplot(data = locsStd, aes(x = locsStd$PC1, y =locsStd$PC2)) +
#   xlim(-5, 5) + ylim(-5, 5) +
#   geom_point(aes(color = rownames(locsStd))) + # TO DO melt and cast on mean Species_max_abundance
#   geom_text(aes(label=rownames(locsStd)),hjust=0, vjust=0) +
#   geom_point(aes(x=predDbFromStd[,1], y=predDbFromStd[,2], color="red")) +
#   geom_text(aes(x=predDbFromStd[,1], y=predDbFromStd[,2], label=rownames(locsDb), color="red"))

# ggplot(data = locsStd, aes(x = locsStd$PC1, y =locsStd$PC2)) +
#   xlim(-5, 5) + ylim(-5, 5) +
#   geom_point(aes(color = rownames(locsStd))) + # TO DO melt and cast on mean Species_max_abundance
#   geom_text(aes(label=rownames(locsStd)),hjust=0, vjust=0) +
#   geom_point(aes(x=predDbFromStd[,1], y=predDbFromStd[,2], color="red")) +
#   geom_text(aes(x=predDbFromStd[,1], y=predDbFromStd[,2], label=rownames(locsDb), color="red")) +
#   #Plot the voronoi lines
#   geom_segment(
#     aes(x = x1, y = y1, xend = x2, yend = y2),
#     size = 2,
#     data = voronoi$dirsgs,
#     linetype = 1,
#     color= "#FFB958")

  

# library(xlsx)
# write.xlsx(stdFreqs, "stdfreqs.xlsx")

