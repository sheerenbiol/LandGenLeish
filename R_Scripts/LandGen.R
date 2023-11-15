#########################################################################
##### Diversity and dissemination of viruses in pathogenic protozoa #####
#-----------------------------------------------------------------------#
#####                   Landscape genomic analyses                  #####
#########################################################################
### Summary:
## - 0. Set up
## - 1. Importing & Filtering data
## - 2. Calculating distances
## - 3. Isolation by distance
## - 4. RDA-based variable selection
## - 5. RDA variation partitioning
## - 6. GDM modelling

#####################-----------------------------------------------------------
##### 0. SET UP #####
#####################

### Set working directory:
setwd('PATH/TO/LandGenLeish/')

### Libraries:
##############
library(data.table); library(readxl); library(geodist); library(adegenet); library(vegan)
library(raster); library(UpSetR); library(ggplot2); library(ggpubr); library(gdm)

### Functions:
##############
###  Read in genotype file:
read.geno <- function(file) {
  # INPUT: file = 012 output file from vcftools
  geno <- fread(file, data.table = F, header = F)[,-1]
  rownames(geno) <- as.character(read.table(paste(file, 'indv', sep='.'))[,1]) #indv names as rownames
  genopos <- read.table(paste(file, 'pos', sep='.'))
  colnames(geno) <- as.character(paste(genopos[,1], genopos[,2], sep=';')) #Chr.pos as colnames
  return(as.data.frame(geno))
}

#########################################---------------------------------------
##### 1. IMPORTING & Filtering DATA #####
#########################################
### Meta data:
info.samples <- read_excel('./Input_files/Metadata/Table1_LRV_isolates_SH.xlsx', sheet=1)
info.samples <- info.samples[order(info.samples$`Analyses Code`),]

### Genotype data:
geno <- read.geno('./Input_files/Genotype_data/LBR.77.GENO.SNP.GATK.QUAL250.CWS.fDP10.fGQ25.PASS.NoNA.IND76.012')

### Removing Clonal isolates (one isolate per clonal group remains in dataset):
exclude.list <- c('CUM97A1', 'CUM180A1','CUM637A1', 'CUM42A1','CUM49-CH55A1',
                  'CUM181A1', 'CUM700A1','CUM152A1','LC1409A1','CUM29A1', 
                  'LC2368A1','PER010A1')
info.samples.NoClonal <- as.data.frame(info.samples)
rownames(info.samples.NoClonal) <- info.samples.NoClonal$`Analyses Code`
for (i in exclude.list) {
  info.samples.NoClonal <- info.samples.NoClonal[rownames(info.samples.NoClonal) != i,]
  geno <- geno[rownames(geno) != i,]
}

### Filter out admixed isolates as inferred by ADMIXTURE, fineSTRUCTURE and PCAdmix:
info.3 <- subset(info.samples.NoClonal, info.samples.NoClonal$Pop != 'ADM') ## ADM contains here ADM and UNC
info.3 <- subset(info.3, info.3$Pop != 'STC') ## Also removing CEN002 and CEN007 --> also showing mixed ancestry

### Filter out isolates without coordinate data:
info.3 <- info.3[which(info.3$Latitude != 'NA'),]

### Retain the same isolates in the genotype data:
rownames(info.3) <- info.3$`Analyses Code` ## To match rownames between the two data objects
geno.3 <- geno[rownames(geno) %in% rownames(info.3),]


#####################################-------------------------------------------
##### 2. CALCULATING DISTANCES  #####
#####################################

##### A. Geographical distance:
###############################
### Collecting coordinates from meta-data:
Coords <- info.3[,c(10,9)]; Coords <- as.data.frame(Coords)
Coords$Longitude <- as.numeric(unlist(Coords$Longitude)); Coords$Latitude <- as.numeric(unlist(Coords$Latitude))
rownames(Coords) <- info.3$`Analyses Code`

### Creating geographical distance matrix | method of distance: https://www.r-bloggers.com/2010/11/great-circle-distance-calculations-in-r/
geograph.mtrx <- geodist(Coords, measure = 'haversine')
rownames(geograph.mtrx) <- rownames(Coords); colnames(geograph.mtrx) <- rownames(Coords)
geograph.mtrx <- geograph.mtrx/1000 #convert to km's


##### B. Genetic distance: 
##########################
### Conversion to genind object:
geno.3.genind <- df2genind(geno.3, ploidy = 2, ind.names = rownames(geno.3), sep = '\t'); geno.3.genind

### Using the inter-individual genetic distance metric (Bray-Curtis Dissimilarity)
### =  proportions of alleles that are different between pairs of individuals
gendist.BCD <- vegdist(geno.3.genind, method = 'bray')
gendist.BCD.mtrx <- as.matrix(gendist.BCD)

##### C. Bioclimatic variables - WorldClim:
##### PRESENT-DAY CLIMATE DATA:
###########################################
### Exctraction from WorldClim Db
w <- getData('worldclim', var='bio', res=0.5, lon=-69, lat=-12). ## Will be saved as folder with name 'wc0.5'
bioclim_vars <- info.3$`Analyses Code`
bioclim_vars<- as.data.frame(bioclim_vars)
bioclim_vars$Sample <- bioclim_vars$bioclim_vars
for (i in w@layers) {
  newcol <- extract(x= i, y=cbind(as.numeric(info.3$Longitude),as.numeric(info.3$Latitude)))
  bioclim_vars[,ncol(bioclim_vars) + 1] <- newcol
  colnames(bioclim_vars)[ncol(bioclim_vars)] <- paste0(i)
}
bioclim_vars <- bioclim_vars[,-1]
colnames(bioclim_vars) <- c('Sample', 'bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')

### Creating a distance matrix for each bioclim variable
for (bio in colnames(bioclim_vars)[2:20]) {
  print(bio)
  tmp <- bioclim_vars[,bio]
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- bioclim_vars$Sample
  tmp.dist <- dist(tmp, method = 'euclidean')
  assign(paste0(bio,'.dist'), tmp.dist)
  assign(paste0(bio,'.mtrx'), as.matrix(tmp.dist))
}

#####################################-------------------------------------------
##### 3. ISOLATION BY DISTANCE  #####
#####################################
### Combining genetic and geographic distance into 1 data frame:
distances <- data.frame(geogr=geograph.mtrx[lower.tri(geograph.mtrx)],
                        genet=gendist.BCD.mtrx[lower.tri(gendist.BCD.mtrx)])

#### A. Regression of the interdeme great-circle distance and W.C.'s Fst of the three ancestral populations:
############################################################################################################
### Fst Matrix - based on the Weir-Cockerham's Fst as calculated by vcftools:
genpops <- matrix(nrow = 3, ncol = 3) ## PAU INP HUP
genpops[,1] <- c(0,0.091, 0.103)
genpops[,2] <- c(0.091,0,0.093)
genpops[,3] <- c(0.103,0.093,0)

### average geographic coordinates per population
### PAU
median(as.numeric(info.3[13:30,10]))
median(as.numeric(info.3[13:30,9]))
### INP
median(as.numeric(info.3[1:11,10]))
median(as.numeric(info.3[1:11,9]))
###HUP
median(as.numeric(info.3[c(11:12,31:38),10]))
median(as.numeric(info.3[c(11:12,31:38),9]))

geopop <- c('PAU','INP','HUP'); geopop <- as.data.frame(geopop)
geopop$Lat <- c(-12.90953,-16.49511, -9.098772 )
geopop$Lon <- c(-71.40449,-65.82669, -74.96517)

## Creating geographical distance matrix | method of distance: https://www.r-bloggers.com/2010/11/great-circle-distance-calculations-in-r/
geopop.dist <- geodist(geopop, measure = 'haversine')
geopop.mtrx <- as.matrix(geopop.dist)
geopop.mtrx <- geopop.mtrx/1000 #convert to km's

distances.pops <- data.frame(geogr=geopop.mtrx[lower.tri(geopop.mtrx)],gen=genpops[lower.tri(genpops)])
##INP-PAU, HUP-PAU, INP-HUP
distances.pops$col <- c('maroon4', 'maroon3', 'maroon')

ggplot(distances.pops, aes(x=geogr, y=gen,colour=factor(col))) + 
  geom_point(cex=3)+scale_colour_identity()+ 
  geom_smooth(method = 'lm', formula = y~x, col='black', se=F) + xlim(500,1300)+ylim(0.0,0.15)+
  stat_cor(data=distances.pops, method = "pearson", label.x = 600, label.y = 0.12, colour='black')+
  xlab('Geographic distance (km)') + ylab('Genetic distance (Fst)')+  #haversine great circle distance & Bray curtis genetic dissimilarity
  theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
  #stat_summary_bin(fun.data='mean_cl_boot', bins = 30,colour='steelblue', size=0.75)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=25))

#### B. Inter-individual regression of great-circle distance with pairwise genetic distance:
############################################################################################
## Loess regression (Loess value = 1.2):
ggplot(distances, aes(x=geogr, y=genet)) + 
  geom_smooth(method = 'loess', formula = y~x, col='cornflowerblue', span=1.2) + xlab('Geographic distance class (10km)') + ylab('Genetic distance (bray)')+  #haversine great circle distance & Bray curtis genetic dissimilarity
  theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
  ggtitle("loess span = 1.2") +
  #stat_cor(method = "pearson", label.x = 3, label.y = 0.2)+
  stat_summary_bin(fun.data='mean_cl_boot', bins = 150,colour='steelblue', size=0.75)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=25))

## Linear regression of intra- and inter population genetic vs. geographic distance:
N <- nrow(gendist.BCD.mtrx) ## Number of rows
## Re-shape matrices into data frames:
gendist.BCD.df<- reshape2::melt(gendist.BCD.mtrx)[seq(from=1, to=N^2,by=1),]
geograph.df<- reshape2::melt(geograph.mtrx)[seq(from=1, to=N^2,by=1),]

## Merge data frames:
distances.df <- geograph.df; distances.df <- cbind(geograph.df, gen=gendist.BCD.df$value)
colnames(distances.df) <- c('IND1', 'IND2', 'geo', 'gen')

## add populations to dataframe
distances.df$POP1 <- c(rep(info.3$Pop2, 38))
distances.df$POP2 <- c(rep('INP',38), rep('INP',38),rep('INP',38),rep('INP',38),rep('INP',38),rep('INP',38),
                       rep('INP',38),rep('INP',38),rep('INP',38),rep('INP',38),rep('INP',38),rep('HUP',38),rep('PAU',38),
                       rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),
                       rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('PAU',38),
                       rep('PAU',38),rep('PAU',38),rep('PAU',38),rep('HUP',38),rep('HUP',38),rep('HUP',38),rep('HUP',38),
                       rep('HUP',38),rep('HUP',38),rep('HUP',38),rep('HUP',38))
## add same/different population...
distances.df$same <- c(rep(c(rep('INP',11),rep('INP-HUP',1), rep('INP-PAU',18), rep('INP-HUP',8)),11), 
                       rep(c(rep('INP-HUP',11),rep('HUP',1), rep('HUP-PAU',18), rep('HUP',8)),1),
                       rep(c(rep('INP-PAU',11),rep('HUP-PAU',1), rep('PAU',18), rep('HUP-PAU',8)),18),
                       rep(c(rep('INP-HUP',11),rep('HUP',1), rep('HUP-PAU',18), rep('HUP',8)),8))
distances.df$same <- as.factor(distances.df$same)
##  Color code            ##
##  INP -> #D95F02        ##
##  PAU -> #1B9E77        ##
##  HUP -> #7570B3        ##
##  INP-HUP -> "maroon"   ##
##  HUP-PAU -> "maroon3"  ##
##  INP-PAU -> "maroon4"  ##
## add colours:
distances.df$sameCol <- c(rep(c(rep('#D95F02',11),rep("maroon",1), rep( "maroon4",18), rep("maroon",8)),11), 
                          rep(c(rep("maroon",11),rep('#7570B3',1), rep("maroon3",18), rep('#7570B3',8)),1),
                          rep(c(rep( "maroon4",11),rep("maroon3",1), rep('#1B9E77',18), rep("maroon3",8)),18),
                          rep(c(rep("maroon",11),rep('#7570B3',1), rep("maroon3",18), rep('#7570B3',8)),8))

## subset - within population:
tmp <- subset(distances.df, distances.df$same == 'INP')
within.pop <- tmp
tmp <- subset(distances.df, distances.df$same == 'PAU')
within.pop <- rbind(within.pop, tmp)
tmp <- subset(distances.df, distances.df$same == 'HUP')
within.pop <- rbind(within.pop, tmp)

within.pop <- within.pop[within.pop$gen != 0,]

## subset - between populations
tmp <- subset(distances.df, distances.df$same == 'INP-HUP')
between.pop <- tmp
tmp <- subset(distances.df, distances.df$same == 'HUP-PAU')
between.pop <- rbind(between.pop, tmp)
tmp <- subset(distances.df, distances.df$same == 'INP-PAU')
between.pop <- rbind(between.pop, tmp)

between.pop <- between.pop[between.pop$gen != 0,]

### Make plot with distance classes:
INP <- subset(distances.df, distances.df$same == 'INP');INP <- INP[INP$gen != 0,]
PAU <- subset(distances.df, distances.df$same == 'PAU');PAU <- PAU[PAU$gen != 0,]
HUP <- subset(distances.df, distances.df$same == 'HUP');HUP <- HUP[HUP$gen != 0,]
INPHUP <- subset(distances.df, distances.df$same == 'INP-HUP');INPHUP <- INPHUP[INPHUP$gen != 0,]
INPPAU <- subset(distances.df, distances.df$same == 'INP-PAU');INPPAU <- INPPAU[INPPAU$gen != 0,]
HUPPAU <- subset(distances.df, distances.df$same == 'HUP-PAU');HUPPAU <- HUPPAU[HUPPAU$gen != 0,]

ggplot(data=distances.df,aes(x=geo, y=gen, colour=factor(sameCol))) + 
  #geom_point(data = as.data.frame(INP)) + 
  #geom_point(data = as.data.frame(PAU))+
  #geom_point(data = as.data.frame(HUP)) + scale_colour_identity()+
  stat_summary_bin(data=INP, fun.data='mean_cl_boot', bins = 150,colour='#D95F02', size=0.75)+
  stat_summary_bin(data=PAU, fun.data='mean_cl_boot', bins = 150,colour='#1B9E77', size=0.75)+
  stat_summary_bin(data=HUP, fun.data='mean_cl_boot', bins = 150,colour='#7570B3', size=0.75)+
  geom_smooth(data = within.pop, method = 'lm', formula = y~x, col='darkgray',se=F)+
  stat_summary_bin(data=INPHUP, fun.data='mean_cl_boot', bins = 150,colour='maroon', size=0.75)+
  stat_summary_bin(data=HUPPAU, fun.data='mean_cl_boot', bins = 150,colour='maroon3', size=0.75)+
  stat_summary_bin(data=INPPAU, fun.data='mean_cl_boot', bins = 150,colour='maroon4', size=0.75)+
  geom_smooth(data=between.pop, method = 'lm', formula = y~x, col='darkgray',se=F)+
  xlab('Geographic distance class (10km)') + ylab('Genetic distance (bray)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=25))

### cor testen in de ggplot:
ggplot(data=distances.df,aes(x=geo, y=gen))+
  geom_smooth(data = within.pop, method = 'lm', formula = y~x, col='black',se=F)+
  geom_smooth(data=between.pop, method = 'lm', formula = y~x, col='black',se=F)+
  stat_cor(data=within.pop,method = "pearson", label.x = 3, label.y = 0.14)+
  stat_cor(data=between.pop,method = "pearson", label.x = 500, label.y = 0.17)


############################################------------------------------------
##### 4. RDA-BASED VARIABLE SELECTION  #####
############################################

#### Response variable:
########################
### Make the sure all loci are properly ordered
geno.3.order <- geno.3[,order(colnames(geno.3))] 
## Perform a Hellinger transformation on genetic data:
geno.3.order.hel <- decostand(geno.3.order, "hellinger")

#### Explanatory variables:
###########################
### Merging all variables into one object:
Num_vars <- cbind(bioclim_vars, Lat=info.3$Latitude, Lon=info.3$Longitude)
rownames(Num_vars) <- Num_vars$Sample; Num_vars <- Num_vars[,-1]
### Make sure longitude & latitude are numeric:
Num_vars$Lon <- as.numeric(Num_vars$Lon); Num_vars$Lat <- as.numeric(Num_vars$Lat)

### Standardization of bioclimatic variables:
Num_vars[,c(1:19)] <- decostand(Num_vars[,c(1:19)], 'standardize') 

#### Defining NULL and FULL model:
##################################
### NULL model:
RDA0.num <- rda(geno.3.order.hel ~ 1, Num_vars)

### FULL model:
RDAfull.num <- rda(geno.3.order.hel ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11+ bio12+ bio13+bio14+bio15+bio16+bio17+bio18+bio19, Num_vars)

#### AUTOMATIC stepwise variable selection:
###########################################
mod_a_full <- ordiR2step(RDA0.num, RDAfull.num, Pin = 0.05, R2scope = T, permutations = 1000, R2permutations = 1000, trace = T)
mod_a_full; vif.cca(mod_a_full); RsquareAdj(mod_a_full); extractAIC(mod_a_full)
## vif not OK; removing bio12, bio16, bio19 | RETAINING: bio3 & bio14
mod_a_full <- rda(geno.3.order.hel ~ bio3 + bio14, Num_vars); vif.cca(mod_All2.R2)
mod_a_full; RsquareAdj(mod_a_full); extractAIC(mod_a_full)

#### MANUAL stepwise variable selection:
########################################
## mod_1 <- rda(geno.3.order.hel ~ bioxx , Num_vars); mod_1; RsquareAdj(mod_1); vif.cca(mod_1); extractAIC(mod_1)
## Compare all options + choose best variable based on Rsq.adj + vif<10 + significance
## mod_2 <- rda(geno.3.order.hel ~ bioxx + bioyy , Num_vars); mod_2; RsquareAdj(mod_2); vif.cca(mod_2); extractAIC(mod_2)
## Compare all options + choose best variable based on Rsq.adj + vif<10 + significance
## ...
## !! anova.cca takes a long time to run with this data !!

## Uncomment to start testing all different combinations ##
#mod_1 <- rda(geno.3.order.hel ~ bio3 , Num_vars)
#mod_1; RsquareAdj(mod_1); vif.cca(mod_1)
#a_mod_1 <- anova.cca(mod_1); a_mod_1

#mod_2 <- rda(geno.3.order.hel ~ bio3 + bio14 , Num_vars)
#mod_2; RsquareAdj(mod_2); vif.cca(mod_2)
#a_mod_2 <- anova.cca(mod_2); a_mod_2

#mod_3 <- rda(geno.3.order.hel ~ bio3 + bio14 + bio18 , Num_vars)
#mod_3; RsquareAdj(mod_3); vif.cca(mod_3)
#a_mod_3 <- anova.cca(mod_3); a_mod_3

#mod_4 <- rda(geno.3.order.hel ~ bio3 + bio14 + bio18 + bio15 , Num_vars)
#mod_4; RsquareAdj(mod_4); vif.cca(mod_4)
#a_mod_4 <- anova.cca(mod_4); a_mod_4

#mod_5 <- rda(geno.3.order.hel ~ bio3 + bio14 + bio18 + bio15 + bio2, Num_vars)
#mod_5; RsquareAdj(mod_5); vif.cca(mod_5)
#a_mod_5 <- anova.cca(mod_5); a_mod_5

#mod_6 <- rda(geno.3.order.hel ~ bio3 + bio14 + bio18 + bio15 + bio2 + bio19, Num_vars)
#mod_6; RsquareAdj(mod_6); vif.cca(mod_6)
#a_mod_6 <- anova.cca(mod_6); a_mod_6

### ==> Conclusion: mod_5 with bio3, bio14, bio18, bio15, bio2: 
mod_m_full <- rda(geno.3.order.hel ~ bio3 + bio14 + bio18 + bio15 + bio2, Num_vars)
mod_m_full; RsquareAdj(mod_m_full); vif.cca(mod_m_full); extractAIC(mod_m_full)
a_mod_m_full <- anova.cca(mod_m_full); a_mod_m_full


##########################################--------------------------------------
##### 5. RDA VARIATION PARTITIONING  #####
##########################################

#### A. AUTOMATIC SELECTION:
############################

### Full model:
pRDA_auto <- rda(geno.3.order.hel ~ Lon + Lat + bio3 + bio14, Num_vars); pRDA_auto; RsquareAdj(pRDA_auto)
## Rsq-adj = 0.1849 | constrained = 27.3%
## anova.cca --> takes a very long time to run!!
a.full.auto <- anova.cca(pRDA_auto); a.full.auto

### Variable partitioning:
## Geography:
pRDA_auto.geo <- rda(geno.3.order.hel ~ Lon + Lat + Condition(bio3 + bio14), Num_vars);pRDA_auto.geo
## --> constrained: 7.5%
a.geo.auto <- anova.cca(pRDA_auto.geo); a.geo.auto

## Bioclimate (bio3 + bio14):
pRDA_auto.clim <- rda(geno.3.order.hel ~ bio3 + bio14 + Condition(Lon + Lat), Num_vars);pRDA_auto.clim
## --> constrained: 10.2%
a.clim.auto <- anova.cca(pRDA_auto.clim); a.clim.auto

## Bio3 only:
pRDA_auto.bio3 <- rda(geno.3.order.hel ~ bio3  + Condition(Lon + Lat+ bio14), Num_vars);pRDA_auto.bio3
## --> constrained: 4.5%
a.bio3.auto <- anova.cca(pRDA_auto.bio3); a.bio3.auto

## Bio14 only:
pRDA_auto.bio14 <- rda(geno.3.order.hel ~ bio14  + Condition(Lon + Lat+ bio3), Num_vars);pRDA_auto.bio14
## --> constrained: 2.4%
a.bio14.auto <- anova.cca(pRDA_auto.bio14); a.bio14.auto

### UpSet plot: (numbers are the constrained proportions of the pRDA outputs)
expressionInput <- c(`geo&bio3&bio14` = 273, geo = 75,`bio3&bio14` = 102, bio3 = 45, bio14 = 24)
upset(fromExpression(expressionInput), order.by = 'freq')

#### B. MANUAL SELECTION:
#########################

### Full model:
pRDA_man <- rda(geno.3.order.hel ~ Lon + Lat + bio3 + bio14 + bio18 + bio15 + bio2, Num_vars); pRDA_man; RsquareAdj(pRDA_man)
## Rsq-adj = 0.1973 | constrained = 34.9%
## anova.cca --> takes a very long time to run!!
a.full.man <- anova.cca(pRDA_man); a.full.man

### Variable partitioning:
### Manual model ###
## Geography:
pRDA_man.geo <- rda(geno.3.order.hel ~ Lon + Lat + Condition(bio3 + bio14+ bio18 + bio15 + bio2), Num_vars);pRDA_man.geo
## --> constrained: 5.3%
a.geo.man <- anova.cca(pRDA_man.geo); a.geo.man

## Bioclimate (bio3, bio14, bio18, bio15, bio2):
pRDA_man.clim <- rda(geno.3.order.hel ~ bio3 + bio14 + bio18 + bio15 + bio2 + Condition(Lon + Lat), Num_vars);pRDA_man.clim
## --> constrained: 17.8%
a.clim.man <- anova.cca(pRDA_man.clim); a.clim.man

## Bio3 only:
pRDA_man.bio3 <- rda(geno.3.order.hel ~ bio3  + Condition(Lon + Lat+ bio14+ bio18 + bio15 + bio2), Num_vars);pRDA_man.bio3
## --> constrained: 2.3%
anova.cca(pRDA_man.bio3)

## Bio14 only:
pRDA_man.bio14 <- rda(geno.3.order.hel ~ bio14  + Condition(Lon + Lat+ bio3+ bio18 + bio15 + bio2), Num_vars);pRDA_man.bio14
## --> constrained: 2.3%
anova.cca(pRDA_man.bio14)

## Bio18 only:
pRDA_man.bio18 <- rda(geno.3.order.hel ~ bio18  + Condition(Lon + Lat+ bio14+ bio3 + bio15 + bio2), Num_vars);pRDA_man.bio18
## --> constrained: 2.4%
anova.cca(pRDA_man.bio18)

## Bio15 only:
pRDA_man.bio15 <- rda(geno.3.order.hel ~ bio15  + Condition(Lon + Lat+ bio3+ bio18 + bio14 + bio2), Num_vars);pRDA_man.bio15
## --> constrained: 2.5%
anova.cca(pRDA_man.bio15)

## Bio2 only:
pRDA_man.bio2 <- rda(geno.3.order.hel ~ bio2  + Condition(Lon + Lat+ bio3+ bio18 + bio15 + bio14), Num_vars);pRDA_man.bio2
## --> constrained: 2.8%
anova.cca(pRDA_man.bio2)

### UpSet plot: (numbers are the constrained proportions of the pRDA outputs)
expressionInput <- c( `geo&bio3&bio14&bio18&bio15&bio2`=349,geo=54, bio3=23, bio14=23, bio18=24, bio15=25, bio2=28,`bio3&bio14&bio18&bio15&bio2`=178)
upset(fromExpression(expressionInput),nsets = 9, order.by = 'freq')

#####################################################---------------------------
##### 6. GDM - Genetic Dissimilarity Modelling  #####
#####################################################

#### Response variable:
########################
### Convert hellinger transformed geno data to genind and create Bray-Curtis dissimilarity matrix:
geno.3.genind <- df2genind(geno.3.order.hel, ploidy = 2, ind.names = rownames(geno.3.order.hel), sep = '\t'); geno.3.genind
gendist.BCD <- vegdist(geno.3.genind, method = 'bray')
gendist.BCD.mtrx <- as.matrix(gendist.BCD)

### Prepare for GDM analysis:
GDM.gen <- gendist.BCD.mtrx
dim(GDM.gen)
GDM.gen[1:4,1:4]
## addin site ID to column:
`Analyses Code` <- rownames(info.3)
GDM.gen.site <- cbind(site=`Analyses Code`, GDM.gen)

#### Explanatory variables:
###########################
### Automatic variable selection:
envTab.aut <- cbind(site=`Analyses Code`, Num_vars[,c(3,14,20,21)])
### Manual variable selection:
envTab.man <- cbind(site=`Analyses Code`, Num_vars[,c(3,14,2,15,18,20,21)])

#### A.1. AUTOMATIC SELECTION - WITH Geography:
###############################################
### Create formatSitePair for GDM:
gdmTab.dis.aut <- formatsitepair(bioData=GDM.gen.site, 
                                 bioFormat=3, #diss matrix 
                                 XColumn="Lon", 
                                 YColumn="Lat", 
                                 predData=envTab.aut, 
                                 siteColumn="site")
summary(gdmTab.dis.aut)
### Fit GDM:
gdm.aut <- gdm(data=gdmTab.dis.aut, geo=T)
summary(gdm.aut)
### Plot GDM:
length(gdm.aut$predictors) # get ideal of number of panels
plot(gdm.aut, plot.layout=c(3,3))
### Permutation of GDM:
set.seed(42)
plotUncertainty(gdmTab.dis.aut,sampleSites = 0.7, geo = T,plot.layout=c(1,1), bsIters = 1000, cores = 3 )

#### A.1. AUTOMATIC SELECTION - WITHOUT Geography:
###############################################
### Create formatSitePair for GDM:
gdmTab.dis.aut <- formatsitepair(bioData=GDM.gen.site, 
                                 bioFormat=3, #diss matrix 
                                 XColumn="Lon", 
                                 YColumn="Lat", 
                                 predData=envTab.aut, 
                                 siteColumn="site")
summary(gdmTab.dis.aut)
### Fit GDM:
gdm.aut <- gdm(data=gdmTab.dis.aut, geo=F)
summary(gdm.aut)
### Plot GDM:
length(gdm.aut$predictors) # get ideal of number of panels
plot(gdm.aut, plot.layout=c(3,3))
### Permutation of GDM:
set.seed(42)
plotUncertainty(gdmTab.dis.aut,sampleSites = 0.7, geo = F,plot.layout=c(1,1), bsIters = 1000, cores = 3 )


#### B.1. MANUAL SELECTION - WITH Geography:
############################################
### Create formatSitePair for GDM:
gdmTab.dis.man <- formatsitepair(bioData=GDM.gen.site, 
                                 bioFormat=3, #diss matrix 
                                 XColumn="Lon", 
                                 YColumn="Lat", 
                                 predData=envTab.man, 
                                 siteColumn="site")
summary(gdmTab.dis.man)

### Fit GDM:
gdm.man <- gdm(data=gdmTab.dis.man, geo=T)
summary(gdm.man)
### Plot GDM:
length(gdm.man$predictors) # get ideal of number of panels
plot(gdm.man, plot.layout=c(3,3))
plot(gdm.man$predicted, gdm.man$observed, ylim=c(0,0.5))
### Permutation of GDM:
set.seed(42)
plotUncertainty(gdmTab.dis.man,sampleSites = 0.7, geo = T,plot.layout=c(1,1), bsIters = 1000, cores = 3 )

#### B.2. MANUAL SELECTION - WITHOUT Geography:
############################################
### Create formatSitePair for GDM:
gdmTab.dis.man <- formatsitepair(bioData=GDM.gen.site, 
                                 bioFormat=3, #diss matrix 
                                 XColumn="Lon", 
                                 YColumn="Lat", 
                                 predData=envTab.man, 
                                 siteColumn="site")
summary(gdmTab.dis.man)

### Fit GDM:
gdm.man <- gdm(data=gdmTab.dis.man, geo=F)
summary(gdm.man)
### Plot GDM:
length(gdm.man$predictors) # get ideal of number of panels
plot(gdm.man, plot.layout=c(3,3))
plot(gdm.man$predicted, gdm.man$observed, ylim=c(0,0.5))
### Permutation of GDM:
set.seed(42)
plotUncertainty(gdmTab.dis.man,sampleSites = 0.7, geo = F,plot.layout=c(1,1), bsIters = 1000, cores = 3 )

### EOF ###---------------------------------------------------------------------

