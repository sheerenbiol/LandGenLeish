#########################################################################
##### Diversity and dissemination of viruses in pathogenic protozoa #####
#-----------------------------------------------------------------------#
#####                  Environmental Niche Modeling                 #####
#########################################################################
### Summary:
## - 0. Set up
## - 1. Importing and filtering data
## - 2. Calculating distances
## - 3. Variable Selection
## - 4. Environmental Niche modeling - Present day data
## - 5. Projection on LGM data
## - 6. Projection on LIG data
## - 7. Additional visualization

### Java Parameters:
options(java.parameters = "-Xmx1g" )

### Set working directory:
setwd('PATH/TO/LandGenLeish/')

### Libraries:
##############
library(data.table); library(readxl);library(sf); library(tibble); library(sfheaders)
library(geodist); library(adegenet); library(vegan); library(ggplot2); library(raster)
library(rJava);library(rworldmap)


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

### Run Maxent in R:
Max.mod <- function(stack,pts, nbgs, rm, replis, outpt){
# stack = raster stack
# pts = Coordinate data
# nbgs = number of background points / pseudo-absence points - recommend 10,000
# rm = set betamultiplier 
# replis = number of replicates
# outpt = Path for output files
dismo::maxent(x=stack, p=pts, a=nbgs, path=outpt,removeDuplicates=T,
              args=c("responsecurves",'jackknife', paste0('replicates=',replis), 'replicatetype=crossvalidate', 
                     'writebackgroundpredictions=true',
                     'linear','quadratic','product','threshold','hinge', paste0('betamultiplier=',rm)))
}

AUC <- function(p,ind) {
  pres <- p[ind]
  combined <- c(pres, bb)
  label <- c(rep(1,length(pres)),rep(0,length(bb)))
  predic <- prediction(combined, label)
  return(performance(predic, "auc")@y.values[[1]])
}

Maxent_eval <- function(P, B, M){
  # ***
  # P = PATH to presence sample predictions --> .._samplePredictions.csv
  # B = PATH to background sample predictions --> .._backgroundPredictions.csv
  # M = MODEL
  # ***
  ## Required library:
  library(ROCR); library(boot)
  ## Read in data:
  pres <- read.csv(P)
  backg <- read.csv(B)
  
  ## Set Presence points + split train from test points
  pp <- pres$Cloglog.prediction
  pp_train <- pp[pres$Test.or.train=='train']
  pp_test <- pp[pres$Test.or.train=='test']
  ## Set background points
  bb <- backg$Cloglog
  ## Combine test points + background
  combi <- c(pp_test, bb)
  
  ## ROC analysis: Generate ROC curve + calculate the AUC
  label <- c(rep(1,length(pp_test)),rep(0,length(bb)))
  pred <- prediction(combi, label) 
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, colorize=TRUE)
  print("AUC:")
  print(performance(pred, "auc")@y.values[[1]]) #AUC
  
  ## Bootstrap of AUC:
  print(boot(pp_test, AUC, 100))
  
  ##dismo evaluate:
  e <- evaluate(pp_test, bb, M)
  plot(e, 'ROC')
  
  ##Return Performance
  return(performance(pred, "auc")@y.values[[1]])
  
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
#### Adding random Jitter to duplicate coordinates: -- To account for uncertainty of original loactions
###################################################
Coords <- info.3[,c(10,9)]; Coords <- as.data.frame(Coords)
Coords$Longitude <- as.numeric(unlist(Coords$Longitude)); Coords$Latitude <- as.numeric(unlist(Coords$Latitude))
rownames(Coords) <- info.3$`Analyses Code`
Coords.sf <- st_as_sf(Coords, coords = c('Longitude', 'Latitude'), crs = 4326); rownames(Coords.sf) <- rownames(Coords)
## Adding a jitter 
set.seed(42)
duplicated(Coords.sf$geometry)
tmp.unique <- Coords.sf[!duplicated(Coords.sf$geometry),]
tmp.dupl <- Coords.sf[duplicated(Coords.sf$geometry),]
tmp.dupl.jit <- tmp.dupl; tmp.dupl.jit <- st_jitter(st_geometry(tmp.dupl.jit), factor = 0.01)
tmp.dupl.jit <- do.call(rbind, st_geometry(tmp.dupl.jit)) %>% 
  as_tibble() %>% setNames(c("lon","lat")); tmp.dupl.jit <- st_as_sf(tmp.dupl.jit, coords = c('lon', 'lat'), crs = 4326)
#Back to sf object + adding rownames
tmp.dupl.jit <- st_as_sf(tmp.dupl.jit, coords = 'geometry', crs = 4326)
class(tmp.dupl.jit) <- class(tmp.dupl.jit)[c(1,4)] #removing tibble class - as this does not allow for adding rownames!!
rownames(tmp.dupl.jit) <- rownames(Coords[c(3,6,7,10,11,14:30,35,37,38),]) 
rownames(tmp.unique) <- rownames(Coords[-c(3,6,7,10,11,14:30,35,37,38),]) 
#Combine everything back together:
Coords.jitter.sf <- rbind(tmp.unique, tmp.dupl.jit); rownames(Coords.jitter.sf) <- rownames(Coords[c(1,2,4,5,8,9,12,13,31,32,33,34,36,3,6,7,10,11,14:30,35,37,38),])
Coords.jitter <-  sf_to_df(Coords.jitter.sf); rownames(Coords.jitter) <- rownames(Coords[c(1,2,4,5,8,9,12,13,31,32,33,34,36,3,6,7,10,11,14:30,35,37,38),])
Coords.jitter <- Coords.jitter[,-c(1:2)]; colnames(Coords.jitter) <- c('Longitude','Latitude')
#### Adding random Jitter to duplicate coordinates: -- To account for uncertainty of original loactions
###################################################----------------------------

### Creating geographical distance matrix | method of distance: https://www.r-bloggers.com/2010/11/great-circle-distance-calculations-in-r/
geograph.mtrx <- geodist(Coords.jitter, measure = 'haversine')
rownames(geograph.mtrx) <- rownames(Coords.jitter); colnames(geograph.mtrx) <- rownames(Coords.jitter)
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
w <- getData('worldclim', var='bio', res=0.5, lon=-69, lat=-12) ## Will be saved as folder with name 'wc0.5'
### Jittered Version:
bioclim_vars.jit <- rownames(Coords.jitter)
bioclim_vars.jit<- as.data.frame(bioclim_vars.jit)
bioclim_vars.jit$Sample <- bioclim_vars.jit$bioclim_vars.jit
for (i in w@layers) {
  newcol <- extract(x= i, y=cbind(as.numeric(Coords.jitter$Longitude),as.numeric(Coords.jitter$Latitude)))
  bioclim_vars.jit[,ncol(bioclim_vars.jit) + 1] <- newcol
  colnames(bioclim_vars.jit)[ncol(bioclim_vars.jit)] <- paste0(i)
}
bioclim_vars.jit <- bioclim_vars.jit[,-1]
colnames(bioclim_vars.jit) <- c('Sample', 'bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')

## Creating a distance matrix for each bioclim variable
for (bio in colnames(bioclim_vars.jit)[2:20]) {
  print(bio)
  tmp <- bioclim_vars.jit[,bio]
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- bioclim_vars.jit$Sample
  tmp.dist <- dist(tmp, method = 'euclidean')
  assign(paste0(bio,'.jit.dist'), tmp.dist)
  assign(paste0(bio,'.jit.mtrx'), as.matrix(tmp.dist))
}

##### D. Adding raster layer objects for SDM: - Same raster data, but different extent - preferred for visualization - 
#############################################
w.per <- geodata::worldclim_country('Peru', var='bio', path = './', res=0.5)   ## Will be saved as folder with name 'wc2.1_country'
w.bol <- geodata::worldclim_country('Bolivia', var='bio', path = './', res=0.5)
bio1 <- merge(w.per[[1]], w.bol[[1]]); plot(bio1);bio1<- bio1 %>% raster()
bio2 <- merge(w.per[[2]], w.bol[[2]]); plot(bio2);bio2<- bio2 %>% raster()
bio3 <- merge(w.per[[3]], w.bol[[3]]); plot(bio3);bio3<- bio3 %>% raster()
bio4 <- merge(w.per[[4]], w.bol[[4]]); plot(bio4);bio4<- bio4 %>% raster()
bio5 <- merge(w.per[[5]], w.bol[[5]]); plot(bio5);bio5<- bio5 %>% raster()
bio6 <- merge(w.per[[6]], w.bol[[6]]); plot(bio6);bio6<- bio6 %>% raster()
bio7 <- merge(w.per[[7]], w.bol[[7]]); plot(bio7);bio7<- bio7 %>% raster()
bio8 <- merge(w.per[[8]], w.bol[[8]]); plot(bio8);bio8<- bio8 %>% raster()
bio9 <- merge(w.per[[9]], w.bol[[9]]); plot(bio9);bio9<- bio9 %>% raster()
bio10 <- merge(w.per[[10]], w.bol[[10]]); plot(bio10);bio10<- bio10 %>% raster()
bio11 <- merge(w.per[[11]], w.bol[[11]]); plot(bio11);bio11<- bio11 %>% raster()
bio12 <- merge(w.per[[12]], w.bol[[12]]); plot(bio12);bio12<- bio12 %>% raster()
bio13 <- merge(w.per[[13]], w.bol[[13]]); plot(bio13);bio13<- bio13 %>% raster()
bio14 <- merge(w.per[[14]], w.bol[[14]]); plot(bio14);bio14<- bio14 %>% raster()
bio15 <- merge(w.per[[15]], w.bol[[15]]); plot(bio15);bio15<- bio15 %>% raster()
bio16 <- merge(w.per[[16]], w.bol[[16]]); plot(bio16);bio16<- bio16 %>% raster()
bio17 <- merge(w.per[[17]], w.bol[[17]]); plot(bio17);bio17<- bio17 %>% raster()
bio18 <- merge(w.per[[18]], w.bol[[18]]); plot(bio18);bio18<- bio18 %>% raster()
bio19 <- merge(w.per[[19]], w.bol[[19]]); plot(bio19);bio19<- bio19 %>% raster()

##################################----------------------------------------------
##### 3. VARIABLE SELECTION  #####
##################################
# We'll start from the variables selected for the RDA & GDM models and then re-check multicolinearity,...
# --> These are: bio3, bio14, bio2, bio15, bio18

## Check correlation of these variables:
psych::pairs.panels(bioclim_vars.jit[,c(4,15,3,16,19)])
## Check VIF:
usdm::vif(bioclim_vars.jit[,c(4,15,3,16,19)])
## --> all VIFs are </=10 (is OK) but we'll remove bio15 as it has a correlation coeff > |0.7| with bio3
usdm::vif(bioclim_vars.jit[,c(4,15,3,19)])
## --> all VIFs are now </= to 5

## Proceeding with bio3, bio14, bio2, bio18 ##


##############################################################------------------
##### 4. Environmental Niche Modeling - PRESENT-DAY data  #####
##############################################################
### Stacking the environmental layers:
######################################
All <- c(bio3, bio14, bio2, bio18); All <- raster::stack(All)

## Jittered coordinates (duplicates are jittered)
Coords.jitter

set.seed(42)
Coords.jitter_split <- rsample::initial_split(Coords.jitter, prop = 0.8)
Coords.jitter_train <- rsample::training(Coords.jitter_split)
Coords.jitter_test <- rsample::testing(Coords.jitter_split)


### Set Background point (Pseudo-absences):
###########################################
set.seed(42)
bg <- sampleRandom(x=All[[1]],size = 10000, na.rm=T, sp=T)
plot(All[[1]]); plot(bg, add=T); plot(SpatialPoints(Coords.jitter) , col='red', add=T)

### MAXENT Modelling:
#####################
mod.coords.jit.rm1 <- Max.mod(stack = All, pts=Coords.jitter_train, nbgs = bg, rm=1, replis = 10, outpt = './EcolNicheMod/rm1')
mod.coords.jit.rm1.5 <- Max.mod(stack = All, pts=Coords.jitter_train, nbgs = bg, rm=1.5, replis = 10, outpt = './EcolNicheMod/rm1.5')
mod.coords.jit.rm2 <- Max.mod(stack = All, pts=Coords.jitter_train, nbgs = bg, rm=2, replis = 10, outpt = './EcolNicheMod/rm2')

### Predict/ project to our raster layer:
#########################################
mod.coords.jit.rm1.pred <- predict(mod.coords.jit.rm1, All, args="threads=4")
mod.coords.jit.rm1.5.pred <- predict(mod.coords.jit.rm1.5, All, args="threads=4")
mod.coords.jit.rm2.pred <- predict(mod.coords.jit.rm2, All, args="threads=4")

### Plot replicates AND Average prediction:
##########################################
mod.rm1.summary <- read.csv('./EcolNicheMod/rm1/maxentResults.csv')
plot(mod.coords.jit.rm1.pred); mod.coords.jit.rm1.pred.mean <- mean(mod.coords.jit.rm1.pred[[1]], mod.coords.jit.rm1.pred[[2]], mod.coords.jit.rm1.pred[[3]],
                                                                    mod.coords.jit.rm1.pred[[4]], mod.coords.jit.rm1.pred[[5]],mod.coords.jit.rm1.pred[[6]],mod.coords.jit.rm1.pred[[7]],
                                                                    mod.coords.jit.rm1.pred[[8]], mod.coords.jit.rm1.pred[[9]], mod.coords.jit.rm1.pred[[10]]); plot(mod.coords.jit.rm1.pred.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

mod.rm1.5.summary <- read.csv('./maxentResults.csv')
plot(mod.coords.jit.rm1.5.pred); mod.coords.jit.rm1.5.pred.mean <- mean(mod.coords.jit.rm1.5.pred[[1]], mod.coords.jit.rm1.5.pred[[2]], mod.coords.jit.rm1.5.pred[[3]],
                                                                        mod.coords.jit.rm1.5.pred[[4]], mod.coords.jit.rm1.5.pred[[5]],mod.coords.jit.rm1.5.pred[[6]],mod.coords.jit.rm1.5.pred[[7]],
                                                                        mod.coords.jit.rm1.5.pred[[8]], mod.coords.jit.rm1.5.pred[[9]], mod.coords.jit.rm1.5.pred[[10]]); plot(mod.coords.jit.rm1.5.pred.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

mod.rm2.summary <- read.csv('./maxentResults.csv')
plot(mod.coords.jit.rm2.pred); mod.coords.jit.rm2.pred.mean <- mean(mod.coords.jit.rm2.pred[[1]], mod.coords.jit.rm2.pred[[2]], mod.coords.jit.rm2.pred[[3]],
                                                                    mod.coords.jit.rm2.pred[[4]], mod.coords.jit.rm2.pred[[5]],mod.coords.jit.rm2.pred[[6]],mod.coords.jit.rm2.pred[[7]],
                                                                    mod.coords.jit.rm2.pred[[8]], mod.coords.jit.rm2.pred[[9]], mod.coords.jit.rm2.pred[[10]]); plot(mod.coords.jit.rm2.pred.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

### Model Evaluation:
#####################
#AUCS for the 10 replications:
aucs.rm1 <- vector()
for (i in 0:9){tmp <- Maxent_eval(P=paste0('./EcolNicheMod/rm1/species_',i,'_samplePredictions.csv'), B=paste0('./EcolNicheMod/rm1/species_',i,'_backgroundPredictions.csv'), M=mod.coords.jit.rm1)
aucs.rm1 <- c(aucs.rm1, tmp)
print(tmp)}
aucs.rm1; mean(aucs.rm1); sd(aucs.rm1)

aucs.rm1.5 <- vector()
for (i in 0:9){tmp <- Maxent_eval(P=paste0('./2_output/2_Maxent_out_coords_jitter_rm1.5/species_',i,'_samplePredictions.csv'), paste0('./2_output/2_Maxent_out_coords_jitter_rm1.5/species_',i,'_backgroundPredictions.csv'))
aucs.rm1.5 <- c(aucs.rm1.5, tmp)}
aucs.rm1.5; mean(aucs.rm1.5); sd(aucs.rm1.5)

aucs.rm2 <- vector()
for (i in 0:9){tmp <- Maxent_eval(P=paste0('./2_output/3_Maxent_out_coords_jitter_rm2/species_',i,'_samplePredictions.csv'), paste0('./2_output/3_Maxent_out_coords_jitter_rm2/species_',i,'_backgroundPredictions.csv'))
aucs.rm2 <- c(aucs.rm2, tmp)}
aucs.rm2; mean(aucs.rm2); sd(aucs.rm2)

### Nice SDM plots:
###################
### Countour layers of Peru and Bolivia:
map_world <- getMap(); projection <- "+proj=longlat +datum=WGS84"
map_SA <- map_world[which(map_world@data$continent == 'South America'),]
per0 <- map_world[which(map_world@data$ADMIN %in% c("Peru")),]
bol0 <- map_world[which(map_world@data$ADMIN %in% c("Bolivia")),]
per0@proj4string <- projection; bol0@proj4string <- projection
crs(per0) <- "+proj=longlat +datum=WGS84"; crs(bol0) <- "+proj=longlat +datum=WGS84"

## Masking SDM rasters
sdm1.rm1 <- raster::mask(mod.coords.jit.rm1.pred.mean, per0) # for Peru
plot(sdm1.rm1)
sdm2.rm1 <- raster::mask(mod.coords.jit.rm1.pred.mean, bol0) # for Bolivia
plot(sdm2.rm1)

sdm1.rm1.5 <- raster::mask(mod.coords.jit.rm1.5.pred.mean, per0) # for Peru
plot(sdm1.rm1.5)
sdm2.rm1.5 <- raster::mask(mod.coords.jit.rm1.5.pred.mean, bol0) # for Bolivia
plot(sdm2.rm1.5)

sdm1.rm2 <- raster::mask(mod.coords.jit.rm2.pred.mean, per0) # for Peru
plot(sdm1.rm2)
sdm2.rm2 <- raster::mask(mod.coords.jit.rm2.pred.mean, bol0) # for Bolivia
plot(sdm2.rm2)

SDM.mod.coords.jit.rm1 <- merge(sdm1.rm1, sdm2.rm1)
SDM.mod.coords.jit.rm1.5 <- merge(sdm1.rm1.5, sdm2.rm1.5)
SDM.mod.coords.jit.rm2 <- merge(sdm1.rm2, sdm2.rm2)

par(mfrow=c(1,3))
SDM.mod.coords.jit.rm1; plot(SDM.mod.coords.jit.rm1, xaxt='n',yaxt='n', ann=FALSE, ylim=c(-24,0), main= 'coords jit cv10 rm1'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
SDM.mod.coords.jit.rm1.5; plot(SDM.mod.coords.jit.rm1.5, xaxt='n',yaxt='n', ann=FALSE, ylim=c(-24,0), main= 'coords jit cv10 rm1.5'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
SDM.mod.coords.jit.rm2; plot(SDM.mod.coords.jit.rm2, ylim=c(-24,0), xaxt='n',yaxt='n', ann=FALSE, main= 'coords jit cv10 rm2'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
par(mfrow=c(1,1))


############################################################--------------------
##### 5. Environmental Nich Modeling - LGM Projection  #####
############################################################
## Bioclim variables: ---> use the same variables: bio3, bio14, bio2, bio18 in models + different regularization values

### Past - LGM - Climate data:
##############################
#download.file('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio02_-69_V1.0.tif', destfile = './Input_files/LGM_rasters/LGM_bio2.tif')
lgm.bio2 <- raster(x = './Input_files/LGM_rasters/LGM_bio2.tif'); plot(lgm.bio2)
r2 <- raster::mask(lgm.bio2, per0);r3 <- raster::mask(lgm.bio2, bol0); lgm.bio2.m <- merge(r2,r3); plot(lgm.bio2.m, xlim=c(-85,-55), ylim=c(-25,0))
#download.file('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio03_-69_V1.0.tif', destfile = './Input_files/LGM_rasters/LGM_bio3.tif')
lgm.bio3 <- raster(x = './Input_files/LGM_rasters/LGM_bio3.tif'); plot(lgm.bio3)
r2 <- raster::mask(lgm.bio3, per0);r3 <- raster::mask(lgm.bio3, bol0); lgm.bio3.m <- merge(r2,r3); plot(lgm.bio3.m, xlim=c(-85,-55), ylim=c(-25,0))
#download.file('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio14_-69_V1.0.tif', destfile = './Input_files/LGM_rasters/LGM_bio14.tif')
lgm.bio14 <- raster(x = './Input_files/LGM_rasters/LGM_bio14.tif'); plot(lgm.bio14)
r2 <- raster::mask(lgm.bio14, per0);r3 <- raster::mask(lgm.bio14, bol0); lgm.bio14.m <- merge(r2,r3); plot(lgm.bio14.m, xlim=c(-85,-55), ylim=c(-25,0))
#download.file('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio18_-69_V1.0.tif', destfile = './Input_files/LGM_rasters/LGM_bio18.tif')
lgm.bio18 <- raster(x = './Input_files/LGM_rasters/LGM_bio18.tif'); plot(lgm.bio18)
r2 <- raster::mask(lgm.bio18, per0);r3 <- raster::mask(lgm.bio18, bol0); lgm.bio18.m <- merge(r2,r3); plot(lgm.bio18.m, xlim=c(-85,-55), ylim=c(-25,0))

## Changing the rasters' extent:
lgm.bio3.crop <- crop(x=lgm.bio3.m, y=w[[1]]); lgm.bio14.crop <- crop(x=lgm.bio14.m, y=w[[1]])
lgm.bio2.crop <- crop(x=lgm.bio2.m, y=w[[1]]); lgm.bio18.crop <- crop(x=lgm.bio18.m, y=w[[1]]) #altitude from gis_maps.R script

## Stack the layers together:
All.LGM <- c(lgm.bio3.crop, lgm.bio14.crop, lgm.bio2.crop, lgm.bio18.crop); All.LGM <- stack(All.LGM)
names(All.LGM) <- names(All) #Necessary for projection with the SDM model!

### Predicting the presence from model using LGM climate data:
##############################################################
LGM.pred.coords.jit.rm1 <- dismo::predict(object = mod.coords.jit.rm1, x = All.LGM, args="threads=4" )
LGM.pred.coords.jit.rm1.5 <- dismo::predict(object = mod.coords.jit.rm1.5, x = All.LGM, args="threads=4" )
LGM.pred.coords.jit.rm2 <- dismo::predict(object = mod.coords.jit.rm2, x = All.LGM, args="threads=4" )

### Plot the LGM projections:
#############################
plot(LGM.pred.coords.jit.rm1); LGM.pred.coords.jit.rm1.mean <- mean(LGM.pred.coords.jit.rm1[[1]], LGM.pred.coords.jit.rm1[[2]], LGM.pred.coords.jit.rm1[[3]],
                                                                    LGM.pred.coords.jit.rm1[[4]], LGM.pred.coords.jit.rm1[[5]],LGM.pred.coords.jit.rm1[[6]],LGM.pred.coords.jit.rm1[[7]],
                                                                    LGM.pred.coords.jit.rm1[[8]], LGM.pred.coords.jit.rm1[[9]], LGM.pred.coords.jit.rm1[[10]]); plot(LGM.pred.coords.jit.rm1.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

plot(LGM.pred.coords.jit.rm1.5); LGM.pred.coords.jit.rm1.5.mean <- mean(LGM.pred.coords.jit.rm1.5[[1]], LGM.pred.coords.jit.rm1.5[[2]], LGM.pred.coords.jit.rm1.5[[3]],
                                                                        LGM.pred.coords.jit.rm1.5[[4]], LGM.pred.coords.jit.rm1.5[[5]],LGM.pred.coords.jit.rm1.5[[6]],LGM.pred.coords.jit.rm1.5[[7]],
                                                                        LGM.pred.coords.jit.rm1.5[[8]], LGM.pred.coords.jit.rm1.5[[9]], LGM.pred.coords.jit.rm1.5[[10]]); plot(LGM.pred.coords.jit.rm1.5.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

plot(LGM.pred.coords.jit.rm2); LGM.pred.coords.jit.rm2.mean <- mean(LGM.pred.coords.jit.rm2[[1]], LGM.pred.coords.jit.rm2[[2]], LGM.pred.coords.jit.rm2[[3]],
                                                                    LGM.pred.coords.jit.rm2[[4]], LGM.pred.coords.jit.rm2[[5]],LGM.pred.coords.jit.rm2[[6]],LGM.pred.coords.jit.rm2[[7]],
                                                                    LGM.pred.coords.jit.rm2[[8]],LGM.pred.coords.jit.rm2[[9]], LGM.pred.coords.jit.rm2[[10]]); plot(LGM.pred.coords.jit.rm2.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

### Nice SDM plots:
###################
## Masking SDM rasters
sdm1.rm1 <- raster::mask(LGM.pred.coords.jit.rm1.mean, per0) # for Peru
plot(sdm1.rm1)
sdm2.rm1 <- raster::mask(LGM.pred.coords.jit.rm1.mean, bol0) # for Bolivia
plot(sdm2.rm1)

sdm1.rm1.5 <- raster::mask(LGM.pred.coords.jit.rm1.5.mean, per0) # for Peru
plot(sdm1.rm1.5)
sdm2.rm1.5 <- raster::mask(LGM.pred.coords.jit.rm1.5.mean, bol0) # for Bolivia
plot(sdm2.rm1.5)

sdm1.rm2 <- raster::mask(LGM.pred.coords.jit.rm2.mean, per0) # for Peru
plot(sdm1.rm2)
sdm2.rm2 <- raster::mask(LGM.pred.coords.jit.rm2.mean, bol0) # for Bolivia
plot(sdm2.rm2)

LGM.mod.coords.jit.rm1 <- merge(sdm1.rm1, sdm2.rm1)
LGM.mod.coords.jit.rm1.5 <- merge(sdm1.rm1.5, sdm2.rm1.5)
LGM.mod.coords.jit.rm2 <- merge(sdm1.rm2, sdm2.rm2)

par(mfrow=c(1,3))
par(mar=c(4,2,4,2), oma=c(1,3,2,4))
LGM.mod.coords.jit.rm1; plot(LGM.mod.coords.jit.rm1, xlim=c(-81,-58), xaxt='n',yaxt='n', ann=FALSE, ylim=c(-24,0), main= 'coords jit cv10 rm1'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
LGM.mod.coords.jit.rm1.5; plot(LGM.mod.coords.jit.rm1.5, xlim=c(-81,-58), xaxt='n',yaxt='n', ann=FALSE, ylim=c(-24,0), main= 'coords jit cv10 rm1.5'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
LGM.mod.coords.jit.rm2; plot(LGM.mod.coords.jit.rm2, xlim=c(-81,-58), ylim=c(-24,0), xaxt='n',yaxt='n', ann=FALSE, main= 'coords jit cv10 rm2'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))




############################################################--------------------
##### 6. Environmental Nich Modeling - LIG Projection  #####
############################################################
## Bioclim variables: ---> use the same variables: bio3, bio14, bio2, bio18 in models + different regularizations values

### Past - LIG - Climate data:
##############################
## Download from: http://sdmtoolbox.org/paleoclim.org/data/LIG/LIG_v1_2_5m.zip
## Save tiffs in: './Input_files/LIG_rasters/'
lig.bio2 <- raster(x = './Input_files/LIG_rasters/LIG_bio_2.tif'); plot(lig.bio2)
r2 <- raster::mask(lig.bio2, per0);r3 <- raster::mask(lig.bio2, bol0); lig.bio2.m <- merge(r2,r3); plot(lig.bio2.m, xlim=c(-85,-55), ylim=c(-25,0))
lig.bio3 <- raster(x = './Input_files/LIG_rasters/LIG_bio_3.tif'); plot(lig.bio3)
r2 <- raster::mask(lig.bio3, per0);r3 <- raster::mask(lig.bio3, bol0); lig.bio3.m <- merge(r2,r3); plot(lig.bio3.m, xlim=c(-85,-55), ylim=c(-25,0))
lig.bio14 <- raster(x = './Input_files/LIG_rasters/LIG_bio_14.tif'); plot(lig.bio14)
r2 <- raster::mask(lig.bio14, per0);r3 <- raster::mask(lig.bio14, bol0); lig.bio14.m <- merge(r2,r3); plot(lig.bio14.m, xlim=c(-85,-55), ylim=c(-25,0))
lig.bio18 <- raster(x = './Input_files/LIG_rasters/LIG_bio_18.tif'); plot(lig.bio18)
r2 <- raster::mask(lig.bio18, per0);r3 <- raster::mask(lig.bio18, bol0); lig.bio18.m <- merge(r2,r3); plot(lig.bio18.m, xlim=c(-85,-55), ylim=c(-25,0))

## Changing the raster's extent:
lig.bio3.crop <- crop(x=lig.bio3.m, y=w[[1]]); lig.bio14.crop <- crop(x=lig.bio14.m, y=w[[1]])
lig.bio2.crop <- crop(x=lig.bio2.m, y=w[[1]]); lig.bio18.crop <- crop(x=lig.bio18.m, y=w[[1]]) #altitude from gis_maps.R script

## Stack the layers together:
All.LIG <- c(lig.bio3.crop, lig.bio14.crop, lig.bio2.crop, lig.bio18.crop); All.LIG <- stack(All.LIG)
names(All.LIG) <- names(All) #Necessary for projection with the SDM model!

### Predicting the presence from model using LIG climate data:
##############################################################
LIG.pred.coords.jit.rm1 <- dismo::predict(object = mod.coords.jit.rm1, x = All.LIG, args="threads=4" )
LIG.pred.coords.jit.rm1.5 <- dismo::predict(object = mod.coords.jit.rm1.5, x = All.LIG, args="threads=4" )
LIG.pred.coords.jit.rm2 <- dismo::predict(object = mod.coords.jit.rm2, x = All.LIG, args="threads=4" )

### Plot the LIG projections:
#############################
plot(LIG.pred.coords.jit.rm1); LIG.pred.coords.jit.rm1.mean <- mean(LIG.pred.coords.jit.rm1[[1]], LIG.pred.coords.jit.rm1[[2]], LIG.pred.coords.jit.rm1[[3]],
                                                                    LIG.pred.coords.jit.rm1[[4]], LIG.pred.coords.jit.rm1[[5]],LIG.pred.coords.jit.rm1[[6]],LIG.pred.coords.jit.rm1[[7]],
                                                                    LIG.pred.coords.jit.rm1[[8]], LIG.pred.coords.jit.rm1[[9]], LIG.pred.coords.jit.rm1[[10]]); plot(LIG.pred.coords.jit.rm1.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

plot(LIG.pred.coords.jit.rm1.5); LIG.pred.coords.jit.rm1.5.mean <- mean(LIG.pred.coords.jit.rm1.5[[1]], LIG.pred.coords.jit.rm1.5[[2]], LIG.pred.coords.jit.rm1.5[[3]],
                                                                        LIG.pred.coords.jit.rm1.5[[4]], LIG.pred.coords.jit.rm1.5[[5]],LIG.pred.coords.jit.rm1.5[[6]],LIG.pred.coords.jit.rm1.5[[7]],
                                                                        LIG.pred.coords.jit.rm1.5[[8]], LIG.pred.coords.jit.rm1.5[[9]], LIG.pred.coords.jit.rm1.5[[10]]); plot(LIG.pred.coords.jit.rm1.5.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

plot(LIG.pred.coords.jit.rm2); LIG.pred.coords.jit.rm2.mean <- mean(LIG.pred.coords.jit.rm2[[1]], LIG.pred.coords.jit.rm2[[2]], LIG.pred.coords.jit.rm2[[3]],
                                                                    LIG.pred.coords.jit.rm2[[4]], LIG.pred.coords.jit.rm2[[5]],LIG.pred.coords.jit.rm2[[6]],LIG.pred.coords.jit.rm2[[7]],
                                                                    LIG.pred.coords.jit.rm2[[8]],LIG.pred.coords.jit.rm2[[9]], LIG.pred.coords.jit.rm2[[10]]); plot(LIG.pred.coords.jit.rm2.mean); plot(SpatialPoints(Coords.jitter), add=T, pch=21)

### Nice SDM plots:
###################
## Masking SDM rasters
sdm1 <- raster::mask(LIG.pred.coords.jit.rm1.5.mean, per0) # for Peru
plot(sdm1)
sdm2 <- raster::mask(LIG.pred.coords.jit.rm1.5.mean, bol0) # for Bolivia
plot(sdm2)

LIG.mod.coords.jit.rm1 <- merge(sdm1, sdm2)
LIG.mod.coords.jit.rm1.5 <- merge(sdm1, sdm2)
LIG.mod.coords.jit.rm2 <- merge(sdm1, sdm2)

par(mfrow=c(1,3))
par(mar=c(4,2,4,2), oma=c(1,3,2,4))
LIG.mod.coords.jit.rm1; plot(LIG.mod.coords.jit.rm1, xlim=c(-81,-58), xaxt='n',yaxt='n', ann=FALSE, ylim=c(-24,0), main= 'coords jit cv10 rm1'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
LIG.mod.coords.jit.rm1.5; plot(LIG.mod.coords.jit.rm1.5, xlim=c(-81,-58), xaxt='n',yaxt='n', ann=FALSE, ylim=c(-24,0), main= 'coords jit cv10 rm1.5'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
LIG.mod.coords.jit.rm2; plot(LIG.mod.coords.jit.rm2, xlim=c(-81,-58), ylim=c(-24,0), xaxt='n',yaxt='n', ann=FALSE, main= 'coords jit cv10 rm2'); scalebar(500)
plot(per0, add=T); plot(bol0, add=T)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))


#################################################################---------------
##### 7. Environmental Niche Modeling - PLOTTING EVERYTHING  #####
#################################################################
### Main figure panel: --> Present day SDM RM = 1.5
######################
## Set the model to Spatial Pixel Data Frame + then to Data Frame:
SDM.mod.coords.jit.rm1.5.spdf <- as(SDM.mod.coords.jit.rm1.5, "SpatialPixelsDataFrame")
SDM.mod.coords.jit.rm1.5.df <- as.data.frame(SDM.mod.coords.jit.rm1.5.spdf); colnames(SDM.mod.coords.jit.rm1.5.df) <- c("value","x","y")

## ggplot:
Coords.unique <- sf_to_df(tmp.unique); rownames(Coords.unique) <- rownames(tmp.unique); Coords.unique <- Coords.unique[,-c(1:2)]
colnames(Coords.unique) <- c('Longitude','Latitude')
info.3.unique <- info.3[which(rownames(info.3) %in% rownames(Coords.unique)),]
info.3.unique$Loc <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)


par(bty='n')
plot(SDM.mod.coords.jit.rm1.5 , zlim=c(0,1),main= 'SDM Present cv10 rm1.5 (bio3, bio14, bio2, bio18)',bty='n', box=F)
plot(per1, add=T); plot(bol1, add=T) #--> Departments
plot(per0, lwd=2.5, add=T); plot(bol0,lwd=2.5, add=T) #--> Country
## Add samples:
popcols<- unique(data.frame(as.character(info.3.unique$Pop2),as.character(info.3.unique$PopCol2),stringsAsFactors = F))
info.3.unique <- info.3.unique[order(info.3.unique$Loc),]  
for (i in 1:length(unique(info.3.unique$Loc))) {
  tmp <- info.3.unique[info.3.unique$Loc == i,]
  tmp2 <- table(as.character(tmp[,25]))
  add.pie(x = as.numeric(tmp[1,16]), y=as.numeric(tmp[1,15]), z = as.numeric(tmp2),col = unique(tmp$PopCol2), radius = 0.4, border = F, labels = NA)
  points(as.numeric(tmp[,16]),as.numeric(tmp[,15]), pch = 21, cex = 2)
  text(as.numeric(tmp[1,16]),as.numeric(tmp[1,15]),labels=nrow(tmp), col = 'white', cex = .7)
}
legend('topright', legend = c('INP','HUP','PAU'), fill = popcols$as.character.info.3.unique.PopCol2., cex = 1.5, bty = 'n')
scalebar(d = 500) # in kilometers



### Suppl Figure: --> All SDMs 
#################

## Set All models to spatial pixel dataframes and then to dataframes:
SDM.mod.coords.jit.rm1.spdf <- as(SDM.mod.coords.jit.rm1, "SpatialPixelsDataFrame")
SDM.mod.coords.jit.rm1.df <- as.data.frame(SDM.mod.coords.jit.rm1.spdf); colnames(SDM.mod.coords.jit.rm1.df) <- c("value","x","y")
SDM.mod.coords.jit.rm1.5.spdf <- as(SDM.mod.coords.jit.rm1.5, "SpatialPixelsDataFrame")
SDM.mod.coords.jit.rm1.5.df <- as.data.frame(SDM.mod.coords.jit.rm1.5.spdf); colnames(SDM.mod.coords.jit.rm1.5.df) <- c("value","x","y")
SDM.mod.coords.jit.rm2.spdf <- as(SDM.mod.coords.jit.rm2, "SpatialPixelsDataFrame")
SDM.mod.coords.jit.rm2.df <- as.data.frame(SDM.mod.coords.jit.rm2.spdf); colnames(SDM.mod.coords.jit.rm2.df) <- c("value","x","y")
LGM.mod.coords.jit.rm1.spdf <- as(LGM.mod.coords.jit.rm1, "SpatialPixelsDataFrame")
LGM.mod.coords.jit.rm1.df <- as.data.frame(LGM.mod.coords.jit.rm1.spdf); colnames(LGM.mod.coords.jit.rm1.df) <- c("value","x","y")
LGM.mod.coords.jit.rm1.5.spdf <- as(LGM.mod.coords.jit.rm1.5, "SpatialPixelsDataFrame")
LGM.mod.coords.jit.rm1.5.df <- as.data.frame(LGM.mod.coords.jit.rm1.5.spdf); colnames(LGM.mod.coords.jit.rm1.5.df) <- c("value","x","y")
LGM.mod.coords.jit.rm2.spdf <- as(LGM.mod.coords.jit.rm2, "SpatialPixelsDataFrame")
LGM.mod.coords.jit.rm2.df <- as.data.frame(LGM.mod.coords.jit.rm2.spdf); colnames(LGM.mod.coords.jit.rm2.df) <- c("value","x","y")
LIG.mod.coords.jit.rm1.spdf <- as(LIG.mod.coords.jit.rm1, "SpatialPixelsDataFrame")
LIG.mod.coords.jit.rm1.df <- as.data.frame(LIG.mod.coords.jit.rm1.spdf); colnames(LIG.mod.coords.jit.rm1.df) <- c("value","x","y")
LIG.mod.coords.jit.rm1.5.spdf <- as(LIG.mod.coords.jit.rm1.5, "SpatialPixelsDataFrame")
LIG.mod.coords.jit.rm1.5.df <- as.data.frame(LIG.mod.coords.jit.rm1.5.spdf); colnames(LIG.mod.coords.jit.rm1.5.df) <- c("value","x","y")
LIG.mod.coords.jit.rm2.spdf <- as(LIG.mod.coords.jit.rm2, "SpatialPixelsDataFrame")
LIG.mod.coords.jit.rm2.df <- as.data.frame(LIG.mod.coords.jit.rm2.spdf); colnames(LIG.mod.coords.jit.rm2.df) <- c("value","x","y")

## Bring raster.df + Country outline to ggplot:
pl1 <- SDM.ggplot(dat=SDM.mod.coords.jit.rm1.df)
pl2 <- SDM.ggplot(dat=SDM.mod.coords.jit.rm1.5.df)
pl3 <- SDM.ggplot(dat=SDM.mod.coords.jit.rm2.df)
pl4 <- SDM.ggplot(dat=LGM.mod.coords.jit.rm1.df)
pl5 <- SDM.ggplot(dat=LGM.mod.coords.jit.rm1.5.df)
pl6 <- SDM.ggplot(dat=LGM.mod.coords.jit.rm2.df)
pl7 <- SDM.ggplot(dat=LIG.mod.coords.jit.rm1.df)
pl8 <- SDM.ggplot(dat=LIG.mod.coords.jit.rm1.5.df)
pl9 <- SDM.ggplot(dat=LIG.mod.coords.jit.rm2.df)

ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9,
          ncol = 3, nrow = 3)


### EOF ###---------------------------------------------------------------------


