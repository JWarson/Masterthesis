#install.packages('psych')
library(vegan)
library(psych)

setwd('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan')

################
##GENETIC DATA##
################
########
#alcon#
#######
#load arlequin allele frequencies
alcon <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Arlequin/AlconAlleleFreqs.txt', h = F)
alcon <- subset(alcon[,1:(ncol(alcon)-1)]) #last column are all NAs
popnamesA <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Arlequin/Alcon4group.res/PopNamesA.txt', h =F)
popnamesA$V3 <- 1:nrow(popnamesA) #numeric popnames, for RDA later
colnames(alcon) <- c("Locus", "Allele", popnamesA[,2]) #for data handling

#remove '?' allele
alconMAF <- subset(alcon, alcon$Allele != "?") #removes 108 '?' allele (zero anyway)

#remove 'allele' column and 'ker' location (last column = ncol(alconMAF))
alconMAF <- subset(alconMAF, select= -c(Allele, Ker))

#remove major allele
library(data.table)
alconMAF$sumAF <- rowSums(alconMAF[,2:ncol(alconMAF)]) #creates column with sum of the allele frequencies at each row
alconMAF <- setDT(alconMAF)[, .SD[which.min(sumAF)], by = Locus] #selects the row with smalles sumAF

#remove sumAF column
alconMAF <- subset(alconMAF, select = -c(sumAF))
#turn around rows and column, 'locus' is header
AlconAF <- as.data.frame(t(alconMAF))
names(AlconAF)<- AlconAF[1,]
AlconAF <- AlconAF[2:nrow(AlconAF),]
str(AlconAF)
AlconAF <- type.convert(AlconAF, as.is =T) #convert to numerical

#add numeric PopID (using 'popnamesA' file to be able to trace back easy later)
#cbind to add to front
#AlconAF <- cbind(PopID = popnamesA[1:(nrow(popnamesA)-1),3], AlconAF) #nrow - 1 because last, 'ker' location was already removed

#########
#gentian#
#########
gentian <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Arlequin/GentianAlleleFreqs.txt', h = F)
gentian <- subset(gentian[,1:(ncol(gentian)-1)]) #last column NA
popnamesG <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Arlequin/Gentian.res/PopNamesG.txt', h =F)
popnamesG$V3 <- 1:nrow(popnamesG) #numeric popnames
colnames(gentian) <- c("Locus", "Allele", popnamesG[,2]) #colnames for data handling

#remove '?' allele
gentianMAF <- subset(gentian, gentian$Allele != "?") #removes 108 '?' alleles
#remove 'allele' column
gentianMAF <- subset(gentianMAF, select= -c(Allele))

#remove major allele
gentianMAF$sumAF <- rowSums(gentianMAF[,2:ncol(gentianMAF)])
gentianMAF <- setDT(gentianMAF)[, .SD[which.min(sumAF)], by = Locus]

#remove sumAF 
gentianMAF <- subset(gentianMAF, select = -c(sumAF))
#turn around rows and column, 'locus' is header
GentianAF <- as.data.frame(t(gentianMAF))
names(GentianAF)<- GentianAF[1,]
GentianAF <- GentianAF[2:nrow(GentianAF),]
GentianAF <- type.convert(GentianAF, as.is =T)
str(GentianAF)

rm(alcon, alconMAF, gentianMAF, gentian) #remove unneeded data from environment

####################
#ENVIRONMENTAL DATA#
####################
#landscape connectivity, patch size, habitat suitability and altitude are to be kept
env <- read.csv('environment', h = T, sep = '\t')
env <- subset(env[2:28,]) #remove 'ker' location
str(env)

#remove locations for which there is no genetic data
setdiff(env$PopID,rownames(AlconAF))    #alcon: SW4, SE2
setdiff(env$PopID, rownames(GentianAF)) #gentian: NW18

env.alc <- env[!(env$PopID=="SE2") & !(env$PopID=="SW4"),]
env.gen <- env[!(env$PopID=="NW18"),]

#check correlations
pairs.panels(env[, 7:13], scale = T) #total and regional diversity correlated 
                                    #(will be tested seperatly), 
#patchsize also with plants and geosuitability (we choose to remove plants 
#because patchsize is of particular intrest -> fragmentation and plants contains NAs)
#Altitude correlated with total conn. -> look for unique contributions
env.alc.reg <- subset(env.alc, select = c(Altitude, Patch.Size, 
                                 RegionalConnectivity, Suitability_Geology))
env.gen.reg <- subset(env.gen, select = c(Altitude, Patch.Size, 
                                          RegionalConnectivity, Suitability_Geology))
#keep PopID for total connectivity because NA locations have to be removed
env.alc.tot <- subset(env.alc, select = c(PopID, Altitude, Patch.Size, 
                                 TotalConnectivity, Suitability_Geology))
env.gen.tot <- subset(env.gen, select = c(PopID,Altitude, Patch.Size, 
                                          TotalConnectivity, Suitability_Geology))

pairs.panels(env.alc.tot, scale = T) #geology and patch size (plants and patchsize OK)
pairs.panels(env.gen.tot, scale = T) #similar
pairs.panels(env.alc.reg, scale = T) #altitude and regional connectivity still highly correlated
pairs.panels(env.gen.reg, scale = T) #similar

#patchsize and suitability correlated (0.81)

sum(is.na(env.alc.tot)) #5 NAs in total connectivity
sum(is.na(env.alc.reg))

##########
###RDA####
##########
#totalconnectivity has NAs

#hellinger transform data (more weight to rare alleles)
library(vegan)
AlconAF.hel <- decostand(AlconAF, 'hellinger')
GentianAF.hel <- decostand(GentianAF, 'hellinger')
setwd('/home/jonas/Documents/Masterthesis/data/Analyses/co-inertia')
write.table(AlconAF, 'alconAF', col.names = T, row.names = T)
write.table(GentianAF, 'gentianAF', col.names = T, row.names = T)

#1 #ALCON + REG. CONN.
rda.AR <- rda(AlconAF.hel~. , data = env.alc.reg, scale = T)
rda.AR
plot(rda.AR)
RsquareAdj(rda.AR) #adjR² = 16.1% -> quite OK

#check if transformation was useful
rda.AR.untr <- rda(AlconAF~., data = env.alc.reg, scale = T)
#par(mfrow=c(1,1))
plot(rda.AR.untr) #snps little less spreaded
RsquareAdj(rda.AR.untr) #15.6% -> less
#model assesment
anova.cca(rda.AR, parallel=getOption("mc.cores")) # ***0.001
anova.cca(rda.AR, by ="axis", parallel=getOption("mc.cores")) #axis 1 and 2 significant

rda.AR.r.sq <- anova.cca(rda.AR, by = "term", parallel=getOption("mc.cores"))
rda.AR.r.sq$Rsq <- lapply(rda.AR.r.sq$Variance, function(x){x/sum(rda.AR.r.sq$Variance)})
p=1
n = 20
rda.AR.r.sq$RsqAdj <- lapply(rda.AR.r.sq$Rsq, function(x){1 - (1 - x) * ((n - 1)/(n-p-1))})

vif.cca(rda.AR) #all below 5

#2 #GENTIAN + REG. CONN.
rda.GR <- rda(GentianAF.hel~., data = env.gen.reg, scale = T)
rda.GR
plot(rda.GR)
#model assesment
RsquareAdj(rda.GR) #12.5%
anova.cca(rda.GR, parallel=getOption("mc.cores")) # **0.002
anova.cca(rda.GR, by = "axis", parallel=getOption("mc.cores")) #only RDA1

rda.GR.r.sq <- anova.cca(rda.GR, by = "term", parallel=getOption("mc.cores"))
rda.GR.r.sq$Rsq <- lapply(rda.GR.r.sq$Variance, function(x){x/sum(rda.GR.r.sq$Variance)})
p=1
n = 22
rda.GR.r.sq$RsqAdj <- lapply(rda.GR.r.sq$Rsq, function(x){1 - (1 - x) * ((n - 1)/(n-p-1))})

vif.cca(rda.GR) #all below 5

##ALCON + TOTAL CONN.
#remove NA locations for calculations of total connectivity
ra <- env$PopID[is.na(env$TotalConnectivity)] #"NW18" "NW5"  "NW9"  "SE5"  "SW6" 
env.alc.tot <- env.alc.tot[!(env.alc.tot$PopID %in% ra),]
AlconAF.hel.tot <- AlconAF.hel[!(rownames(AlconAF.hel) %in% ra),]
env.alc.tot <- subset(env.alc.tot, select = -c(PopID)) #PopID can be removed
#RDA
rda.AT <- rda(AlconAF.hel.tot ~ . , data = env.alc.tot, scale = T)
rda.AT
plot(rda.AT)
#model assesment
RsquareAdj(rda.AT) #13.3% -> regional conn explains more variation, but less populations here
anova.cca(rda.AT, parallel=getOption("mc.cores")) # **0.04
anova.cca(rda.AT, by = "axis", parallel=getOption("mc.cores")) #only RDA1 sign

rda.AT.r.sq <- anova.cca(rda.AT, by = "term", parallel=getOption("mc.cores"))
rda.AT.r.sq$Rsq <- lapply(rda.AT.r.sq$Variance, function(x){x/sum(rda.AT.r.sq$Variance)})
p=1
n = 20
rda.AT.r.sq$RsqAdj <- lapply(rda.AT.r.sq$Rsq, function(x){1 - (1 - x) * ((n - 1)/(n-p-1))})

vif.cca(rda.AT) #all below 2

## GENTIAN + TOT. CONN.
#remove NA localities
rg <- ra[-1] #NW18 already deleted in gentian
env.gen.tot <- env.gen.tot[!(env.gen.tot$PopID %in% rg),]
GentianAF.hel.tot <- GentianAF.hel[!(rownames(GentianAF.hel) %in% rg),]
env.gen.tot <- subset(env.gen.tot, select = -c(PopID))
#RDA
rda.GT <- rda(GentianAF.hel.tot ~. , data= env.gen.tot, scale = T)
rda.GT
plot(rda.GT) # some SNPs seem correlated with total connectivity and negatively with suitability 
#model assesment
RsquareAdj(rda.GT) #12.6 -> Gentian better explained by total connectivity, even with less populations!
anova.cca(rda.GT, parallel=getOption("mc.cores")) # ** 0.006
anova.cca(rda.GT, by = "axis", parallel=getOption("mc.cores")) #only 1

rda.GT.r.sq <- anova.cca(rda.GT, by = "term", parallel=getOption("mc.cores"))
rda.GT.r.sq$Rsq <- lapply(rda.GT.r.sq$Variance, function(x){x/sum(rda.GT.r.sq$Variance)})
p=1
n = 22
rda.GT.r.sq$RsqAdj <- lapply(rda.GT.r.sq$Rsq, function(x){1 - (1 - x) * ((n - 1)/(n-p-1))})

vif.cca(rda.GT) #all below 2
anova.cca(rda.GT)


##VISUALISATION##
#################
par(mfrow=c(2,2))
plot(rda.GT)
plot(rda.AT)
plot(rda.GR)
plot(rda.AR)

#GGPLOT#

#install.packages('devtools')
#library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(ggplot2)

#gentian
popnamesG <- popnamesG[!(popnamesG$V2 %in% rg),]
regionsG <- gsub('[0-9]+', '', popnamesG$V2) #removes number from name-> NE, NW, etc.
ggord(rda.GT, regionsG, grp_title = "Region", alpha_el = 0.3, 
      cols =c('plum3','gold2','steelblue2','springgreen4'),
      addsize = 0.4, addcol = "grey", size = 2,
      vec_lab = list(Suitability_Geology = 'Suitability', Patch.Size = 'PatchSize', TotalConnectivity = ' Connectivity', Altitude = 'Altitude'),
      ext = 0.90, vec_ext = 1, labcol ='grey11', arrow=0.2, txt = 2.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Gentian')
?ggord
#alcon
popnamesA <- popnamesA[1:25,] #removes ker location
regionsA <- gsub('[0-9]+', '', popnamesA$V2)
ggord(rda.AR, regionsA, grp_title = "Region", alpha_el = 0.3,
      cols =c('plum3','gold2','steelblue2','springgreen4'),
      addsize = 0.4, addcol = "grey", size = 2,
      vec_lab = list(Suitability_Geology = 'Suitability', Patch.Size = 'PatchSize', RegionalConnectivity = 'Connectivity', Altitude = 'Altitude'),
      ext = c(0.9,1.1), vec_ext = 1, labcol ='grey11', arrow=0.2, txt = 2.75) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none')+
  ggtitle('Alcon')

##SEARCH CANDIDATE SNPs##
#########################
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

outliertest <- function(RDA){
  load <- scores(RDA, choices = c(1:3), display = "species")
  candA1 <- outliers(load[,1],2.5)
  candA2 <- outliers(load[,2],2.5)
  candA3 <- outliers(load[,3],2.5)
  return(c(candA1, candA2, candA3))
}

outliertest2 <- function(RDA){
  load <- scores(RDA, choices = c(1:3), display = "species")
  candA1 <- outliers(load[,1],2)
  candA2 <- outliers(load[,2],2)
  candA3 <- outliers(load[,3],2)
  return(c(candA1, candA2, candA3))
}

#partial RDA gentian-total connectivity
#all combinations
envfactors <- c('Altitude', 'Suitability_Geology', 'TotalConnectivity', 'Patch.Size')

#Alcon -> regional connectivity and 2 axes significant
load.rda.AR <- scores(rda.AR, choices = c(1,2,3), display = "species")
hist(load.rda.AR[,1], main="Loadings on RDA1")
hist(load.rda.AR[,2], main="Loadings on RDA2")
hist(load.rda.AR[,3], main="Loadings on RDA3") 

cand.rda1.AR <- outliers(load.rda.AR[,1],3) #no candidates (even at 2.5 SD cutoff)
cand.rda2.AR <- outliers(load.rda.AR[,2],3) #none
cand.rda3.AR <- outliers(load.rda.AR[,3],3)

as.numeric(gsub("Locus_", "", names(outliertest(rda.AR))))
as.numeric(gsub("Locus_", "", names(outliertest2(rda.AR))))

#Gentian -> total connectivity
load.rda.GT <- scores(rda.GT, choices = c(1,2), display = "species")
hist(load.rda.GT[,1], main="Loadings on RDA1")
hist(load.rda.GT[,2], main="Loadings on RDA2")

cand.rda1.GT <- outliers(load.rda.GT[,1],3)
cand.rda2.GT <- outliers(load.rda.GT[,2],3) #not siginificant! (scale)

as.numeric(gsub("Locus_", "", names(outliertest(rda.GT))))
as.numeric(gsub("Locus_", "", names(outliertest2(rda.GT))))

#
##
###
####
###
##
#
#############
#PARTIAL RDA#
#############

##ALCON##
###########
#RDA Function
#######
#function to return r²adj and pvalue from rda
p.rda.Alc <- function(AF, test){
  c <- c(1:4)
  if (length(test)==1){
    t <- grep(test, colnames(env.alc.reg)) #retrieve index of env var to test
    c <- c[-t] #all others are conditional
    a <- rda(AF ~ env.alc.reg[,t] + Condition(env.alc.reg[,c[1]],
                                              env.alc.reg[,c[2]],
                                              env.alc.reg[,c[3]]), scale = T)
  }
  if (length(test)==2){
    t <- c(which(colnames(env.alc.reg)==test[1]), which(colnames(env.alc.reg)==test[2])) #retrieve index of env var to test
    c <- c[-t[1:2]] #all others are conditional
    a <- rda(AF ~ env.alc.reg[,t[1]] + env.alc.reg[,t[2]] + 
               Condition(env.alc.reg[,c[1]], env.alc.reg[,c[2]]), scale = T)
  }
  if (length(test)==3){
    t <- which(colnames(env.alc.reg)==test[1] | colnames(env.alc.reg)==test[2] |
                 colnames(env.alc.reg)==test[3] ) #retrieve index of env var to test
    c <- c[-t[1:3]] #all others are conditional
    a <- rda(AF ~ env.alc.reg[,t[1]] + env.alc.reg[,t[2]] + env.alc.reg[,t[3]] +
               Condition(env.alc.reg[,c]), scale = T)
  }
  return(a)
}
?rda
##########
#1 test variable, 3 conditions
envfactors2 <- c("Altitude","Suitability_Geology","RegionalConnectivity","Patch.Size")
rda.res1 <- data.frame(row.names = envfactors2)
for (i in 1:4){
  b <- p.rda.Alc(AF = AlconAF.hel, test = envfactors2[i])
  r1 <- anova.cca(b, step = 50000)
  r2 <- RsquareAdj(b)
  rda.res1[i,1] <- round(r1['Model',4], digits = 4)
  rda.res1[i,2] <- round(r2$adj.r.squared, digits = 4)
}
colnames(rda.res1) <- c('P', 'R²Adj')

#2 test variables, 2 conditions
vars2 <- data.frame(t(combn(envfactors2, m = 2, simplify = T)))
rda.res2 <- data.frame(row.names= c(paste(vars2[,1],vars2[,2], sep = "-")))
for (i in 1:6){
  b <- p.rda.Alc(AF = AlconAF.hel, test = c(vars2[i,1], vars2[i,2]))
  r1 <- anova.cca(b, step = 50000)
  r2 <- RsquareAdj(b)
  rda.res2[i,1] <- round(r1['Model',4], digits = 4)
  rda.res2[i,2] <- round(r2$adj.r.squared, digits = 4)
}
colnames(rda.res2) <- c('P', 'R²Adj')

#3 test, 1 conditions
varss2 <- data.frame(t(combn(envfactors2, m = 3, simplify = T)))
rda.res3 <- data.frame(row.names= c(paste(varss2[,1],varss2[,2],varss2[,3], sep = "-")))
for (i in 1:4){
  b <- p.rda.Alc(AF = AlconAF.hel, test = c(varss2[i,1], varss2[i,2], varss2[i,3]))
  r1 <- anova.cca(b, step = 50000)
  r2 <- RsquareAdj(b)
  rda.res3[i,1] <- round(r1['Model',4], digits = 4)
  rda.res3[i,2] <- round(r2$adj.r.squared, digits = 4)
}
colnames(rda.res3) <- c('P', 'R²Adj')

rda.res.A <- rbind(rda.res1, rda.res2, rda.res3)

write.table(rda.res.A, 'pRDAResultsAlc')

#which significant models have outliers
n <- rownames(rda.res.A[(rda.res.A[,'P'] < 0.05), ]) #these are significant (ALL!)
k <- strsplit(n[1:length(n)],"-")                        #make ready to calculate outliers

outlA <- list()
for (i in 1:length(k)){
  a <- p.rda.Alc(AF = AlconAF.hel, test = c(k[[i]]))
  o <- outliertest(a)
  if (length(o) != 0){
    outlA[[paste(i)]] <- gsub("Locus_","",names(o))
  }
}
names(outlA) <- n[as.numeric(names(outlA))]

outlA2 <- list()
for (i in 1:length(k)){
  a <- p.rda.Alc(AF = AlconAF.hel, test = c(k[[i]]))
  o <- outliertest2(a)
  if (length(o) != 0){
    outlA2[[paste(i)]] <- gsub("Locus_","",names(o))
  }
}

names(outlA2) <- n

save(outlA2, file = 'outlA2.RData')

##GENTIAN##
##############
#RDA Function
#######
#function to return r²adj and pvalue from rda
p.rda.Gen <- function(AF, test){
  c <- c(1:4)
  if (length(test)==1){
    t <- grep(test, colnames(env.gen.tot)) #retrieve index of env var to test
    c <- c[-t] #all others are conditional
    a <- rda(AF ~ env.gen.tot[,t] + Condition(env.gen.tot[,c[1]],
                                              env.gen.tot[,c[2]],
                                              env.gen.tot[,c[3]]), scale = T)
  }
  if (length(test)==2){
    t <- c(which(colnames(env.gen.tot)==test[1]), which(colnames(env.gen.tot)==test[2])) #retrieve index of env var to test
    c <- c[-t[1:2]] #all others are conditional
    a <- rda(AF ~ env.gen.tot[,t[1]] + env.gen.tot[,t[2]] + 
               Condition(env.gen.tot[,c[1]], env.gen.tot[,c[2]]), scale = T)
  }
  if (length(test)==3){
    t <- which(colnames(env.gen.tot)==test[1] | colnames(env.gen.tot)==test[2] |
                 colnames(env.gen.tot)==test[3] ) #retrieve index of env var to test
    c <- c[-t[1:3]] #all others are conditional
    a <- rda(AF ~ env.gen.tot[,t[1]] + env.gen.tot[,t[2]] + env.gen.tot[,t[3]] +
               Condition(env.gen.tot[,c]), scale = T)
  }
  return(a)
}
##########
#1 test variable, 3 conditions
rda.res1 <- data.frame(row.names = envfactors)
for (i in 1:4){
  b <- p.rda.Gen(AF = GentianAF.hel.tot, test = envfactors[i])
  r1 <- anova.cca(b, step = 50000)
  r2 <- RsquareAdj(b)
  rda.res1[i,1] <- round(r1['Model',4], digits = 3)
  rda.res1[i,2] <- round(r2$adj.r.squared, digits = 3)
}
colnames(rda.res1) <- c('P', 'R²Adj')

#2 test variables, 2 conditions
vars <- data.frame(t(combn(envfactors, m = 2, simplify = T)))
rda.res2 <- data.frame(row.names= c(paste(vars[,1],vars[,2], sep = "-")))
for (i in 1:6){
  b <- p.rda.Gen(AF = GentianAF.hel.tot, test = c(vars[i,1], vars[i,2]))
  r1 <- anova.cca(b, step = 50000)
  r2 <- RsquareAdj(b)
  rda.res2[i,1] <- round(r1['Model',4], digits = 3)
  rda.res2[i,2] <- round(r2$adj.r.squared, digits = 3)
}
colnames(rda.res2) <- c('P', 'R²Adj')

#3 test, 1 conditions
varss <- data.frame(t(combn(envfactors, m = 3, simplify = T)))
rda.res3 <- data.frame(row.names= c(paste(varss[,1],varss[,2],varss[,3], sep = "-")))
for (i in 1:4){
  b <- p.rda.Gen(AF = GentianAF.hel.tot, test = c(varss[i,1], varss[i,2], varss[i,3]))
  r1 <- anova.cca(b, step = 50000)
  r2 <- RsquareAdj(b)
  rda.res3[i,1] <- round(r1['Model',4], digits = 3)
  rda.res3[i,2] <- round(r2$adj.r.squared, digits = 3)
}
colnames(rda.res3) <- c('P', 'R²Adj')

rda.res.G <- rbind(rda.res1, rda.res2, rda.res3)

write.table(rda.res.G, 'pRDAResultsGen')

#which significant models have outliers
n <- rownames(rda.res.G[(rda.res.G[,'P'] < 0.05), ]) #these are significant
k <- strsplit(n[1:length(n)],"-")                        #make ready to calculate outliers

outlG <- list()
for (i in 1:length(k)){
  a <- p.rda.Gen(AF = GentianAF.hel.tot, test = c(k[[i]]))
  o <- outliertest(a)
  if (length(o) != 0){
    outlG[[paste(i)]] <- gsub("Locus_","",names(o))
  }
}
names(outlG) <- n[as.numeric(names(outlG))]

outlG2 <- list()
for (i in 1:length(k)){
  a <- p.rda.Gen(AF = GentianAF.hel.tot, test = c(k[[i]]))
  o <- outliertest2(a)
  if (length(o) != 0){
    outlG2[[paste(i)]] <- gsub("Locus_","",names(o))
  }
}
names(outlG2) <- n

save(outlG2, file = 'outlG2.RData')

for (i in 1:9){
  print(outlG2[i])
}

for(i in 1:length(outlA2)){
  print(outlA2[i])
}

combo <- list()
for (i in 1:length(outlG2)){
  for (j in 1:length(combn(1:length(outlG2), m =i)[1,])){
    if (length(Reduce(intersect, outlG2[combn(1:length(outlG2), m = i)[,j]]))!= 0){
      combo[[paste(i,j, sep = "-")]] <- Reduce(intersect, outlG2[combn(1:length(outlG2), m = i)[,j]])
    }
  }
}
strsplit(names(combo), "-")
for (i in 1:length(combo)){
  
}

length(combn(1:9, m =1)[1,])+
length(combn(1:9, m =2)[1,])+
length(combn(1:9, m =3)[1,])+
length(combn(1:9, m =4)[1,])+
length(combn(1:9, m =5)[1,])+
length(combn(1:9, m =6)[1,])+
length(combn(1:9, m =7)[1,])+
length(combn(1:9, m =8)[1,])+
length(combn(1:9, m =9)[1,])


names(combo[i+i*j]) <- paste(names(outlG2[combn(1:length(outlG2), m = i)[,j]]), sep ="-")
##Venn diagram
rda.res.A$promille <- rda.res.A$`R²Adj`*1000
rda.res.A$promille <- round(rda.res.A$promille)

nmb <- list('1' = 1:89)
for (i in 2:14){
  if (rda.res.A[i,3] != 0){
    nmb[[paste(i)]] <- c(as.integer(sum(rda.res.A[1:(i-1),3])+1):as.integer(sum(rda.res.A[1:i,3])))
  }
  else nmb[[paste(i)]] <- NULL
}
nmb$`14`[length(nmb$`14`)] == sum(rda.res.A$promille) #confirms are operations were correct

#so now we got unique 'promille' index for each partial rda model, now we can compile a list for each environmental variable

alcprop <- list(Altitude = c(nmb$`1`, nmb$`5`,nmb$`6`,nmb$`7`,nmb$`11`,nmb$`12`,nmb$`13`),
                Connectivity = c(nmb$`3`, nmb$`6`,nmb$`8`,nmb$`10`,nmb$`11`,nmb$`13`,nmb$`14`),
                Suitablility = c(nmb$`2`, nmb$`5`,nmb$`8`,nmb$`11`,nmb$`12`,nmb$`14`),
                PatchSize = c(nmb$`4`,nmb$`7`,nmb$`10`,nmb$`12`,nmb$`13`,nmb$`14`))

library(ggVennDiagram)
library(ggplot2)
ggVennDiagram(alcprop, label = "percent", label_alpha = 0, label_size = 2.5, 
              set_size = 3, label_percent_digit = 1,
              category.names= c('Altitude', 'Connectivity', 'Suitability', 'Patch Size'))+
  ggplot2::scale_fill_gradient(low="white",high = "#116E8A")+
  ggtitle('Alcon')+
  theme(legend.position = 'none',  plot.margin=unit(c(0,0.3,0,0),"cm"))



#Gentian
########
rda.res.G$promille <- rda.res.G$`R²Adj`*1000
rda.res.G[rda.res.G[,2]<0,2:3] <- c(0,0) #negative value to zero

nmbG <- list('1' = 1:60)
for (i in 2:14){
  if (rda.res.G[i,3] != 0){
    nmbG[[paste(i)]] <- c(as.integer(sum(rda.res.G[1:(i-1),3])+1):as.integer(sum(rda.res.G[1:i,3])))
  }
  else nmbG[[paste(i)]] <- NULL
}
nmbG$`14`[length(nmbG$`14`)] == sum(rda.res.G$promille) #confirms are operations were correct
sum(rda.res.G$`R.Adj`)

genprop <- list(Altitude = c(nmbG$`1`, nmbG$`5`,nmbG$`6`,nmbG$`7`,nmbG$`11`,nmbG$`12`,nmbG$`13`),
                  Connectivity = c(nmbG$`3`, nmbG$`6`,nmbG$`8`,nmbG$`10`,nmbG$`11`,nmbG$`13`,nmbG$`14`),
                  Suitablility = c(nmbG$`5`,nmbG$`8`,nmbG$`11`,nmbG$`12`,nmbG$`14`),
                  PatchSize = c(nmbG$`7`,nmbG$`10`,nmbG$`12`,nmbG$`13`,nmbG$`14`))

library(ggVennDiagram)
ggVennDiagram(genprop, label = "percent", label_alpha = 0, label_size = 2.5, 
              set_size = 3, label_percent_digit = 1,
              category.names= c('Altitude', 'Connectivity', 'Suitability', 'Patch Size'))+
  ggplot2::scale_fill_gradient(low="white",high = "#116E8A")+
  ggtitle('Gentian')+
  theme(legend.position = 'none',  plot.margin=unit(c(0,0.3,0,0),"cm"))

rownames(rda.res.A[rda.res.A$P<0.05,])
rownames(rda.res.G[rda.res.G$P<0.05,])



#for (i in 1:length(outlA2)){
  outlA2[[i]] <- as.numeric(outlA2[[i]])
}
names(outlA2) <- c('Alt', 'Conn','Alt-Sui', 'Alt-Conn','Alt-PS','Suit-Conn','Conn-PS','Alt-Sui-PS','Alt-Conn-PS','Suit-Conn-PS','Full')

#for (i in 1:length(outlG2)){
  outlG2[[i]] <- as.numeric(outlG2[[i]])
}
names(outlG2) <- c('Alt','Conn','Alt-Sui','Alt-Conn','Alt-PS','Sui-Conn', 'Conn-PS','Alt-Sui-Conn','Alt-Sui-PS','Alt-Conn-PS','Sui-Conn-PS','Full')



###############

#altitude 
p.rda.GT.A <- rda(GentianAF.hel.tot ~ Altitude + Condition(Suitability_Geology, 
                                                           TotalConnectivity, 
                                                           Patch.Size),
                  data = env.gen.tot, scale =T)
ano<- anova.cca(p.rda.GT.A, by = "axis", step = 1000) #p(altitude ) = 0.014
RsquareAdj(p.rda.GT.A) #6% of variance
outliertest(p.rda.GT.A) #empty

#total connectivity
p.rda.GT.C <- rda(GentianAF.hel.tot ~ TotalConnectivity + Condition(Suitability_Geology, 
                                                           Altitude, 
                                                           Patch.Size),
                  data = env.gen.tot, scale =T)
anova.cca(p.rda.GT.C, by = "axis", step = 1000) #p(total connectivity) = 0.001
RsquareAdj(p.rda.GT.C) #10.8% of var
outliertest(p.rda.GT.C) #empty

#suitability of geology
p.rda.GT.S <- rda(GentianAF.hel.tot ~ Suitability_Geology + Condition(TotalConnectivity , 
                                                                    Altitude, 
                                                                    Patch.Size),
                  data = env.gen.tot, scale =T)
anova.cca(p.rda.GT.S, by = "axis", step = 1000) #insiginificant
RsquareAdj(p.rda.GT.S) #0

#patch size
p.rda.GT.P <- rda(GentianAF.hel.tot ~ Patch.Size + Condition(TotalConnectivity , 
                                                                      Altitude, 
                                                                      Suitability_Geology ),
                  data = env.gen.tot, scale =T)
anova.cca(p.rda.GT.P, by = "axis", step = 1000) #insiginificant
RsquareAdj(p.rda.GT.P) #-1.2%

###############
#alcon - regional conn.
#altitude
plot(rda.AR)
p.rda.AR.A <- rda(AlconAF.hel ~ Altitude + Condition(RegionalConnectivity,
                                                     Patch.Size,
                                                     Suitability_Geology),
                  data = env.alc.reg, scale = T)
anova.cca(p.rda.AR.A, by = "axis", step = 1000) #0.006
RsquareAdj(p.rda.AR.A) #5.1% of variance
outliertest(p.rda.AR.A) #empty

#regional connectivity
p.rda.AR.C <- rda(AlconAF.hel ~ RegionalConnectivity + Condition(Altitude,
                                                                 Patch.Size,
                                                                 Suitability_Geology),
                  data = env.alc.reg, scale = T)
anova.cca(p.rda.AR.C, by = "axis", step = 1000) #0.002
RsquareAdj(p.rda.AR.C) #5.4% of variance
outliertest(p.rda.AR.C) #empty

#patch size
p.rda.AR.P <- rda(AlconAF.hel ~ Patch.Size + Condition(Altitude,
                                                       RegionalConnectivity,
                                                       Suitability_Geology),
                  data = env.alc.reg, scale = T)
anova.cca(p.rda.AR.P) # 0.038 patch size signicantly correlated with alcon allele frequencies but not with gentian allele frequencies
RsquareAdj(p.rda.AR.P) #2.8% -> only small part of the variance correlated with patch size
outliertest(p.rda.AR.P) #empty

#suitability of geology
p.rda.AR.S <- rda(AlconAF.hel ~ Suitability_Geology + Condition(Altitude,
                                                       RegionalConnectivity,
                                                       Patch.Size),
                  data = env.alc.reg, scale = T)
anova.cca(p.rda.AR.S) #0.016 -> significant! 
RsquareAdj(p.rda.AR.S) #4.5%
outliertest(p.rda.AR.S) #empty
