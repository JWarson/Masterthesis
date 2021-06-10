setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setwd to directory where script is located

#index of divergent and balancing SNPs 
index.div.alc <- c(98,99,109,  #bayescan
                   36,         #lfmm altitude (98 was double)
                   84,81,79,126,67,41,47, 26,23,14, 65,28,40,51,69)

index.bal.alc <- c(29,25,116,108,7)

index.div.gen <- c(5,23,29,59, #bayescan
                   39,82)      #lfmm (23 was double)
                  
index.bal.gen <- c(21,53,84,72,77,73,80,87,88,32,97)

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("made4")
library(made4)

AF.gen <- read.table('gentianAF', h = T , sep= " ")
AF.alc <- read.table('alconAF', h = T, sep = " ")

x <- rownames(AF.alc[rownames(AF.alc)%in%rownames(AF.gen),])
#rownames(AF.gen[rownames(AF.gen)%in%rownames(AF.alc),])

#data of 24 populations is available for both species
#select only these and transpose
?cia #genes as rows, pops as column

AF.gen <- data.frame(t(AF.gen[x,]))
AF.alc <- data.frame(t(AF.alc[x,]))

cia(AF.gen, AF.alc) #RV 0.74 (total correlation quite strong as expected)

#balancing only
cia.ABal.GBal <- cia(AF.gen[index.bal.gen,], AF.alc[index.bal.alc,])
cia.ABal.GBal$coinertia$RV #0.270

#Balancing - Neutral
cia.Abal.Gneu <- cia(AF.alc[index.bal.alc,], 
                     AF.gen[-(c(index.bal.gen,index.div.gen)),])
cia.Abal.Gneu$coinertia$RV #0.256

cia.Aneu.Gbal <- cia(AF.alc[-(c(index.bal.alc,index.div.alc)),],
                     AF.gen[index.bal.gen,])
cia.Aneu.Gbal$coinertia$RV #0.419 (highest correlation)
#selecting retained axes based on scree plot (cia.scan = T) gives same results

#divergent only
cia.Adiv.Gdiv <- cia(AF.gen[index.div.gen,], AF.alc[index.div.alc,])
cia.Adiv.Gdiv$coinertia$RV #0.506
#divergent - neutral
cia.Adiv.Gneu <- cia(AF.alc[index.div.alc,], 
                     AF.gen[-(c(index.bal.gen,index.div.gen)),])
cia.Adiv.Gneu$coinertia$RV #0.636

cia.Aneu.Gdiv <- cia(AF.alc[-(c(index.bal.alc,index.div.alc)),],
                     AF.gen[index.div.gen,])
cia.Aneu.Gdiv$coinertia$RV #0.562

#divergent - balancing
cia.Adiv.Gbal <- cia(AF.alc[index.div.alc,], AF.gen[index.bal.gen,])
cia.Adiv.Gbal$coinertia$RV
cia.Abal.Gdiv <- cia(AF.alc[index.bal.alc,], AF.gen[index.div.gen,])
cia.Abal.Gdiv$coinertia$RV

#neutral - neutral
cia.Aneu.Gneu <- cia(AF.alc[-(c(index.bal.alc,index.div.alc)),],
                     AF.gen[-(c(index.bal.gen,index.div.gen)),])
cia.Aneu.Gneu$coinertia$RV

##Sample Sizes-> equal(11)
#load bayescan results
bayes.alc <- read.delim(file = '/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/alcon/alcon_fst.txt', h = T,sep = ' ')
bayes.alc <- bayes.alc[,3:7]
colnames(bayes.alc) <- c('prob','log10(PO)','qval','alpha','fst')

bayes.gen <- read.delim(file = '/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/gentian/gentian_fst.txt', h = T,sep = ' ')
bayes.gen <- bayes.gen[,3:7]
colnames(bayes.gen) <- c('prob','log10(PO)','qval','alpha','fst')

#select 11 most 'balanced' genes (lowest alpha)
rownames(bayes.alc[order(bayes.alc$alpha),])[1:11]
index.bal.alc
iba2 <- order(bayes.alc$alpha)[1:11]
ina2 <- setdiff(c(1:127), c(iba2,index.div.alc)) #too sample random neutral
ina3 <- order(abs(bayes.alc$alpha))[1:11] #11 snps with alpha closest to 0

ing2 <- setdiff(c(1:105),(c(index.bal.gen,index.div.gen)))
ing3 <- order(abs(bayes.gen$alpha))[1:11]

#11 balancing inertia
cia.ABal.GBal2 <- cia(AF.gen[index.bal.gen,], AF.alc[iba2,])
cia.ABal.GBal2$coinertia$RV #0.240

#random sample of 11 neutral
#balancing(alcon) - neutral(gentian)
cia.Abal.Gneu2 <- cia(AF.alc[iba2,], 
                     AF.gen[sample(ing2, size=11),])
cia.Abal.Gneu2$coinertia$RV #0.28

res <- c()
for (i in 1:1000){
  cia <- cia(AF.alc[iba2,], 
             AF.gen[sample(ing2, size=11),])
  res[i] <- cia$coinertia$RV
}
c(mean(res), sd(res)) #0.297 +/- 0.054

#neutral(alcon) - balancing(gentian)
cia.Aneu.Gbal2 <- cia(AF.alc[sample(ina2, size =11),],
                     AF.gen[index.bal.gen,])
cia.Aneu.Gbal2$coinertia$RV #0.33 (highest correlation)

res1 <- c()
for (i in 1:1000){
  cia <- cia(AF.alc[sample(ina2, size =11),],
             AF.gen[index.bal.gen,])
  res1[i] <- cia$coinertia$RV
}
c(mean(res1), sd(res1)) #0.309 +/- 0.052

#11 'most neutral' snps
cia.Abal.Gneu3 <- cia(AF.alc[iba2,], 
                      AF.gen[ing3,])
cia.Abal.Gneu3$coinertia$RV #0.360

cia.Aneu.Gbal3 <- cia(AF.alc[ina3,],
                      AF.gen[index.bal.gen,])
cia.Aneu.Gbal3$coinertia$RV #0.355 (highest correlation)
?cia

####################
#instead of 6 extra alcon, only take 5 balanced gentian genes
ibg4 <- order(bayes.gen$alpha)[1:5]
cia.ABal.GBal4 <- cia(AF.alc[index.bal.alc,], AF.gen[ibg4,])
cia.ABal.GBal4$coinertia$RV #0.152

res3 <- c()
for (i in 1:500){
  cia <- cia(AF.alc[sample(ina2, size =5),],
             AF.gen[ibg4,])
  res3[i] <- cia$coinertia$RV
}
c(mean(res3), sd(res3)) #0.123 +/- 0.046

res4 <- c()
for (i in 1:500){
  cia <- cia(AF.alc[index.bal.alc,],
             AF.gen[sample(ing2, size = 5),])
  res4[i] <- cia$coinertia$RV
}
c(mean(res4), sd(res4)) #0.151 +/- 0.046
#similar results

####################
#divergent snps
##
idg <- order(bayes.gen$alpha)[(105-18):105]
#div - div 
cia.Adiv.Gdiv <- cia(AF.alc[index.div.alc,], AF.gen[idg,])
cia.Adiv.Gdiv$coinertia$RV  ##0.484

#Aneu-Gdiv
res5 <- c()
for (i in 1:500){
  cia <- cia(AF.alc[sample(ina2, size =19),],
             AF.gen[idg,])
  res5[i] <- cia$coinertia$RV
}
c(mean(res5), sd(res5)) #0.520 +/- 0.061

#Adiv-Gneu
res6 <- c()
for (i in 1:500){
  cia <- cia(AF.alc[index.div.alc,],
             AF.gen[sample(ing2, size = 19),])
  res6[i] <- cia$coinertia$RV
}
c(mean(res6), sd(res6)) #0.545 +/- 0.038


##via dudi.pca and coinertia() -> not working because of unequal rows
dudi.div.alc <- dudi.pca(AF.alc[index.div.alc,], scannf = FALSE, nf = 2) 
dudi.bal.alc <- dudi.pca(AF.alc[index.bal.alc,], scannf = FALSE, nf = 2)
dudi.neu.alc <- dudi.pca(AF.alc[-(c(index.bal.alc,index.div.alc)),],scannf = FALSE, nf = 1)

dudi.div.gen <- dudi.pca(AF.gen[index.div.gen,], scannf = FALSE, nf = 3) #3 AF.genes was best
dudi.bal.gen <- dudi.pca(AF.gen[index.bal.gen,], scannf = FALSE, nf = 1) #also 3, although less clear
dudi.neu.gen <- dudi.pca(AF.gen[-(c(index.bal.gen,index.div.gen)),],
                         scannf = FALSE, nf = 1) #only 2 seems best
coinertia(dudi.bal.alc, dudi.bal.gen)
