library(hierfstat)
library(readxl)
alcon <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/alconSNP', h = T, sep = "\t")
gentian <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/gentianSNP', h = T, sep = "\t")

#conversion to hierfstat format (two alleles in 1 column)

##Alcon
alc <- data.frame(matrix(nrow = 661, ncol = 127))
colnA <- vector()
for (i in 1:ncol(alc)){
  alc[,i] <- paste(alcon[,2*i+1], alcon[,2*i+2], sep="")
  colnA <- append(colnA, colnames(alcon[2*i+1]))
}
alc[] <- lapply(alc, function(x) as.numeric(as.character(x)))
alc <- cbind(alcon$PopID, alc)
colnames(alc) <- c('PopID',colnA)
alc$PopID <- as.integer(as.factor(alc$PopID))

##Gentian
gen <- data.frame(matrix(nrow = 652, ncol = 105))
colnG <- vector()
for (i in 1:ncol(gen)){
  gen[,i] <- paste(gentian[,2*i+1], gentian[,2*i+2], sep="")
  colnG <- append(colnG, colnames(gentian[2*i+1]))
}
gen[] <- lapply(gen, function(x) as.numeric(as.character(x)))
gen <- cbind(gentian$PopID, gen)
colnames(gen) <- c('PopID',colnG)
gen$PopID <- as.integer(as.factor(gen$PopID))

#calculate rarefied allele counts
#Alcon
str(alc)
alcAR <- allelic.richness(alc)
alcAR <- alcAR$Ar
colnames(alcAR)<-unique(alcon$PopID)
alcAR$mean <- rowMeans(alcAR)

#Gentian
str(gen)
genAR <- allelic.richness(gen)
genAR <- genAR$Ar
colnames(genAR) <- unique(gentian$PopID)
genAR$mean <- rowMeans(genAR)

#index of divergent and balancing SNPs 
index.div.alc <- c(98,99,109,  #bayescan
                   36,         #lfmm altitude (98 was double)
                   84,81,79,126,67,41,47, 26,23,14,65,28,40,51,69)
index.bal.alc <- c(29,25,116,108,7)

index.div.gen <- c(5,23,29,59, #bayescan
                   39,82)      #lfmm (23 was double)
index.bal.gen <- c(21,53,84,72,77,73,80,87,88,32,97)

#calculate average diversity per population for balancing and divergent SNPs
#alcon
alc.avAR <- data.frame(row.names = c('Divergent', 'Balancing', 'Neutral'))
for (i in 1:length(alcAR)){
  alc.avAR[1,i] <- mean(alcAR[index.div.alc, i])
  alc.avAR[2,i] <- mean(alcAR[index.bal.alc, i])
  alc.avAR[3,i] <- mean(alcAR[-c(index.div.alc, index.bal.alc), i])
}
colnames(alc.avAR) <- colnames(alcAR)

#gentian
gen.avAR <- data.frame(row.names = c('Divergent', 'Balancing', 'Neutral'))
for (i in 1:length(genAR)){
  gen.avAR[1,i] <- mean(genAR[index.div.gen, i])
  gen.avAR[2,i] <- mean(genAR[index.bal.gen, i])
  gen.avAR[3,i] <- mean(genAR[-c(index.div.gen, index.bal.gen),i])
}
colnames(gen.avAR) <- colnames(genAR)

##correlate adaptive potential to connectivity
env <- read.csv('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/environment', h = T, sep = '\t')

#alcon
setdiff(env$PopID, colnames(alc.avAR))
setdiff(colnames(alc.avAR),env$PopID)
match(env$PopID, colnames(alc.avAR))
match(colnames(alc.avAR),env$PopID)

#divergent
x.a = env$TotalConnectivity[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(alc.avAR))]
y.da = as.numeric(alc.avAR['Divergent', colnames(alc.avAR) %in% env$PopID[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(alc.avAR))]])
cor.test(x.a, y.da) #insignificant

#balancing
y.ba = as.numeric(alc.avAR['Balancing', colnames(alc.avAR) %in% env$PopID[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(alc.avAR))]])
cor.test(x.a, y.ba) #insignificant

#neutral
y.na = as.numeric(alc.avAR['Neutral', colnames(alc.avAR) %in% env$PopID[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(alc.avAR))]])
cor.test(x.a, y.na) #neutral alcon significant!

library(ggplot2)
ggplot(data = NULL, aes(x=x.a, y = y.na))+
  geom_point()+
  ylab('Rarefied Allelic Richness')+ xlab('Connectivity')+
  theme_classic()+
  geom_smooth(method ='lm', color='#116E8A', alpha = 0)+
  ggtitle('Alcon - Neutral')

ggplot(data = NULL, aes(x=x.a, y = y.ba))+
  geom_point()+
  ylab('Rarefied Allelic Richness')+ xlab('Connectivity')+
  theme_classic()+
  geom_smooth(method ='lm', color='#116E8A', alpha = 0, lty = 2)+
  ggtitle('Alcon - Balancing')



#gentian
setdiff(env$PopID, colnames(gen.avAR))
setdiff(colnames(gen.avAR),env$PopID)

#divergent
x.g = env$TotalConnectivity[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(gen.avAR))]
y.dg = as.numeric(gen.avAR['Divergent', colnames(gen.avAR) %in% env$PopID[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(gen.avAR))]])
cor.test(x.g, y.dg) #insignificant but close (0.08), corr = 38%

#balancing
y.bg = as.numeric(gen.avAR['Balancing', colnames(gen.avAR) %in% env$PopID[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(gen.avAR))]])
corr <- cor.test(x.g, y.bg) #significant! positive correlation between connectivity and rAR

#neutral
y.ng <- as.numeric(gen.avAR['Neutral', colnames(gen.avAR) %in% env$PopID[!is.na(env$TotalConnectivity)&(env$PopID %in% colnames(gen.avAR))]])
cor.test(x.g, y.ng)

dfg <- data.frame(Type = rep(c('Divergent', 'Balancing'), each = length(x.g)),
                  Connectivity = c(x.g,x.g),
                  AR = c(y.dg, y.bg))

ggplot(data = NULL, aes(x=x.g, y = y.bg))+
  geom_point()+
  ylab('Rarefied Allelic Richness')+ xlab('Connectivity')+
  theme_classic()+
  geom_smooth(method ='lm', color='#116E8A', alpha = 0) +
  ggtitle('Gentian - Balancing')
annotate(geom="text", x=8.4, y=2.5, family = 'Arial',
         label=paste("P =",round(corr$p.value, digits = 3))) 

ggplot(data = NULL, aes(x=x.g, y = y.ng))+
  geom_point()+
  ylab('Rarefied Allelic Richness')+ xlab('Connectivity')+
  theme_classic()+
  geom_smooth(method ='lm', color='#116E8A', alpha = 0, lty = 2)+
  ggtitle('Gentian - Neutral')

#both divergent and balancing together
ggplot(data = dfg, aes(x= Connectivity, y = AR, color = Type, fill = Type))+
  geom_point()+
  ylab('Rarefied Allelic Richness')+ xlab('Connectivity')+
  theme_classic()+
  geom_smooth(method = 'lm', alpha = 0)+
  scale_color_manual(values = c('#E5302D','#116E8A'))
annotate(geom="text", x=8.4, y=2.5, family = 'Arial',
         label=paste("P =",round(corr$p.value, digits = 3))) 


#correlate to bayescan probabilty
#Alcon
bayes.alc <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/alcon/alcon_fst.txt', 
                        h = T, sep = ' ')
bayes.alc <- bayes.alc[,3:7]
colnames(bayes.alc) <- c('prob','log10(PO)','qval','alpha','fst')
ARprob <- data.frame(cbind(alcAR$mean, bayes.alc$prob))
colnames(ARprob) <- c('meanAR', 'probability')

library(ggplot2)
library(ggthemes)
library(ggrepel)
library(dplyr)
ggalc <- ggplot(data = ARprob, aes(y = meanAR, x = probability))+
  geom_point()+
  #geom_smooth(method = 'lm')+
  xlab('Q-value') + 
  ylab('Mean Rarefied AR')+
  ggtitle('Alcon')+
  theme_classic()+
  geom_label_repel(data = ARprob %>% filter(probability > 0.70), 
                   aes(label=rownames(ARprob[ARprob$probability>0.7,])),
                   color = '#116E8A')
cor.test(ARprob$probability, ARprob$meanAR) #p-value = 0.1973
ggalc

#Gentian
bayes.gen <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/gentian/gentian_fst.txt', 
                        h = T, sep = ' ')
bayes.gen <- bayes.gen[,3:7]
colnames(bayes.gen) <- c('prob','log10(PO)','qval','alpha','fst')
ARprobG <- data.frame(cbind(genAR$mean, bayes.gen$prob))
colnames(ARprobG) <- c('meanAR', 'probability')

gggen <- ggplot(data = ARprobG, aes(y = meanAR, x = probability))+
  geom_point()+
  #geom_smooth(method = 'lm')+
  xlab('Q-value') + 
  ylab('Mean Rarefied AR')+
  ggtitle('Gentian')+
  theme_classic()+
  geom_label_repel(data = ARprobG %>% filter(probability > 0.70), 
                   aes(label=rownames(ARprobG[ARprobG$probability>0.7,])), 
                   max.overlaps = 11, color = '#116E8A')
gggen
cor.test(ARprobG$probability, ARprobG$meanAR) #p-value = 0.3776
?geom_label_repel
