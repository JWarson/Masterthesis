library(hierfstat)
library(readxl)
alcon <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Vegan/alconSNP', h = T, sep = "\t")
gentian <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Vegan/gentianSNP', h = T, sep = "\t")

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

#Gentian
str(gen)
genAR <- allelic.richness(gen)
genAR <- genAR$Ar
colnames(genAR) <- unique(gentian$PopID)

#select seperate connected and unconnected populations
#alcon: SE (18:21) and Ker (26)
unique(alcon$PopID)
alcAR.conn <- alcAR[,-26] #removes ker
isol.pops <- c(1,5,6,11,14,18,19,20,21,24,25)
alcAR.conn <- alcAR.conn[,-isol.pops] #removes isolated pops
alcAR.unconn <- alcAR[,isol.pops] #selects isolated pops

#gentian (based on dbMEM figure)
colnames(genAR[colnames(genAR) %in% rownames(gen.mems.sampled[gen.mems.sampled$MEM1>1,])])
#note that population NW18 is missing from gentian data set
genAR.conn <- genAR[colnames(genAR) %in% rownames(gen.mems.sampled[gen.mems.sampled$MEM1>1,])]
genAR.unconn <- genAR[!(colnames(genAR) %in% rownames(gen.mems.sampled[gen.mems.sampled$MEM1>1,]))]


#index of divergent and balancing SNPs 
index.div.alc <- c(98,99,109,  #bayescan
                   36,         #lfmm altitude (98 was double)
                   84,81,79,126,67,41,47, 26,23,14,65,28,40,51,69)

index.bal.alc <- c(108,7) #see from bayescanAlcon figure

index.div.gen <- c(5,23,29,59, #bayescan
                   39,82)      #lfmm (23 was double)
index.bal.gen <- as.numeric(rownames(bayes.gen[bayes.gen$balancing==1,])) #bayes.gen dataframe comes from bayescanPlot2.R script

#ALCON
#select the balancing SNPs from connected alcon populations
alcAR.bal <- alcAR[index.bal.alc,-26] #remove also ker location

alc.data <- data.frame(matrix(nrow = 2*25, ncol=4))
colnames(alc.data) <- c('AR', 'SNP', 'POP', 'CONN')

for (i in 1:length(alcAR.bal[1,])){
  alc.data$AR[(1:(i*2))] <- append(na.omit(alc.data$AR), alcAR.bal[,i])
  alc.data$SNP[(1:(i*2))] <- append(na.omit(alc.data$SNP), rownames(alcAR.bal))
  alc.data$POP[(1:(i*2))] <- append(na.omit(alc.data$POP), rep(colnames(alcAR.bal[i]), 
                                                      times = length(alcAR.bal[,i])))
}

alc.data$CONN <- as.factor(ifelse(test = alc.data[,3] %in% c('SE1','SE3','SE4','SE5'), 'Unconnected', 'Connected'))

alc.data$SNP <- as.factor(alc.data$SNP)
alc.data$POP <- as.factor(alc.data$POP)

library(lme4)
library(lmerTest)
str(alc.data)
fit.alc <- lmer(AR ~ CONN + (1|SNP) + (1|POP), data = alc.data, REML = T)

summary(fit.alc)
car::Anova(fit.alc, type = 'III')
anova(fit.alc)
lmerTest::ranova(fit.alc)
plot(density(resid(fit.alc, type = 'deviance')), 
     main = 'Residual Density of Alcon model') #much better after incorporation of random effects
boxplot(AR ~ CONN, data = alc.data)


##Gentian
genAR.bal <- genAR[index.bal.gen,]

gen.data <- data.frame(matrix(nrow = 260, ncol=4))
colnames(gen.data) <- c('AR', 'SNP', 'POP', 'CONN')

for (i in 1:length(genAR.bal[1,])){
  gen.data$AR[(1:(i*10))] <- append(na.omit(gen.data$AR), genAR.bal[,i])
  gen.data$SNP[(1:(i*10))] <- append(na.omit(gen.data$SNP), rownames(genAR.bal))
  gen.data$POP[(1:(i*10))] <- append(na.omit(gen.data$POP), rep(colnames(genAR.bal[i]), 
                                                               times = length(genAR.bal[,i])))
}

gen.data$CONN <- as.factor(ifelse(test = gen.data[,3] %in% rownames(gen.mems.sampled[gen.mems.sampled$MEM1>1,])
                                  , 'Connected', 'Unconnected')) #gen.mems.sampled comes from dbMEM.R file

gen.data$SNP <- as.factor(gen.data$SNP)
gen.data$POP <- as.factor(gen.data$POP)

library(lme4)
library(lmerTest)
str(gen.data)
fit.gen <- lmer(AR ~ CONN + (1|SNP) + (1|POP), data = gen.data, REML = T) #singular fit
ranova(fit.gen) #there is no effect of population
fit.gen2 <- lmer(AR ~ CONN + (1|SNP) , data = gen.data, REML = T) #removing POP resolves it

summary(fit.gen2)
car::Anova(fit.gen2, type = 'III')
anova(fit.gen2)
lmerTest::ranova(fit.gen2)
plot(density(resid(fit.gen, type = 'deviance')), 
     main = 'Residual Density of Gentian model') #not good, two peaks


#test whether plant number affects AR
env <- read.csv('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Vegan/environment', h = T, sep = '\t')
envgen <- env[env$PopID %in% colnames(genAR),]
#all SNPs
cor.test(colMeans(genAR), envgen$Plants)
plot(colMeans(genAR), envgen$Plants)
#balanced SNPs
cor.test(colMeans(genAR.bal), envgen$Plants)
plot(colMeans(genAR.bal), envgen$Plants)


#shapiro.test(alcAR.conn.val)  #no normal distribution
#hist(alcAR.conn.val, main = 'Alcon - Connected (W = 0.621)')
#shapiro.test(alcAR.unconn.val) #no normal distribution
#hist(alcAR.unconn.val, main = 'Alcon - Unconnected (W = 0.748)')

#library(car)
#c(alcAR.conn.val, alcAR.unconn.val)
#data <- data.frame(AR = c(alcAR.conn.val, alcAR.unconn.val), Group = c(rep(1, times = 189), rep(0, times = 36)))
#leveneTest(AR~as.factor(Group), data = data) #no homogeneity of variances
#fligner.test(AR~Group, data = data) #no homogeneity of variances
#t.test(alcAR.conn.val, alcAR.unconn.val) #significantly higher AR of balancing 
#genes in connected populations compared to the AR of these genes in unconnected pops

#but assumptions are violated so we will do wilcoxon rank test
?wilcox.test #default is two sided
wtest.alc <- wilcox.test(alcAR.conn.val, alcAR.unconn.val) #significant differance

#make 1 data frame of mean and SD
datAR <- data.frame(Connected = c('Connected', 'Unconnected'), 
                    Mean = c(mean(alcAR.conn.val), mean(alcAR.unconn.val)), 
                    SD = c(sd(alcAR.conn.val), sd(alcAR.unconn.val)))

#gg boxplot
library(ggplot2)
library(ggsignif)
?geom_signif
ggplot(data = datAR, aes(x=Connected, y = Mean))+
  geom_bar(stat='identity', width = 0.5, fill = '#116E8A')+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width = 0.1)+
  geom_signif(comparisons = list(c('Connected', 'Unconnected')), annotations=c(paste('P = ', round(wtest.alc$p.value, digits = 3))),
              textsize = 4, y_position = c(2.3), vjust = 0.1)+
  theme_classic()+
  ylab('Rarefied Allelic Richness') +
  ggtitle('Alcon') +
  theme(axis.title.x = element_blank())

#Gentian
#select the balancing SNPs from connected alcon populations
genAR.conn[index.bal.gen,]

genAR.conn.val <- vector()
for (i in 1:length(genAR.conn[1,])){
  genAR.conn.val <- append(genAR.conn.val, genAR.conn[index.bal.gen,i])
}

genAR.unconn.val <- vector()
for (i in 1:length(genAR.unconn[1,])){
  genAR.unconn.val <- append(genAR.unconn.val, genAR.unconn[index.bal.gen,i])
}

shapiro.test(genAR.conn.val)
hist(genAR.conn.val, main = 'Gentian - Connected (W = 0.734)')
shapiro.test(genAR.unconn.val) #no homogeneity of variances
hist(genAR.unconn.val, main = 'Gentian - Unconnected (W = 0.781)')

leveneTest(c(genAR.conn.val, genAR.unconn.val)~c(rep('connected', times = length(genAR.conn.val)), 
                                                  rep('unconnected', times = length(genAR.unconn.val))))
#data is homogeneous, yet not normally distributed
#we'll therefore perform a rank test (wilcoxon)
wilcox.test(genAR.conn.val, genAR.unconn.val)

#make 1 data frame of mean and SD
datgAR <- data.frame(Connected = c('Connected', 'Unconnected'), 
                    Mean = c(mean(genAR.conn.val), mean(genAR.unconn.val)), 
                    SD = c(sd(genAR.conn.val), sd(genAR.unconn.val)))

#gg boxplot
ggplot(data = datgAR, aes(x=Connected, y = Mean))+
  geom_bar(stat='identity', width = 0.5, fill = '#116E8A')+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width = 0.1)+
  geom_signif(comparisons = list(c('Connected', 'Unconnected')), annotations=c('NS'),
              textsize = 4, y_position = c(2.8), vjust = 0.05)+
  theme_classic()+
  ylab('Rarefied Allelic Richness') +
  ggtitle('Gentian') +
  theme(axis.title.x = element_blank())



########
ggplot(data = NULL, aes(x=x.a, y = y.ba))+
  geom_point()+
  ylab('Rarefied Allelic Richness')+ xlab = ('Connectivity')+
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
